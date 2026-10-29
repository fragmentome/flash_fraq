"""Admixture design: coverage harvest, feasibility table, mixture plan.

The background fixture is the *real* 44-sample Snyder coverage distribution
(`SNYDER_NATIVE_COVERAGE`, harvested from the mosdepth summaries on the
workstation) rather than a comfortable synthetic one, because the whole design
turns on that distribution's 0.75x floor: a fixture built from round numbers
would never exercise the constraint the module exists to respect.
"""

from pathlib import Path

import pytest
from pandas import read_csv

from flash_fraq.admixture import AdmixtureBuilder, AdmixtureConfig
from flash_fraq.admixture_design import (
    MATCHED_PAIR_SPLIT_GROUP,
    REGIME_MATCHED_PAIR,
    REGIME_REAL_BACKGROUND,
    AdmixtureDesignConfig,
    AdmixtureDesigner,
    AdmixtureFeasibilityAnalyzer,
    MixturePlanner,
    NativeCoverageConfig,
    NativeCoverageHarvester,
)
from flash_fraq.feature_matrix import SPLIT_GROUP_ANY
from flash_fraq.input_manifest import InputManifestBuilder

# Real genome-wide means from /mnt/8tb_hdd/flash_fraq_derived/*/*/coverage, harvested
# on the workstation. Snyder is the full 44; Cristiano is an evenly-spaced 44-point
# subsample of its 503 so the two sources are the same size here while keeping the real
# shape -- including the 0.67x floor and the 13x tail.
SNYDER_NATIVE_COVERAGE = [
    0.75, 0.91, 0.92, 0.92, 0.93, 0.94, 0.96, 0.96, 0.97, 0.97, 0.98, 0.99,
    1.02, 1.05, 1.06, 1.07, 1.07, 1.10, 1.10, 1.14, 1.14, 1.17, 1.18, 1.19,
    1.19, 1.21, 1.24, 1.24, 1.24, 1.29, 1.33, 1.34, 1.35, 1.36, 1.36, 1.37,
    1.41, 1.41, 1.42, 1.42, 1.43, 1.43, 1.58, 1.58,
]
# Cristiano sits a full octave higher (median 2.26x over the full 503), which is what
# makes the two sources drop out of the feasible pool at different target coverages.
CRISTIANO_NATIVE_COVERAGE = [
    0.67, 1.33, 1.53, 1.67, 1.74, 1.78, 1.82, 1.88, 1.91, 1.94, 1.97, 2.00,
    2.04, 2.06, 2.08, 2.11, 2.13, 2.15, 2.17, 2.19, 2.23, 2.25, 2.28, 2.29,
    2.31, 2.35, 2.37, 2.40, 2.43, 2.45, 2.48, 2.56, 2.63, 2.68, 2.73, 2.76,
    2.82, 2.94, 3.12, 3.39, 4.30, 7.95, 10.01, 13.00,
]

FRACTION_LADDER = (0.0, 0.01, 0.05, 0.1, 0.2, 0.5)
CANDIDATE_COVERAGES = (0.1, 0.25, 0.5, 0.75, 1.0)
SPLIT_FRACTIONS = {"train": 0.7, "validation": 0.15, "test": 0.15}


def _mosdepth_summary(coverage: float) -> str:
    return (
        "chrom\tlength\tbases\tmean\tmin\tmax\n"
        f"chr1\t249250621\t240725769\t{coverage}\t0\t19217\n"
        f"total\t3137161264\t3029116188\t{coverage}\t0\t39949\n"
        f"total_region\t3137161264\t3029116188\t{coverage}\t0\t39949\n"
    )


def _write_derived_tree(derived_root: Path, datasets: dict[str, list[float]], incomplete: set[str] = frozenset()):
    for dataset, coverages in datasets.items():
        for index, coverage in enumerate(coverages, start=1):
            sample_id = f"{dataset[:2].upper()}{index:03d}"
            sample_dir = derived_root / dataset / sample_id
            (sample_dir / "coverage").mkdir(parents=True)
            (sample_dir / "coverage" / f"{sample_id}.mosdepth.summary.txt").write_text(_mosdepth_summary(coverage))
            if sample_id not in incomplete:
                (sample_dir / "localfeat_sample_complete.flag").touch()


DATASET_ARGUMENTS = [
    {
        "name": "Snyder_Shendure_WPS",
        "platform": "HiSeq2000",
        "source_group": "snyder",
        "tumor_status": "healthy",
    },
    {
        "name": "Cristiano_Velculescu_DELFI",
        "platform": "HiSeq2500",
        "source_group": "cristiano_delfi",
        "tumor_status": "tumor",
    },
]


@pytest.fixture
def coverage_manifest(tmp_path: Path) -> Path:
    derived_root = tmp_path / "derived"
    _write_derived_tree(
        derived_root,
        {
            "Snyder_Shendure_WPS": SNYDER_NATIVE_COVERAGE,
            "Cristiano_Velculescu_DELFI": CRISTIANO_NATIVE_COVERAGE,
        },
    )
    harvester = NativeCoverageHarvester(
        NativeCoverageConfig.from_sources(
            derived_root=derived_root,
            cram_root=tmp_path / "crams",
            datasets=DATASET_ARGUMENTS,
            output_path=tmp_path / "native_coverage.csv",
        )
    )
    return harvester.write(harvester.harvest())


def _design_config(coverage_manifest: Path, output_root: Path, **overrides) -> AdmixtureDesignConfig:
    arguments = {
        "coverage_manifest_path": coverage_manifest,
        "output_root": output_root,
        "tumor": {
            "sample_id": "HCC1395",
            "platform": "HiSeq",
            "native_coverage": 35.0,
            "bam_path": "/nonexistent/HCC1395.cram",
        },
        "matched_pair_background": {
            "sample_id": "HCC1395BL",
            "platform": "HiSeq",
            "native_coverage": 34.0,
            "bam_path": "/nonexistent/HCC1395BL.cram",
        },
        "fraction_ladder": FRACTION_LADDER,
        "candidate_target_coverages": CANDIDATE_COVERAGES,
        "split_fractions": SPLIT_FRACTIONS,
        "mixtures_per_fraction": 20,
        "matched_pair_replicates_per_fraction": 1,
        "min_background_samples_per_source": 20,
        "min_background_retention_per_source": 0.9,
        "seed": 42,
    }
    return AdmixtureDesignConfig.from_sources(**(arguments | overrides))


def test_harvester_reads_total_row_and_gates_on_the_complete_flag(tmp_path: Path):
    derived_root = tmp_path / "derived"
    _write_derived_tree(
        derived_root,
        {"Snyder_Shendure_WPS": SNYDER_NATIVE_COVERAGE},
        incomplete={"SN001", "SN044"},
    )
    harvester = NativeCoverageHarvester(
        NativeCoverageConfig.from_sources(
            derived_root=derived_root,
            cram_root=tmp_path / "crams",
            datasets=DATASET_ARGUMENTS,
            output_path=tmp_path / "native_coverage.csv",
        )
    )
    entries = harvester.harvest()

    assert len(entries) == len(SNYDER_NATIVE_COVERAGE) - 2
    assert {entry.native_coverage for entry in entries} == set(SNYDER_NATIVE_COVERAGE[1:-1])
    assert all(entry.source_group == "snyder" and entry.platform == "HiSeq2000" for entry in entries)
    assert entries[0].bam_path == tmp_path / "crams" / "Snyder_Shendure_WPS" / "CRAM" / "SN002" / "SN002.sorted.cram"


def _cristiano_phenotype_harvester(tmp_path: Path, phenotype_rows: dict[str, str]) -> NativeCoverageHarvester:
    derived_root = tmp_path / "derived"
    _write_derived_tree(derived_root, {"Cristiano_Velculescu_DELFI": CRISTIANO_NATIVE_COVERAGE})
    phenotype_path = tmp_path / "phenotype.csv"
    phenotype_path.write_text(
        "sample_id,cohort,tumor_status,source\n"
        + "".join(f"{sample_id},x,{status},test\n" for sample_id, status in phenotype_rows.items())
    )
    return NativeCoverageHarvester(
        NativeCoverageConfig.from_sources(
            derived_root=derived_root,
            cram_root=tmp_path / "crams",
            datasets=[DATASET_ARGUMENTS[1] | {"phenotype_path": phenotype_path}],
            output_path=tmp_path / "native_coverage.csv",
        )
    )


def test_background_tumor_status_keeps_only_per_sample_healthy_rows(tmp_path: Path):
    # Odd-numbered Cristiano samples healthy, even ones cancer -- the dataset-level
    # "tumor" must be overridden per sample, then filtered on.
    phenotype = {
        f"CR{index:03d}": "healthy" if index % 2 else "tumor" for index in range(1, len(CRISTIANO_NATIVE_COVERAGE) + 1)
    }
    harvester = _cristiano_phenotype_harvester(tmp_path, phenotype)
    manifest = harvester.write(harvester.harvest())

    config = _design_config(
        manifest,
        tmp_path / "out",
        background_source_groups=["cristiano_delfi"],
        background_tumor_status="healthy",
    )
    background = AdmixtureFeasibilityAnalyzer(config).load_background()

    assert {entry.sample_id for entry in background} == {s for s, status in phenotype.items() if status == "healthy"}
    assert all(entry.tumor_status == "healthy" for entry in background)


def test_harvester_raises_when_the_phenotype_join_leaves_a_sample_unmatched(tmp_path: Path):
    phenotype = {f"CR{index:03d}": "healthy" for index in range(1, len(CRISTIANO_NATIVE_COVERAGE))}  # drops the last
    harvester = _cristiano_phenotype_harvester(tmp_path, phenotype)

    with pytest.raises(ValueError, match="CR044.*missing from"):
        harvester.harvest()


def test_harvested_manifest_round_trips_through_input_manifest(coverage_manifest: Path):
    """The harvest populates ManifestEntry.native_coverage, not a parallel schema."""
    entries = InputManifestBuilder.load(coverage_manifest)

    assert len(entries) == len(SNYDER_NATIVE_COVERAGE) + len(CRISTIANO_NATIVE_COVERAGE)
    assert all(entry.native_coverage is not None for entry in entries)
    assert min(entry.native_coverage for entry in entries) == 0.67


def test_coverage_summary_matches_the_real_snyder_distribution(coverage_manifest: Path):
    entries = [entry for entry in InputManifestBuilder.load(coverage_manifest) if entry.source_group == "snyder"]
    snyder = NativeCoverageHarvester.summarise(entries)["coverage_by_source_group"]["snyder"]

    # Nearest-rank, so every number is a depth some library actually has. The p25/p75
    # match the ticket's interpolated table; the median reads 1.18 rather than its 1.17
    # because 44 is even and nearest-rank does not average the two middle samples.
    assert snyder == {"n": 44, "min": 0.75, "p25": 0.99, "median": 1.18, "p75": 1.36, "max": 1.58, "below_1x": 12}


def test_feasibility_marks_the_low_coverage_sample_infeasible_at_low_f_only(coverage_manifest: Path, tmp_path: Path):
    """The 0.75x floor sample is the one C_target actually has to be chosen around."""
    analyzer = AdmixtureFeasibilityAnalyzer(_design_config(coverage_manifest, tmp_path / "out"))
    table = analyzer.build_table(analyzer.load_background())

    floor_sample = table[(table["healthy_native_coverage"] == 0.75) & (table["target_coverage"] == 1.0)]
    verdicts = floor_sample.set_index("tumor_fraction")["feasible"].to_dict()

    # (1-f)*1.0/0.75 <= 1 first holds at f = 0.25.
    assert verdicts[0.0] is False
    assert verdicts[0.1] is False
    assert verdicts[0.5] is True
    # The tumor side never binds: 0.5*1.0/35 is far below 1 at every rung.
    assert (floor_sample["tumor_downsample_fraction"] <= 1.0).all()


def test_chosen_target_coverage_is_the_deepest_that_keeps_every_source(coverage_manifest: Path, tmp_path: Path):
    analyzer = AdmixtureFeasibilityAnalyzer(_design_config(coverage_manifest, tmp_path / "out"))
    table = analyzer.build_table(analyzer.load_background())
    chosen = analyzer.choose_target_coverage(table)

    # 0.75x is exactly Snyder's floor, so it keeps all 44 at f=0 and loses only
    # Cristiano's single 0.67x sample. The next candidate up, 1.0x, drops all 12 Snyder
    # samples below 1x -- 32/44 = 73%, under the 90% floor -- so 0.75x is the deepest.
    assert chosen == 0.75
    retention = analyzer.retention_by_source(table, chosen)
    assert retention["snyder"] == {"eligible": 44, "total": 44}
    assert retention["cristiano_delfi"] == {"eligible": 43, "total": 44}

    deeper = analyzer.retention_by_source(table, 1.0)
    assert deeper["snyder"]["eligible"] == 32


def test_choose_target_coverage_raises_when_no_candidate_clears_the_floors(coverage_manifest: Path, tmp_path: Path):
    config = _design_config(
        coverage_manifest, tmp_path / "out", candidate_target_coverages=[2.0, 3.0], min_background_samples_per_source=5
    )
    analyzer = AdmixtureFeasibilityAnalyzer(config)
    with pytest.raises(ValueError, match="no candidate target coverage"):
        analyzer.choose_target_coverage(analyzer.build_table(analyzer.load_background()))


def test_no_subject_spans_two_split_groups(coverage_manifest: Path, tmp_path: Path):
    config = _design_config(coverage_manifest, tmp_path / "out")
    planner = MixturePlanner(config)
    background = AdmixtureFeasibilityAnalyzer(config).load_background()
    splits = planner.assign_split_groups(background)

    assert list(splits.columns) == ["subject_id", "split_group"]
    assert splits["subject_id"].is_unique
    assert set(splits["split_group"]) == {"train", "validation", "test", SPLIT_GROUP_ANY}
    # Both sources reach all three splits -- a single-source split would turn a
    # validation score into a cross-cohort generalisation measurement.
    subject_to_source = {entry.subject_id: entry.source_group for entry in background}
    for split_group in ("train", "validation", "test"):
        members = splits.loc[splits["split_group"] == split_group, "subject_id"]
        assert {subject_to_source[subject] for subject in members} == {"snyder", "cristiano_delfi"}


def test_plan_rows_all_pass_the_builders_own_feasibility_check(coverage_manifest: Path, tmp_path: Path):
    summary = AdmixtureDesigner(_design_config(coverage_manifest, tmp_path / "out")).run()
    plan = read_csv(tmp_path / "out" / "admixture_plan.csv")

    # real_background: 20 per rung; matched_pair: 1 replicate per rung (test config).
    assert len(plan) == summary["n_mixtures"] == len(FRACTION_LADDER) * 20 + len(FRACTION_LADDER) * 1
    for row in plan.to_dict("records"):
        builder = AdmixtureBuilder(
            AdmixtureConfig.from_sources(
                tumor_bam=Path(row["tumor_bam_path"]),
                healthy_bam=Path(row["healthy_bam_path"]),
                reference_fasta=Path("unused"),
                snv_vaf_paths={"unused": Path("unused")},
                output_dir=Path("unused"),
                tumor_fraction=row["tumor_fraction"],
                target_coverage=row["target_coverage"],
                tumor_native_coverage=row["tumor_native_coverage"],
                healthy_native_coverage=row["healthy_native_coverage"],
            )
        )
        assert builder.compute_downsample_fractions() == pytest.approx(
            (row["tumor_downsample_fraction"], row["healthy_downsample_fraction"])
        )


def test_plan_carries_every_field_admixture_config_needs(coverage_manifest: Path, tmp_path: Path):
    AdmixtureDesigner(_design_config(coverage_manifest, tmp_path / "out")).run()
    plan = read_csv(tmp_path / "out" / "admixture_plan.csv")

    required = {
        "sample_name", "tumor_source_id", "healthy_source_id", "tumor_subject_id", "healthy_subject_id",
        "tumor_fraction", "target_coverage", "tumor_native_coverage", "healthy_native_coverage", "seed",
        "platform_tumor", "platform_healthy",
    }
    assert required <= set(plan.columns)
    assert plan["sample_name"].is_unique
    assert plan["seed"].is_unique
    assert plan["tumor_source_id"].unique().tolist() == ["HCC1395"]


def test_source_proportions_are_orthogonal_to_the_mixing_ratio(coverage_manifest: Path, tmp_path: Path):
    """The §6 constraint: no f may over-represent either cohort's backgrounds.

    Computed over real_background rows only -- matched_pair has one background
    subject by construction, so it is evaluated separately (never pooled), per
    docs/admixture-design.md's reopened decision.
    """
    summary = AdmixtureDesigner(_design_config(coverage_manifest, tmp_path / "out")).run()
    plan = read_csv(tmp_path / "out" / "admixture_plan.csv")
    real_background = plan[plan["regime"] == REGIME_REAL_BACKGROUND]

    pooled = real_background["healthy_source_group"].value_counts(normalize=True)
    for _, group in real_background.groupby("tumor_fraction"):
        per_fraction = group["healthy_source_group"].value_counts(normalize=True).reindex(pooled.index).fillna(0.0)
        assert (per_fraction - pooled).abs().max() == pytest.approx(0.0, abs=1e-9)
    assert summary["max_source_proportion_deviation"] == pytest.approx(0.0, abs=1e-9)

    # Same must hold per split group, or a split could be all-Snyder at low f.
    for _, group in real_background.groupby("split_group"):
        assert group["tumor_fraction"].nunique() == len(FRACTION_LADDER)


def test_matched_pair_regime_is_unsplit_and_replicated_across_the_ladder(coverage_manifest: Path, tmp_path: Path):
    config = _design_config(coverage_manifest, tmp_path / "out")
    in_memory = MixturePlanner(config).build_matched_pair_plan(1.0)
    assert set(in_memory["split_group"]) == {MATCHED_PAIR_SPLIT_GROUP}

    AdmixtureDesigner(config).run()
    plan = read_csv(tmp_path / "out" / "admixture_plan.csv")
    matched_pair = plan[plan["regime"] == REGIME_MATCHED_PAIR]

    assert len(matched_pair) == len(FRACTION_LADDER)
    assert set(matched_pair["tumor_fraction"]) == set(FRACTION_LADDER)
    assert matched_pair["healthy_source_id"].unique().tolist() == ["HCC1395BL"]
    assert matched_pair["tumor_source_id"].unique().tolist() == ["HCC1395"]
    # Unsplit: no split_group at all, not even the SPLIT_GROUP_ANY sentinel. Written
    # as "" (MATCHED_PAIR_SPLIT_GROUP) and round-trips through CSV as NaN.
    assert matched_pair["split_group"].isna().all()


def test_matched_pair_plan_raises_when_the_pair_is_infeasible(coverage_manifest: Path, tmp_path: Path):
    config = _design_config(
        coverage_manifest,
        tmp_path / "out",
        matched_pair_background={
            "sample_id": "HCC1395BL",
            "platform": "HiSeq",
            "native_coverage": 0.5,
            "bam_path": "/nonexistent/HCC1395BL.cram",
        },
        candidate_target_coverages=[1.0],
    )
    planner = MixturePlanner(config)
    with pytest.raises(ValueError, match="matched-pair background"):
        planner.build_matched_pair_plan(1.0)


def test_regime_column_partitions_the_plan(coverage_manifest: Path, tmp_path: Path):
    AdmixtureDesigner(_design_config(coverage_manifest, tmp_path / "out")).run()
    plan = read_csv(tmp_path / "out" / "admixture_plan.csv")

    assert set(plan["regime"]) == {REGIME_MATCHED_PAIR, REGIME_REAL_BACKGROUND}
    assert plan["seed"].is_unique


def test_plan_is_byte_identical_under_a_fixed_seed(coverage_manifest: Path, tmp_path: Path):
    outputs = []
    for run_index in range(2):
        output_root = tmp_path / f"run{run_index}"
        AdmixtureDesigner(_design_config(coverage_manifest, output_root)).run()
        outputs.append(
            [(output_root / name).read_bytes() for name in ("admixture_plan.csv", "split_assignment.csv")]
        )
    assert outputs[0] == outputs[1]

    other_seed = tmp_path / "other"
    AdmixtureDesigner(_design_config(coverage_manifest, other_seed, seed=7)).run()
    assert (other_seed / "admixture_plan.csv").read_bytes() != outputs[0][0]


def test_split_assignment_feeds_the_feature_matrix_schema(coverage_manifest: Path, tmp_path: Path):
    """split_assignment.csv is consumed by feature_matrix.assign_split_groups verbatim.

    Scoped to real_background: matched_pair's background subject (HCC1395BL) is
    deliberately absent from split_assignment.csv, since the regime is unsplit. Making
    that combination resolve sensibly through assign_split_groups is MAT-17's job, not
    this design's -- see the module's "Reopened decision" note.
    """
    from flash_fraq.feature_matrix import FeatureMatrixAssembler, FeatureMatrixConfig

    AdmixtureDesigner(_design_config(coverage_manifest, tmp_path / "out")).run()
    plan = read_csv(tmp_path / "out" / "admixture_plan.csv")
    real_background = plan[plan["regime"] == REGIME_REAL_BACKGROUND]
    splits = read_csv(tmp_path / "out" / "split_assignment.csv")

    metadata = real_background.rename(columns={"sample_name": "sample_id"})[
        ["sample_id", "tumor_subject_id", "healthy_subject_id"]
    ]
    assembler = FeatureMatrixAssembler(FeatureMatrixConfig.from_sources(tmp_path, tmp_path / 'matrix.csv'))
    resolved = assembler.assign_split_groups(metadata, splits)

    assert resolved["split_group"].tolist() == real_background["split_group"].tolist()


def test_zero_matched_pair_replicates_produces_real_background_only(coverage_manifest: Path, tmp_path: Path):
    """`matched_pair_replicates_per_fraction=0` is how the reshaped-tumor rebuild
    (docs/todo/todo/reshaped-mixture-downstream-rerun.md) builds *only*
    `real_background` rows against the reshaped pooled NV_1+NV_3 tumor source, while
    reusing the already-built, deliberately-unreshaped `matched_pair` mixtures from
    the original build -- fragment-size-reshaping.md's finding that reshaping does
    not apply to a sheared-vs-sheared regime. This must not raise (matched-pair
    feasibility is still checked once at the binding fraction, per
    `build_matched_pair_plan`'s docstring) and must not otherwise disturb the
    `real_background` rows or the plan's schema.
    """
    config = _design_config(coverage_manifest, tmp_path / "out", matched_pair_replicates_per_fraction=0)
    summary = AdmixtureDesigner(config).run()
    plan = read_csv(tmp_path / "out" / "admixture_plan.csv")

    assert set(plan["regime"]) == {REGIME_REAL_BACKGROUND}
    assert len(plan) == summary["n_mixtures"] == len(FRACTION_LADDER) * 20
    assert summary["n_mixtures_by_regime"].get(REGIME_MATCHED_PAIR, 0) == 0
