"""MAT-23 Block C: lane -> library CRAM merge and the merged-library coverage
report.

The merge tests build real (tiny) CRAMs with samtools against a 2-contig
synthetic FASTA and then merge them for real, rather than stubbing subprocess:
the properties under test -- that the merged header advertises one sample, that
every lane's reads survive, that a mismatched reference is refused -- are
properties of samtools' behaviour, and a stub would only assert that the
arguments were assembled the way the test author imagined.
"""

import csv
import subprocess
from pathlib import Path

import pytest

from flash_fraq.hcc1395_merge import (
    ARM_NORMAL,
    ARM_TUMOR,
    CramLibraryMerger,
    HCC1395MergeConfig,
    LibraryMergePlanner,
    MergedCoverageReporter,
    MergeUnit,
    merged_sample_id,
)

PAIRS_HEADER = "run_accession,library_name,prep,site,arm,variant,read_count,base_count\n"

PAIRS_ROWS = (
    # WGS_IL_T_2: two top-up lanes of one library -- the multi-lane merge case.
    "SRR0001,WGS_IL_T_2,WGS,IL,T,2,100,30000\n"
    "SRR0002,WGS_IL_T_2,WGS,IL,T,2,100,30000\n"
    # Its matched normal, one lane.
    "SRR0003,WGS_IL_N_2,WGS,IL,N,2,100,30000\n"
    # Same site, different replicate: must NOT be merged into WGS_IL_T_2.
    "SRR0004,WGS_IL_T_3,WGS,IL,T,3,100,30000\n"
    # Different site: must NOT be merged into any IL library.
    "SRR0005,WGS_NV_T_1,WGS,NV,T,1,100,30000\n"
    # A prep family outside `preps`, and a pre-pooled mixture with no arm.
    "SRR0006,CHR_EA_N_1,CHR,EA,N,1,100,30000\n"
    "SRR0007,SPP_GT_1-0_1,SPP,GT,,1-0_1,100,30000\n"
)

CONTIGS = (("chr1", 2000), ("chr2", 1000))


def _write_pairs_manifest(tmp_path: Path, rows: str = PAIRS_ROWS) -> Path:
    path = tmp_path / "hcc1395_library_pairs.csv"
    path.write_text(PAIRS_HEADER + rows)
    return path


def _write_reference(tmp_path: Path, contigs=CONTIGS, name: str = "ref.fa") -> Path:
    """A synthetic FASTA plus its .fai, so CRAM has a reference to encode against."""
    fasta = tmp_path / name
    with open(fasta, "w") as handle:
        for contig, length in contigs:
            handle.write(f">{contig}\n")
            sequence = ("ACGT" * (length // 4 + 1))[:length]
            for start in range(0, length, 60):
                handle.write(sequence[start : start + 60] + "\n")
    subprocess.run(["samtools", "faidx", str(fasta)], check=True)
    return fasta


def _write_lane_cram(
    cram_root: Path,
    dataset: str,
    run_accession: str,
    reference: Path,
    contigs=CONTIGS,
    n_reads: int = 4,
    start: int = 100,
) -> Path:
    """One tiny coordinate-sorted CRAM + .crai, laid out exactly like the mirror."""
    sam_lines = ["@HD\tVN:1.6\tSO:coordinate"]
    sam_lines += [f"@SQ\tSN:{contig}\tLN:{length}" for contig, length in contigs]
    sam_lines.append(f"@RG\tID:{run_accession}\tSM:{run_accession}\tPL:ILLUMINA")
    for index in range(n_reads):
        position = start + index * 10
        sam_lines.append(
            f"{run_accession}.read{index}\t0\t{contigs[0][0]}\t{position}\t60\t8M\t*\t0\t0\t"
            f"ACGTACGT\tIIIIIIII\tRG:Z:{run_accession}"
        )
    sam_path = cram_root / f"{run_accession}.sam"
    sam_path.parent.mkdir(parents=True, exist_ok=True)
    sam_path.write_text("\n".join(sam_lines) + "\n")

    cram_dir = cram_root / dataset / "CRAM" / run_accession
    cram_dir.mkdir(parents=True, exist_ok=True)
    cram_path = cram_dir / f"{run_accession}.sorted.cram"
    subprocess.run(
        ["samtools", "view", "-C", "-T", str(reference), "-o", str(cram_path), str(sam_path)],
        check=True,
    )
    subprocess.run(["samtools", "index", str(cram_path)], check=True)
    sam_path.unlink()
    return cram_path


def _config(tmp_path: Path, reference: Path, **overrides) -> HCC1395MergeConfig:
    # The default manifest is written lazily: writing it eagerly would clobber a
    # caller-supplied one, since both land at the same tmp_path filename.
    if "pairs_manifest" not in overrides:
        overrides["pairs_manifest"] = _write_pairs_manifest(tmp_path)
    return HCC1395MergeConfig.from_sources(
        pairs_manifest=overrides.pop("pairs_manifest"),
        cram_root=overrides.pop("cram_root", tmp_path / "crams"),
        reference_fasta=reference,
        **overrides,
    )


def _header(cram_path: Path, reference: Path) -> list[str]:
    result = subprocess.run(
        ["samtools", "view", "-H", "--reference", str(reference), str(cram_path)],
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout.splitlines()


# --- planning ---------------------------------------------------------------


def test_plan_groups_lanes_by_library_not_by_site_or_cell_line(tmp_path):
    """The merge-unit decision, asserted: lanes of one library merge; different
    replicates, different sites and different arms stay apart."""
    planner = LibraryMergePlanner(_config(tmp_path, _write_reference(tmp_path)))
    units = {unit.sample_id: unit for unit in planner.plan()}

    assert set(units) == {
        "HCC1395_WGS_IL_2",
        "HCC1395_WGS_IL_3",
        "HCC1395_WGS_NV_1",
        "HCC1395BL_WGS_IL_2",
    }
    assert units["HCC1395_WGS_IL_2"].run_accessions == ("SRR0001", "SRR0002")
    # Same site, other replicate -- its own unit, not folded into IL_2.
    assert units["HCC1395_WGS_IL_3"].run_accessions == ("SRR0004",)
    # Different site -- never merged with IL.
    assert units["HCC1395_WGS_NV_1"].run_accessions == ("SRR0005",)


def test_plan_excludes_other_preps_and_armless_spp_rows(tmp_path):
    planner = LibraryMergePlanner(_config(tmp_path, _write_reference(tmp_path)))
    accessions = {run for unit in planner.plan() for run in unit.run_accessions}
    # CHR (prep outside `preps`) and SPP (no T/N arm) are both dropped.
    assert "SRR0006" not in accessions
    assert "SRR0007" not in accessions


def test_plan_carries_arm_metadata_and_subject_ids(tmp_path):
    planner = LibraryMergePlanner(_config(tmp_path, _write_reference(tmp_path)))
    units = {unit.sample_id: unit for unit in planner.plan()}

    tumor = units["HCC1395_WGS_IL_2"]
    normal = units["HCC1395BL_WGS_IL_2"]
    assert (tumor.subject_id, tumor.source_group, tumor.tumor_status) == ("HCC1395", "hcc1395", "tumor")
    assert (normal.subject_id, normal.source_group, normal.tumor_status) == (
        "HCC1395BL",
        "hcc1395bl",
        "healthy",
    )
    assert tumor.arm == ARM_TUMOR and normal.arm == ARM_NORMAL
    assert tumor.library_name == "WGS_IL_T_2"
    assert merged_sample_id("HCC1395", "WGS", "IL", "2") == "HCC1395_WGS_IL_2"


def test_plan_available_keeps_only_lanes_with_cram_and_index(tmp_path):
    """MAT-23 closed with 11 runs unretrieved, so the manifest is a superset of
    what is banked. A lane whose .crai has not landed is still in flight."""
    reference = _write_reference(tmp_path)
    cram_root = tmp_path / "crams"
    _write_lane_cram(cram_root, "HCC1395_SEQC2", "SRR0001", reference)
    _write_lane_cram(cram_root, "HCC1395_SEQC2", "SRR0002", reference)
    _write_lane_cram(cram_root, "HCC1395_SEQC2", "SRR0003", reference)
    # SRR0003's index removed -> the lane, and with it its single-lane library,
    # drops out. SRR0004/5 were never retrieved at all.
    Path(f"{cram_root / 'HCC1395_SEQC2' / 'CRAM' / 'SRR0003' / 'SRR0003.sorted.cram'}.crai").unlink()

    planner = LibraryMergePlanner(_config(tmp_path, reference, cram_root=cram_root))
    available = {unit.sample_id: unit for unit in planner.plan_available()}

    assert set(available) == {"HCC1395_WGS_IL_2"}
    assert available["HCC1395_WGS_IL_2"].run_accessions == ("SRR0001", "SRR0002")


def test_plan_paired_available_drops_libraries_with_no_partner_arm(tmp_path):
    """WGS_IL_T_3 and WGS_NV_T_1 are on disk but have no matched normal in
    PAIRS_ROWS -- the report script must not expect coverage for them any
    more than `rule hcc1395_merge_library` (pairs_only) merges them."""
    reference = _write_reference(tmp_path)
    cram_root = tmp_path / "crams"
    for run_accession in ("SRR0001", "SRR0002", "SRR0003", "SRR0004", "SRR0005"):
        _write_lane_cram(cram_root, "HCC1395_SEQC2", run_accession, reference)

    planner = LibraryMergePlanner(_config(tmp_path, reference, cram_root=cram_root))
    assert {unit.sample_id for unit in planner.plan_available()} == {
        "HCC1395_WGS_IL_2",
        "HCC1395BL_WGS_IL_2",
        "HCC1395_WGS_IL_3",
        "HCC1395_WGS_NV_1",
    }
    assert {unit.sample_id for unit in planner.plan_paired_available()} == {
        "HCC1395_WGS_IL_2",
        "HCC1395BL_WGS_IL_2",
    }


def test_matched_pairs_requires_the_same_replicate(tmp_path):
    """The NS defect MAT-23 documented (WGS_NS_N_8 x WGS_NS_T_7) must not be
    silently paired: replicate is part of the pairing key."""
    manifest = _write_pairs_manifest(
        tmp_path,
        "SRR1001,WGS_NS_T_7,WGS,NS,T,7,100,30000\n" "SRR1002,WGS_NS_N_8,WGS,NS,N,8,100,30000\n",
    )
    planner = LibraryMergePlanner(_config(tmp_path, _write_reference(tmp_path), pairs_manifest=manifest))
    assert planner.matched_pairs(planner.plan()) == []


def test_matched_pairs_returns_same_site_same_replicate_arms(tmp_path):
    planner = LibraryMergePlanner(_config(tmp_path, _write_reference(tmp_path)))
    pairs = planner.matched_pairs(planner.plan())
    assert [(tumor.sample_id, normal.sample_id) for tumor, normal in pairs] == [
        ("HCC1395_WGS_IL_2", "HCC1395BL_WGS_IL_2")
    ]


def test_plan_rejects_a_manifest_missing_a_column(tmp_path):
    manifest = tmp_path / "bad.csv"
    manifest.write_text("run_accession,library_name\nSRR0001,WGS_IL_T_2\n")
    planner = LibraryMergePlanner(_config(tmp_path, _write_reference(tmp_path), pairs_manifest=manifest))
    with pytest.raises(ValueError, match="missing column"):
        planner.plan()


def test_merge_unit_paths_follow_the_mirror_layout(tmp_path):
    planner = LibraryMergePlanner(_config(tmp_path, _write_reference(tmp_path)))
    unit = next(unit for unit in planner.plan() if unit.sample_id == "HCC1395_WGS_IL_2")
    root = tmp_path / "crams"

    assert unit.merged_cram_path(root, "HCC1395_SEQC2") == (
        root / "HCC1395_SEQC2" / "CRAM" / "HCC1395_WGS_IL_2" / "HCC1395_WGS_IL_2.sorted.cram"
    )
    assert unit.lane_cram_paths(root, "HCC1395_SEQC2")[0] == (
        root / "HCC1395_SEQC2" / "CRAM" / "SRR0001" / "SRR0001.sorted.cram"
    )
    assert unit.base_count == 60000
    assert unit.approximate_coverage(genome_size_bp=30000) == pytest.approx(2.0)


# --- merging ----------------------------------------------------------------


def _merged(tmp_path: Path, **overrides) -> tuple[Path, Path, MergeUnit, HCC1395MergeConfig]:
    reference = _write_reference(tmp_path)
    cram_root = tmp_path / "crams"
    for run_accession in ("SRR0001", "SRR0002"):
        _write_lane_cram(cram_root, "HCC1395_SEQC2", run_accession, reference)
    config = _config(tmp_path, reference, cram_root=cram_root, threads=1, **overrides)
    unit = next(
        unit
        for unit in LibraryMergePlanner(config).plan_available()
        if unit.sample_id == "HCC1395_WGS_IL_2"
    )
    return CramLibraryMerger(config).merge(unit), reference, unit, config


def test_merge_produces_an_indexed_quickcheck_clean_cram(tmp_path):
    merged_cram, reference, _, config = _merged(tmp_path)

    assert merged_cram.exists()
    assert Path(f"{merged_cram}.crai").exists()
    assert subprocess.run(["samtools", "quickcheck", str(merged_cram)]).returncode == 0
    assert CramLibraryMerger(config).quickcheck(merged_cram)


def test_merge_keeps_every_lane_read(tmp_path):
    merged_cram, reference, unit, config = _merged(tmp_path)
    merger = CramLibraryMerger(config)

    lane_total = sum(
        merger.flagstat_total_reads(path)
        for path in unit.lane_cram_paths(config.cram_root, config.dataset)
    )
    assert lane_total == 8
    assert merger.flagstat_total_reads(merged_cram) == lane_total


def test_merge_header_declares_one_sample_and_keeps_lane_read_groups(tmp_path):
    """Each lane was aligned with ID:<run> SM:<run>, so an unrewritten merge
    would advertise two samples in one file."""
    merged_cram, reference, _, _ = _merged(tmp_path)
    read_groups = [line for line in _header(merged_cram, reference) if line.startswith("@RG")]

    assert len(read_groups) == 2
    assert {field for line in read_groups for field in line.split("\t") if field.startswith("ID:")} == {
        "ID:SRR0001",
        "ID:SRR0002",
    }
    assert {field for line in read_groups for field in line.split("\t") if field.startswith("SM:")} == {
        "SM:HCC1395_WGS_IL_2"
    }


def test_merged_reads_still_resolve_their_read_group(tmp_path):
    """SM was rewritten but ID was not, so every read's RG:Z tag must still name
    a declared group -- otherwise lanes stop being separable."""
    merged_cram, reference, _, _ = _merged(tmp_path)
    result = subprocess.run(
        ["samtools", "view", "--reference", str(reference), str(merged_cram)],
        check=True,
        capture_output=True,
        text=True,
    )
    tags = {
        field
        for line in result.stdout.splitlines()
        for field in line.split("\t")
        if field.startswith("RG:Z:")
    }
    assert tags == {"RG:Z:SRR0001", "RG:Z:SRR0002"}


def test_merge_verifies_read_counts_when_asked(tmp_path):
    merged_cram, _, _, _ = _merged(tmp_path, verify_read_counts=True)
    assert merged_cram.exists()


def test_merge_of_a_single_lane_library_still_indexes_and_verifies(tmp_path):
    """Most banked libraries are one lane; the merge path must not need two."""
    reference = _write_reference(tmp_path)
    cram_root = tmp_path / "crams"
    _write_lane_cram(cram_root, "HCC1395_SEQC2", "SRR0003", reference)
    config = _config(tmp_path, reference, cram_root=cram_root, threads=1, verify_read_counts=True)
    unit = next(
        unit for unit in LibraryMergePlanner(config).plan_available() if unit.arm == ARM_NORMAL
    )

    merged_cram = CramLibraryMerger(config).merge(unit)
    assert merged_cram.name == "HCC1395BL_WGS_IL_2.sorted.cram"
    assert CramLibraryMerger(config).flagstat_total_reads(merged_cram) == 4


def test_merge_refuses_lanes_aligned_to_a_different_reference(tmp_path):
    """samtools merge would concatenate the two @SQ blocks into a structurally
    valid, biologically meaningless file."""
    reference = _write_reference(tmp_path)
    other_contigs = (("chr1", 2000), ("chr3", 500))
    other_reference = _write_reference(tmp_path, contigs=other_contigs, name="other.fa")
    cram_root = tmp_path / "crams"
    _write_lane_cram(cram_root, "HCC1395_SEQC2", "SRR0001", reference)
    _write_lane_cram(cram_root, "HCC1395_SEQC2", "SRR0002", other_reference, contigs=other_contigs)

    config = _config(tmp_path, reference, cram_root=cram_root, threads=1)
    unit = next(
        unit
        for unit in LibraryMergePlanner(config).plan_available()
        if unit.sample_id == "HCC1395_WGS_IL_2"
    )
    with pytest.raises(ValueError, match="not aligned to the same reference"):
        CramLibraryMerger(config).merge(unit)


def test_merge_raises_on_an_absent_lane_cram(tmp_path):
    reference = _write_reference(tmp_path)
    config = _config(tmp_path, reference, cram_root=tmp_path / "crams", threads=1)
    unit = next(
        unit for unit in LibraryMergePlanner(config).plan() if unit.sample_id == "HCC1395_WGS_IL_2"
    )
    with pytest.raises(FileNotFoundError):
        CramLibraryMerger(config).merge(unit)


def test_rewrite_sample_adds_sm_when_the_lane_had_none(tmp_path):
    rewrite = CramLibraryMerger._rewrite_sample
    assert rewrite("@RG\tID:SRR0001\tPL:ILLUMINA", "S") == "@RG\tID:SRR0001\tPL:ILLUMINA\tSM:S"
    assert rewrite("@RG\tID:SRR0001\tSM:SRR0001", "S") == "@RG\tID:SRR0001\tSM:S"


# --- merged-library coverage report -----------------------------------------


def _write_mosdepth_summary(path: Path, total_mean: float, region_mean: float) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "chrom\tlength\tbases\tmean\tmin\tmax\n"
        f"chr1\t2000\t100\t{total_mean}\t0\t99\n"
        f"total\t3000\t150\t{total_mean}\t0\t99\n"
        f"total_region\t3000\t150\t{region_mean}\t0\t99\n"
    )


def test_coverage_report_reads_the_total_row_not_total_region(tmp_path):
    """Same field as NativeCoverageHarvester.parse_mosdepth_summary, which is
    what makes these numbers comparable with native_coverage.csv's 547 rows."""
    reference = _write_reference(tmp_path)
    cram_root = tmp_path / "crams"
    derived_root = tmp_path / "derived"
    for run_accession in ("SRR0001", "SRR0002", "SRR0003"):
        _write_lane_cram(cram_root, "HCC1395_SEQC2", run_accession, reference)

    config = _config(tmp_path, reference, cram_root=cram_root)
    units = LibraryMergePlanner(config).plan_available()
    reporter = MergedCoverageReporter(config)
    _write_mosdepth_summary(
        reporter.summary_path(derived_root, config.dataset, "HCC1395_WGS_IL_2"), 41.75, 0.0
    )
    _write_mosdepth_summary(
        reporter.summary_path(derived_root, config.dataset, "HCC1395BL_WGS_IL_2"), 37.5, 0.0
    )

    rows = {row["sample_id"]: row for row in reporter.rows(units, derived_root)}
    assert rows["HCC1395_WGS_IL_2"]["native_coverage"] == pytest.approx(41.75)
    assert rows["HCC1395BL_WGS_IL_2"]["native_coverage"] == pytest.approx(37.5)
    assert rows["HCC1395_WGS_IL_2"]["n_lanes"] == 2
    assert rows["HCC1395_WGS_IL_2"]["run_accessions"] == "SRR0001;SRR0002"
    assert rows["HCC1395BL_WGS_IL_2"]["tumor_status"] == "healthy"


def test_coverage_report_summary_path_matches_the_localfeat_layout(tmp_path):
    """NativeCoverageHarvester walks <derived>/<dataset>/<sample>/coverage/; a
    different layout here would make the merged libraries unharvestable."""
    config = _config(tmp_path, _write_reference(tmp_path))
    path = MergedCoverageReporter(config).summary_path(Path("/derived"), "HCC1395_SEQC2", "HCC1395_WGS_IL_2")
    assert path == Path(
        "/derived/HCC1395_SEQC2/HCC1395_WGS_IL_2/coverage/HCC1395_WGS_IL_2.mosdepth.summary.txt"
    )


def test_coverage_report_writes_a_csv_and_a_pasteable_yaml_block(tmp_path):
    reference = _write_reference(tmp_path)
    cram_root = tmp_path / "crams"
    derived_root = tmp_path / "derived"
    _write_lane_cram(cram_root, "HCC1395_SEQC2", "SRR0003", reference)

    config = _config(tmp_path, reference, cram_root=cram_root)
    units = LibraryMergePlanner(config).plan_available()
    reporter = MergedCoverageReporter(config)
    _write_mosdepth_summary(
        reporter.summary_path(derived_root, config.dataset, "HCC1395BL_WGS_IL_2"), 37.5, 0.0
    )

    rows = reporter.rows(units, derived_root)
    output_path = reporter.write(rows, tmp_path / "out" / "coverage.csv")
    with open(output_path, newline="") as handle:
        written = list(csv.DictReader(handle))

    assert [row["sample_id"] for row in written] == ["HCC1395BL_WGS_IL_2"]
    assert written[0]["native_coverage"] == "37.5"
    assert written[0]["cram_path"].endswith("HCC1395BL_WGS_IL_2/HCC1395BL_WGS_IL_2.sorted.cram")
    assert "native_coverage: 37.50" in reporter.config_yaml_block(rows)


def test_coverage_report_raises_on_a_summary_with_no_total_row(tmp_path):
    reference = _write_reference(tmp_path)
    cram_root = tmp_path / "crams"
    derived_root = tmp_path / "derived"
    _write_lane_cram(cram_root, "HCC1395_SEQC2", "SRR0003", reference)

    config = _config(tmp_path, reference, cram_root=cram_root)
    units = LibraryMergePlanner(config).plan_available()
    reporter = MergedCoverageReporter(config)
    summary = reporter.summary_path(derived_root, config.dataset, "HCC1395BL_WGS_IL_2")
    summary.parent.mkdir(parents=True, exist_ok=True)
    summary.write_text("chrom\tlength\tbases\tmean\tmin\tmax\nchr1\t2000\t100\t1.0\t0\t9\n")

    with pytest.raises(ValueError, match="no 'total' row"):
        reporter.rows(units, derived_root)


# --- the real manifest ------------------------------------------------------


def test_real_library_pairs_manifest_plans_the_expected_sites(tmp_path):
    """Guards the decomposition against a manifest edit: every WGS library in
    manifests/ must land in exactly one merge unit, tagged by its own site."""
    manifest = Path(__file__).resolve().parents[1] / "manifests" / "hcc1395" / "hcc1395_library_pairs.csv"
    if not manifest.exists():
        pytest.skip("manifests/hcc1395/hcc1395_library_pairs.csv not in this checkout")

    config = _config(tmp_path, _write_reference(tmp_path), pairs_manifest=manifest)
    units = LibraryMergePlanner(config).plan()

    assert {unit.site for unit in units} == {"EA", "FD", "IL", "LL", "NC", "NS", "NV"}
    assert {unit.subject_id for unit in units} == {"HCC1395", "HCC1395BL"}
    # No run may be claimed by two libraries -- that would double-count reads.
    accessions = [run for unit in units for run in unit.run_accessions]
    assert len(accessions) == len(set(accessions))
    # Every unit's lanes share one library_name, which is the merge-unit claim.
    assert all(len({lane.library_name for lane in unit.lanes}) == 1 for unit in units)
