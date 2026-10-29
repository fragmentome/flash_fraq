"""Admixture *design*: native-coverage harvest, feasibility table, mixture plan.

`admixture.py` knows how to build one mixed sample once someone tells it which two
libraries to mix, at what `f`, and to what depth. This module is the part that
decides *what to build*, and it is codeable before HCC1395 exists because none of it
touches the tumor reads -- the tumor enters as one sample id and one native-coverage
number (phase2 §7: 35-47x), which stays a declared placeholder until MAT-23 lands.

Three stages, three files:

1. `NativeCoverageHarvester` -> `native_coverage.csv`. Parses the `total` row of every
   `mosdepth.summary.txt` the local_features workflow already wrote and emits it as an
   `input_manifest` CSV, populating `ManifestEntry.native_coverage` rather than
   inventing a parallel schema. `AdmixtureConfig` requires both native coverages and
   nothing else produces either.
2. `AdmixtureFeasibilityAnalyzer` -> `admixture_feasibility.csv`. One row per
   (background sample, `f`, candidate `C_target`) with both phase2 §7 step-1
   downsample fractions and whether the pair clears the `<= 1` rule. The arithmetic is
   imported from `admixture.downsample_fractions`, not restated here.
3. `MixturePlanner` -> `admixture_plan.csv` + `split_assignment.csv`. Chooses
   `C_target` off the feasibility table, assigns `split_group` at source-subject level,
   and pairs each `f` with background partners under the orthogonality constraint below.

The plan is two regimes, tagged by the `regime` column and never pooled (see
"Reopened decision" in docs/admixture-design.md): `matched_pair` mixes HCC1395 against
its own matched normal HCC1395BL for a clean, unconfounded `f` ladder and carries no
`split_group`, since one background subject makes subject-level splitting impossible;
`real_background` is everything below -- the multi-subject pool, its orthogonality
constraint, and its source-subject-level splits.

**The binding constraint.** The healthy-side downsample fraction is
`(1-f)*C_target/C_healthy`, which is largest at `f=0`; so `C_target` is capped by the
*lowest-coverage* background sample the design wants to keep. Measured over the 547
samples on disk, Snyder runs 0.75-1.58x (median 1.17x) and Cristiano 0.67-13.0x
(median 2.26x), so `C_target` is set by the background pool, never by the tumor.

**The orthogonality constraint.** The phase2 §6 gate run in
`docs/cohort-qc-baseline.md` separates Snyder from Cristiano at accuracy 1.000 against
a 0.890 baseline: source membership is unambiguously legible in the features. If the
background pool composition drifted with `f`, a model could hit the regression target
by reading batch identity instead of tumor signal. The plan therefore draws every `f`
from **one pool -- the samples feasible at every `f` in the ladder** -- and allocates
each `f` across (source, split) cells in the same proportions. Per-`f` source
proportions are then equal to the pooled proportions by construction, not by luck.

Note on "healthy": unfiltered, the background pool is every completed cohort sample,
and many of those libraries are cancer-patient plasma (39 of the 44 Snyder samples are
the IC0x lung-cancer series; 261 of the 503 Cristiano DELFI samples are from cancer
patients, per `manifests/cristiano_delfi_phenotype.csv`). A cancer-patient background
carries its own ctDNA, so the true fraction is `f + (1-f)*TF_bg`. The HCC1395-SNV label
cannot see `TF_bg`, but the CNV and fragment features can. `background_tumor_status`
restricts the pool to samples whose per-sample `tumor_status` matches (`healthy`);
the harvester fills that column from the phenotype table.
"""

from __future__ import annotations

import json
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path
from random import Random

from pandas import DataFrame, concat, read_csv

from flash_fraq.admixture import MAX_DOWNSAMPLE_FRACTION, downsample_fractions
from flash_fraq.feature_matrix import SPLIT_GROUP_ANY
from flash_fraq.input_manifest import InputManifestBuilder, InputManifestConfig, ManifestEntry

MOSDEPTH_TOTAL_ROW: str = "total"
MOSDEPTH_MEAN_COLUMN: int = 3

SPLIT_ASSIGNMENT_COLUMNS: tuple[str, ...] = ("subject_id", "split_group")

# The two-regime design (docs/admixture-design.md, "Reopened decision 2026-09-01"):
# matched-pair mixtures give a clean, unconfounded f ladder against the tumor's own
# matched normal; real-background mixtures keep the multi-subject pool for background
# diversity. `regime` in the plan is a hard partition key -- downstream consumers must
# never pool rows across it.
REGIME_MATCHED_PAIR: str = "matched_pair"
REGIME_REAL_BACKGROUND: str = "real_background"

# matched_pair has exactly one background subject, so a subject-level split is
# impossible by construction and no non-subject-level substitute prevents leakage
# (the same pair appears in every split regardless of how rows are sliced). The
# regime is therefore unsplit -- not even the SPLIT_GROUP_ANY sentinel, which still
# implies participation in the leakage check.
MATCHED_PAIR_SPLIT_GROUP: str = ""

# Every field AdmixtureConfig.from_sources needs, so the eventual MAT-13 Snakemake
# rule is a loop over these columns and nothing has to be re-derived at build time.
PLAN_COLUMNS: tuple[str, ...] = (
    "sample_name",
    "regime",
    "tumor_source_id",
    "healthy_source_id",
    "tumor_subject_id",
    "healthy_subject_id",
    "tumor_fraction",
    "target_coverage",
    "tumor_native_coverage",
    "healthy_native_coverage",
    "seed",
    "platform_tumor",
    "platform_healthy",
    "reference_build",
    "tumor_bam_path",
    "healthy_bam_path",
    "healthy_source_group",
    "split_group",
    "tumor_downsample_fraction",
    "healthy_downsample_fraction",
)


@dataclass(frozen=True)
class CoverageDataset:
    """One derived-tree dataset directory and the source metadata its samples carry.

    Platform, source group and tumor status are properties of the *cohort*, not of
    anything recoverable from a mosdepth summary, so they are declared per dataset --
    same reasoning as `cohort_matrix.CohortDataset`, which this deliberately mirrors.

    `phenotype_path`, when set, replaces the dataset-level `tumor_status` with a
    per-sample one: the DELFI release mixes cancer patients and controls, so one
    value for all 503 Cristiano samples is wrong for about half of them. The table is
    `manifests/cristiano_delfi_phenotype.csv` (scripts/extract_delfi_phenotype.py).
    """

    name: str
    platform: str
    source_group: str
    tumor_status: str
    phenotype_path: Path | None = None


def load_phenotype(phenotype_path: Path) -> dict[str, str]:
    """`sample_id -> tumor_status` from a phenotype CSV."""
    table = read_csv(phenotype_path, dtype=str)
    if not table["sample_id"].is_unique:
        raise ValueError(f"{phenotype_path}: sample_id is not unique")
    return dict(zip(table["sample_id"], table["tumor_status"]))


@dataclass(frozen=True)
class NativeCoverageConfig:
    derived_root: Path
    cram_root: Path
    datasets: tuple[CoverageDataset, ...]
    output_path: Path
    reference_build: str = "hg19"
    require_complete_flag: bool = True

    @classmethod
    def from_sources(
        cls,
        derived_root: Path,
        cram_root: Path,
        datasets: Iterable[dict[str, str]],
        output_path: Path,
        reference_build: str = "hg19",
        require_complete_flag: bool = True,
    ) -> "NativeCoverageConfig":
        return cls(
            derived_root=Path(derived_root),
            cram_root=Path(cram_root),
            datasets=tuple(
                CoverageDataset(
                    name=dataset["name"],
                    platform=dataset["platform"],
                    source_group=dataset["source_group"],
                    tumor_status=dataset["tumor_status"],
                    phenotype_path=Path(dataset["phenotype_path"]) if dataset.get("phenotype_path") else None,
                )
                for dataset in datasets
            ),
            output_path=Path(output_path),
            reference_build=reference_build,
            require_complete_flag=require_complete_flag,
        )


class NativeCoverageHarvester:
    """Genome-wide mean coverage per sample, harvested from the mosdepth summaries."""

    def __init__(self, config: NativeCoverageConfig) -> None:
        self.config = config

    @staticmethod
    def parse_mosdepth_summary(summary_path: Path) -> float:
        """The `mean` column of the `total` row -- genome-wide mean depth.

        `total` and not `total_region`: with no `--by`, mosdepth writes the two
        identically, but the region row is empty when the workflow is ever run with a
        BED, and silently reading 0.0x there would make every sample look infeasible.
        """
        for line in summary_path.read_text().splitlines():
            fields = line.split("\t")
            if fields[0] == MOSDEPTH_TOTAL_ROW:
                return float(fields[MOSDEPTH_MEAN_COLUMN])
        raise ValueError(f"{summary_path}: no '{MOSDEPTH_TOTAL_ROW}' row")

    def cram_path(self, dataset: CoverageDataset, sample_id: str) -> Path:
        return self.config.cram_root / dataset.name / "CRAM" / sample_id / f"{sample_id}.sorted.cram"

    def harvest(self) -> list[ManifestEntry]:
        """One entry per completed sample, gated on `localfeat_sample_complete.flag`.

        Same gate `cohort_matrix.discover_samples` uses and for the same reason: the
        derived tree is written incrementally, so a sample whose run is still in flight
        has a readable but partial `coverage/` directory.
        """
        entries: list[ManifestEntry] = []
        seen: dict[str, str] = {}
        for dataset in self.config.datasets:
            dataset_dir = self.config.derived_root / dataset.name
            if not dataset_dir.is_dir():
                continue
            phenotype = load_phenotype(dataset.phenotype_path) if dataset.phenotype_path is not None else None
            for sample_dir in sorted(path for path in dataset_dir.iterdir() if path.is_dir()):
                sample_id = sample_dir.name
                if self.config.require_complete_flag and not (sample_dir / "localfeat_sample_complete.flag").exists():
                    continue
                summary_path = sample_dir / "coverage" / f"{sample_id}.mosdepth.summary.txt"
                if not summary_path.exists():
                    continue
                if sample_id in seen:
                    raise ValueError(
                        f"sample_id {sample_id!r} appears in both {seen[sample_id]!r} and {dataset.name!r}; "
                        f"sample_id keys the coverage manifest and must be unique"
                    )
                seen[sample_id] = dataset.name
                if phenotype is not None and sample_id not in phenotype:
                    # An unmatched sample would silently keep the dataset-level value,
                    # which is exactly the mislabel the phenotype table exists to fix.
                    raise ValueError(f"sample_id {sample_id!r} of {dataset.name!r} is missing from {dataset.phenotype_path}")
                entries.append(
                    ManifestEntry(
                        sample_id=sample_id,
                        # No sample-to-subject mapping is published for either cohort;
                        # both are one library per donor, so sample is subject here.
                        # Split assignment reads this column, so if a multi-library
                        # subject ever appears it must be fixed *here*, not downstream.
                        subject_id=sample_id,
                        source_group=dataset.source_group,
                        tumor_status=phenotype[sample_id] if phenotype is not None else dataset.tumor_status,
                        platform=dataset.platform,
                        # These are CRAMs; `bam_path` is input_manifest's aligned-reads
                        # column and every samtools call in admixture.py passes
                        # --reference, so CRAM needs no separate file_type.
                        file_type="bam",
                        reference_build=self.config.reference_build,
                        reference_build_verified=True,
                        native_coverage=self.parse_mosdepth_summary(summary_path),
                        bam_path=self.cram_path(dataset, sample_id),
                    )
                )
        return entries

    def write(self, entries: Sequence[ManifestEntry]) -> Path:
        builder = InputManifestBuilder(InputManifestConfig.from_sources(output_path=self.config.output_path))
        for entry in entries:
            builder.add_entry(entry)
        return builder.write()

    @staticmethod
    def summarise(entries: Sequence[ManifestEntry]) -> dict[str, object]:
        """Per-source coverage distribution -- the table the C_target choice rests on."""
        by_source: dict[str, list[float]] = {}
        for entry in entries:
            by_source.setdefault(entry.source_group, []).append(float(entry.native_coverage))
        return {
            "n_samples": len(entries),
            "coverage_by_source_group": {
                source: _describe_coverage(sorted(values)) for source, values in sorted(by_source.items())
            },
        }


def _describe_coverage(sorted_values: list[float]) -> dict[str, float | int]:
    def quantile(q: float) -> float:
        # Nearest-rank, so every reported number is a coverage a real sample actually
        # has -- an interpolated p25 is not a value any library can be downsampled to.
        return sorted_values[min(len(sorted_values) - 1, int(q * len(sorted_values)))]

    return {
        "n": len(sorted_values),
        "min": sorted_values[0],
        "p25": quantile(0.25),
        "median": quantile(0.50),
        "p75": quantile(0.75),
        "max": sorted_values[-1],
        "below_1x": sum(1 for value in sorted_values if value < 1.0),
    }


@dataclass(frozen=True)
class TumorSource:
    """The tumor side of every mixture. One sample, declared rather than discovered.

    HCC1395 is not aligned yet (MAT-23), so `native_coverage` is the phase2 §7 figure
    and `bam_path` is a path that does not exist. Neither is load-bearing for the
    design: the tumor-side downsample fraction is `f*C_target/35x <= 0.03` across the
    whole ladder, so it never binds and a later correction to the real coverage cannot
    change which mixtures are feasible.
    """

    sample_id: str
    subject_id: str
    platform: str
    native_coverage: float
    bam_path: Path


@dataclass(frozen=True)
class AdmixtureDesignConfig:
    coverage_manifest_path: Path
    output_root: Path
    tumor: TumorSource
    matched_pair_background: ManifestEntry
    fraction_ladder: tuple[float, ...]
    candidate_target_coverages: tuple[float, ...]
    split_fractions: tuple[tuple[str, float], ...]
    background_source_groups: tuple[str, ...] = ()
    background_tumor_status: str | None = None
    mixtures_per_fraction: int = 60
    matched_pair_replicates_per_fraction: int = 1
    min_background_samples_per_source: int = 20
    min_background_retention_per_source: float = 0.75
    reference_build: str = "hg19"
    seed: int = 42

    @property
    def feasibility_path(self) -> Path:
        return self.output_root / "admixture_feasibility.csv"

    @property
    def plan_path(self) -> Path:
        return self.output_root / "admixture_plan.csv"

    @property
    def split_assignment_path(self) -> Path:
        return self.output_root / "split_assignment.csv"

    @classmethod
    def from_sources(
        cls,
        coverage_manifest_path: Path,
        output_root: Path,
        tumor: dict[str, object],
        matched_pair_background: dict[str, object],
        fraction_ladder: Iterable[float],
        candidate_target_coverages: Iterable[float],
        split_fractions: dict[str, float],
        background_source_groups: Iterable[str] = (),
        background_tumor_status: str | None = None,
        mixtures_per_fraction: int = 60,
        matched_pair_replicates_per_fraction: int = 1,
        min_background_samples_per_source: int = 20,
        min_background_retention_per_source: float = 0.75,
        reference_build: str = "hg19",
        seed: int = 42,
    ) -> "AdmixtureDesignConfig":
        return cls(
            coverage_manifest_path=Path(coverage_manifest_path),
            output_root=Path(output_root),
            tumor=TumorSource(
                sample_id=str(tumor["sample_id"]),
                subject_id=str(tumor.get("subject_id", tumor["sample_id"])),
                platform=str(tumor["platform"]),
                native_coverage=float(tumor["native_coverage"]),
                bam_path=Path(str(tumor["bam_path"])),
            ),
            matched_pair_background=ManifestEntry(
                sample_id=str(matched_pair_background["sample_id"]),
                subject_id=str(matched_pair_background.get("subject_id", matched_pair_background["sample_id"])),
                source_group=str(matched_pair_background.get("source_group", "hcc1395bl")),
                tumor_status="healthy",
                platform=str(matched_pair_background["platform"]),
                file_type="bam",
                reference_build=reference_build,
                reference_build_verified=True,
                native_coverage=float(matched_pair_background["native_coverage"]),
                bam_path=Path(str(matched_pair_background["bam_path"])),
            ),
            fraction_ladder=tuple(sorted(float(value) for value in fraction_ladder)),
            candidate_target_coverages=tuple(sorted(float(value) for value in candidate_target_coverages)),
            split_fractions=tuple((name, float(weight)) for name, weight in split_fractions.items()),
            background_source_groups=tuple(background_source_groups),
            background_tumor_status=background_tumor_status or None,
            mixtures_per_fraction=int(mixtures_per_fraction),
            matched_pair_replicates_per_fraction=int(matched_pair_replicates_per_fraction),
            min_background_samples_per_source=int(min_background_samples_per_source),
            min_background_retention_per_source=float(min_background_retention_per_source),
            reference_build=reference_build,
            seed=int(seed),
        )


class AdmixtureFeasibilityAnalyzer:
    """The (background sample x `f` x `C_target`) grid and its `<= 1` verdicts."""

    def __init__(self, config: AdmixtureDesignConfig) -> None:
        self.config = config

    def load_background(self) -> list[ManifestEntry]:
        entries = [
            entry
            for entry in InputManifestBuilder.load(self.config.coverage_manifest_path)
            if entry.native_coverage is not None
            and (not self.config.background_source_groups or entry.source_group in self.config.background_source_groups)
            and (self.config.background_tumor_status is None or entry.tumor_status == self.config.background_tumor_status)
        ]
        if not entries:
            raise ValueError(
                f"no background samples with a native coverage in {self.config.coverage_manifest_path} "
                f"(source groups: {self.config.background_source_groups or 'all'}, "
                f"tumor status: {self.config.background_tumor_status or 'any'})"
            )
        return sorted(entries, key=lambda entry: entry.sample_id)

    def build_table(self, background: Sequence[ManifestEntry]) -> DataFrame:
        tumor = self.config.tumor
        rows = []
        for entry in background:
            for tumor_fraction in self.config.fraction_ladder:
                for target_coverage in self.config.candidate_target_coverages:
                    tumor_downsample, healthy_downsample = downsample_fractions(
                        tumor_fraction=tumor_fraction,
                        target_coverage=target_coverage,
                        tumor_native_coverage=tumor.native_coverage,
                        healthy_native_coverage=entry.native_coverage,
                    )
                    rows.append(
                        {
                            "healthy_source_id": entry.sample_id,
                            "healthy_source_group": entry.source_group,
                            "healthy_native_coverage": entry.native_coverage,
                            "tumor_fraction": tumor_fraction,
                            "target_coverage": target_coverage,
                            "tumor_downsample_fraction": tumor_downsample,
                            "healthy_downsample_fraction": healthy_downsample,
                            "feasible": (
                                tumor_downsample <= MAX_DOWNSAMPLE_FRACTION
                                and healthy_downsample <= MAX_DOWNSAMPLE_FRACTION
                            ),
                        }
                    )
        return DataFrame(rows)

    def retention_by_source(self, feasibility: DataFrame, target_coverage: float) -> dict[str, dict[str, int]]:
        """How many samples per source survive at the *binding* fraction.

        The healthy fraction falls monotonically in `f`, so the lowest rung of the
        ladder is the worst case and a sample feasible there is feasible everywhere.
        """
        binding_fraction = min(self.config.fraction_ladder)
        at_binding = feasibility[
            (feasibility["target_coverage"] == target_coverage)
            & (feasibility["tumor_fraction"] == binding_fraction)
        ]
        return {
            str(source): {
                "eligible": int(group["feasible"].sum()),
                "total": int(len(group)),
            }
            for source, group in at_binding.groupby("healthy_source_group")
        }

    def choose_target_coverage(self, feasibility: DataFrame) -> float:
        """The deepest candidate `C_target` that still keeps enough of every source.

        Deeper is better -- more reads per mixture means less counting noise in the
        short:long ratios -- so this walks the candidates from the top and takes the
        first that clears both retention floors for *every* source. Retention is
        measured at the binding fraction, which is what makes the choice honest: a
        `C_target` picked off the median would quietly drop the low-coverage tail at
        exactly the low-`f` end where the regression label matters most.
        """
        for target_coverage in sorted(self.config.candidate_target_coverages, reverse=True):
            retention = self.retention_by_source(feasibility, target_coverage)
            if not retention:
                continue
            if all(
                counts["eligible"] >= self.config.min_background_samples_per_source
                and counts["eligible"] >= self.config.min_background_retention_per_source * counts["total"]
                for counts in retention.values()
            ):
                return target_coverage
        raise ValueError(
            f"no candidate target coverage in {self.config.candidate_target_coverages} keeps "
            f">= {self.config.min_background_samples_per_source} samples and "
            f">= {self.config.min_background_retention_per_source:.0%} of every source at "
            f"f={min(self.config.fraction_ladder)}; lower the floors or add deeper background libraries"
        )


@dataclass
class _Cycle:
    """Deterministic round-robin over one (source, split) cell's background samples.

    Round-robin rather than sampling with replacement so usage is as even as the counts
    allow: a background library that appeared at many `f` values and another that
    appeared at none would put library identity back into correlation with `f`, which
    is the thing this whole design is arranged to avoid.
    """

    members: list[ManifestEntry]
    position: int = 0

    def take(self, count: int) -> list[ManifestEntry]:
        drawn = [self.members[(self.position + offset) % len(self.members)] for offset in range(count)]
        self.position = (self.position + count) % len(self.members)
        return drawn


class MixturePlanner:
    """Turns the feasibility table into a split assignment and a mixture plan."""

    def __init__(self, config: AdmixtureDesignConfig) -> None:
        self.config = config

    def assign_split_groups(self, background: Sequence[ManifestEntry]) -> DataFrame:
        """Source-subject-level train/validation/test, stratified within each source.

        Assigned over *every* background subject in the manifest rather than over the
        eligible pool only, so re-running with a different `C_target` cannot move a
        subject between splits. Stratifying within source keeps both cohorts in all
        three splits -- a split that was purely one source would make the held-out
        score a cross-source generalisation measurement rather than a validation score.

        The tumor subject gets `SPLIT_GROUP_ANY`: v1 has exactly one tumor source and
        it is in every mixture, so binding it to one split would collapse the matrix
        into that split. See the sentinel's note in `feature_matrix`.
        """
        subjects_by_source: dict[str, list[str]] = {}
        for entry in sorted(background, key=lambda entry: entry.sample_id):
            subjects_by_source.setdefault(entry.source_group, []).append(entry.subject_id)

        rows = [{"subject_id": self.config.tumor.subject_id, "split_group": SPLIT_GROUP_ANY}]
        for source, subjects in sorted(subjects_by_source.items()):
            unique_subjects = sorted(set(subjects))
            shuffled = list(unique_subjects)
            Random(f"{self.config.seed}:split:{source}").shuffle(shuffled)
            for subject, split_group in zip(shuffled, self._split_labels(len(shuffled))):
                rows.append({"subject_id": subject, "split_group": split_group})

        return DataFrame(rows, columns=list(SPLIT_ASSIGNMENT_COLUMNS)).sort_values("subject_id", ignore_index=True)

    def _split_labels(self, count: int) -> list[str]:
        names = [name for name, _ in self.config.split_fractions]
        weights = [weight for _, weight in self.config.split_fractions]
        return [label for name, size in zip(names, _largest_remainder(count, weights)) for label in [name] * size]

    def eligible_background(self, background: Sequence[ManifestEntry], target_coverage: float) -> list[ManifestEntry]:
        """The one pool every `f` draws from: samples feasible at *every* rung.

        Restricting to the intersection rather than to each rung's own feasible set is
        what makes source proportions constant across `f`. Per-rung pools would grow as
        `f` rises (the healthy downsample fraction shrinks), and they would grow
        unevenly -- Snyder's 0.75x floor drops out long before Cristiano's tail does --
        so the pool composition would itself encode `f`.
        """
        binding_fraction = min(self.config.fraction_ladder)
        eligible = []
        for entry in background:
            _, healthy_downsample = downsample_fractions(
                tumor_fraction=binding_fraction,
                target_coverage=target_coverage,
                tumor_native_coverage=self.config.tumor.native_coverage,
                healthy_native_coverage=entry.native_coverage,
            )
            if healthy_downsample <= MAX_DOWNSAMPLE_FRACTION:
                eligible.append(entry)
        if not eligible:
            raise ValueError(f"no background sample is feasible at f={binding_fraction}, C_target={target_coverage}")
        return eligible

    def build_plan(
        self, background: Sequence[ManifestEntry], target_coverage: float, splits: DataFrame, row_offset: int = 0
    ) -> DataFrame:
        eligible = self.eligible_background(background, target_coverage)
        subject_to_split = splits.set_index("subject_id")["split_group"]

        cells: dict[tuple[str, str], _Cycle] = {}
        for entry in sorted(eligible, key=lambda entry: entry.sample_id):
            key = (entry.source_group, str(subject_to_split[entry.subject_id]))
            cells.setdefault(key, _Cycle(members=[])).members.append(entry)
        for key, cycle in cells.items():
            Random(f"{self.config.seed}:cell:{key}").shuffle(cycle.members)

        ordered_keys = sorted(cells)
        allocation = _largest_remainder(
            self.config.mixtures_per_fraction, [len(cells[key].members) for key in ordered_keys]
        )

        rows: list[dict[str, object]] = []
        used_names: dict[str, int] = {}
        for tumor_fraction in self.config.fraction_ladder:
            for key, count in zip(ordered_keys, allocation):
                for entry in cells[key].take(count):
                    rows.append(
                        self._plan_row(
                            entry,
                            tumor_fraction,
                            target_coverage,
                            key[1],
                            REGIME_REAL_BACKGROUND,
                            row_offset + len(rows),
                            used_names,
                        )
                    )
        return DataFrame(rows, columns=list(PLAN_COLUMNS))

    def build_matched_pair_plan(self, target_coverage: float, row_offset: int = 0) -> DataFrame:
        """The matched-pair regime: HCC1395 x HCC1395BL across the full ladder.

        One background subject by construction, so none of `build_plan`'s
        orthogonality/split-cell machinery applies -- every row draws the same pair.
        Replicated `matched_pair_replicates_per_fraction` times per rung, under
        distinct seeds, so the regime is more than a single sample per rung.
        Feasibility is checked once at the binding fraction (lowest `f`, same logic
        as `eligible_background`) so an infeasible pair fails loudly here rather
        than being caught later by `validate_plan`.
        """
        entry = self.config.matched_pair_background
        binding_fraction = min(self.config.fraction_ladder)
        _, healthy_downsample = downsample_fractions(
            tumor_fraction=binding_fraction,
            target_coverage=target_coverage,
            tumor_native_coverage=self.config.tumor.native_coverage,
            healthy_native_coverage=entry.native_coverage,
        )
        if healthy_downsample > MAX_DOWNSAMPLE_FRACTION:
            raise ValueError(
                f"matched-pair background {entry.sample_id!r} is infeasible at "
                f"f={binding_fraction}, C_target={target_coverage} "
                f"(healthy downsample fraction {healthy_downsample:.4f} > {MAX_DOWNSAMPLE_FRACTION})"
            )

        rows: list[dict[str, object]] = []
        used_names: dict[str, int] = {}
        for tumor_fraction in self.config.fraction_ladder:
            for _ in range(self.config.matched_pair_replicates_per_fraction):
                rows.append(
                    self._plan_row(
                        entry,
                        tumor_fraction,
                        target_coverage,
                        MATCHED_PAIR_SPLIT_GROUP,
                        REGIME_MATCHED_PAIR,
                        row_offset + len(rows),
                        used_names,
                    )
                )
        return DataFrame(rows, columns=list(PLAN_COLUMNS))

    def _plan_row(
        self,
        entry: ManifestEntry,
        tumor_fraction: float,
        target_coverage: float,
        split_group: str,
        regime: str,
        row_index: int,
        used_names: dict[str, int],
    ) -> dict[str, object]:
        tumor = self.config.tumor
        tumor_downsample, healthy_downsample = downsample_fractions(
            tumor_fraction=tumor_fraction,
            target_coverage=target_coverage,
            tumor_native_coverage=tumor.native_coverage,
            healthy_native_coverage=entry.native_coverage,
        )
        base_name = f"mix_{entry.sample_id}_f{tumor_fraction:.3f}".replace(".", "p")
        # A cell smaller than its allocation makes the round-robin wrap within one `f`,
        # which is a legitimate replicate (different seed, same pair) but needs a
        # distinct sample_name -- sample_name keys the whole §9 matrix.
        replicate = used_names.get(base_name, 0) + 1
        used_names[base_name] = replicate
        return {
            "sample_name": base_name if replicate == 1 else f"{base_name}_rep{replicate}",
            "regime": regime,
            "tumor_source_id": tumor.sample_id,
            "healthy_source_id": entry.sample_id,
            "tumor_subject_id": tumor.subject_id,
            "healthy_subject_id": entry.subject_id,
            "tumor_fraction": tumor_fraction,
            "target_coverage": target_coverage,
            "tumor_native_coverage": tumor.native_coverage,
            "healthy_native_coverage": entry.native_coverage,
            # Per-row rather than one global seed: samtools' -s is `seed.fraction`, so
            # two mixtures sharing a seed and a fraction would retain the very same
            # reads from the same background library.
            "seed": self.config.seed + row_index,
            "platform_tumor": tumor.platform,
            "platform_healthy": entry.platform,
            "reference_build": self.config.reference_build,
            "tumor_bam_path": str(tumor.bam_path),
            "healthy_bam_path": str(entry.bam_path),
            "healthy_source_group": entry.source_group,
            "split_group": split_group,
            "tumor_downsample_fraction": tumor_downsample,
            "healthy_downsample_fraction": healthy_downsample,
        }

    def validate_plan(self, plan: DataFrame) -> None:
        """Assert the plan is feasible *by construction*, through the builder itself.

        `AdmixtureBuilder.compute_downsample_fractions` is the only authority on what
        MAT-13 will accept, so the plan is checked by asking it rather than by
        re-testing the `<= 1` rule against a second copy of the arithmetic.
        """
        from flash_fraq.admixture import AdmixtureBuilder, AdmixtureConfig

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
                    sample_name=row["sample_name"],
                )
            )
            builder.compute_downsample_fractions()

        duplicated = plan["sample_name"][plan["sample_name"].duplicated()]
        if not duplicated.empty:
            raise ValueError(f"duplicate sample_name in plan: {sorted(set(duplicated))}")

    def summarise(self, plan: DataFrame, splits: DataFrame, target_coverage: float) -> dict[str, object]:
        """Everything a reviewer needs to check the two constraints held.

        `max_source_proportion_deviation` is the orthogonality claim as a number: the
        largest gap between any single `f`'s source mix and the pooled source mix. The
        design makes it exactly 0; reporting it means a future change that breaks the
        property shows up in the run log instead of silently in a model score.

        Computed over `real_background` rows only: `matched_pair` has one background
        subject by construction, so folding it in would make the metric meaningless
        rather than merely small, and the two regimes are never pooled anyway.
        """
        real_background = plan[plan["regime"] == REGIME_REAL_BACKGROUND]
        pooled = real_background["healthy_source_group"].value_counts(normalize=True)
        deviation = 0.0
        per_fraction = {}
        for tumor_fraction, group in real_background.groupby("tumor_fraction"):
            proportions = group["healthy_source_group"].value_counts(normalize=True)
            aligned = proportions.reindex(pooled.index).fillna(0.0)
            deviation = max(deviation, float((aligned - pooled).abs().max()))
            per_fraction[f"{tumor_fraction:g}"] = {str(k): round(float(v), 6) for k, v in aligned.items()}
        return {
            "target_coverage": target_coverage,
            "fraction_ladder": list(self.config.fraction_ladder),
            "n_mixtures": int(len(plan)),
            "n_mixtures_by_regime": {str(k): int(v) for k, v in plan["regime"].value_counts().items()},
            "n_background_libraries_used": int(real_background["healthy_source_id"].nunique()),
            "mixtures_by_split_group": {
                str(k): int(v) for k, v in real_background["split_group"].value_counts().items()
            },
            "subjects_by_split_group": {str(k): int(v) for k, v in splits["split_group"].value_counts().items()},
            "source_proportions_by_fraction": per_fraction,
            "max_source_proportion_deviation": round(deviation, 9),
        }


def _largest_remainder(total: int, weights: Sequence[float]) -> list[int]:
    """Apportion `total` over `weights` so the parts sum to exactly `total`.

    Plain rounding leaves the counts short or over by a few, which would make the
    per-`f` allocation drift and break the orthogonality property it exists to hold.
    Ties resolve by descending remainder then ascending index, so it is deterministic.
    """
    weight_total = sum(weights)
    if weight_total <= 0:
        return [0] * len(weights)
    exact = [total * weight / weight_total for weight in weights]
    counts = [int(value) for value in exact]
    remaining = total - sum(counts)
    order = sorted(range(len(weights)), key=lambda index: (-(exact[index] - counts[index]), index))
    for index in order[:remaining]:
        counts[index] += 1
    return counts


class AdmixtureDesigner:
    """Runs the three stages end to end and writes the design's four output files."""

    def __init__(self, config: AdmixtureDesignConfig) -> None:
        self.config = config
        self.analyzer = AdmixtureFeasibilityAnalyzer(config)
        self.planner = MixturePlanner(config)

    def run(self) -> dict[str, object]:
        background = self.analyzer.load_background()
        feasibility = self.analyzer.build_table(background)
        target_coverage = self.analyzer.choose_target_coverage(feasibility)

        splits = self.planner.assign_split_groups(background)
        matched_pair_plan = self.planner.build_matched_pair_plan(target_coverage)
        real_background_plan = self.planner.build_plan(
            background, target_coverage, splits, row_offset=len(matched_pair_plan)
        )
        plan = concat([matched_pair_plan, real_background_plan], ignore_index=True)
        self.planner.validate_plan(plan)

        self.config.output_root.mkdir(parents=True, exist_ok=True)
        feasibility.to_csv(self.config.feasibility_path, index=False)
        plan.to_csv(self.config.plan_path, index=False)
        splits.to_csv(self.config.split_assignment_path, index=False)

        summary = self.planner.summarise(plan, splits, target_coverage)
        summary["retention_by_source_at_binding_fraction"] = self.analyzer.retention_by_source(
            feasibility, target_coverage
        )
        summary["candidate_target_coverages"] = list(self.config.candidate_target_coverages)
        return summary

    def write_summary(self, summary: dict[str, object], summary_path: Path) -> Path:
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        summary_path.write_text(json.dumps(summary, indent=2))
        return summary_path
