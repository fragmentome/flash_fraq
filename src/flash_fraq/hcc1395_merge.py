"""Lane -> library merge for the HCC1395/HCC1395BL (SEQC2) CRAM archive (MAT-23
Block C).

MAT-23 closed with 54+ verified *per-run* CRAMs
(`<cram_root>/HCC1395_SEQC2/CRAM/<RUN>/<RUN>.sorted.cram`), while
`workflow/config/admixture_design.yaml` names *sample-level* paths that nothing
produces. This module is the missing step.

The merge unit is the **library**, not the cell line
=====================================================

A library here is one `(prep, site, arm, variant)` cell of
`manifests/hcc1395/hcc1395_library_pairs.csv` -- e.g. `WGS_IL_T_2`, whose four
banked runs `SRR7890863`/`SRR7910043`/`44`/`45` are top-up lanes of a single
prepped library. Merging those four is loss-free: same molecules, same prep,
same site; splitting them across runs was a sequencing-run artefact, not
biology.

Merging any *wider* than that is not loss-free, and the design does not want it:

* **Across sites.** The 2026-09-04 finding in the parent ticket is that the
  sequencing site lives only in ENA's `library_name` (`center_name` is `FDA` on
  all 118 rows), and that tumor and normal must come from the *same* site or the
  batch confound that demoting Snyder was meant to remove comes back silently.
  A single `HCC1395` CRAM pooling all six sites' tumor lanes destroys exactly
  that property: there is then no site a matched normal could be matched *to*.
* **Across replicates within a site.** `WGS_IL_N_1/2/3` are three separate
  library preps of the same cell line. This project's entire feature set is
  fragment-length structure (DELFI short:long ratios), which is a property of
  the library prep -- pooling three preps yields a blended fragment
  distribution that no real library has, and it breaks the replicate pairing
  (IL has normal replicates 1/2/3 banked but tumor 2/3, so a "per site/arm"
  merge would pair 3 normal preps against 2 tumor preps).
* **Nothing consumes the extra depth.** `admixture_design.yaml` asks for 35x
  tumor / 34x normal, and MAT-13 mixes down to `C_target <= 1x`. Every banked
  library already clears 35x on its own lanes (the weakest, `WGS_NS_N_8`, is
  ~36x). Pooling to 250x buys depth the ladder immediately downsamples away,
  at the cost of the two properties above.

So: one merged CRAM per library, named `<subject>_<prep>_<site>_<variant>` --
`HCC1395_WGS_IL_2` (tumor) and `HCC1395BL_WGS_IL_2` (its matched normal). The
subject id stays `HCC1395`/`HCC1395BL` so subject-level splits and the
matched-pair regime still see one tumor subject and one background subject.

`SPP_*` runs are excluded: they are pre-pooled genotype mixtures with no T/N
arm at all, so they are not a tumor library or a background library and must
not be merged into either.
"""

from __future__ import annotations

import csv
from collections.abc import Iterable, Sequence
from dataclasses import dataclass, replace
from pathlib import Path
from subprocess import run
from tempfile import TemporaryDirectory

# manifests/hcc1395/hcc1395_library_pairs.csv columns.
PAIRS_MANIFEST_COLUMNS: tuple[str, ...] = (
    "run_accession",
    "library_name",
    "prep",
    "site",
    "arm",
    "variant",
    "read_count",
    "base_count",
)

ARM_TUMOR: str = "T"
ARM_NORMAL: str = "N"

# One entry per T/N arm: (subject_id, source_group, tumor_status). These are
# input_manifest's vocabulary and match hcc1395_manifest.SraRun's derivation, so
# a merged library lands in the cohort tables tagged exactly like its lanes.
ARM_METADATA: dict[str, tuple[str, str, str]] = {
    ARM_TUMOR: ("HCC1395", "hcc1395", "tumor"),
    ARM_NORMAL: ("HCC1395BL", "hcc1395bl", "healthy"),
}

# Preps that carry a real tumor/normal arm. CHR/LBP are other prep families in
# SRP162370; SPP is the pre-pooled genotype titration and has no arm at all.
DEFAULT_MERGE_PREPS: tuple[str, ...] = ("WGS",)

# The dataset directory the per-run CRAMs were mirrored into, and that the
# merged libraries are written back into, so that
# NativeCoverageHarvester.cram_path and local_features' discover_local_crams
# both resolve them without a special case.
HCC1395_DATASET: str = "HCC1395_SEQC2"

# Matches a merged library id and never a bare run accession (SRR#######), so
# local_features.yaml can register this dataset without also enqueueing the
# per-run lane archive that lives beside it.
MERGED_SAMPLE_PATTERN: str = r"HCC1395(BL)?_[A-Z]+_[A-Z]+_[0-9]+"


@dataclass(frozen=True)
class LibraryLane:
    """One sequencing run: a row of hcc1395_library_pairs.csv."""

    run_accession: str
    library_name: str
    prep: str
    site: str
    arm: str
    variant: str
    read_count: int
    base_count: int

    @classmethod
    def from_row(cls, row: dict[str, str]) -> "LibraryLane":
        return cls(
            run_accession=row["run_accession"],
            library_name=row["library_name"],
            prep=row["prep"],
            site=row["site"],
            arm=row["arm"],
            variant=row["variant"],
            read_count=int(row["read_count"] or 0),
            base_count=int(row["base_count"] or 0),
        )

    def cram_path(self, cram_root: Path, dataset: str = HCC1395_DATASET) -> Path:
        return cram_root / dataset / "CRAM" / self.run_accession / f"{self.run_accession}.sorted.cram"

    def crai_path(self, cram_root: Path, dataset: str = HCC1395_DATASET) -> Path:
        return Path(f"{self.cram_path(cram_root, dataset)}.crai")


@dataclass(frozen=True)
class MergeUnit:
    """One library's worth of lanes, and the sample-level CRAM they merge into."""

    sample_id: str
    subject_id: str
    source_group: str
    tumor_status: str
    prep: str
    site: str
    arm: str
    variant: str
    lanes: tuple[LibraryLane, ...]

    @property
    def library_name(self) -> str:
        return f"{self.prep}_{self.site}_{self.arm}_{self.variant}"

    @property
    def run_accessions(self) -> tuple[str, ...]:
        return tuple(lane.run_accession for lane in self.lanes)

    @property
    def base_count(self) -> int:
        return sum(lane.base_count for lane in self.lanes)

    def approximate_coverage(self, genome_size_bp: int = 3_100_000_000) -> float:
        """Bases/genome -- the arithmetic the close-out table used. A planning
        number only: the real one comes from mosdepth on the merged CRAM."""
        return self.base_count / genome_size_bp

    def lane_cram_paths(self, cram_root: Path, dataset: str = HCC1395_DATASET) -> list[Path]:
        return [lane.cram_path(cram_root, dataset) for lane in self.lanes]

    def lane_crai_paths(self, cram_root: Path, dataset: str = HCC1395_DATASET) -> list[Path]:
        return [lane.crai_path(cram_root, dataset) for lane in self.lanes]

    def merged_cram_path(self, cram_root: Path, dataset: str = HCC1395_DATASET) -> Path:
        return cram_root / dataset / "CRAM" / self.sample_id / f"{self.sample_id}.sorted.cram"


def merged_sample_id(subject_id: str, prep: str, site: str, variant: str) -> str:
    return f"{subject_id}_{prep}_{site}_{variant}"


@dataclass(frozen=True)
class HCC1395MergeConfig:
    pairs_manifest: Path
    cram_root: Path
    reference_fasta: Path
    dataset: str = HCC1395_DATASET
    preps: tuple[str, ...] = DEFAULT_MERGE_PREPS
    threads: int = 8
    # Off by default: it is a second full pass over a ~200 GB CRAM plus one
    # over every lane, i.e. hours of IO to re-derive a number `samtools merge`
    # does not lose. quickcheck (structural integrity + EOF container) is the
    # gate that catches the failure that actually happens -- a truncated
    # output. Turn it on for a one-off audit, or in tests, where it is free.
    verify_read_counts: bool = False

    @classmethod
    def from_sources(
        cls,
        pairs_manifest: Path,
        cram_root: Path,
        reference_fasta: Path,
        dataset: str = HCC1395_DATASET,
        preps: Iterable[str] = DEFAULT_MERGE_PREPS,
        threads: int = 8,
        verify_read_counts: bool = False,
    ) -> "HCC1395MergeConfig":
        return cls(
            pairs_manifest=Path(pairs_manifest),
            cram_root=Path(cram_root),
            reference_fasta=Path(reference_fasta),
            dataset=dataset,
            preps=tuple(preps),
            threads=threads,
            verify_read_counts=verify_read_counts,
        )


class LibraryMergePlanner:
    """Groups hcc1395_library_pairs.csv rows into the merge units above."""

    def __init__(self, config: HCC1395MergeConfig) -> None:
        self.config = config

    def load_lanes(self) -> list[LibraryLane]:
        with open(self.config.pairs_manifest, newline="") as handle:
            reader = csv.DictReader(handle)
            missing = set(PAIRS_MANIFEST_COLUMNS) - set(reader.fieldnames or ())
            if missing:
                raise ValueError(
                    f"{self.config.pairs_manifest}: missing column(s) {sorted(missing)}; "
                    f"expected {list(PAIRS_MANIFEST_COLUMNS)}"
                )
            return [LibraryLane.from_row(row) for row in reader]

    def plan(self) -> list[MergeUnit]:
        """Every library the manifest declares, whether or not its lanes are on
        disk. `plan_available` is the on-disk-filtered view."""
        grouped: dict[tuple[str, str, str, str], list[LibraryLane]] = {}
        for lane in self.load_lanes():
            if lane.prep not in self.config.preps:
                continue
            # SPP rows carry an empty arm: a pre-pooled genotype mixture is
            # neither a tumor library nor a background one.
            if lane.arm not in ARM_METADATA:
                continue
            grouped.setdefault((lane.prep, lane.site, lane.arm, lane.variant), []).append(lane)

        units: list[MergeUnit] = []
        for (prep, site, arm, variant), lanes in sorted(grouped.items()):
            subject_id, source_group, tumor_status = ARM_METADATA[arm]
            units.append(
                MergeUnit(
                    sample_id=merged_sample_id(subject_id, prep, site, variant),
                    subject_id=subject_id,
                    source_group=source_group,
                    tumor_status=tumor_status,
                    prep=prep,
                    site=site,
                    arm=arm,
                    variant=variant,
                    lanes=tuple(sorted(lanes, key=lambda lane: lane.run_accession)),
                )
            )
        return units

    def plan_available(self) -> list[MergeUnit]:
        """Libraries restricted to the lanes whose CRAM *and* index are banked.

        Requiring the `.crai` is the same gate `discover_local_crams` uses: the
        mirror is written by rsync --partial, so a CRAM at its final filename
        may still be in flight, while the index only appears once its own
        transfer finishes. A library with no banked lane is dropped entirely --
        11 of MAT-23's 65 runs were never retrieved, and a merge unit that
        silently produced a partial library would be worse than an absent one.
        """
        available: list[MergeUnit] = []
        for unit in self.plan():
            lanes = tuple(
                lane
                for lane in unit.lanes
                if lane.cram_path(self.config.cram_root, self.config.dataset).exists()
                and lane.crai_path(self.config.cram_root, self.config.dataset).exists()
            )
            if lanes:
                available.append(replace(unit, lanes=lanes))
        return available

    def matched_pairs(self, units: Sequence[MergeUnit]) -> list[tuple[MergeUnit, MergeUnit]]:
        """(tumor, normal) for every (prep, site, variant) with both arms.

        The replicate is part of the key on purpose: `WGS_NS_N_8` x
        `WGS_NS_T_7` is the one banked pair that is replicate-mismatched, and
        MAT-23's close-out asks for that to stay visible rather than be papered
        over, so this returns nothing for NS and the caller must choose.
        """
        by_key = {(unit.prep, unit.site, unit.variant, unit.arm): unit for unit in units}
        pairs = []
        for (prep, site, variant, arm), unit in sorted(by_key.items()):
            if arm != ARM_TUMOR:
                continue
            normal = by_key.get((prep, site, variant, ARM_NORMAL))
            if normal is not None:
                pairs.append((unit, normal))
        return pairs

    def plan_paired_available(self) -> list[MergeUnit]:
        """`plan_available` restricted to libraries in a complete T/N pair.

        The one filter both `hcc1395_merge_units` (workflow/rules/common.smk)
        and `run_hcc1395_merged_coverage.py` need: a library whose partner arm
        was never sequenced is never merged, so it must not show up as an
        expected coverage row either -- see common.smk's `hcc1395_merge_units`
        for the full "why pairs_only" rationale.
        """
        units = self.plan_available()
        paired = {unit.sample_id for pair in self.matched_pairs(units) for unit in pair}
        return [unit for unit in units if unit.sample_id in paired]


class CramLibraryMerger:
    """`samtools merge` over one library's lanes, indexed and quickchecked."""

    def __init__(self, config: HCC1395MergeConfig) -> None:
        self.config = config

    def header_lines(self, cram_path: Path) -> list[str]:
        result = run(
            ["samtools", "view", "-H", "--reference", str(self.config.reference_fasta), str(cram_path)],
            check=True,
            capture_output=True,
            text=True,
        )
        return result.stdout.splitlines()

    @staticmethod
    def sequence_lines(header: Sequence[str]) -> list[str]:
        return [line for line in header if line.startswith("@SQ")]

    def build_merged_header(self, unit: MergeUnit, lane_crams: Sequence[Path]) -> str:
        """@HD + the lanes' shared @SQ block + one @RG per lane, SM rewritten.

        `samtools merge` carries every lane's @RG through untouched, and each
        lane was aligned with `ID:<run> SM:<run>` (alignment.FastqAligner), so
        an unrewritten merge advertises four different samples in one file.
        Rewriting SM -- and only SM -- keeps each read's RG:Z tag resolvable
        (the IDs are unchanged, so lanes stay separable) while making the file
        answer "which sample is this?" with the merged library.

        Also the place the @SQ equality check happens: `samtools merge` will
        combine mismatched references by concatenating their @SQ blocks, which
        would produce a file that is structurally valid and biologically
        meaningless.
        """
        headers = [self.header_lines(path) for path in lane_crams]
        reference_sq = self.sequence_lines(headers[0])
        if not reference_sq:
            raise ValueError(f"{lane_crams[0]}: no @SQ lines in header")
        for path, header in zip(lane_crams[1:], headers[1:]):
            if self.sequence_lines(header) != reference_sq:
                raise ValueError(
                    f"{path}: @SQ header differs from {lane_crams[0]}; "
                    f"lanes of {unit.sample_id} were not aligned to the same reference"
                )

        lines = ["@HD\tVN:1.6\tSO:coordinate", *reference_sq]
        for path, header in zip(lane_crams, headers):
            for line in header:
                if line.startswith("@RG"):
                    lines.append(self._rewrite_sample(line, unit.sample_id))
        lines.append(
            f"@CO\tflash_fraq hcc1395_merge: library {unit.library_name} "
            f"from lanes {','.join(unit.run_accessions)}"
        )
        return "\n".join(lines) + "\n"

    @staticmethod
    def _rewrite_sample(read_group_line: str, sample_id: str) -> str:
        fields = read_group_line.split("\t")
        rewritten = [field if not field.startswith("SM:") else f"SM:{sample_id}" for field in fields]
        if not any(field.startswith("SM:") for field in fields):
            rewritten.append(f"SM:{sample_id}")
        return "\t".join(rewritten)

    def flagstat_total_reads(self, alignment_path: Path) -> int:
        result = run(
            [
                "samtools", "flagstat", "-@", str(self.config.threads),
                "--input-fmt-option", f"reference={self.config.reference_fasta}",
                str(alignment_path),
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        return int(result.stdout.splitlines()[0].split()[0])

    def quickcheck(self, cram_path: Path) -> bool:
        return run(["samtools", "quickcheck", str(cram_path)]).returncode == 0

    def merge(self, unit: MergeUnit) -> Path:
        """Merge `unit`'s lanes into one indexed, verified sample-level CRAM.

        Raises (leaving the partial output in place for Snakemake to delete)
        rather than returning a file that failed quickcheck.
        """
        lane_crams = unit.lane_cram_paths(self.config.cram_root, self.config.dataset)
        missing = [path for path in lane_crams if not path.exists()]
        if missing:
            raise FileNotFoundError(f"{unit.sample_id}: lane CRAM(s) absent: {[str(path) for path in missing]}")

        merged_cram = unit.merged_cram_path(self.config.cram_root, self.config.dataset)
        merged_cram.parent.mkdir(parents=True, exist_ok=True)

        with TemporaryDirectory() as temporary_dir:
            header_path = Path(temporary_dir) / f"{unit.sample_id}.header.sam"
            header_path.write_text(self.build_merged_header(unit, lane_crams))
            run(
                [
                    "samtools", "merge",
                    "-f",
                    "-@", str(self.config.threads),
                    "-O", "CRAM",
                    "--reference", str(self.config.reference_fasta),
                    # -c/-p: the lanes' @RG and @PG IDs are already unique (the
                    # run accession), so combining rather than uniquifying them
                    # keeps every read's RG:Z tag pointing at a declared group.
                    "-c", "-p",
                    "-h", str(header_path),
                    "-o", str(merged_cram),
                    *[str(path) for path in lane_crams],
                ],
                check=True,
            )

        if not self.quickcheck(merged_cram):
            raise RuntimeError(f"samtools quickcheck failed for {merged_cram}")

        run(["samtools", "index", "-@", str(self.config.threads), str(merged_cram)], check=True)

        if self.config.verify_read_counts:
            expected = sum(self.flagstat_total_reads(path) for path in lane_crams)
            observed = self.flagstat_total_reads(merged_cram)
            if expected != observed:
                raise RuntimeError(
                    f"{merged_cram}: merged read count {observed} != sum of lanes {expected}"
                )
        return merged_cram


class MergedCoverageReporter:
    """mosdepth summaries for the merged libraries -> the two numbers
    `workflow/config/admixture_design.yaml` still declares as placeholders.

    Reads the `mean` column of the `total` row via
    `NativeCoverageHarvester.parse_mosdepth_summary` -- the same parser, and so
    the same field, as `rule localfeat_coverage` feeds for the 547 cohort rows
    already in `native_coverage.csv`. `total`, not `total_region`: with
    `--by 1000000` the two differ, and the region row is the wrong one.
    """

    COLUMNS: tuple[str, ...] = (
        "sample_id",
        "subject_id",
        "source_group",
        "tumor_status",
        "library_name",
        "site",
        "arm",
        "variant",
        "n_lanes",
        "run_accessions",
        "native_coverage",
        "cram_path",
    )

    def __init__(self, config: HCC1395MergeConfig) -> None:
        self.config = config

    @staticmethod
    def summary_path(derived_root: Path, dataset: str, sample_id: str) -> Path:
        """Exactly where `rule localfeat_coverage` would have put it, so the
        merged libraries are harvestable by NativeCoverageHarvester unchanged if
        HCC1395_SEQC2 is ever added to admixture_design's `datasets:`."""
        return derived_root / dataset / sample_id / "coverage" / f"{sample_id}.mosdepth.summary.txt"

    def rows(self, units: Sequence[MergeUnit], derived_root: Path) -> list[dict[str, object]]:
        # Imported here rather than at module scope: admixture_design pulls in
        # pandas and this module is also imported by workflow/rules/common.smk
        # at DAG-build time, where the only thing needed is the planner.
        from flash_fraq.admixture_design import NativeCoverageHarvester

        rows: list[dict[str, object]] = []
        for unit in units:
            summary = self.summary_path(derived_root, self.config.dataset, unit.sample_id)
            rows.append(
                {
                    "sample_id": unit.sample_id,
                    "subject_id": unit.subject_id,
                    "source_group": unit.source_group,
                    "tumor_status": unit.tumor_status,
                    "library_name": unit.library_name,
                    "site": unit.site,
                    "arm": unit.arm,
                    "variant": unit.variant,
                    "n_lanes": len(unit.lanes),
                    "run_accessions": ";".join(unit.run_accessions),
                    "native_coverage": NativeCoverageHarvester.parse_mosdepth_summary(summary),
                    "cram_path": str(unit.merged_cram_path(self.config.cram_root, self.config.dataset)),
                }
            )
        return rows

    def write(self, rows: Sequence[dict[str, object]], output_path: Path) -> Path:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(self.COLUMNS))
            writer.writeheader()
            for row in sorted(rows, key=lambda row: str(row["sample_id"])):
                writer.writerow(row)
        return output_path

    @staticmethod
    def config_yaml_block(rows: Sequence[dict[str, object]]) -> str:
        """The `native_coverage:`/`bam_path:` lines to paste into
        admixture_design.yaml, so replacing the 35.0/34.0 placeholders is a
        copy rather than a re-derivation."""
        lines = ["# measured by rule hcc1395_merged_coverage (mosdepth total-row mean)"]
        for row in sorted(rows, key=lambda row: str(row["sample_id"])):
            lines.append(
                f"# {row['sample_id']}: native_coverage: {float(row['native_coverage']):.2f}  "
                f"bam_path: {row['cram_path']}"
            )
        return "\n".join(lines)
