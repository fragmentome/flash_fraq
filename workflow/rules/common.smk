"""Helper functions shared by the rule files included from workflow/Snakefile.

Kept out of the individual rule files so those stay rules-only (Snakemake's
--lint flags mixing rules and functions in one snakefile). Globals these
close over (e.g. HCC1395_OUTPUT_ROOT) are defined in the rule file that owns
them; include: shares one Python namespace, and these bodies resolve at call
time, so the owning file need not be included before this one.
"""

import csv
import re
from pathlib import Path


def _read_run_list(run_list_path: str) -> list[dict[str, str]]:
    with open(run_list_path, newline="") as handle:
        return list(csv.DictReader(handle))


def _run_metadata(run_accession: str) -> dict[str, str]:
    run_list_path = checkpoints.resolve_hcc1395_runs.get().output.run_list
    return next(row for row in _read_run_list(run_list_path) if row["run_accession"] == run_accession)


def all_run_accessions(wildcards) -> list[str]:
    run_list_path = checkpoints.resolve_hcc1395_runs.get(**wildcards).output.run_list
    return [row["run_accession"] for row in _read_run_list(run_list_path)]


def all_aligned_bams(wildcards) -> list[str]:
    return expand(str(HCC1395_OUTPUT_ROOT / "bam" / "{run}.sorted.bam"), run=all_run_accessions(wildcards))


def dataset_name(dataset) -> str:
    """The directory name of a local_features `datasets:` entry.

    Entries are either a bare string (`"Snyder_Shendure_WPS"`) or a mapping
    carrying a `sample_pattern` alongside the name. The mapping form exists for
    HCC1395_SEQC2, whose CRAM/ directory holds both the 56 per-run lane
    archives and the merged libraries built from them -- without a filter,
    registering that dataset would enqueue the full feature chain on every lane
    as well as on the library it was merged into.
    """
    return dataset if isinstance(dataset, str) else dataset["name"]


def dataset_sample_pattern(dataset) -> str | None:
    return None if isinstance(dataset, str) else dataset.get("sample_pattern")


def dataset_names(datasets: list) -> list[str]:
    return [dataset_name(dataset) for dataset in datasets]


def discover_local_crams(cram_root: Path, datasets: list) -> list[tuple[str, str]]:
    """(dataset, sample_id) for every mirrored CRAM that has its index alongside.

    Requiring the .crai is what keeps a still-transferring sample out of the
    DAG. sync_crams_to_hdd.sh runs rsync with --partial and no --partial-dir,
    so an in-flight CRAM sits at its *final* filename holding partial bytes --
    globbing '*.sorted.cram' alone would happily enqueue a truncated file. The
    index only appears once its own transfer finishes, and rule verify_cram
    runs `samtools quickcheck` on top of that to catch the remaining case
    where the .crai landed first (rsync orders by size, not by pair).

    A dataset may additionally declare a `sample_pattern`, which the sample id
    must match in full; see dataset_name() for why HCC1395_SEQC2 needs one.
    """
    discovered: list[tuple[str, str]] = []
    for dataset in datasets:
        name = dataset_name(dataset)
        pattern = dataset_sample_pattern(dataset)
        compiled = re.compile(pattern) if pattern else None
        for cram in sorted((cram_root / name / "CRAM").glob("*/*.sorted.cram")):
            sample_id = cram.name.removesuffix(".sorted.cram")
            if compiled is not None and compiled.fullmatch(sample_id) is None:
                continue
            if cram.with_suffix(".cram.crai").exists():
                discovered.append((name, sample_id))
    return discovered


def local_cram_size_mb(cram_root: Path, dataset: str, sample_id: str, multiplier: float = 1.0) -> int:
    """Size of the *HDD-mirrored* CRAM, in MB, scaled by multiplier.

    Deliberately stats the mirror rather than the staged copy. Rules that
    consume a staged file would otherwise size their ssd_mb resource against a
    path that does not exist when the resource callable is evaluated, silently
    yielding 0 and letting unlimited jobs stage at once -- which is exactly the
    case the budget exists to prevent.
    """
    cram = cram_root / dataset / "CRAM" / sample_id / f"{sample_id}.sorted.cram"
    return max(1, int(cram.stat().st_size / 1e6 * multiplier))


def all_aligned_crams(wildcards) -> list[str]:
    return expand(str(HCC1395_OUTPUT_ROOT / "bam" / "{run}.sorted.cram"), run=all_run_accessions(wildcards))


def _admixture_plan_rows() -> list[dict[str, str]]:
    plan_path = checkpoints.admixture_design.get().output.plan
    return _read_run_list(str(plan_path))


def _admixture_plan_row(sample_name: str) -> dict[str, str]:
    return next(row for row in _admixture_plan_rows() if row["sample_name"] == sample_name)


def all_admixture_sample_names(wildcards) -> list[str]:
    checkpoints.admixture_design.get(**wildcards)
    return [row["sample_name"] for row in _admixture_plan_rows()]


def all_admixture_sample_metadata(wildcards) -> list[str]:
    return expand(
        str(ADMIX_OUTPUT_ROOT / "{sample_name}" / "sample_metadata.csv"),
        sample_name=all_admixture_sample_names(wildcards),
    )


def _tumor_sub_batch_dir(wildcards) -> Path:
    """The one-pass tumor-downsample batch's output directory (admixture.smk's
    checkpoint admixture_tumor_downsample_batch), re-evaluated per call like
    the other checkpoint-backed helpers above -- see tumor_batch_sampling.py
    for why every admixture_plan.csv row's tumor-side downsample is built
    there in a single shared-CRAM pass rather than per row.

    A completed batch (tumor_total_reads.txt is written last) is used as-is:
    on a resumed run Snakemake neither resolves checkpoints.get() for the
    existing temp() directory nor reschedules the checkpoint, so every
    admixture_mixture job would otherwise get a literal `<TBD>` input."""
    batch_dir = ADMIX_OUTPUT_ROOT / "tumor_sub_batch"
    if (batch_dir / "tumor_total_reads.txt").is_file():
        return batch_dir
    return Path(checkpoints.admixture_tumor_downsample_batch.get(**wildcards).output.tumor_sub_dir)


def admixture_tumor_sub_bam(wildcards) -> str:
    return str(_tumor_sub_batch_dir(wildcards) / f"{wildcards.sample_name}.tumor_sub.bam")


def admixture_tumor_total_reads(wildcards) -> str:
    total_reads_path = _tumor_sub_batch_dir(wildcards) / "tumor_total_reads.txt"
    return total_reads_path.read_text().strip()


def _row_is_pure_background(wildcards) -> bool:
    """True for a plan row with tumor_downsample_fraction == 0 -- a pure-background
    negative control (nominal admixture fraction f=0). tumor_batch_sampling.py's
    _load_targets skips building a tumor_sub.bam for these rows entirely (there is
    no tumor read to keep), so rule admixture_mixture must not request one for them
    -- see admixture_tumor_sub_input/admixture_tumor_sub_cli_args below."""
    return float(_admixture_plan_row(wildcards.sample_name)["tumor_downsample_fraction"]) == 0


def admixture_tumor_sub_input(wildcards) -> list[str]:
    """rule admixture_mixture's `tumor_sub` input. Empty for a pure-background row,
    so that row does not need to wait on checkpoint admixture_tumor_downsample_batch
    at all -- there is nothing for it to consume from that batch."""
    if _row_is_pure_background(wildcards):
        return []
    return [admixture_tumor_sub_bam(wildcards)]


def admixture_tumor_sub_cli_args(wildcards) -> str:
    """The `--tumor-sub-bam`/`--tumor-total-reads` flags for run_admixture.py's
    command line, or an empty string for a pure-background row --
    AdmixtureBuilder.run skips tumor-side downsampling entirely when the row's
    tumor fraction is 0, so passing neither flag is exactly what it expects
    (both default to None)."""
    if _row_is_pure_background(wildcards):
        return ""
    return (
        f"--tumor-sub-bam {admixture_tumor_sub_bam(wildcards)} "
        f"--tumor-total-reads {admixture_tumor_total_reads(wildcards)}"
    )


def hcc1395_merge_units():
    """The HCC1395 libraries that can be merged from what is banked on disk.

    Resolved once at DAG construction, like discover_local_crams and for the
    same reason: the CRAM mirror is append-only while a sync runs, so a
    per-job re-plan would grow the DAG underneath a running workflow.

    Imported here rather than at module scope so that a checkout without the
    HCC1395 mirror mounted still parses the workflow -- the planner only
    touches the filesystem when this is called.
    """
    from flash_fraq.hcc1395_merge import HCC1395MergeConfig, LibraryMergePlanner

    config = HCC1395MergeConfig.from_sources(
        pairs_manifest=HCC1395MERGE_PAIRS_MANIFEST,
        cram_root=HCC1395MERGE_CRAM_ROOT,
        reference_fasta=HCC1395MERGE_REFERENCE_FASTA,
        dataset=HCC1395MERGE_DATASET,
        preps=HCC1395MERGE_PREPS,
    )
    planner = LibraryMergePlanner(config)

    # Merge only libraries that form a complete tumor/normal pair, unless the
    # config says otherwise. Decided 2026-09-16: the merged libraries exist to
    # feed MAT-13's matched_pair regime, which consumes a T and an N from the
    # same site and replicate. A library whose partner arm was never sequenced
    # -- there are six, e.g. HCC1395_WGS_NS_4 and HCC1395BL_WGS_IL_1 -- cannot
    # enter that regime, is not a background for any other tumor (mixing across
    # sites rebuilds the batch confound the 2026-08-31 demotion removed), and
    # so has no consumer. Skipping them saves ~250 GB of the ~1.07 TB the full
    # 24 would cost, and costs nothing: set pairs_only: false and they build.
    units = planner.plan_paired_available() if HCC1395MERGE_PAIRS_ONLY else planner.plan_available()

    return {unit.sample_id: unit for unit in units}
