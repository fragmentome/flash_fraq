"""In silico tumor/healthy admixture mixture build (phase2 §7 step 2, MAT-13's
wiring half -- see docs/todo/completed/admixture-mixture-snakemake-wiring.md).

One rule per admixture_plan.csv row: downsample the row's tumor and healthy
sources, merge, dedup once (never per-source), and compute both SNV-VAF
empirical tumor-fraction estimators against the merged dedup BAM. Consumes
flash_fraq.admixture.AdmixtureBuilder via scripts/run_admixture.py; nothing
here re-derives the mixing arithmetic.

admixture_plan.csv already carries `tumor_bam_path`/`healthy_bam_path`
(AdmixtureDesigner._plan_row) -- unlike the HCC1395/cohort CRAM trees this
workflow otherwise globs (discover_local_crams), the plan resolves its own
source paths at design time, so this file only reads the plan, it does not
re-derive a path from `regime` + a source id.

admixture_plan.csv is itself checkpoint admixture_design's output, so its rows
-- and their count -- are unknown at DAG-parse time. common.smk's
_admixture_plan_row/all_admixture_sample_names re-read it through
checkpoints.admixture_design.get(...) on every call, the same pattern
all_run_accessions/_run_metadata use for hcc1395's run-time-resolved SRA
accessions -- kept there rather than here per --lint's "mixed rules and
functions" rule.

Globals are prefixed ADMIX_ (see admixture_design.smk's ADMIXDES_ for why:
include: shares one Python namespace, so a generic name would collide).
"""

from pathlib import Path

ADMIX_CONFIG = config["admixture"]
ADMIX_OUTPUT_ROOT = Path(ADMIX_CONFIG["output_root"])
ADMIX_REFERENCE_FASTA = ADMIX_CONFIG["reference_fasta"]
ADMIX_SNV_VAF_PATHS = ADMIX_CONFIG["snv_vaf_paths"]
ADMIX_THREADS = ADMIX_CONFIG.get("threads", 8)
ADMIX_MEM_MB = ADMIX_CONFIG.get("mem_mb", 8192)
ADMIX_RUNTIME = ADMIX_CONFIG.get("runtime", 480)

ADMIX_SNV_VAF_ARGUMENTS = " ".join(f"--snv-vaf-path {label}={path}" for label, path in ADMIX_SNV_VAF_PATHS.items())
ADMIX_ESTIMATOR_LABELS = sorted(ADMIX_SNV_VAF_PATHS)


checkpoint admixture_tumor_downsample_batch:
    """Downsample the shared tumor CRAM for every admixture_plan.csv row, in
    one pass.

    admixture_design.py's `matched_pair` and `real_background` regimes both
    downsample from the *same* single global `TumorSource` CRAM (its
    `tumor_bam_path` column is identical across every row of either regime).
    `samtools view -s` -- what rule admixture_mixture used to call directly,
    once per row -- has no index-assisted subsampling path: it always decodes
    the full source linearly. Doing that once per plan row meant decoding the
    ~88GB tumor CRAM once per row (~490 times for the current plan) to
    produce outputs that are each a small fraction of it.
    flash_fraq.tumor_batch_sampling.TumorBatchSampler instead decodes it
    exactly once here and writes every row's tumor_sub.bam from that one
    pass -- see its docstring for the sampling scheme and why it does not
    reproduce samtools -s's exact read selection (reproducibility and
    per-row distinctness are preserved; bit-for-bit agreement with an
    already-completed row built the old way is not, and is not needed).

    A checkpoint, like admixture_design: the row count -- and so the count of
    tumor_sub.bam files this produces -- is only known once admixture_design
    has run, and the batch script itself decides internally how many writers
    to open; Snakemake only needs to track the directory as a unit. temp()
    on the directory: once every admixture_mixture job that reads from it has
    finished, the whole directory (all per-row tumor_sub.bam plus the shared
    tumor_total_reads.txt) is removed together, same "don't multiply
    peak-scratch across the run" intent as admixture_mixture's own temp()
    outputs -- just coarser-grained, because a one-pass producer can't know
    when an individual row's consumer is done before every row's is.
    """
    input:
        plan=lambda wc: checkpoints.admixture_design.get(**wc).output.plan,
    output:
        tumor_sub_dir=temp(directory(str(ADMIX_OUTPUT_ROOT / "tumor_sub_batch"))),
    log:
        str(ADMIX_OUTPUT_ROOT / "logs" / "admixture_tumor_downsample_batch.log"),
    conda:
        "flash_fraq"
    threads: ADMIX_THREADS
    resources:
        mem_mb=ADMIX_MEM_MB,
        # One job now does the tumor-side decode work that used to be spread
        # (redundantly) across every admixture_mixture job -- budget it like
        # several of those combined rather than like a single row.
        runtime=ADMIX_RUNTIME * 6,
    shell:
        f"python {SCRIPTS_DIR / 'run_tumor_batch_sampling.py'} "
        "--plan {input.plan} "
        f"--reference-fasta {ADMIX_REFERENCE_FASTA} "
        "--output-dir {output.tumor_sub_dir} "
        "--threads {threads} "
        "> {log} 2>&1"


rule admixture_mixture:
    """Merge + dedup + both empirical tumor-fraction estimators for one
    admixture_plan.csv row.

    The tumor-side downsample is consumed as input, not produced here -- see
    checkpoint admixture_tumor_downsample_batch above. A pure-background row
    (tumor_downsample_fraction == 0, the fraction_ladder's f=0 negative
    controls that tf_eval.compute_detection_floor's null pool needs) has no
    tumor_sub.bam at all -- admixture_tumor_sub_input/admixture_tumor_sub_cli_args
    (common.smk) leave that input empty and drop the --tumor-sub-bam/
    --tumor-total-reads flags for such a row, and AdmixtureBuilder.run
    (admixture.py) builds its mixture from the healthy downsample alone. Only
    the dedup BAM, its
    index, the per-estimator fraction files and sample_metadata.csv persist
    from this rule -- the healthy-side downsample and every intermediate
    merge/sort stage are temp() so a plan of hundreds of mixtures does not
    multiply the ~100GB-per-run peak-scratch cost documented in
    hcc1395-local-smoke-benchmark.md across the whole run.
    """
    input:
        plan=lambda wc: checkpoints.admixture_design.get(**wc).output.plan,
        tumor_sub=admixture_tumor_sub_input,
    output:
        dedup_bam=str(ADMIX_OUTPUT_ROOT / "{sample_name}" / "{sample_name}.dedup.bam"),
        dedup_bai=str(ADMIX_OUTPUT_ROOT / "{sample_name}" / "{sample_name}.dedup.bam.bai"),
        fractions=expand(
            str(ADMIX_OUTPUT_ROOT / "{{sample_name}}" / "{{sample_name}}.empirical_tumor_fraction.{label}.txt"),
            label=ADMIX_ESTIMATOR_LABELS,
        ),
        metadata=str(ADMIX_OUTPUT_ROOT / "{sample_name}" / "sample_metadata.csv"),
        healthy_sub=temp(str(ADMIX_OUTPUT_ROOT / "{sample_name}" / "{sample_name}.healthy_sub.bam")),
        merged=temp(str(ADMIX_OUTPUT_ROOT / "{sample_name}" / "{sample_name}.merged.bam")),
        namesort=temp(str(ADMIX_OUTPUT_ROOT / "{sample_name}" / "{sample_name}.namesort.bam")),
        fixmate=temp(str(ADMIX_OUTPUT_ROOT / "{sample_name}" / "{sample_name}.fixmate.bam")),
        sorted_bam=temp(str(ADMIX_OUTPUT_ROOT / "{sample_name}" / "{sample_name}.sorted.bam")),
    params:
        tumor_bam_path=lambda wc: _admixture_plan_row(wc.sample_name)["tumor_bam_path"],
        healthy_bam_path=lambda wc: _admixture_plan_row(wc.sample_name)["healthy_bam_path"],
        tumor_fraction=lambda wc: _admixture_plan_row(wc.sample_name)["tumor_fraction"],
        target_coverage=lambda wc: _admixture_plan_row(wc.sample_name)["target_coverage"],
        tumor_native_coverage=lambda wc: _admixture_plan_row(wc.sample_name)["tumor_native_coverage"],
        healthy_native_coverage=lambda wc: _admixture_plan_row(wc.sample_name)["healthy_native_coverage"],
        seed=lambda wc: _admixture_plan_row(wc.sample_name)["seed"],
        tumor_source_id=lambda wc: _admixture_plan_row(wc.sample_name)["tumor_source_id"],
        healthy_source_id=lambda wc: _admixture_plan_row(wc.sample_name)["healthy_source_id"],
        tumor_subject_id=lambda wc: _admixture_plan_row(wc.sample_name)["tumor_subject_id"],
        healthy_subject_id=lambda wc: _admixture_plan_row(wc.sample_name)["healthy_subject_id"],
        platform_tumor=lambda wc: _admixture_plan_row(wc.sample_name)["platform_tumor"],
        platform_healthy=lambda wc: _admixture_plan_row(wc.sample_name)["platform_healthy"],
        reference_build=lambda wc: _admixture_plan_row(wc.sample_name)["reference_build"],
        tumor_sub_args=admixture_tumor_sub_cli_args,
        output_dir=lambda wc, output: str(Path(output.dedup_bam).parent),
    log:
        str(ADMIX_OUTPUT_ROOT / "logs" / "admixture_mixture" / "{sample_name}.log"),
    conda:
        "flash_fraq"
    threads: ADMIX_THREADS
    resources:
        mem_mb=ADMIX_MEM_MB,
        runtime=ADMIX_RUNTIME,
    shell:
        f"python {SCRIPTS_DIR / 'run_admixture.py'} "
        "--tumor-bam {params.tumor_bam_path} "
        "{params.tumor_sub_args} "
        "--healthy-bam {params.healthy_bam_path} "
        f"--reference-fasta {ADMIX_REFERENCE_FASTA} "
        f"{ADMIX_SNV_VAF_ARGUMENTS} "
        "--output-dir {params.output_dir} "
        "--tumor-fraction {params.tumor_fraction} "
        "--target-coverage {params.target_coverage} "
        "--tumor-native-coverage {params.tumor_native_coverage} "
        "--healthy-native-coverage {params.healthy_native_coverage} "
        "--seed {params.seed} "
        "--sample-name {wildcards.sample_name} "
        "--tumor-source-id {params.tumor_source_id} "
        "--healthy-source-id {params.healthy_source_id} "
        "--tumor-subject-id {params.tumor_subject_id} "
        "--healthy-subject-id {params.healthy_subject_id} "
        "--platform-tumor {params.platform_tumor} "
        "--platform-healthy {params.platform_healthy} "
        "--reference-build {params.reference_build} "
        "--threads {threads} "
        "> {log} 2>&1"


rule admixture_mixtures_complete:
    input:
        all_admixture_sample_metadata,
    output:
        touch(str(ADMIX_OUTPUT_ROOT / "admixture_mixtures_complete.flag")),
    conda:
        "flash_fraq"
