"""Alignment of Cristiano/DELFI FASTQ pairs to hg19, with @SQ header
verification (docs/phase2-preprocessing-plan.md §1b; MAT-19).

Sample IDs are read straight off disk via InputManifestBuilder's existing
discover_cristiano_fastq() glob (docs/phase2-preprocessing-plan.md §1b) --
unlike MAT-23's SRA-resolved runs, every Cristiano FASTQ pair already exists
locally, so the sample list is known at DAG-construction time and doesn't
need a checkpoint.

Globals below are prefixed MAT19_ (rather than e.g. OUTPUT_ROOT) because
include: shares one Python namespace across all rule files pulled into
workflow/Snakefile -- a generic name would be silently reassigned by
whichever included file loads last, and any lambda that closes over it
(e.g. a params: callable) would then resolve to the wrong workflow's value
at execution time instead of the value in scope when this file was parsed.
"""

from pathlib import Path

from flash_fraq.input_manifest import InputManifestBuilder, InputManifestConfig

MAT19_CONFIG = config["mat19"]
MAT19_CRISTIANO_FASTQ_ROOT = Path(MAT19_CONFIG["cristiano_fastq_root"])
MAT19_REFERENCE_FASTA = MAT19_CONFIG["reference_fasta"]
MAT19_OUTPUT_ROOT = Path(MAT19_CONFIG["output_root"])
MAT19_THREADS = MAT19_CONFIG.get("threads", 8)

MAT19_SAMPLE_IDS = [
    entry.sample_id
    for entry in InputManifestBuilder(
        InputManifestConfig.from_sources(
            output_path=MAT19_OUTPUT_ROOT / "manifest" / "cristiano_manifest.csv",
            cristiano_fastq_root=MAT19_CRISTIANO_FASTQ_ROOT,
        )
    ).discover_cristiano_fastq()
]


rule align_cristiano:
    output:
        bam=str(MAT19_OUTPUT_ROOT / "{sample_id}" / "{sample_id}.sorted.bam"),
    params:
        output_dir=lambda wc: str(MAT19_OUTPUT_ROOT / wc.sample_id),
    log:
        str(MAT19_OUTPUT_ROOT / "logs" / "align_cristiano" / "{sample_id}.log"),
    conda:
        "flash_fraq"
    threads: MAT19_THREADS
    resources:
        mem_mb=163840,
        runtime=2880,
    shell:
        f"python {SCRIPTS_DIR / 'run_mat19_cristiano_alignment.py'} "
        "--sample-id {wildcards.sample_id} "
        f"--cristiano-fastq-root {MAT19_CRISTIANO_FASTQ_ROOT} "
        f"--reference-fasta {MAT19_REFERENCE_FASTA} "
        "--output-dir {params.output_dir} "
        "--threads {threads} "
        "> {log} 2>&1"


rule cram_cristiano:
    input:
        bam=str(MAT19_OUTPUT_ROOT / "{sample_id}" / "{sample_id}.sorted.bam"),
    output:
        cram=str(MAT19_OUTPUT_ROOT / "{sample_id}" / "{sample_id}.sorted.cram"),
    log:
        str(MAT19_OUTPUT_ROOT / "logs" / "cram_cristiano" / "{sample_id}.log"),
    conda:
        "flash_fraq"
    # 8, not 1: samtools view -C is the conversion's slow half and threads
    # roughly 2x at eight (19 -> 40 MB/s measured), and the two verification
    # flagstat passes thread ~8x. Until this finishes, the raw BAM is still on
    # the fileset holding a bam_align slot.
    threads: 8
    resources:
        mem_mb=16384,
        runtime=240,
    shell:
        f"python {SCRIPTS_DIR / 'run_cram_conversion.py'} "
        "--bam-path {input.bam} "
        f"--reference-fasta {MAT19_REFERENCE_FASTA} "
        "--threads {threads} "
        "> {log} 2>&1"


rule mat19_alignment_complete:
    input:
        expand(str(MAT19_OUTPUT_ROOT / "{sample_id}" / "{sample_id}.sorted.cram"), sample_id=MAT19_SAMPLE_IDS),
    output:
        touch(str(MAT19_OUTPUT_ROOT / "mat19_alignment_complete.flag")),
    conda:
        "flash_fraq"
