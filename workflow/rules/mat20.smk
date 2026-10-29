"""Download and alignment of GSE71378 (Snyder) FASTQ pairs to hg19, with @SQ
header verification (docs/phase2-preprocessing-plan.md §1a; MAT-20).

Sample titles are read straight from install_snyder_data.SNYDER_SAMPLES --
the fixed 60-sample GSE71378 title-to-SRX mapping -- so, like MAT-19, the
sample list is known at DAG-construction time and doesn't need a checkpoint.

Globals below are prefixed MAT20_ (rather than e.g. OUTPUT_ROOT) because
include: shares one Python namespace across all rule files pulled into
workflow/Snakefile -- a generic name would be silently reassigned by
whichever included file loads last, and any lambda that closes over it
(e.g. a params: callable) would then resolve to the wrong workflow's value
at execution time instead of the value in scope when this file was parsed.
"""

from pathlib import Path

from flash_fraq.install_snyder_data import SNYDER_SAMPLES

MAT20_CONFIG = config["mat20"]
MAT20_SNYDER_DOWNLOAD_DIR = MAT20_CONFIG["snyder_download_dir"]
MAT20_REFERENCE_FASTA = MAT20_CONFIG["reference_fasta"]
MAT20_OUTPUT_ROOT = Path(MAT20_CONFIG["output_root"])
MAT20_THREADS = MAT20_CONFIG.get("threads", 8)
MAT20_STAGING_DIR = MAT20_CONFIG["staging_dir"]
MAT20_STAGING_NODE = MAT20_CONFIG["staging_node"]

MAT20_SAMPLE_TITLES = list(SNYDER_SAMPLES)


rule install_align_snyder:
    output:
        bam=str(MAT20_OUTPUT_ROOT / "{sample_title}" / "{sample_title}.sorted.bam"),
    params:
        output_dir=lambda wc: str(MAT20_OUTPUT_ROOT / wc.sample_title),
    log:
        str(MAT20_OUTPUT_ROOT / "logs" / "install_align_snyder" / "{sample_title}.log"),
    conda:
        "flash_fraq"
    threads: MAT20_THREADS
    resources:
        mem_mb=196608,
        runtime=2880,
        # Pinned to the node with the large local scratch disk (staging_dir
        # below) -- concurrent samples were saturating isilon NFS throughput,
        # driving bwa mem alone past 24-44h real time on samples that should
        # take a few hours.
        slurm_extra=f"--nodelist={MAT20_STAGING_NODE}",
        # Caps concurrent ENA FASTQ downloads (see the ena_download pool in
        # workflow/profiles/*/config.yaml) -- up to `jobs: 30` of these were
        # hitting ENA's FTP at once, which reads as the cause of the
        # repeated identical wget failures across resubmissions.
        ena_download=1,
    shell:
        f"python {SCRIPTS_DIR / 'run_mat20_snyder_install.py'} "
        "--sample-title {wildcards.sample_title} "
        f"--download-dir {MAT20_SNYDER_DOWNLOAD_DIR} "
        f"--reference-fasta {MAT20_REFERENCE_FASTA} "
        "--output-dir {params.output_dir} "
        f"--staging-dir {MAT20_STAGING_DIR} "
        "--threads {threads} "
        "--sort-memory-per-thread 8G "
        "> {log} 2>&1"


rule cram_snyder:
    input:
        bam=str(MAT20_OUTPUT_ROOT / "{sample_title}" / "{sample_title}.sorted.bam"),
    output:
        cram=str(MAT20_OUTPUT_ROOT / "{sample_title}" / "{sample_title}.sorted.cram"),
    log:
        str(MAT20_OUTPUT_ROOT / "logs" / "cram_snyder" / "{sample_title}.log"),
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
        f"--reference-fasta {MAT20_REFERENCE_FASTA} "
        "--threads {threads} "
        "> {log} 2>&1"


rule mat20_install_complete:
    input:
        expand(str(MAT20_OUTPUT_ROOT / "{sample_title}" / "{sample_title}.sorted.cram"), sample_title=MAT20_SAMPLE_TITLES),
    output:
        touch(str(MAT20_OUTPUT_ROOT / "mat20_install_complete.flag")),
    conda:
        "flash_fraq"
