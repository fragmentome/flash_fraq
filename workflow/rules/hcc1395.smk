"""Retrieval and hg19 alignment of HCC1395/HCC1395BL (SEQC2) from SRA study
SRP162370 (docs/phase1-sourcing.md, docs/phase2-preprocessing-plan.md §0/§1b; MAT-23).

Per-sample SRR accessions aren't known until `resolve_hcc1395_runs` runs, so
downstream fetch/align targets are expanded dynamically off its output via a
checkpoint rather than a static sample list, per Snakemake's recommended
pattern for run-time-determined targets.

Globals below are prefixed HCC1395_ (rather than e.g. OUTPUT_ROOT) because
include: shares one Python namespace across all rule files pulled into
workflow/Snakefile -- a generic name would be silently reassigned by
whichever included file loads last, and any lambda that closes over it
(e.g. a params: callable) would then resolve to the wrong workflow's value
at execution time instead of the value in scope when this file was parsed.
"""

from pathlib import Path

HCC1395_CONFIG = config["hcc1395"]
HCC1395_OUTPUT_ROOT = Path(HCC1395_CONFIG["output_root"])
HCC1395_REFERENCE_FASTA = HCC1395_CONFIG["reference_fasta"]
HCC1395_THREADS = HCC1395_CONFIG.get("threads", 8)
# Separate from the align threads: a fetch job spends them on pigz (and
# fasterq-dump) only on the SRA-toolkit fallback path, so asking for an
# alignment-sized allocation would just idle cores on ENA-served runs.
HCC1395_FETCH_THREADS = HCC1395_CONFIG.get("fetch_threads", 8)
# fasterq-dump's memory use scales with --threads; see the config key's
# comment for the OOM this is fixing.
HCC1395_FETCH_MEM_MB_PER_THREAD = HCC1395_CONFIG.get("fetch_mem_mb_per_thread", 4000)
# Per-thread samtools sort cap; see the flag's comment in
# scripts/run_hcc1395_alignment.py for why this is deployment-configurable.
HCC1395_SORT_MEMORY_PER_THREAD = HCC1395_CONFIG.get("sort_memory_per_thread", "16G")


# --- Lane -> library merge + native coverage (MAT-23 Block C) ---------------
#
# Workstation-local, unlike everything above it: the inputs are the per-run
# CRAMs that scripts/sync_crams_to_hdd.sh mirrored onto the 8 TB HDD, not the
# isilon tree HCC1395_OUTPUT_ROOT points at. Kept in this file because it is
# the same chain (it is what turns this workflow's per-run products into the
# sample-level CRAMs admixture_design.yaml names), with its own config section
# so the cluster and workstation paths never share a key.
HCC1395MERGE_CONFIG = config["hcc1395_merge"]
HCC1395MERGE_CRAM_ROOT = Path(HCC1395MERGE_CONFIG["cram_root"])
HCC1395MERGE_DATASET = HCC1395MERGE_CONFIG["dataset"]
HCC1395MERGE_DERIVED_ROOT = Path(HCC1395MERGE_CONFIG["derived_root"])
HCC1395MERGE_REFERENCE_FASTA = HCC1395MERGE_CONFIG["reference_fasta"]
HCC1395MERGE_PREPS = tuple(HCC1395MERGE_CONFIG.get("preps", ("WGS",)))
HCC1395MERGE_PAIRS_ONLY = HCC1395MERGE_CONFIG.get("pairs_only", True)
HCC1395MERGE_MERGE_THREADS = HCC1395MERGE_CONFIG.get("merge_threads", 8)
# 4 to match rule localfeat_coverage; see the config key's comment.
HCC1395MERGE_MOSDEPTH_THREADS = HCC1395MERGE_CONFIG.get("mosdepth_threads", 4)

# Resolved against REPO_ROOT, never against the invocation directory: the
# drivers run Snakemake from snakemake_workdirs/<name>, where a relative
# "manifests/..." does not exist (CLAUDE.md, "Never write an
# invocation-relative path into the workflow").
HCC1395MERGE_PAIRS_MANIFEST = (
    Path(HCC1395MERGE_CONFIG["pairs_manifest"])
    if Path(HCC1395MERGE_CONFIG["pairs_manifest"]).is_absolute()
    else REPO_ROOT / HCC1395MERGE_CONFIG["pairs_manifest"]
)

# {sample_id: MergeUnit} for every library with at least one banked lane.
# Resolved once here rather than per job, for the reason common.smk documents.
HCC1395MERGE_UNITS = hcc1395_merge_units()

# A merged library's CRAM sits in the same CRAM/<id>/<id>.sorted.cram layout as
# the per-run lanes it is built from, so without this constraint Snakemake also
# offers rule merge_hcc1395_library as a way to *produce* a lane CRAM and then
# fails resolving that lane against HCC1395MERGE_UNITS. Pinned to the planned
# ids rather than to a shape, so a library with no banked lane is not offered
# either; "(?!)" is a never-matching regex for the no-mirror case, which is
# what a checkout without the HDD sees.
HCC1395MERGE_SAMPLE_CONSTRAINT = "|".join(HCC1395MERGE_UNITS) or "(?!)"

HCC1395MERGE_COVERAGE_DIR = HCC1395MERGE_DERIVED_ROOT / HCC1395MERGE_DATASET / "{sample_id}" / "coverage"
HCC1395MERGE_COVERAGE_REPORT = (
    HCC1395MERGE_DERIVED_ROOT / HCC1395MERGE_DATASET / "hcc1395_merged_native_coverage.csv"
)



checkpoint resolve_hcc1395_runs:
    output:
        run_list=str(HCC1395_OUTPUT_ROOT / "manifest" / "hcc1395_runs.csv"),
    params:
        sra_study_accession=HCC1395_CONFIG["sra_study_accession"],
    log:
        str(HCC1395_OUTPUT_ROOT / "logs" / "resolve_hcc1395_runs.log"),
    conda:
        "flash_fraq"
    shell:
        f"python {SCRIPTS_DIR / 'resolve_hcc1395_runs.py'} "
        "--sra-study-accession {params.sra_study_accession} "
        "--output {output.run_list} "
        "> {log} 2>&1"


rule fetch_hcc1395_fastq:
    output:
        r1=str(HCC1395_OUTPUT_ROOT / "fastq" / "{run}" / "{run}_1.fastq.gz"),
        r2=str(HCC1395_OUTPUT_ROOT / "fastq" / "{run}" / "{run}_2.fastq.gz"),
    log:
        str(HCC1395_OUTPUT_ROOT / "logs" / "fetch_hcc1395_fastq" / "{run}.log"),
    conda:
        "flash_fraq"
    threads: HCC1395_FETCH_THREADS
    resources:
        mem_mb=HCC1395_FETCH_THREADS * HCC1395_FETCH_MEM_MB_PER_THREAD,
        runtime=2880,
        # Caps concurrent fetches against the sra_fetch pool in
        # workflow/profiles/*/config.yaml -- a disk limit, not a bandwidth
        # one; see the comment there.
        sra_fetch=1,
    shell:
        f"python {SCRIPTS_DIR / 'fetch_hcc1395_fastq.py'} "
        "--run-accession {wildcards.run} "
        f"--output-dir {HCC1395_OUTPUT_ROOT / 'fastq'} "
        "--threads {threads} "
        "> {log} 2>&1"


rule align_hcc1395:
    input:
        r1=str(HCC1395_OUTPUT_ROOT / "fastq" / "{run}" / "{run}_1.fastq.gz"),
        r2=str(HCC1395_OUTPUT_ROOT / "fastq" / "{run}" / "{run}_2.fastq.gz"),
        # params: below read this manifest via _run_metadata(), so Snakemake
        # requires the checkpoint's output be declared an input here too --
        # otherwise requesting one <run>.sorted.bam/.cram directly fails to
        # build a DAG ("Rule parameter depends on checkpoint but checkpoint
        # output is not defined as input file for the rule"). Reaching these
        # outputs via hcc1395_alignment_complete masked that, since its input
        # function resolves the checkpoint first.
        run_list=str(HCC1395_OUTPUT_ROOT / "manifest" / "hcc1395_runs.csv"),
    output:
        bam=str(HCC1395_OUTPUT_ROOT / "bam" / "{run}.sorted.bam"),
    params:
        source_group=lambda wc: _run_metadata(wc.run)["source_group"],
        tumor_status=lambda wc: _run_metadata(wc.run)["tumor_status"],
        platform=lambda wc: _run_metadata(wc.run)["platform"],
    log:
        str(HCC1395_OUTPUT_ROOT / "logs" / "align_hcc1395" / "{run}.log"),
    conda:
        "flash_fraq"
    threads: HCC1395_THREADS
    resources:
        # 48G of samtools sort buffer (24 threads x 2G) + ~10G for bwa's hg19
        # index + headroom. The previous 160G was sized for the old implicit
        # 8 x 16G sort and made this rule a scheduling problem in its own
        # right: a 500G node fits 3 such jobs no matter how many CPUs are
        # free, so jobs queued on memory the workflow was never going to use.
        mem_mb=65536,
        runtime=2880,
        # Caps concurrent aligns against the bam_align pool in
        # workflow/profiles/*/config.yaml -- the sorted BAM stays on the
        # fileset until cram_hcc1395 replaces it; see the comment there.
        bam_align=1,
    shell:
        f"python {SCRIPTS_DIR / 'run_hcc1395_alignment.py'} "
        "--run-accession {wildcards.run} "
        "--source-group {params.source_group} "
        "--tumor-status {params.tumor_status} "
        "--platform {params.platform} "
        "--fastq-r1 {input.r1} "
        "--fastq-r2 {input.r2} "
        f"--reference-fasta {HCC1395_REFERENCE_FASTA} "
        f"--output-dir {HCC1395_OUTPUT_ROOT / 'bam'} "
        "--threads {threads} "
        f"--sort-memory-per-thread {HCC1395_SORT_MEMORY_PER_THREAD} "
        "> {log} 2>&1"


rule cram_hcc1395:
    input:
        bam=str(HCC1395_OUTPUT_ROOT / "bam" / "{run}.sorted.bam"),
    output:
        cram=str(HCC1395_OUTPUT_ROOT / "bam" / "{run}.sorted.cram"),
    log:
        str(HCC1395_OUTPUT_ROOT / "logs" / "cram_hcc1395" / "{run}.log"),
    conda:
        "flash_fraq"
    # 8, not 1: samtools view -C is the conversion's slow half and threads
    # roughly 2x at eight (19 -> 40 MB/s measured), and the two verification
    # flagstat passes thread ~8x. Until this finishes, the raw BAM is still on
    # the fileset holding a bam_align slot.
    threads: 8
    resources:
        mem_mb=16384,
        # 8h, not 4: at one thread a 350 GB HCC1395 BAM needed ~5.8h for the
        # encode alone, so the largest runs in the cohort were being killed on
        # walltime after hours of work. Threaded they land near 3h, and this
        # leaves margin for a contended fileset.
        runtime=480,
    shell:
        f"python {SCRIPTS_DIR / 'run_cram_conversion.py'} "
        "--bam-path {input.bam} "
        f"--reference-fasta {HCC1395_REFERENCE_FASTA} "
        "--threads {threads} "
        "> {log} 2>&1"


rule hcc1395_alignment_complete:
    input:
        all_aligned_crams,
    output:
        touch(str(HCC1395_OUTPUT_ROOT / "hcc1395_alignment_complete.flag")),
    conda:
        "flash_fraq"


rule merge_hcc1395_library:
    """Merge one SEQC2 library's top-up lanes into a sample-level CRAM.

    The merge unit is the library -- one (prep, site, arm, replicate) cell of
    manifests/hcc1395/hcc1395_library_pairs.csv -- not the cell line and not
    the site/arm. flash_fraq.hcc1395_merge's module docstring carries the
    justification; the short form is that merging across sites destroys the
    same-site tumor/normal pairing the 2026-09-04 finding made a hard
    requirement, merging across replicates blends distinct library preps into a
    fragment-length distribution no real library has, and neither buys depth
    that anything consumes (admixture_design asks for 35x/34x against
    C_target <= 1x, and every banked library clears that on its own lanes).
    """
    wildcard_constraints:
        sample_id=HCC1395MERGE_SAMPLE_CONSTRAINT,
    input:
        crams=lambda wc: [
            str(path)
            for path in HCC1395MERGE_UNITS[wc.sample_id].lane_cram_paths(
                HCC1395MERGE_CRAM_ROOT, HCC1395MERGE_DATASET
            )
        ],
        crais=lambda wc: [
            str(path)
            for path in HCC1395MERGE_UNITS[wc.sample_id].lane_crai_paths(
                HCC1395MERGE_CRAM_ROOT, HCC1395MERGE_DATASET
            )
        ],
    output:
        cram=str(HCC1395MERGE_CRAM_ROOT / HCC1395MERGE_DATASET / "CRAM" / "{sample_id}" / "{sample_id}.sorted.cram"),
        crai=str(
            HCC1395MERGE_CRAM_ROOT / HCC1395MERGE_DATASET / "CRAM" / "{sample_id}" / "{sample_id}.sorted.cram.crai"
        ),
    log:
        str(HCC1395MERGE_DERIVED_ROOT / "logs" / "merge_hcc1395_library" / "{sample_id}.log"),
    conda:
        "flash_fraq"
    threads: HCC1395MERGE_MERGE_THREADS
    resources:
        # samtools merge streams; the memory is the CRAM encoder's block
        # buffers per thread, not a sort buffer, so this does not scale with
        # the file size.
        mem_mb=8192,
        runtime=480,
    shell:
        f"python {SCRIPTS_DIR / 'run_hcc1395_merge.py'} "
        "--sample-id {wildcards.sample_id} "
        f"--pairs-manifest {HCC1395MERGE_PAIRS_MANIFEST} "
        f"--cram-root {HCC1395MERGE_CRAM_ROOT} "
        f"--dataset {HCC1395MERGE_DATASET} "
        f"--reference-fasta {HCC1395MERGE_REFERENCE_FASTA} "
        + "".join(f"--prep {prep} " for prep in HCC1395MERGE_PREPS)
        + "--threads {threads} "
        "> {log} 2>&1"


rule hcc1395_merged_coverage:
    """Native coverage of a merged library.

    The mosdepth invocation is character-for-character rule
    localfeat_coverage's -- same --by 1000000, same --no-per-base, same
    --threads, same --fasta -- and the output lands at the same
    <dataset>/<sample_id>/coverage/<sample_id>.mosdepth.summary.txt path in the
    derived tree. That is the entire point: the number this produces has to be
    comparable with the 547 rows already in native_coverage.csv, and it has to
    be readable by NativeCoverageHarvester without a special case.

    No SSD staging step, unlike local_features: a merged library is read once
    here, so the HDD->SSD copy would cost a second full pass to save nothing.
    """
    wildcard_constraints:
        sample_id=HCC1395MERGE_SAMPLE_CONSTRAINT,
    input:
        cram=str(HCC1395MERGE_CRAM_ROOT / HCC1395MERGE_DATASET / "CRAM" / "{sample_id}" / "{sample_id}.sorted.cram"),
        crai=str(
            HCC1395MERGE_CRAM_ROOT / HCC1395MERGE_DATASET / "CRAM" / "{sample_id}" / "{sample_id}.sorted.cram.crai"
        ),
    output:
        summary=str(HCC1395MERGE_COVERAGE_DIR / "{sample_id}.mosdepth.summary.txt"),
        # All five of mosdepth's outputs are declared for the reason
        # localfeat_coverage documents: Snakemake only cleans up what it knows
        # about, so an undeclared file survives a failed job and trips its retry.
        global_dist=str(HCC1395MERGE_COVERAGE_DIR / "{sample_id}.mosdepth.global.dist.txt"),
        region_dist=str(HCC1395MERGE_COVERAGE_DIR / "{sample_id}.mosdepth.region.dist.txt"),
        regions=str(HCC1395MERGE_COVERAGE_DIR / "{sample_id}.regions.bed.gz"),
        regions_index=str(HCC1395MERGE_COVERAGE_DIR / "{sample_id}.regions.bed.gz.csi"),
    params:
        prefix=lambda wc, output: str(Path(output.summary).parent / wc.sample_id),
    log:
        str(HCC1395MERGE_DERIVED_ROOT / "logs" / "hcc1395_merged_coverage" / "{sample_id}.log"),
    conda:
        # mosdepth is not in the flash_fraq env -- same env split as
        # rule localfeat_coverage.
        "mosdepth"
    threads: HCC1395MERGE_MOSDEPTH_THREADS
    resources:
        mem_mb=8192,
        runtime=300,
    shell:
        "mkdir -p $(dirname {params.prefix}) && "
        "mosdepth --threads {threads} --by 1000000 --no-per-base "
        f"--fasta {HCC1395MERGE_REFERENCE_FASTA} "
        "{params.prefix} {input.cram} > {log} 2>&1"


rule hcc1395_merged_coverage_report:
    """One CSV over every merged library's total-row mean, plus the exact lines
    to paste over admixture_design.yaml's 35.0/34.0 placeholders."""
    input:
        summaries=[
            str(HCC1395MERGE_DERIVED_ROOT / HCC1395MERGE_DATASET / sample_id / "coverage" / f"{sample_id}.mosdepth.summary.txt")
            for sample_id in HCC1395MERGE_UNITS
        ],
    output:
        report=str(HCC1395MERGE_COVERAGE_REPORT),
    log:
        str(HCC1395MERGE_DERIVED_ROOT / "logs" / "hcc1395_merged_coverage_report.log"),
    conda:
        "flash_fraq"
    threads: 1
    resources:
        mem_mb=2048,
        runtime=15,
    shell:
        f"python {SCRIPTS_DIR / 'run_hcc1395_merged_coverage.py'} "
        f"--pairs-manifest {HCC1395MERGE_PAIRS_MANIFEST} "
        f"--cram-root {HCC1395MERGE_CRAM_ROOT} "
        f"--dataset {HCC1395MERGE_DATASET} "
        f"--derived-root {HCC1395MERGE_DERIVED_ROOT} "
        f"--reference-fasta {HCC1395MERGE_REFERENCE_FASTA} "
        + "".join(f"--prep {prep} " for prep in HCC1395MERGE_PREPS)
        + "--output-path {output.report} "
        "> {log} 2>&1"


rule hcc1395_merge_complete:
    input:
        [
            str(HCC1395MERGE_CRAM_ROOT / HCC1395MERGE_DATASET / "CRAM" / sample_id / f"{sample_id}.sorted.cram")
            for sample_id in HCC1395MERGE_UNITS
        ],
    output:
        touch(str(HCC1395MERGE_CRAM_ROOT / HCC1395MERGE_DATASET / "hcc1395_merge_complete.flag")),
    conda:
        "flash_fraq"
