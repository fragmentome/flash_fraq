"""MAT-14/15/16: workstation-local feature extraction from hg19 CRAMs
(docs/phase2-preprocessing-plan.md §8b.3, CRAM edition).

Consumes the CRAM mirror that scripts/sync_crams_to_hdd.sh writes to the 8 TB
HDD and produces, per sample, the MAT-14/15/16 tables that feature_matrix.py
joins: fragment_stats.csv, window_ratios.csv, arm_zscores.csv, cnv_bins.csv,
ichor_estimates.csv.

**SSD staging is the point of this file.** §8b.4 measured the HDD at ~150 MB/s
saturated and ~27.5 MB/s while a transfer runs, which is why the plan insisted
stages 3/5/6 share a single sequential pass -- three passes over a spinning
disk cost ~29 h instead of ~10 h for the full cohort. Staging each CRAM onto
the NVMe SSD first inverts that trade: one HDD-sequential read per sample, and
every subsequent pass (coverage, fragments, readCounter) hits the SSD, where
random-ish re-reads are effectively free. So the rules below deliberately do
NOT fuse into one pass -- they stay separate, restartable, and independently
cacheable, which the single-pass design could not be.

Both staged copies are temp(): Snakemake deletes a sample's SSD footprint as
soon as its last consumer finishes, so peak SSD usage tracks the number of
concurrent samples, not the cohort. Cap that with --resources ssd_mb=... (see
workflow/profiles/local/config.yaml).

Globals are prefixed LOCALFEAT_ because include: shares one Python namespace
across every rule file pulled into workflow/Snakefile -- an unprefixed name
would be silently reassigned by whichever file is included last.
"""

from pathlib import Path

from flash_fraq.genome_mask import DEFAULT_MASK_MARGIN_BP

LOCALFEAT_CONFIG = config["local_features"]
LOCALFEAT_CRAM_ROOT = Path(LOCALFEAT_CONFIG["cram_root"])
# Entries may be a bare name or a {name, sample_pattern} mapping; see
# common.smk:dataset_name(). LOCALFEAT_DATASET_NAMES is the directory-name
# view the wildcard constraint needs.
LOCALFEAT_DATASETS = LOCALFEAT_CONFIG["datasets"]
LOCALFEAT_DATASET_NAMES = dataset_names(LOCALFEAT_DATASETS)
LOCALFEAT_OUTPUT_ROOT = Path(LOCALFEAT_CONFIG["output_root"])
LOCALFEAT_STAGING_ROOT = Path(LOCALFEAT_CONFIG["staging_root"])
LOCALFEAT_REFERENCE_FASTA = LOCALFEAT_CONFIG["reference_fasta"]
LOCALFEAT_PRIMARY_CONTIGS = LOCALFEAT_CONFIG["primary_contigs"]
LOCALFEAT_DELFI_WINDOW_BP = LOCALFEAT_CONFIG["delfi_window_bp"]
LOCALFEAT_ICHOR_BIN_BP = LOCALFEAT_CONFIG["ichor_bin_bp"]
LOCALFEAT_EXTDATA = Path(LOCALFEAT_CONFIG["ichorcna_extdata"])
LOCALFEAT_ICHORCNA_SCRIPT = LOCALFEAT_CONFIG["ichorcna_script"]
LOCALFEAT_SAMTOOLS_THREADS = LOCALFEAT_CONFIG.get("samtools_threads", 4)
LOCALFEAT_MOSDEPTH_THREADS = LOCALFEAT_CONFIG.get("mosdepth_threads", 4)
# Centromere/chromosome-end mask (flash_fraq.genome_mask); 0 disables. Not required in
# local_features.yaml -- falls back to genome_mask.DEFAULT_MASK_MARGIN_BP (5 Mb) so an
# existing config file need not be touched to pick this up.
LOCALFEAT_MASK_MARGIN_BP = LOCALFEAT_CONFIG.get("mask_margin_bp", DEFAULT_MASK_MARGIN_BP)

LOCALFEAT_WINDOWS_DIR = LOCALFEAT_OUTPUT_ROOT / "windows"
LOCALFEAT_DELFI_BED = LOCALFEAT_WINDOWS_DIR / f"delfi_windows_{LOCALFEAT_DELFI_WINDOW_BP}.bed"
LOCALFEAT_ICHOR_BED = LOCALFEAT_WINDOWS_DIR / f"ichor_bins_{LOCALFEAT_ICHOR_BIN_BP}.bed"

# Resolved once at DAG construction, not per-job: the mirror is append-only
# while a sync runs, so a mid-run re-glob would grow the DAG underneath us.
LOCALFEAT_SAMPLES = discover_local_crams(LOCALFEAT_CRAM_ROOT, LOCALFEAT_DATASETS)

LOCALFEAT_DERIVED = LOCALFEAT_OUTPUT_ROOT / "{dataset}" / "{sample_id}"
LOCALFEAT_STAGED = LOCALFEAT_STAGING_ROOT / "{dataset}" / "{sample_id}"


wildcard_constraints:
    # Sample IDs are bare accessions/titles (IA01, EGAF00002727132); pinning
    # this keeps the two-level {dataset}/{sample_id} paths from being split at
    # the wrong slash.
    dataset="|".join(LOCALFEAT_DATASET_NAMES),
    sample_id=r"[A-Za-z0-9._-]+",


rule localfeat_make_windows:
    """DELFI 5 Mb windows + ichorCNA 1 Mb bins, primary contigs only (§3, §4).

    Neither BED existed on this machine -- fragment_extraction.py has always
    taken them as inputs, and only the test suite ever synthesised them.
    """
    input:
        fai=f"{LOCALFEAT_REFERENCE_FASTA}.fai",
    output:
        delfi=str(LOCALFEAT_DELFI_BED),
        ichor=str(LOCALFEAT_ICHOR_BED),
    log:
        str(LOCALFEAT_OUTPUT_ROOT / "logs" / "localfeat_make_windows.log"),
    conda:
        "flash_fraq"
    threads: 1
    resources:
        mem_mb=2048,
        runtime=10,
    shell:
        f"python {SCRIPTS_DIR / 'make_genome_windows.py'} "
        "--fai {input.fai} "
        f"--delfi-window-bp {LOCALFEAT_DELFI_WINDOW_BP} "
        f"--ichor-bin-bp {LOCALFEAT_ICHOR_BIN_BP} "
        "--delfi-bed {output.delfi} "
        "--ichor-bed {output.ichor} "
        f"--contigs {' '.join(LOCALFEAT_PRIMARY_CONTIGS)} "
        "> {log} 2>&1"


rule localfeat_stage_cram:
    """HDD -> SSD. The one sequential HDD read per sample; everything downstream
    reads the SSD copy.

    ssd_mb is a plain Snakemake resource, not a real reservation -- it exists
    so `--resources ssd_mb=N` caps how many samples can be staged at once and
    keeps the 611 GB SSD from filling when many jobs run in parallel.
    """
    input:
        cram=str(LOCALFEAT_CRAM_ROOT / "{dataset}" / "CRAM" / "{sample_id}" / "{sample_id}.sorted.cram"),
        crai=str(LOCALFEAT_CRAM_ROOT / "{dataset}" / "CRAM" / "{sample_id}" / "{sample_id}.sorted.cram.crai"),
    output:
        cram=temp(str(LOCALFEAT_STAGED / "{sample_id}.sorted.cram")),
        crai=temp(str(LOCALFEAT_STAGED / "{sample_id}.sorted.cram.crai")),
    log:
        str(LOCALFEAT_OUTPUT_ROOT / "logs" / "localfeat_stage_cram" / "{dataset}.{sample_id}.log"),
    conda:
        "flash_fraq"
    threads: 1
    resources:
        mem_mb=2048,
        runtime=120,
        ssd_mb=lambda wc: local_cram_size_mb(LOCALFEAT_CRAM_ROOT, wc.dataset, wc.sample_id),
    shell:
        # cp, not rsync: the source is a local mount and a partial staged copy
        # must never be resumed -- Snakemake deletes and re-runs on failure.
        "cp --no-preserve=mode {input.cram} {output.cram} > {log} 2>&1 && "
        "cp --no-preserve=mode {input.crai} {output.crai} >> {log} 2>&1 && "
        # rsync preserved the cluster's mtimes, so a mirrored .crai is usually
        # older than its CRAM and samtools warns about a stale index on every
        # read (§8b.3 step 2). Stamping the staged copy keeps that warning
        # meaningful instead of constant.
        "touch {output.crai} >> {log} 2>&1"


rule localfeat_verify_cram:
    """§8b.3 step 1 + the truncation gate.

    quickcheck rejects a CRAM whose EOF container is missing (an interrupted
    rsync --partial); the @SQ comparison closes the phase2 "reference FASTA
    variant" open item by proving these CRAMs were aligned against the exact
    contig set in hg19.fa.fai -- which HCC1395's later alignment has to match
    for `samtools merge` to work at admixture time.
    """
    input:
        cram=str(LOCALFEAT_STAGED / "{sample_id}.sorted.cram"),
        crai=str(LOCALFEAT_STAGED / "{sample_id}.sorted.cram.crai"),
        fai=f"{LOCALFEAT_REFERENCE_FASTA}.fai",
    output:
        report=str(LOCALFEAT_DERIVED / "cram_header_verification.json"),
    log:
        str(LOCALFEAT_OUTPUT_ROOT / "logs" / "localfeat_verify_cram" / "{dataset}.{sample_id}.log"),
    conda:
        "flash_fraq"
    threads: 1
    resources:
        mem_mb=4096,
        runtime=15,
    shell:
        f"python {SCRIPTS_DIR / 'run_verify_cram_header.py'} "
        "--cram-path {input.cram} "
        "--reference-fai {input.fai} "
        f"--reference-fasta {LOCALFEAT_REFERENCE_FASTA} "
        "--sample-id {wildcards.sample_id} "
        "--output-path {output.report} "
        "> {log} 2>&1"


rule localfeat_coverage:
    """§8b.3 step 3 -- native coverage, which feeds the §7.1 downsampling
    fractions and closes the open BH01-coverage item."""
    input:
        cram=str(LOCALFEAT_STAGED / "{sample_id}.sorted.cram"),
        crai=str(LOCALFEAT_STAGED / "{sample_id}.sorted.cram.crai"),
        report=str(LOCALFEAT_DERIVED / "cram_header_verification.json"),
    output:
        summary=str(LOCALFEAT_DERIVED / "coverage" / "{sample_id}.mosdepth.summary.txt"),
        # mosdepth writes five files off one prefix. Declaring only the summary
        # leaves the other four undeclared, so a job that dies partway through
        # leaves them behind for the retry to trip over -- Snakemake only
        # cleans outputs it knows about.
        global_dist=str(LOCALFEAT_DERIVED / "coverage" / "{sample_id}.mosdepth.global.dist.txt"),
        region_dist=str(LOCALFEAT_DERIVED / "coverage" / "{sample_id}.mosdepth.region.dist.txt"),
        regions=str(LOCALFEAT_DERIVED / "coverage" / "{sample_id}.regions.bed.gz"),
        regions_index=str(LOCALFEAT_DERIVED / "coverage" / "{sample_id}.regions.bed.gz.csi"),
    params:
        prefix=lambda wc, output: str(Path(output.summary).parent / wc.sample_id),
    log:
        str(LOCALFEAT_OUTPUT_ROOT / "logs" / "localfeat_coverage" / "{dataset}.{sample_id}.log"),
    conda:
        # mosdepth is not in the flash_fraq env (§8b.2 lists it in its own).
        "mosdepth"
    threads: LOCALFEAT_MOSDEPTH_THREADS
    resources:
        mem_mb=8192,
        runtime=180,
    shell:
        "mkdir -p $(dirname {params.prefix}) && "
        "mosdepth --threads {threads} --by 1000000 --no-per-base "
        f"--fasta {LOCALFEAT_REFERENCE_FASTA} "
        "{params.prefix} {input.cram} > {log} 2>&1"


rule localfeat_fragments:
    """§8b.3 step 5 -- the pysam TLEN/midpoint pass (MAT-14).

    §8b.4 measured a single pysam worker at ~56 MB/s, well under the HDD's
    ~150 MB/s, and concluded 3-4 workers were needed to saturate the spindle.
    Reading from the SSD instead makes each worker CPU-bound rather than
    IO-bound, so scale this with `jobs`, not with per-job threads.
    """
    input:
        cram=str(LOCALFEAT_STAGED / "{sample_id}.sorted.cram"),
        crai=str(LOCALFEAT_STAGED / "{sample_id}.sorted.cram.crai"),
        report=str(LOCALFEAT_DERIVED / "cram_header_verification.json"),
        delfi=str(LOCALFEAT_DELFI_BED),
        ichor=str(LOCALFEAT_ICHOR_BED),
    output:
        window_counts=str(LOCALFEAT_DERIVED / "window_fragment_counts.tsv"),
        fragment_stats=str(LOCALFEAT_DERIVED / "fragment_stats.csv"),
        # fragment_midpoints.bed is one line per fragment -- tens of GB for a
        # deep sample, and nothing downstream of this rule reads it.
        midpoints=temp(str(LOCALFEAT_DERIVED / "fragment_midpoints.bed")),
        delfi_intersected=temp(str(LOCALFEAT_DERIVED / "fragments_delfi_5mb.bed")),
        ichor_intersected=temp(str(LOCALFEAT_DERIVED / "fragments_ichor_1mb.bed")),
    params:
        output_dir=lambda wc, output: str(Path(output.window_counts).parent),
    log:
        str(LOCALFEAT_OUTPUT_ROOT / "logs" / "localfeat_fragments" / "{dataset}.{sample_id}.log"),
    conda:
        "flash_fraq"
    threads: 1
    resources:
        mem_mb=16384,
        runtime=600,
    shell:
        f"python {SCRIPTS_DIR / 'run_local_fragment_extraction.py'} "
        "--input-alignment {input.cram} "
        "--output-dir {params.output_dir} "
        "--delfi-windows-bed {input.delfi} "
        "--ichor-bins-bed {input.ichor} "
        f"--reference-fasta {LOCALFEAT_REFERENCE_FASTA} "
        "--sample-id {wildcards.sample_id} "
        "> {log} 2>&1"


rule localfeat_gc_correct:
    """§8b.3 step 6 (DELFI half) -- MAT-15 LOESS GC correction. Pure pandas /
    statsmodels over the small per-window table, so it never touches the CRAM.

    priority=1: the ILP scheduler's objective rewards core utilization and
    temp() release, and this rule offers neither (1 thread, non-temp input) --
    without a priority bump it loses every selection round to 4-thread
    coverage/cram_to_bam jobs whenever enough of those are ready to fill the
    core budget, starving it indefinitely even though it's the cheapest rule
    in the DAG.
    """
    input:
        window_counts=str(LOCALFEAT_DERIVED / "window_fragment_counts.tsv"),
    output:
        window_ratios=str(LOCALFEAT_DERIVED / "window_ratios.csv"),
        arm_zscores=str(LOCALFEAT_DERIVED / "arm_zscores.csv"),
        corrected=str(LOCALFEAT_DERIVED / "window_fragment_counts.gc_corrected.tsv"),
    log:
        str(LOCALFEAT_OUTPUT_ROOT / "logs" / "localfeat_gc_correct" / "{dataset}.{sample_id}.log"),
    conda:
        "flash_fraq"
    threads: 1
    priority: 1
    resources:
        mem_mb=8192,
        runtime=60,
    shell:
        f"python {SCRIPTS_DIR / 'run_local_gc_correction.py'} "
        "--input-windows-path {input.window_counts} "
        "--output-windows-path {output.corrected} "
        "--window-ratios-path {output.window_ratios} "
        "--arm-zscore-path {output.arm_zscores} "
        f"--mask-margin-bp {LOCALFEAT_MASK_MARGIN_BP} "
        "> {log} 2>&1"


rule localfeat_cram_to_bam:
    """readCounter (HMMcopy Suite) is BamTools-based and cannot open a CRAM --
    it fails with 'BgzfStream ERROR: read block failed - invalid block header'.

    Rather than reimplement 1 Mb read counting (§4 forbids that: HMMcopy stays
    an external subprocess), decode one temp BAM on the SSD for readCounter's
    sole benefit. Every other consumer reads the CRAM directly, so this exists
    only for rule localfeat_ichor and is deleted with it.
    """
    input:
        cram=str(LOCALFEAT_STAGED / "{sample_id}.sorted.cram"),
        crai=str(LOCALFEAT_STAGED / "{sample_id}.sorted.cram.crai"),
        report=str(LOCALFEAT_DERIVED / "cram_header_verification.json"),
    output:
        bam=temp(str(LOCALFEAT_STAGED / "{sample_id}.sorted.bam")),
        bai=temp(str(LOCALFEAT_STAGED / "{sample_id}.sorted.bam.bai")),
    log:
        str(LOCALFEAT_OUTPUT_ROOT / "logs" / "localfeat_cram_to_bam" / "{dataset}.{sample_id}.log"),
    conda:
        "flash_fraq"
    threads: LOCALFEAT_SAMTOOLS_THREADS
    resources:
        mem_mb=8192,
        runtime=240,
        # A BAM runs roughly 2.5-3x its CRAM; budget 3x plus headroom.
        ssd_mb=lambda wc: local_cram_size_mb(LOCALFEAT_CRAM_ROOT, wc.dataset, wc.sample_id, multiplier=3.0),
    shell:
        "samtools view -@ {threads} -b "
        f"-T {LOCALFEAT_REFERENCE_FASTA} "
        "-o {output.bam} {input.cram} > {log} 2>&1 && "
        "samtools index -@ {threads} {output.bam} >> {log} 2>&1"


rule localfeat_ichor:
    """§8b.3 step 6 (CN half) -- MAT-16 readCounter + runIchorCNA.R.

    gcCounter/mapCounter are skipped: ichorCNA ships hg19 1 Mb GC and
    mappability WIGs, gcCounter would regenerate a byte-identical GC file for
    every sample, and mapCounter cannot run without a mappability bigWig that
    is not on this machine.
    """
    input:
        bam=str(LOCALFEAT_STAGED / "{sample_id}.sorted.bam"),
        bai=str(LOCALFEAT_STAGED / "{sample_id}.sorted.bam.bai"),
        gc_wig=str(LOCALFEAT_EXTDATA / LOCALFEAT_CONFIG["gc_wig"]),
        map_wig=str(LOCALFEAT_EXTDATA / LOCALFEAT_CONFIG["map_wig"]),
        centromere=str(LOCALFEAT_EXTDATA / LOCALFEAT_CONFIG["centromere"]),
    output:
        cnv_bins=str(LOCALFEAT_DERIVED / "ichor" / "cnv_bins.csv"),
        ichor_estimates=str(LOCALFEAT_DERIVED / "ichor" / "ichor_estimates.csv"),
    params:
        output_dir=lambda wc, output: str(Path(output.cnv_bins).parent),
    log:
        str(LOCALFEAT_OUTPUT_ROOT / "logs" / "localfeat_ichor" / "{dataset}.{sample_id}.log"),
    conda:
        # readCounter comes from TSO500ichorCNA and Rscript+ichorCNA from
        # problem2 (§8b.2 flags the split); the driver resolves both by
        # absolute path, so this env only needs Python.
        "flash_fraq"
    threads: 1
    resources:
        mem_mb=16384,
        runtime=300,
    shell:
        f"python {SCRIPTS_DIR / 'run_local_cnv_binning.py'} "
        "--alignment-path {input.bam} "
        f"--reference-fasta {LOCALFEAT_REFERENCE_FASTA} "
        "--output-dir {params.output_dir} "
        "--sample-id {wildcards.sample_id} "
        "--gc-wig {input.gc_wig} "
        "--map-wig {input.map_wig} "
        "--centromere {input.centromere} "
        f"--ichorcna-script {LOCALFEAT_ICHORCNA_SCRIPT} "
        f"--ichorcna-libdir {LOCALFEAT_CONFIG['ichorcna_libdir']} "
        f"--readcounter-bin {LOCALFEAT_CONFIG['readcounter_bin']} "
        f"--rscript-bin {LOCALFEAT_CONFIG['rscript_bin']} "
        f"--bin-size-bp {LOCALFEAT_ICHOR_BIN_BP} "
        f"--genome-build {LOCALFEAT_CONFIG['genome_build']} "
        f"--genome-style {LOCALFEAT_CONFIG['genome_style']} "
        f"--chromosomes {','.join(LOCALFEAT_PRIMARY_CONTIGS)} "
        # Quoted: these are R-vector literals containing parentheses and commas.
        f"--ploidy '{LOCALFEAT_CONFIG['ichor_ploidy']}' "
        f"--normal '{LOCALFEAT_CONFIG['ichor_normal']}' "
        f"--sc-states '{LOCALFEAT_CONFIG['ichor_sc_states']}' "
        f"--max-cn {LOCALFEAT_CONFIG['ichor_max_cn']} "
        f"--txn-e {LOCALFEAT_CONFIG['ichor_txn_e']} "
        f"--txn-strength {LOCALFEAT_CONFIG['ichor_txn_strength']} "
        + (f"--normal-panel {LOCALFEAT_CONFIG['normal_panel']} " if LOCALFEAT_CONFIG.get("normal_panel") else "")
        + "> {log} 2>&1"


rule localfeat_sample_complete:
    """Per-sample join point. Its existence is what lets Snakemake free the
    sample's staged CRAM/BAM -- temp() outputs survive until every consuming
    job has run, and this is the last of them."""
    input:
        report=str(LOCALFEAT_DERIVED / "cram_header_verification.json"),
        coverage=str(LOCALFEAT_DERIVED / "coverage" / "{sample_id}.mosdepth.summary.txt"),
        fragment_stats=str(LOCALFEAT_DERIVED / "fragment_stats.csv"),
        window_ratios=str(LOCALFEAT_DERIVED / "window_ratios.csv"),
        arm_zscores=str(LOCALFEAT_DERIVED / "arm_zscores.csv"),
        cnv_bins=str(LOCALFEAT_DERIVED / "ichor" / "cnv_bins.csv"),
        ichor_estimates=str(LOCALFEAT_DERIVED / "ichor" / "ichor_estimates.csv"),
    output:
        touch(str(LOCALFEAT_DERIVED / "localfeat_sample_complete.flag")),
    conda:
        "flash_fraq"


rule local_features_complete:
    input:
        [
            str(LOCALFEAT_OUTPUT_ROOT / dataset / sample_id / "localfeat_sample_complete.flag")
            for dataset, sample_id in LOCALFEAT_SAMPLES
        ],
    output:
        touch(str(LOCALFEAT_OUTPUT_ROOT / "local_features_complete.flag")),
    conda:
        "flash_fraq"
