"""End-to-end reconciliation check: runs MAT-13..18's actual writer methods against each
other (not hand-written fixtures) and confirms the real output of each stage is consumable
by the next, all the way through to the assembled feature matrix and the decoy QC gate.
"""

from pathlib import Path
from subprocess import run

import pysam
from pandas import DataFrame

from flash_fraq.admixture import AdmixtureBuilder, AdmixtureConfig
from flash_fraq.cnv_binning import IchorCnaConfig, IchorCnaRunner
from flash_fraq.decoy_qc import DecoyClassifierQc, DecoyClassifierQcConfig
from flash_fraq.feature_matrix import FeatureMatrixAssembler, FeatureMatrixConfig
from flash_fraq.fragment_extraction import FragmentExtractionConfig, FragmentExtractor
from flash_fraq.gc_correction import GcCorrectionConfig, GcCorrector
from test_admixture import _build_reference_fasta, _build_snv_vaf_file

CONTIG_LENGTH = 1000
SNV_POS = 500  # 1-based
READ_LENGTH = 60
FRAGMENT_LEN = 150  # within the DELFI "short" range (100-150bp)
REF_BASE = "A"
ALT_BASE = "G"


def _sam_pair(name: str, i: int, mate1_start: int, fragment_len: int, seq1: str, seq2: str, qual: str) -> list[str]:
    mate2_start = mate1_start + fragment_len - READ_LENGTH
    return [
        f"{name}_{i}\t99\tchr1\t{mate1_start}\t60\t{READ_LENGTH}M\t=\t{mate2_start}\t{fragment_len}\t{seq1}\t{qual}",
        f"{name}_{i}\t147\tchr1\t{mate2_start}\t60\t{READ_LENGTH}M\t=\t{mate1_start}\t{-fragment_len}\t{seq2}\t{qual}",
    ]


def _build_paired_source_bam(tmp_path: Path, name: str, n_pairs: int, alt_ref: str) -> Path:
    """Proper-pair reads (unlike test_admixture's single-end helper) so the merged BAM
    is usable by both MAT-13's mpileup step and MAT-14's fragment extraction filters."""
    sam_path = tmp_path / f"{name}.sam"
    header = f"@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:{CONTIG_LENGTH}\n"
    mate1_start = SNV_POS - 10  # 1-based; SNV falls within read1's span.
    base = ALT_BASE if alt_ref == "alt" else REF_BASE
    qual = "I" * READ_LENGTH
    lines = []
    for i in range(n_pairs):
        seq1 = [REF_BASE] * READ_LENGTH
        seq1[SNV_POS - mate1_start] = base
        lines += _sam_pair(name, i, mate1_start, FRAGMENT_LEN, "".join(seq1), REF_BASE * READ_LENGTH, qual)
    # Filler pairs elsewhere on the contig (long fragments, no SNV overlap) so a second
    # DELFI window has nonzero coverage too -- otherwise GC-correction's LOESS fit and
    # the short:long ratio both degenerate on a single populated window.
    for i in range(n_pairs):
        lines += _sam_pair(f"{name}_filler", i, 100, 220, REF_BASE * READ_LENGTH, REF_BASE * READ_LENGTH, qual)
    sam_path.write_text(header + "\n".join(lines) + "\n")
    bam_path = tmp_path / f"{name}.bam"
    run(["samtools", "sort", "-o", str(bam_path), str(sam_path)], check=True)
    run(["samtools", "index", str(bam_path)], check=True)
    return bam_path


def test_full_pipeline_reconciles_across_all_six_modules(tmp_path: Path, monkeypatch) -> None:
    shared_output_dir = tmp_path / "shared"
    sample_id = "mix_001"

    # --- MAT-13: admixture -> merged/dedup BAM + sample_metadata.csv ---
    tumor_bam = _build_paired_source_bam(tmp_path, "tumor", n_pairs=200, alt_ref="alt")
    healthy_bam = _build_paired_source_bam(tmp_path, "healthy", n_pairs=200, alt_ref="ref")
    reference_fasta = _build_reference_fasta(tmp_path)
    snv_vaf_path = _build_snv_vaf_file(tmp_path)

    admixture_config = AdmixtureConfig.from_sources(
        tumor_bam=tumor_bam,
        healthy_bam=healthy_bam,
        reference_fasta=reference_fasta,
        snv_vaf_paths={"sites": snv_vaf_path},
        output_dir=shared_output_dir,
        tumor_fraction=0.5,
        target_coverage=100,
        tumor_native_coverage=200,
        healthy_native_coverage=200,
        sample_name=sample_id,
    )
    empirical_tumor_fraction = AdmixtureBuilder(admixture_config).run()["sites"]
    dedup_bam = shared_output_dir / f"{sample_id}.dedup.bam"
    assert dedup_bam.exists()

    # --- MAT-14: fragment extraction -> window_fragment_counts.tsv + fragment_stats.csv ---
    delfi_windows_bed = tmp_path / "delfi_windows.bed"
    # Use the same tiny contig length as the synthetic reference built for MAT-13, split
    # into two windows so GC-correction's LOESS fit has more than a single data point.
    contig_length = pysam.AlignmentFile(str(dedup_bam)).header.get_reference_length("chr1")
    half = contig_length // 2
    delfi_windows_bed.write_text(f"chr1\t0\t{half}\twindow1\nchr1\t{half}\t{contig_length}\twindow2\n")
    ichor_bins_bed = tmp_path / "ichor_bins.bed"
    ichor_bins_bed.write_text(f"chr1\t0\t{contig_length}\tbin1\n")

    fragment_config = FragmentExtractionConfig.from_sources(
        input_bam=dedup_bam,
        output_dir=shared_output_dir,
        delfi_windows_bed=delfi_windows_bed,
        ichor_bins_bed=ichor_bins_bed,
        reference_fasta=reference_fasta,
        sample_id=sample_id,
    )
    FragmentExtractor(fragment_config).run_full()
    window_counts_path = shared_output_dir / "window_fragment_counts.tsv"
    assert window_counts_path.exists()

    # --- MAT-15: GC correction -> window_ratios.csv + arm_zscores.csv ---
    gc_config = GcCorrectionConfig.from_sources(
        input_windows_path=window_counts_path,
        output_windows_path=shared_output_dir / f"{sample_id}.windows_corrected.tsv",
        # This test's contig is a synthetic 1000 bp toy sequence, not real hg19 chr1 --
        # genome_mask's default margin (5 Mb) would mask every window on it. Masking
        # itself is covered by tests/test_gc_correction.py and tests/test_genome_mask.py.
        mask_margin_bp=0,
    )
    GcCorrector(gc_config).run_full()
    assert (shared_output_dir / "window_ratios.csv").exists()
    assert (shared_output_dir / "arm_zscores.csv").exists()

    # --- MAT-16: ichorCNA/HMMcopy wrapper (mocked subprocess) -> cnv_bins.csv + ichor_estimates.csv ---
    reference_dir = tmp_path / "ichor_ref"
    reference_dir.mkdir()
    (reference_dir / "hg19.fa").write_text(">chr1\nACGT\n")
    (reference_dir / "mappability_hg19_1000kb.bw").write_bytes(b"")
    cnv_config = IchorCnaConfig.from_sources(
        bam_path=dedup_bam, reference_dir=reference_dir, output_dir=shared_output_dir, sample_id=sample_id,
    )

    def fake_run(cmd: list[str], check: bool, capture_output: bool):
        from subprocess import CompletedProcess

        if cmd[0] == "Rscript":
            cnv_config.cna_seg_path.write_text(
                "ID\tchr\tstart\tend\tcopy.number\tevent\tlogR\n"
                f"{sample_id}\tchr1\t1\t1000000\t2\tNEUT\t0.01\n"
            )
            cnv_config.params_path.write_text("Tumor Fraction: 0.18\nPloidy: 2.0\n")
        return CompletedProcess(cmd, returncode=0, stdout=b"fixedStep\n1.0\n")

    monkeypatch.setattr("flash_fraq.cnv_binning.run", fake_run)
    IchorCnaRunner(cnv_config).run()
    assert (shared_output_dir / "cnv_bins.csv").exists()
    assert (shared_output_dir / "ichor_estimates.csv").exists()

    # split_assignment.csv is a hand-maintained subject-split manifest, not produced by
    # any pipeline stage (assigned before mixing, per phase0 §5.2) -- supplied directly here.
    DataFrame(
        {"subject_id": [admixture_config.tumor_subject_id, admixture_config.healthy_subject_id], "split_group": ["train", "train"]}
    ).to_csv(shared_output_dir / "split_assignment.csv", index=False)

    # --- MAT-17: assemble the real outputs of MAT-13..16 into the §9 feature matrix ---
    feature_matrix_config = FeatureMatrixConfig.from_sources(
        input_dir=shared_output_dir,
        output_path=shared_output_dir / "feature_matrix.csv",
        # Same synthetic-contig reasoning as gc_config above: the lone cnv bin here sits
        # at chr1:1, which the real hg19 chrom-size table would treat as chromosome-end-
        # adjacent and mask out.
        mask_margin_bp=0,
    )
    assembler = FeatureMatrixAssembler(feature_matrix_config)
    matrix = assembler.assemble()
    assembler.write(matrix)

    assert matrix.shape[0] == 1
    assert matrix["sample_id"].iloc[0] == sample_id
    assert matrix["empirical_tumor_fraction"].iloc[0] == empirical_tumor_fraction
    assert matrix["ichorCNA_tumor_fraction_estimate"].iloc[0] == 0.18
    # bin_index is derived from (chr, start), not row order: chr1's first 1 Mb
    # bin is chromosome rank 0, ordinal 0.
    assert matrix["cnv_logr_bin_0"].iloc[0] == 0.01
    assert (shared_output_dir / "feature_matrix.csv").exists()

    # --- MAT-18: decoy classifier QC gate runs against the real assembled matrix ---
    # The QC gate needs >= 2 source classes to fit a classifier; this pipeline run only
    # produced one real mixed sample, so append one synthetic second-platform row with
    # the real matrix's exact column layout to exercise the gate against real column names.
    second_row = matrix.iloc[0].copy()
    second_row["sample_id"] = "mix_002"
    second_row["platform_tumor"] = "other_platform"
    second_row["empirical_tumor_fraction"] = 0.35
    two_sample_matrix = DataFrame([matrix.iloc[0], second_row])
    two_sample_matrix.to_csv(feature_matrix_config.output_path, index=False)

    decoy_config = DecoyClassifierQcConfig.from_sources(
        feature_matrix_path=shared_output_dir / "feature_matrix.csv", output_dir=shared_output_dir,
    )
    report = DecoyClassifierQc(decoy_config).run()
    assert (shared_output_dir / "decoy_qc_report.json").exists()
    assert report.correlation_check is not None
