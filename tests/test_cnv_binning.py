"""Tests for the IchorCNA/HMMcopy subprocess wrapper.

readCounter/gcCounter/mapCounter/runIchorCNA.R are not installed in this local environment
(they are Institute-HPC-only tools per phase2 §8), so `subprocess.run` is mocked at the
flash_fraq.cnv_binning.run boundary. These tests verify command construction, subprocess
failure handling, and output parsing — not the real tool behavior. An HPC test against a real
merged BAM and hg19 reference is a follow-up, not covered here.
"""

from dataclasses import replace
from pathlib import Path
from subprocess import CalledProcessError, CompletedProcess

import pytest
from pandas import DataFrame, read_csv

from flash_fraq.cnv_binning import IchorCnaConfig, IchorCnaRunner


@pytest.fixture
def config(tmp_path: Path) -> IchorCnaConfig:
    bam_path = tmp_path / "merged.bam"
    bam_path.write_bytes(b"")
    reference_dir = tmp_path / "ref"
    reference_dir.mkdir()
    (reference_dir / "hg19.fa").write_text(">chr1\nACGT\n")
    (reference_dir / "mappability_hg19_1000kb.bw").write_bytes(b"")
    output_dir = tmp_path / "out"
    return IchorCnaConfig.from_sources(bam_path=bam_path, reference_dir=reference_dir, output_dir=output_dir)


def test_from_sources_derives_expected_paths(config: IchorCnaConfig) -> None:
    assert config.sample_id == "merged"
    assert config.reference_fasta.name == "hg19.fa"
    assert config.mappability_bigwig.name == "mappability_hg19_1000kb.bw"
    assert config.bin_size_bp == 1_000_000
    assert config.read_wig_path == config.output_dir / "merged.readcounter.wig"
    assert config.cna_seg_path == config.output_dir / "merged.cna.seg"


def test_run_read_counter_builds_expected_command(config: IchorCnaConfig, monkeypatch: pytest.MonkeyPatch) -> None:
    config.output_dir.mkdir()
    captured_cmd: list[str] = []

    def fake_run(cmd: list[str], check: bool, capture_output: bool) -> CompletedProcess[bytes]:
        captured_cmd.extend(cmd)
        return CompletedProcess(cmd, returncode=0, stdout=b"fixedStep chrom=chr1 start=1 step=1000000\n1.0\n")

    monkeypatch.setattr("flash_fraq.cnv_binning.run", fake_run)

    wig_path = IchorCnaRunner(config).run_read_counter()

    assert captured_cmd == [
        "readCounter",
        "--window", "1000000",
        "--quality", "20",
        "--chromosome", ",".join(config.chromosomes),
        str(config.bam_path),
    ]
    assert wig_path == config.read_wig_path
    assert wig_path.read_bytes().startswith(b"fixedStep")


def test_run_ichor_cna_builds_expected_command(config: IchorCnaConfig, monkeypatch: pytest.MonkeyPatch) -> None:
    config.output_dir.mkdir()
    captured_cmd: list[str] = []

    def fake_run(cmd: list[str], check: bool, capture_output: bool) -> CompletedProcess[bytes]:
        captured_cmd.extend(cmd)
        return CompletedProcess(cmd, returncode=0, stdout=b"")

    monkeypatch.setattr("flash_fraq.cnv_binning.run", fake_run)

    cna_seg_path = IchorCnaRunner(config).run_ichor_cna()

    assert captured_cmd == [
        "Rscript",
        "runIchorCNA.R",
        "--id", "merged",
        "--WIG", str(config.read_wig_path),
        "--gcWig", str(config.gc_wig_path),
        "--mapWig", str(config.map_wig_path),
        "--outDir", str(config.output_dir),
        "--genomeBuild", "hg19",
        "--genomeStyle", "UCSC",
        "--chrs", "c(" + ",".join(repr(c) for c in config.chromosomes) + ")",
        "--chrTrain", "c(" + ",".join(repr(c) for c in config.chromosomes if c != "chrX") + ")",
        "--includeHOMD", "FALSE",
        "--ploidy", "c(2,3)",
        "--normal", "c(0.5,0.6,0.7,0.8,0.9)",
        "--scStates", "c(1,3)",
        "--estimateScPrevalence", "TRUE",
        "--maxCN", "5",
        "--txnE", "0.9999",
        "--txnStrength", "10000",
    ]
    assert cna_seg_path == config.cna_seg_path


def test_run_ichor_cna_passes_optional_libdir_centromere_and_normal_panel(
    config: IchorCnaConfig, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """All three are None by default; runIchorCNA.R must not receive a bare flag.

    libdir matters most: without it the script does `library(ichorCNA)` and
    silently binds to whatever version conda installed rather than the sources
    it shipped with.
    """
    config.output_dir.mkdir()
    libdir = tmp_path / "ichorCNA"
    centromere = tmp_path / "centromere.txt"
    normal_panel = tmp_path / "pon.rds"
    with_optionals = replace(config, libdir=libdir, centromere=centromere, normal_panel=normal_panel)
    captured_cmd: list[str] = []

    def fake_run(cmd: list[str], check: bool, capture_output: bool) -> CompletedProcess[bytes]:
        captured_cmd.extend(cmd)
        return CompletedProcess(cmd, returncode=0, stdout=b"")

    monkeypatch.setattr("flash_fraq.cnv_binning.run", fake_run)
    IchorCnaRunner(with_optionals).run_ichor_cna()

    assert captured_cmd[-6:] == [
        "--libdir", str(libdir),
        "--centromere", str(centromere),
        "--normalPanel", str(normal_panel),
    ]


def test_write_cnv_bins_reads_ichorcna_sample_prefixed_logr_column(
    config: IchorCnaConfig, tmp_path: Path
) -> None:
    """Real ichorCNA output names the column "<id>.logR", not "logR" -- a
    fixture using the bare name hid this until the first run against a real
    .cna.seg."""
    config.output_dir.mkdir()
    bins = DataFrame({"chr": ["chr1", "chr1"], "start": [1, 1_000_001], "merged.logR": [0.0146, 0.0215]})

    output_path = IchorCnaRunner(config).write_cnv_bins(bins)

    written = read_csv(output_path)
    assert list(written["cnv_logr"]) == [0.0146, 0.0215]
    assert list(written["bin_index"]) == [0, 1]


def test_write_cnv_bins_errors_clearly_when_no_logr_column(config: IchorCnaConfig) -> None:
    config.output_dir.mkdir()
    with pytest.raises(ValueError, match="no logR column"):
        IchorCnaRunner(config).write_cnv_bins(DataFrame({"chr": ["chr1"], "start": [1]}))


def test_bin_index_keys_on_position_so_filtered_bins_do_not_shift_samples(
    config: IchorCnaConfig,
) -> None:
    """ichorCNA drops bins its filters reject, so the row count differs between
    samples. bin_index has to name a locus, not a row: the chr2 bin below must
    get the same index whether or not the chr1 bins survived."""
    config.output_dir.mkdir()
    full = DataFrame(
        {"chr": ["chr1", "chr1", "chr2"], "start": [1, 1_000_001, 1], "merged.logR": [0.1, 0.2, 0.3]}
    )
    filtered = DataFrame({"chr": ["chr2"], "start": [1], "merged.logR": [0.3]})

    runner = IchorCnaRunner(config)
    assert runner._bin_indices(full)[-1] == runner._bin_indices(filtered)[0]

    # Same loci, 0-based starts: the ordinals must come out identical.
    zero_based = DataFrame(
        {"chr": ["chr1", "chr1", "chr2"], "start": [0, 1_000_000, 0], "merged.logR": [0.1, 0.2, 0.3]}
    )
    assert runner._bin_indices(zero_based) == runner._bin_indices(full)


def test_write_cnv_bins_rejects_a_chromosome_outside_the_configured_set(
    config: IchorCnaConfig,
) -> None:
    config.output_dir.mkdir()
    with pytest.raises(ValueError, match="outside the configured set"):
        IchorCnaRunner(config).write_cnv_bins(
            DataFrame({"chr": ["chrY"], "start": [1], "merged.logR": [0.1]})
        )


def test_prebuilt_wigs_bypass_gc_and_map_counter(
    config: IchorCnaConfig, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The local workstation path supplies ichorCNA's shipped hg19 1 Mb tracks,
    so gcCounter must not run and mapCounter -- which needs a bigWig that is not
    present locally -- must not run either."""
    gc_wig = tmp_path / "gc_hg19_1000kb.wig"
    map_wig = tmp_path / "map_hg19_1000kb.wig"
    with_prebuilt = replace(config, gc_wig_source=gc_wig, map_wig_source=map_wig)

    # gc/map_wig_path is what run_ichor_cna hands to runIchorCNA.R; it must
    # resolve to the prebuilt track, not the per-sample file gcCounter writes.
    assert with_prebuilt.gc_wig_path == gc_wig
    assert with_prebuilt.map_wig_path == map_wig

    invoked: list[str] = []

    def fake_run(cmd: list[str], check: bool, capture_output: bool) -> CompletedProcess[bytes]:
        invoked.append(Path(cmd[0]).name)
        return CompletedProcess(cmd, returncode=0, stdout=b"fixedStep chrom=chr1 start=1 step=1000000\n1.0\n")

    monkeypatch.setattr("flash_fraq.cnv_binning.run", fake_run)
    monkeypatch.setattr(
        IchorCnaRunner,
        "parse_cna_bins",
        lambda self: DataFrame({"chr": ["chr1", "chr1"], "start": [1, 1_000_001], "logR": [0.1, 0.2]}),
    )
    monkeypatch.setattr(IchorCnaRunner, "parse_tumor_fraction_estimate", lambda self: 0.12)

    IchorCnaRunner(with_prebuilt).run()

    assert "gcCounter" not in invoked
    assert "mapCounter" not in invoked
    assert invoked == ["readCounter", "Rscript"]


def test_subprocess_failure_is_surfaced_as_runtime_error(
    config: IchorCnaConfig, monkeypatch: pytest.MonkeyPatch
) -> None:
    def fake_run(cmd: list[str], check: bool, capture_output: bool) -> CompletedProcess[bytes]:
        raise CalledProcessError(returncode=1, cmd=cmd, stderr=b"contig not found")

    monkeypatch.setattr("flash_fraq.cnv_binning.run", fake_run)

    with pytest.raises(RuntimeError, match="readCounter failed \\(exit 1\\): contig not found"):
        IchorCnaRunner(config).run_read_counter()


def test_parse_cna_bins_reads_1mb_bin_output(config: IchorCnaConfig) -> None:
    config.output_dir.mkdir()
    config.cna_seg_path.write_text(
        "ID\tchr\tstart\tend\tcopy.number\tevent\tlogR\n"
        "merged\tchr1\t1\t1000000\t2\tNEUT\t0.01\n"
        "merged\tchr1\t1000001\t2000000\t3\tGAIN\t0.35\n"
    )

    bins = IchorCnaRunner(config).parse_cna_bins()

    assert bins.shape == (2, 7)
    assert list(bins["chr"]) == ["chr1", "chr1"]
    assert list(bins["copy.number"]) == [2, 3]
    assert bins["logR"].iloc[1] == pytest.approx(0.35)


def test_run_orchestrates_full_pipeline_and_returns_parsed_bins(
    config: IchorCnaConfig, monkeypatch: pytest.MonkeyPatch
) -> None:
    def fake_run(cmd: list[str], check: bool, capture_output: bool) -> CompletedProcess[bytes]:
        if cmd[0] == "Rscript":
            config.cna_seg_path.write_text(
                "ID\tchr\tstart\tend\tcopy.number\tevent\tlogR\n" "merged\tchr1\t1\t1000000\t2\tNEUT\t0.0\n"
            )
            config.params_path.write_text("Tumor Fraction: 0.23\nPloidy: 2.1\n")
        return CompletedProcess(cmd, returncode=0, stdout=b"fixedStep\n1.0\n")

    monkeypatch.setattr("flash_fraq.cnv_binning.run", fake_run)

    bins = IchorCnaRunner(config).run()

    assert config.output_dir.is_dir()
    assert config.read_wig_path.exists()
    assert config.gc_wig_path.exists()
    assert config.map_wig_path.exists()
    assert config.gc_wig_path == config.generated_gc_wig_path
    assert bins.shape == (1, 7)

    cnv_bins_written = read_csv(config.cnv_bins_output_path)
    assert list(cnv_bins_written.columns) == ["sample_id", "bin_index", "chrom", "start", "cnv_logr"]
    assert cnv_bins_written["cnv_logr"].iloc[0] == pytest.approx(0.0)

    ichor_estimate_written = read_csv(config.ichor_estimate_output_path)
    assert ichor_estimate_written["ichorCNA_tumor_fraction_estimate"].iloc[0] == pytest.approx(0.23)


def test_parse_tumor_fraction_estimate_reads_params_file(config: IchorCnaConfig) -> None:
    config.output_dir.mkdir()
    config.params_path.write_text("Tumor Fraction: 0.42\nPloidy: 1.9\n")

    assert IchorCnaRunner(config).parse_tumor_fraction_estimate() == pytest.approx(0.42)


def test_run_gc_and_map_counter_build_expected_commands(
    config: IchorCnaConfig, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The other two HMMcopy Suite binaries. run() skips both on the local
    workstation path (prebuilt tracks), so nothing else pins their argv."""
    config.output_dir.mkdir()
    captured: list[list[str]] = []

    def fake_run(cmd: list[str], check: bool, capture_output: bool) -> CompletedProcess[bytes]:
        captured.append(list(cmd))
        return CompletedProcess(cmd, returncode=0, stdout=b"fixedStep chrom=chr1 start=1 step=1000000\n0.4\n")

    monkeypatch.setattr("flash_fraq.cnv_binning.run", fake_run)
    runner = IchorCnaRunner(config)

    gc_wig_path = runner.run_gc_counter()
    map_wig_path = runner.run_map_counter()

    assert captured == [
        ["gcCounter", "--window", "1000000", str(config.reference_fasta)],
        ["mapCounter", "--window", "1000000", str(config.mappability_bigwig)],
    ]
    assert gc_wig_path == config.generated_gc_wig_path
    assert map_wig_path == config.generated_map_wig_path
    assert gc_wig_path.read_bytes().startswith(b"fixedStep")
    assert map_wig_path.read_bytes().startswith(b"fixedStep")


def test_generated_wigs_never_overwrite_a_configured_prebuilt_track(
    config: IchorCnaConfig, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """ichorCNA's shipped extdata WIGs are a shared, read-only asset: calling
    the counters explicitly with prebuilt tracks configured must still write to
    the per-sample paths, not back over the source."""
    config.output_dir.mkdir()
    gc_wig = tmp_path / "gc_hg19_1000kb.wig"
    map_wig = tmp_path / "map_hg19_1000kb.wig"
    gc_wig.write_text("shipped gc track\n")
    map_wig.write_text("shipped map track\n")
    runner = IchorCnaRunner(replace(config, gc_wig_source=gc_wig, map_wig_source=map_wig))

    monkeypatch.setattr(
        "flash_fraq.cnv_binning.run",
        lambda cmd, check, capture_output: CompletedProcess(cmd, returncode=0, stdout=b"regenerated\n"),
    )
    runner.run_gc_counter()
    runner.run_map_counter()

    assert gc_wig.read_text() == "shipped gc track\n"
    assert map_wig.read_text() == "shipped map track\n"
    assert config.generated_gc_wig_path.read_bytes() == b"regenerated\n"
    assert config.generated_map_wig_path.read_bytes() == b"regenerated\n"


def test_parse_tumor_fraction_estimate_errors_when_the_line_is_absent(config: IchorCnaConfig) -> None:
    config.output_dir.mkdir()
    config.params_path.write_text("Ploidy: 1.9\n")

    with pytest.raises(ValueError, match="no 'Tumor Fraction:' line found"):
        IchorCnaRunner(config).parse_tumor_fraction_estimate()
