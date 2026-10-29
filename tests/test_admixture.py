from dataclasses import replace
from pathlib import Path
from subprocess import CompletedProcess, run

import pytest

from flash_fraq.admixture import AdmixtureBuilder, AdmixtureConfig

CONTIG = "chr1"
CONTIG_LENGTH = 1000
SNV_POS = 500  # 1-based
READ_LENGTH = 60
REF_BASE = "A"
ALT_BASE = "G"


def _write_bam(sam_lines: list[str], header: str, out_bam: Path) -> Path:
    sam_path = out_bam.with_suffix(".sam")
    sam_path.write_text(header + "\n".join(sam_lines) + "\n")
    run(["samtools", "sort", "-o", str(out_bam), str(sam_path)], check=True)
    run(["samtools", "index", str(out_bam)], check=True)
    return out_bam


def _build_source_bam(tmp_path: Path, name: str, n_reads: int, alt_ref: str) -> Path:
    """alt_ref: 'alt' -> every read carries ALT_BASE at SNV_POS (tumor, VAF=1.0 pure);
    'ref' -> every read carries REF_BASE (healthy background, VAF=0)."""
    header = f"@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:{CONTIG}\tLN:{CONTIG_LENGTH}\n"
    lo = max(1, SNV_POS - READ_LENGTH + 1)
    hi = SNV_POS
    lines = []
    for i in range(n_reads):
        start = lo + (i % (hi - lo + 1))
        bases = [REF_BASE] * READ_LENGTH
        if alt_ref == "alt":
            bases[SNV_POS - start] = ALT_BASE
        seq = "".join(bases)
        qual = "I" * READ_LENGTH
        lines.append(f"{name}_{i}\t0\t{CONTIG}\t{start}\t60\t{READ_LENGTH}M\t*\t0\t0\t{seq}\t{qual}")
    return _write_bam(lines, header, tmp_path / f"{name}.bam")


def _build_reference_fasta(tmp_path: Path) -> Path:
    fasta_path = tmp_path / "ref.fa"
    sequence = REF_BASE * CONTIG_LENGTH
    with open(fasta_path, "w") as handle:
        handle.write(f">{CONTIG}\n")
        for i in range(0, len(sequence), 70):
            handle.write(sequence[i:i + 70] + "\n")
    run(["samtools", "faidx", str(fasta_path)], check=True)
    return fasta_path


def _build_snv_vaf_file(tmp_path: Path) -> Path:
    vaf_path = tmp_path / "snvs.tsv"
    vaf_path.write_text(f"chrom\tpos\tref\talt\tvaf_tumor_pure\n{CONTIG}\t{SNV_POS}\t{REF_BASE}\t{ALT_BASE}\t1.0\n")
    return vaf_path


@pytest.fixture
def config(tmp_path: Path) -> AdmixtureConfig:
    tumor_bam = _build_source_bam(tmp_path, "tumor", n_reads=200, alt_ref="alt")
    healthy_bam = _build_source_bam(tmp_path, "healthy", n_reads=200, alt_ref="ref")
    reference_fasta = _build_reference_fasta(tmp_path)
    snv_vaf_path = _build_snv_vaf_file(tmp_path)
    output_dir = tmp_path / "out"
    return AdmixtureConfig.from_sources(
        tumor_bam=tumor_bam,
        healthy_bam=healthy_bam,
        reference_fasta=reference_fasta,
        snv_vaf_paths={"sites": snv_vaf_path},
        output_dir=output_dir,
        tumor_fraction=0.5,
        target_coverage=100,
        tumor_native_coverage=200,
        healthy_native_coverage=200,
        seed=42,
        sample_name="test_mix",
    )


def test_compute_downsample_fractions(config: AdmixtureConfig) -> None:
    builder = AdmixtureBuilder(config)
    tumor_fraction, healthy_fraction = builder.compute_downsample_fractions()
    assert tumor_fraction == pytest.approx(0.25)
    assert healthy_fraction == pytest.approx(0.25)
    assert tumor_fraction <= 1
    assert healthy_fraction <= 1


def test_compute_downsample_fractions_raises_when_native_coverage_too_low(config: AdmixtureConfig) -> None:
    from dataclasses import replace

    low_coverage_config = replace(config, tumor_native_coverage=10)
    builder = AdmixtureBuilder(low_coverage_config)
    with pytest.raises(ValueError, match="tumor downsample fraction"):
        builder.compute_downsample_fractions()


def test_full_admixture_pipeline_matches_independent_count(config: AdmixtureConfig, tmp_path: Path) -> None:
    builder = AdmixtureBuilder(config)
    empirical_tumor_fractions = builder.run()
    empirical_tumor_fraction = empirical_tumor_fractions["sites"]

    dedup_bam = config.output_dir / f"{config.sample_name}.dedup.bam"
    assert dedup_bam.exists()

    result_path = config.output_dir / f"{config.sample_name}.empirical_tumor_fraction.sites.txt"
    assert result_path.exists()
    assert float(result_path.read_text().strip()) == pytest.approx(empirical_tumor_fraction)

    # Independent oracle: count ref/alt bases at SNV_POS directly from the
    # non-duplicate reads in the merged+dedup BAM, bypassing bcftools entirely.
    view_result = run(
        ["samtools", "view", "-F", "0x400", str(dedup_bam)],
        check=True, capture_output=True, text=True,
    )
    alt_count = 0
    ref_count = 0
    for line in view_result.stdout.splitlines():
        fields = line.split("\t")
        pos = int(fields[3])
        seq = fields[9]
        if pos <= SNV_POS <= pos + len(seq) - 1:
            base = seq[SNV_POS - pos]
            if base == ALT_BASE:
                alt_count += 1
            elif base == REF_BASE:
                ref_count += 1

    depth = alt_count + ref_count
    assert depth > 0
    expected_empirical_tumor_fraction = alt_count / depth  # vaf_tumor_pure == 1.0

    assert empirical_tumor_fraction == pytest.approx(expected_empirical_tumor_fraction, abs=1e-6)
    # Sanity check against the nominal admixture ratio (large N, generous tolerance).
    assert empirical_tumor_fraction == pytest.approx(0.5, abs=0.2)


def test_run_writes_sample_metadata_matching_feature_matrix_schema(config: AdmixtureConfig) -> None:
    from pandas import read_csv

    builder = AdmixtureBuilder(config)
    empirical_tumor_fractions = builder.run()

    metadata_path = config.output_dir / "sample_metadata.csv"
    assert metadata_path.exists()

    metadata = read_csv(metadata_path)
    assert list(metadata.columns) == [
        "sample_id", "tumor_source_id", "healthy_source_id", "platform_tumor", "platform_healthy",
        "reference_build", "target_coverage", "tumor_subject_id", "healthy_subject_id",
        "nominal_tumor_fraction", "empirical_tumor_fraction", "empirical_tumor_fraction_sites",
        "mtdna_fraction_log10",
    ]
    assert metadata["sample_id"].iloc[0] == config.sample_name
    assert metadata["empirical_tumor_fraction"].iloc[0] == pytest.approx(empirical_tumor_fractions["sites"])
    assert metadata["empirical_tumor_fraction_sites"].iloc[0] == pytest.approx(empirical_tumor_fractions["sites"])
    assert metadata["nominal_tumor_fraction"].iloc[0] == pytest.approx(config.tumor_fraction)


def test_run_computes_both_snv_estimators_when_two_snv_sources_given(
    config: AdmixtureConfig, tmp_path: Path
) -> None:
    from dataclasses import replace

    from pandas import read_csv

    two_source_config = replace(
        config, snv_vaf_paths={"sites": config.snv_vaf_paths["sites"], "clonal": config.snv_vaf_paths["sites"]}
    )
    builder = AdmixtureBuilder(two_source_config)
    empirical_tumor_fractions = builder.run()

    assert set(empirical_tumor_fractions) == {"sites", "clonal"}
    assert empirical_tumor_fractions["sites"] == pytest.approx(empirical_tumor_fractions["clonal"])

    metadata = read_csv(two_source_config.output_dir / "sample_metadata.csv")
    assert "empirical_tumor_fraction_sites" in metadata.columns
    assert "empirical_tumor_fraction_clonal" in metadata.columns
    # "sites" is the first snv_vaf_paths entry, so it is the canonical column
    # feature_matrix.REGRESSION_LABEL_COLUMN reads.
    assert metadata["empirical_tumor_fraction"].iloc[0] == pytest.approx(empirical_tumor_fractions["sites"])
    assert (two_source_config.output_dir / f"{two_source_config.sample_name}.empirical_tumor_fraction.clonal.txt").exists()


def test_load_snvs_dispatches_on_vcf_extension(tmp_path: Path, config: AdmixtureConfig) -> None:
    vcf_path = tmp_path / "snvs.vcf"
    header = (
        '##fileformat=VCFv4.2\n'
        '##INFO=<ID=bwaTVAF,Number=1,Type=Float,Description="tumor VAF from bwa data">\n'
        f"##contig=<ID={CONTIG},length={CONTIG_LENGTH}>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    )
    record = f"{CONTIG}\t{SNV_POS}\t.\t{REF_BASE}\t{ALT_BASE}\t.\t.\tbwaTVAF=1\n"
    vcf_path.write_text(header + record)

    builder = AdmixtureBuilder(config)
    snvs = builder.load_snvs(vcf_path)
    assert len(snvs) == 1
    assert snvs[0].chrom == CONTIG
    assert snvs[0].pos == SNV_POS
    assert snvs[0].ref == REF_BASE
    assert snvs[0].alt == ALT_BASE
    assert snvs[0].vaf_tumor_pure == pytest.approx(1.0)


def test_load_snvs_vcf_skips_a_record_with_a_missing_bwatvaf_instead_of_raising(
    tmp_path: Path, config: AdmixtureConfig
) -> None:
    """bcftools query renders a missing INFO/bwaTVAF as '.' -- confirmed against
    the real SEQC2 sites.vcf.gz that this actually happens (1/39141 sites there),
    crashing every single admixture_mixture row since they all share that one
    VCF. One missing annotation among thousands of sites should not invalidate
    the whole row: _load_snvs_vcf must skip that record and keep parsing the
    rest, not raise ValueError."""
    vcf_path = tmp_path / "snvs_with_missing_vaf.vcf"
    header = (
        '##fileformat=VCFv4.2\n'
        '##INFO=<ID=bwaTVAF,Number=1,Type=Float,Description="tumor VAF from bwa data">\n'
        f"##contig=<ID={CONTIG},length={CONTIG_LENGTH}>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    )
    missing_vaf_pos = SNV_POS - 10
    valid_record = f"{CONTIG}\t{SNV_POS}\t.\t{REF_BASE}\t{ALT_BASE}\t.\t.\tbwaTVAF=0.8\n"
    missing_vaf_record = f"{CONTIG}\t{missing_vaf_pos}\t.\t{REF_BASE}\t{ALT_BASE}\t.\t.\t.\n"
    vcf_path.write_text(header + missing_vaf_record + valid_record)

    builder = AdmixtureBuilder(config)
    snvs = builder.load_snvs(vcf_path)

    # The missing-bwaTVAF site is skipped entirely; the valid one still parses.
    assert len(snvs) == 1
    assert snvs[0].pos == SNV_POS
    assert snvs[0].vaf_tumor_pure == pytest.approx(0.8)


def test_verify_downsample_depth_raises_outside_tolerance(config: AdmixtureConfig, tmp_path: Path) -> None:
    builder = AdmixtureBuilder(config)
    config.output_dir.mkdir(parents=True, exist_ok=True)
    tumor_sub = builder.subsample_bam(
        config.tumor_bam, 0.25, config.output_dir / "tumor_sub.bam"
    )
    with pytest.raises(ValueError, match="deviates from"):
        builder.verify_downsample_depth(config.tumor_bam, tumor_sub, expected_fraction=0.9, tolerance=0.05)


def test_config_defaults_threads_to_eight() -> None:
    assert AdmixtureConfig.__dataclass_fields__["threads"].default == 8


def test_subsample_bam_passes_threads_to_samtools_view(
    config: AdmixtureConfig, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    threaded_config = replace(config, threads=6)
    threaded_config.output_dir.mkdir(parents=True, exist_ok=True)
    captured_cmd: list[str] = []

    def fake_run(cmd: list[str], check: bool) -> CompletedProcess[bytes]:
        captured_cmd.extend(cmd)
        return CompletedProcess(cmd, returncode=0)

    monkeypatch.setattr("flash_fraq.admixture.run", fake_run)

    builder = AdmixtureBuilder(threaded_config)
    builder.subsample_bam(threaded_config.tumor_bam, 0.5, threaded_config.output_dir / "sub.bam")

    assert captured_cmd[:2] == ["samtools", "view"]
    assert "-@" in captured_cmd
    assert captured_cmd[captured_cmd.index("-@") + 1] == "6"


def test_merge_and_dedup_passes_threads_to_every_subprocess_call(
    config: AdmixtureConfig, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    threaded_config = replace(config, threads=5)
    threaded_config.output_dir.mkdir(parents=True, exist_ok=True)
    captured_cmds: list[list[str]] = []

    def fake_run(cmd: list[str], check: bool) -> CompletedProcess[bytes]:
        captured_cmds.append(list(cmd))
        return CompletedProcess(cmd, returncode=0)

    monkeypatch.setattr("flash_fraq.admixture.run", fake_run)

    builder = AdmixtureBuilder(threaded_config)
    fake_bam = threaded_config.output_dir / "sub.bam"
    fake_bam.write_bytes(b"")
    builder.merge_and_dedup([fake_bam, fake_bam])

    # samtools merge, sort (namesort), fixmate, sort (coord), markdup, index --
    # every one of the six calls in merge_and_dedup should carry -@ threads.
    assert len(captured_cmds) == 6
    for cmd in captured_cmds:
        assert "-@" in cmd, cmd
        assert cmd[cmd.index("-@") + 1] == "5"


def test_compute_empirical_tumor_fraction_passes_threads_to_bcftools_mpileup(
    config: AdmixtureConfig, monkeypatch: pytest.MonkeyPatch
) -> None:
    from flash_fraq.admixture import SomaticSnv

    threaded_config = replace(config, threads=3)
    threaded_config.output_dir.mkdir(parents=True, exist_ok=True)
    captured_cmd: list[str] = []

    fake_mpileup_line = "\t".join(
        [CONTIG, str(SNV_POS), REF_BASE, "*", ALT_BASE, ".", ".", ".", "DP:AD", "10:5,5"]
    )

    def fake_run(cmd: list[str], check: bool, capture_output: bool, text: bool) -> CompletedProcess[str]:
        captured_cmd.extend(cmd)
        return CompletedProcess(cmd, returncode=0, stdout=fake_mpileup_line + "\n")

    monkeypatch.setattr("flash_fraq.admixture.run", fake_run)

    builder = AdmixtureBuilder(threaded_config)
    snv = SomaticSnv(chrom=CONTIG, pos=SNV_POS, ref=REF_BASE, alt=ALT_BASE, vaf_tumor_pure=1.0)
    fraction = builder.compute_empirical_tumor_fraction(
        threaded_config.output_dir / "dedup.bam", [snv], label="sites"
    )

    assert fraction == pytest.approx(0.5)
    assert captured_cmd[:2] == ["bcftools", "mpileup"]
    assert "--threads" in captured_cmd
    assert captured_cmd[captured_cmd.index("--threads") + 1] == "3"


def _fake_mpileup_stdout(sites: list[tuple[int, int, int]]) -> str:
    """sites: (pos, depth, alt_count), all on CONTIG with a single ALT_BASE allele
    and a single sample column -- enough to drive `compute_empirical_tumor_fraction`'s
    per-site AD parsing without a real BAM/mpileup round trip."""
    lines = []
    for pos, depth, alt_count in sites:
        ref_count = depth - alt_count
        lines.append(
            "\t".join(
                [CONTIG, str(pos), REF_BASE, "*", ALT_BASE, ".", ".", ".", "DP:AD", f"{depth}:{ref_count},{alt_count}"]
            )
        )
    return "\n".join(lines) + "\n"


def test_compute_empirical_tumor_fraction_estimator_unbiased_and_lower_variance(
    config: AdmixtureConfig, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Synthetic pileup at a known f=0.05, with per-site VAF spread down to 0.04 and
    depth down to 1x -- the low-depth, low-VAF regime docs/results-corruption-audit.md's
    finding 7 flags as the one the old unweighted mean of per-site ratios blows up on.
    Both estimators are unbiased under this exact-VAF binomial model (E[alt/depth] = p
    for a binomial draw, whatever depth is), so the property that actually distinguishes
    them is variance: the new depth-weighted ratio of sums must have noticeably less of
    it, and its mean must still track the true f."""
    import numpy

    from flash_fraq.admixture import SomaticSnv

    f_true = 0.05
    # (pos, depth, vaf_tumor_pure) -- depth 1-2 and vaf as low as 0.04 mirror the
    # audit's bwaTVAF p1/median figures at 1x coverage.
    site_specs = [
        (100, 2, 0.04),
        (200, 1, 0.058),
        (300, 3, 0.10),
        (400, 2, 0.20),
        (500, 5, 0.335),
        (600, 10, 0.50),
    ]
    snvs = [
        SomaticSnv(chrom=CONTIG, pos=pos, ref=REF_BASE, alt=ALT_BASE, vaf_tumor_pure=vaf)
        for pos, _depth, vaf in site_specs
    ]

    config.output_dir.mkdir(parents=True, exist_ok=True)
    rng = numpy.random.default_rng(42)
    n_replicates = 400
    new_estimates = []
    old_estimates = []
    for _ in range(n_replicates):
        alt_counts = [
            int(rng.binomial(depth, min(f_true * vaf, 1.0))) for _pos, depth, vaf in site_specs
        ]
        fake_stdout = _fake_mpileup_stdout(
            [(pos, depth, alt_count) for (pos, depth, _vaf), alt_count in zip(site_specs, alt_counts)]
        )

        def fake_run(cmd: list[str], check: bool, capture_output: bool, text: bool) -> CompletedProcess[str]:
            return CompletedProcess(cmd, returncode=0, stdout=fake_stdout)

        monkeypatch.setattr("flash_fraq.admixture.run", fake_run)

        builder = AdmixtureBuilder(config)
        new_estimates.append(builder.compute_empirical_tumor_fraction(
            config.output_dir / "dedup.bam", snvs, label="sites"
        ))

        old_ratios = [
            (alt_count / depth) / vaf for (_pos, depth, vaf), alt_count in zip(site_specs, alt_counts)
        ]
        old_estimates.append(sum(old_ratios) / len(old_ratios))

    new_mean = sum(new_estimates) / n_replicates
    old_mean = sum(old_estimates) / n_replicates
    new_variance = sum((x - new_mean) ** 2 for x in new_estimates) / n_replicates
    old_variance = sum((x - old_mean) ** 2 for x in old_estimates) / n_replicates

    assert new_mean == pytest.approx(f_true, abs=0.01)
    assert old_mean == pytest.approx(f_true, abs=0.03)
    assert new_variance < old_variance / 5  # far lower variance, not just marginally lower


def test_compute_empirical_tumor_fraction_writes_readable_per_site_dump(
    config: AdmixtureConfig, monkeypatch: pytest.MonkeyPatch
) -> None:
    from pandas import read_csv

    from flash_fraq.admixture import SomaticSnv

    config.output_dir.mkdir(parents=True, exist_ok=True)
    site_specs = [(100, 2, 1), (200, 4, 0)]  # (pos, depth, alt_count)
    snvs = [
        SomaticSnv(chrom=CONTIG, pos=pos, ref=REF_BASE, alt=ALT_BASE, vaf_tumor_pure=0.5)
        for pos, _depth, _alt_count in site_specs
    ]
    fake_stdout = _fake_mpileup_stdout(site_specs)

    def fake_run(cmd: list[str], check: bool, capture_output: bool, text: bool) -> CompletedProcess[str]:
        return CompletedProcess(cmd, returncode=0, stdout=fake_stdout)

    monkeypatch.setattr("flash_fraq.admixture.run", fake_run)

    builder = AdmixtureBuilder(config)
    fraction = builder.compute_empirical_tumor_fraction(
        config.output_dir / "dedup.bam", snvs, label="sites"
    )
    assert fraction == pytest.approx((1 + 0) / (2 * 0.5 + 4 * 0.5))

    sites_path = config.output_dir / f"{config.sample_name}.empirical_tumor_fraction.sites.sites.tsv.gz"
    assert sites_path.exists()
    dumped = read_csv(sites_path, sep="\t")
    assert list(dumped.columns) == ["chrom", "pos", "alt_count", "depth", "vaf_tumor_pure"]
    assert list(dumped["pos"]) == [100, 200]
    assert list(dumped["alt_count"]) == [1, 0]
    assert list(dumped["depth"]) == [2, 4]
    assert list(dumped["vaf_tumor_pure"]) == pytest.approx([0.5, 0.5])


def test_verify_downsample_depth_passes_threads_to_samtools_flagstat(
    config: AdmixtureConfig, monkeypatch: pytest.MonkeyPatch
) -> None:
    threaded_config = replace(config, threads=7)
    captured_cmds: list[list[str]] = []

    def fake_run(cmd: list[str], check: bool, capture_output: bool, text: bool) -> CompletedProcess[str]:
        captured_cmds.append(list(cmd))
        return CompletedProcess(cmd, returncode=0, stdout="100 + 0 in total\n")

    monkeypatch.setattr("flash_fraq.admixture.run", fake_run)

    builder = AdmixtureBuilder(threaded_config)
    builder.verify_downsample_depth(Path("source.bam"), Path("sub.bam"), expected_fraction=1.0, tolerance=1.0)

    assert len(captured_cmds) == 2
    for cmd in captured_cmds:
        assert cmd[:2] == ["samtools", "flagstat"]
        assert "-@" in cmd
        assert cmd[cmd.index("-@") + 1] == "7"


def test_compute_mtdna_fraction_log10_passes_threads_to_samtools_idxstats(
    config: AdmixtureConfig, monkeypatch: pytest.MonkeyPatch
) -> None:
    threaded_config = replace(config, threads=2)
    captured_cmd: list[str] = []

    def fake_run(cmd: list[str], check: bool, capture_output: bool, text: bool) -> CompletedProcess[str]:
        captured_cmd.extend(cmd)
        return CompletedProcess(cmd, returncode=0, stdout=f"{CONTIG}\t{CONTIG_LENGTH}\t100\t0\nchrM\t16569\t5\t0\n")

    monkeypatch.setattr("flash_fraq.admixture.run", fake_run)

    builder = AdmixtureBuilder(threaded_config)
    builder.compute_mtdna_fraction_log10(Path("dedup.bam"))

    assert captured_cmd[:2] == ["samtools", "idxstats"]
    assert "-@" in captured_cmd
    assert captured_cmd[captured_cmd.index("-@") + 1] == "2"


def test_verify_downsample_depth_accepts_a_precomputed_reads_before(
    config: AdmixtureConfig, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A caller (run()'s precomputed_tumor_total_reads, wired from
    tumor_batch_sampling.TumorBatchSampler's single-pass byproduct) can pass
    reads_before directly -- this should skip the source-side flagstat call
    entirely, leaving only the one flagstat against the subsampled BAM."""
    captured_cmds: list[list[str]] = []

    def fake_run(cmd: list[str], check: bool, capture_output: bool, text: bool) -> CompletedProcess[str]:
        captured_cmds.append(list(cmd))
        return CompletedProcess(cmd, returncode=0, stdout="50 + 0 in total\n")

    monkeypatch.setattr("flash_fraq.admixture.run", fake_run)

    builder = AdmixtureBuilder(config)
    actual_fraction = builder.verify_downsample_depth(
        Path("source.bam"), Path("sub.bam"), expected_fraction=0.5, tolerance=0.01, reads_before=100
    )

    assert actual_fraction == pytest.approx(0.5)
    # Only the subsampled-BAM flagstat ran -- not one against "source.bam".
    assert len(captured_cmds) == 1
    assert captured_cmds[0][:2] == ["samtools", "flagstat"]
    assert captured_cmds[0][-1] == "sub.bam"


def test_run_uses_a_precomputed_tumor_sub_bam_instead_of_calling_subsample_bam(
    config: AdmixtureConfig,
) -> None:
    """run()'s precomputed_tumor_sub -- wired from
    tumor_batch_sampling.TumorBatchSampler's one-pass output in the
    Snakemake DAG -- should be used directly as the tumor side of the merge,
    with no independent samtools view -s against the full tumor source."""
    builder = AdmixtureBuilder(config)
    config.output_dir.mkdir(parents=True, exist_ok=True)

    precomputed_tumor_sub = config.output_dir / "precomputed.tumor_sub.bam"
    builder.subsample_bam(config.tumor_bam, 0.25, precomputed_tumor_sub)

    subsample_bam_calls = []
    original_subsample_bam = builder.subsample_bam

    def counting_subsample_bam(source_bam, fraction, output_bam):
        subsample_bam_calls.append(source_bam)
        return original_subsample_bam(source_bam, fraction, output_bam)

    builder.subsample_bam = counting_subsample_bam

    empirical_tumor_fractions = builder.run(precomputed_tumor_sub=precomputed_tumor_sub)

    assert "sites" in empirical_tumor_fractions
    # Only the healthy side went through subsample_bam -- the tumor side used
    # the precomputed BAM directly.
    assert subsample_bam_calls == [config.healthy_bam]
    dedup_bam = config.output_dir / f"{config.sample_name}.dedup.bam"
    assert dedup_bam.exists()


def test_run_passes_precomputed_tumor_total_reads_into_verify_downsample_depth(
    config: AdmixtureConfig, monkeypatch: pytest.MonkeyPatch
) -> None:
    captured_reads_before: list[int | None] = []
    builder = AdmixtureBuilder(config)
    original_verify = builder.verify_downsample_depth

    def capturing_verify(source_bam, subsampled_bam, expected_fraction, tolerance=0.25, reads_before=None):
        captured_reads_before.append(reads_before)
        return original_verify(source_bam, subsampled_bam, expected_fraction, tolerance, reads_before)

    builder.verify_downsample_depth = capturing_verify

    # 200 matches the fixture's real tumor-source read count (n_reads=200), so
    # the fraction check downstream still passes comfortably -- this test is
    # only about what gets threaded into verify_downsample_depth's
    # reads_before, not about exercising the tolerance check itself.
    builder.run(precomputed_tumor_total_reads=200)

    # First call is the tumor side (gets the precomputed total), second is
    # the healthy side (still None -- unaffected by the tumor-side wiring).
    assert captured_reads_before == [200, None]


def test_run_with_zero_tumor_fraction_builds_mixture_from_healthy_downsample_alone(
    config: AdmixtureConfig,
) -> None:
    """A tumor_fraction == 0 row is a legitimate pure-background negative
    control (MAT-13's fraction_ladder f=0 rows, needed by
    tf_eval.compute_detection_floor's null pool) -- it must not crash
    subsample_bam's `0 < fraction <= 1` check, and its dedup BAM should be
    built from the healthy downsample alone."""
    from pandas import read_csv

    zero_config = replace(config, tumor_fraction=0.0)
    builder = AdmixtureBuilder(zero_config)

    subsample_bam_calls = []
    original_subsample_bam = builder.subsample_bam

    def counting_subsample_bam(source_bam, fraction, output_bam):
        subsample_bam_calls.append(source_bam)
        return original_subsample_bam(source_bam, fraction, output_bam)

    builder.subsample_bam = counting_subsample_bam

    empirical_tumor_fractions = builder.run()

    # Only the healthy side was ever downsampled -- no tumor_sub.bam at all.
    assert subsample_bam_calls == [zero_config.healthy_bam]
    assert not (zero_config.output_dir / f"{zero_config.sample_name}.tumor_sub.bam").exists()

    dedup_bam = zero_config.output_dir / f"{zero_config.sample_name}.dedup.bam"
    assert dedup_bam.exists()

    # No tumor reads at all were merged in -- the empirical estimator comes
    # out at (essentially) zero on its own, with no special-casing needed
    # downstream of merge_and_dedup.
    assert empirical_tumor_fractions["sites"] == pytest.approx(0.0, abs=1e-6)

    metadata = read_csv(zero_config.output_dir / "sample_metadata.csv")
    assert metadata["nominal_tumor_fraction"].iloc[0] == pytest.approx(0.0)
    assert metadata["empirical_tumor_fraction"].iloc[0] == pytest.approx(0.0, abs=1e-6)


def test_run_ignores_a_precomputed_tumor_sub_when_tumor_fraction_is_zero(
    config: AdmixtureConfig,
) -> None:
    """Even when a caller (the Snakemake DAG) passes a precomputed tumor_sub
    and total-reads, a tumor_fraction == 0 row must still skip using them --
    there is no tumor side to this mixture at all, so those arguments must be
    ignored rather than fed into verify_downsample_depth or merge_and_dedup."""
    zero_config = replace(config, tumor_fraction=0.0)
    builder = AdmixtureBuilder(zero_config)

    empirical_tumor_fractions = builder.run(
        precomputed_tumor_sub=Path("/nonexistent/should-not-be-touched.bam"),
        precomputed_tumor_total_reads=999,
    )

    assert "sites" in empirical_tumor_fractions
    dedup_bam = zero_config.output_dir / f"{zero_config.sample_name}.dedup.bam"
    assert dedup_bam.exists()
