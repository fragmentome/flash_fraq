from pathlib import Path
from unittest.mock import patch

from flash_fraq.hcc1395_manifest import HCC1395ManifestConfig, HCC1395RunResolver, SraRun

RUNINFO_CSV = (
    "Run,SampleName,CenterName,Platform,LibraryStrategy\n"
    "SRR001,HCC1395,CENTER_A,ILLUMINA,WGS\n"
    "SRR002,HCC1395BL,CENTER_A,ILLUMINA,WGS\n"
    "SRR003,HCC1395,CENTER_B,ION_TORRENT,WGS\n"
    "SRR004,HCC1395_RNA,CENTER_A,ILLUMINA,RNA-Seq\n"
)


def test_parse_runinfo_builds_sra_run_rows() -> None:
    runs = HCC1395RunResolver.parse_runinfo(RUNINFO_CSV)

    assert len(runs) == 4
    assert runs[0] == SraRun(
        run_accession="SRR001",
        sample_name="HCC1395",
        center_name="CENTER_A",
        platform="ILLUMINA",
        library_strategy="WGS",
    )


def test_resolve_wgs_illumina_runs_filters_non_illumina_wgs(tmp_path: Path) -> None:
    resolver = HCC1395RunResolver(HCC1395ManifestConfig.from_sources(output_path=tmp_path / "runs.csv"))

    with patch.object(resolver, "fetch_runinfo_csv", return_value=RUNINFO_CSV):
        runs = resolver.resolve_wgs_illumina_runs()

    assert [r.run_accession for r in runs] == ["SRR001", "SRR002"]


def test_source_group_and_tumor_status_split_tumor_and_normal() -> None:
    tumor = SraRun(
        run_accession="SRR001", sample_name="HCC1395", center_name="C", platform="ILLUMINA", library_strategy="WGS"
    )
    normal = SraRun(
        run_accession="SRR002", sample_name="HCC1395BL", center_name="C", platform="ILLUMINA", library_strategy="WGS"
    )

    assert (tumor.source_group, tumor.tumor_status) == ("hcc1395", "tumor")
    assert (normal.source_group, normal.tumor_status) == ("hcc1395bl", "healthy")


def test_write_run_list_writes_csv(tmp_path: Path) -> None:
    output_path = tmp_path / "manifest" / "hcc1395_runs.csv"
    resolver = HCC1395RunResolver(HCC1395ManifestConfig.from_sources(output_path=output_path))

    with patch.object(resolver, "fetch_runinfo_csv", return_value=RUNINFO_CSV):
        result_path = resolver.write_run_list()

    assert result_path == output_path
    content = output_path.read_text()
    assert "SRR001" in content
    assert "SRR002" in content
    assert "SRR003" not in content
    assert "SRR004" not in content


def test_to_manifest_entries_sets_fastq_paths_under_fastq_root(tmp_path: Path) -> None:
    resolver = HCC1395RunResolver(HCC1395ManifestConfig.from_sources(output_path=tmp_path / "runs.csv"))
    run = SraRun(
        run_accession="SRR001", sample_name="HCC1395", center_name="C", platform="ILLUMINA", library_strategy="WGS"
    )

    entries = resolver.to_manifest_entries([run], fastq_root=tmp_path / "fastq")

    assert len(entries) == 1
    entry = entries[0]
    assert entry.sample_id == "SRR001"
    assert entry.source_group == "hcc1395"
    assert entry.tumor_status == "tumor"
    assert entry.fastq_r1_path == tmp_path / "fastq" / "SRR001" / "SRR001_1.fastq.gz"
    assert entry.fastq_r2_path == tmp_path / "fastq" / "SRR001" / "SRR001_2.fastq.gz"
