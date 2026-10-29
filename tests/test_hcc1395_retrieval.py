import json
from pathlib import Path
from unittest.mock import patch

import pytest
import requests

from flash_fraq.hcc1395_retrieval import (
    RETRIEVAL_SOURCE_ENA,
    RETRIEVAL_SOURCE_SRA,
    EnaSubmittedFastqPair,
    HCC1395FastqRetriever,
    HCC1395RetrievalConfig,
)

ENA_SUBMITTED_TSV = (
    "run_accession\tsubmitted_ftp\tsubmitted_md5\tsubmitted_format\n"
    "SRR001\t"
    "ftp.sra.ebi.ac.uk/vol1/run/SRR001/SRR001_1.fastq.gz;"
    "ftp.sra.ebi.ac.uk/vol1/run/SRR001/SRR001_2.fastq.gz\t"
    "d41d8cd98f00b204e9800998ecf8427e;d41d8cd98f00b204e9800998ecf8427e\t"
    "FASTQ;FASTQ\n"
)
ENA_SUBMITTED_BAM_TSV = (
    "run_accession\tsubmitted_ftp\tsubmitted_md5\tsubmitted_format\n"
    "SRR001\tftp.sra.ebi.ac.uk/vol1/run/SRR001/SRR001.bam\tabc\tBAM\n"
)

# md5 of the empty bytes the curl stub writes, so the checksum check is
# exercised end to end rather than skipped.
EMPTY_MD5 = "d41d8cd98f00b204e9800998ecf8427e"


def _sra_retriever(tmp_path: Path) -> HCC1395FastqRetriever:
    """Retriever pinned to the SRA-toolkit path, for the tests covering it."""
    return HCC1395FastqRetriever(
        HCC1395RetrievalConfig.from_sources(output_dir=tmp_path, prefer_ena=False)
    )


def _fake_run_factory(tmp_path: Path, extracted_names: tuple[str, ...]):
    """Stubs prefetch/fasterq-dump/gzip. extracted_names is what fasterq-dump
    writes, which varies with the basename it's handed (see download_run)."""

    def fake_run(cmd, **kwargs):
        sample_dir = tmp_path / "SRR001"
        if cmd[0] == "prefetch":
            sample_dir.mkdir(parents=True, exist_ok=True)
            (sample_dir / "SRR001.sra").write_bytes(b"")
        elif cmd[0] == "fasterq-dump":
            for name in extracted_names:
                (sample_dir / name).write_bytes(b"")
        elif cmd[0] in ("gzip", "pigz"):
            fastq_path = Path(cmd[-1])
            fastq_path.rename(fastq_path.with_suffix(fastq_path.suffix + ".gz"))
        return None

    return fake_run


def _fake_curl_factory():
    def fake_run(cmd, **kwargs):
        assert cmd[0] == "curl"
        Path(cmd[cmd.index("--output") + 1]).write_bytes(b"")
        return None

    return fake_run


class _FakeResponse:
    def __init__(self, text: str) -> None:
        self.text = text

    def raise_for_status(self) -> None:
        return None


# ------------------------------------------------------------------ ENA


def test_parse_ena_filereport_returns_submitted_pair() -> None:
    pair = HCC1395FastqRetriever.parse_ena_filereport(ENA_SUBMITTED_TSV)

    assert pair == EnaSubmittedFastqPair(
        r1_url="https://ftp.sra.ebi.ac.uk/vol1/run/SRR001/SRR001_1.fastq.gz",
        r2_url="https://ftp.sra.ebi.ac.uk/vol1/run/SRR001/SRR001_2.fastq.gz",
        r1_md5=EMPTY_MD5,
        r2_md5=EMPTY_MD5,
    )


@pytest.mark.parametrize(
    "filereport_tsv",
    [
        # Submitted as BAM -- reading it would need the very conversion round
        # trip this path exists to avoid, so the run belongs on the toolkit.
        ENA_SUBMITTED_BAM_TSV,
        # Header only: ENA knows the run but serves no submitted files.
        "run_accession\tsubmitted_ftp\tsubmitted_md5\tsubmitted_format\n",
        "",
    ],
)
def test_parse_ena_filereport_returns_none_when_no_submitted_pair(filereport_tsv: str) -> None:
    assert HCC1395FastqRetriever.parse_ena_filereport(filereport_tsv) is None


def test_resolve_ena_submitted_fastq_returns_none_when_ena_unreachable(tmp_path: Path) -> None:
    # ENA being down must not fail the fetch -- it just means this run takes
    # the SRA-toolkit path.
    retriever = HCC1395FastqRetriever(HCC1395RetrievalConfig.from_sources(output_dir=tmp_path))

    with patch(
        "flash_fraq.hcc1395_retrieval.requests.get", side_effect=requests.ConnectionError("down")
    ):
        assert retriever.resolve_ena_submitted_fastq("SRR001") is None


def test_download_run_prefers_ena_and_skips_the_sra_toolkit(tmp_path: Path) -> None:
    retriever = HCC1395FastqRetriever(HCC1395RetrievalConfig.from_sources(output_dir=tmp_path))

    with patch(
        "flash_fraq.hcc1395_retrieval.requests.get", return_value=_FakeResponse(ENA_SUBMITTED_TSV)
    ):
        with patch("flash_fraq.hcc1395_retrieval.run", side_effect=_fake_curl_factory()) as mock_run:
            record = retriever.download_run("SRR001")

    assert record.source == RETRIEVAL_SOURCE_ENA
    assert record.r1_path == tmp_path / "SRR001" / "SRR001_1.fastq.gz"
    assert record.r2_path == tmp_path / "SRR001" / "SRR001_2.fastq.gz"
    assert record.r1_path.exists()
    assert record.r2_path.exists()
    # Two curls and nothing else: no prefetch, no .sra staging copy, no
    # decompress/recompress round trip.
    assert all(call.args[0][0] == "curl" for call in mock_run.call_args_list)
    assert mock_run.call_count == 2
    assert not list((tmp_path / "SRR001").glob("*.partial"))


def test_download_run_rejects_ena_pair_failing_its_published_md5(tmp_path: Path) -> None:
    # A truncated transfer must not be left behind looking like a finished
    # output, since Snakemake only checks that the file exists.
    mismatched_tsv = ENA_SUBMITTED_TSV.replace(EMPTY_MD5, "0" * 32)
    retriever = HCC1395FastqRetriever(HCC1395RetrievalConfig.from_sources(output_dir=tmp_path))

    with patch(
        "flash_fraq.hcc1395_retrieval.requests.get", return_value=_FakeResponse(mismatched_tsv)
    ):
        with patch("flash_fraq.hcc1395_retrieval.run", side_effect=_fake_curl_factory()):
            with pytest.raises(RuntimeError, match="md5 mismatch"):
                retriever.download_run("SRR001")

    assert not (tmp_path / "SRR001" / "SRR001_1.fastq.gz").exists()
    assert not list((tmp_path / "SRR001").glob("*.partial"))


def test_download_run_falls_back_to_sra_toolkit_when_ena_has_no_fastq(tmp_path: Path) -> None:
    # Availability is per run, not per cohort: a BAM-submitted run in an
    # otherwise ENA-served study still has to work.
    retriever = HCC1395FastqRetriever(HCC1395RetrievalConfig.from_sources(output_dir=tmp_path))
    fake_run = _fake_run_factory(tmp_path, ("SRR001_1.fastq", "SRR001_2.fastq"))

    with patch(
        "flash_fraq.hcc1395_retrieval.requests.get",
        return_value=_FakeResponse(ENA_SUBMITTED_BAM_TSV),
    ):
        with patch("flash_fraq.hcc1395_retrieval.run", side_effect=fake_run) as mock_run:
            record = retriever.download_run("SRR001")

    assert record.source == RETRIEVAL_SOURCE_SRA
    assert record.r1_path.exists()
    assert record.r2_path.exists()
    assert mock_run.call_args_list[0].args[0][0] == "prefetch"


def test_download_run_writes_a_retrieval_record_naming_the_path_used(tmp_path: Path) -> None:
    retriever = HCC1395FastqRetriever(HCC1395RetrievalConfig.from_sources(output_dir=tmp_path))

    with patch(
        "flash_fraq.hcc1395_retrieval.requests.get", return_value=_FakeResponse(ENA_SUBMITTED_TSV)
    ):
        with patch("flash_fraq.hcc1395_retrieval.run", side_effect=_fake_curl_factory()):
            retriever.download_run("SRR001")

    record = json.loads((tmp_path / "SRR001" / "SRR001_retrieval.json").read_text())
    assert record["run_accession"] == "SRR001"
    assert record["source"] == RETRIEVAL_SOURCE_ENA
    assert record["source_urls"] == [
        "https://ftp.sra.ebi.ac.uk/vol1/run/SRR001/SRR001_1.fastq.gz",
        "https://ftp.sra.ebi.ac.uk/vol1/run/SRR001/SRR001_2.fastq.gz",
    ]
    assert record["r1_path"] == str(tmp_path / "SRR001" / "SRR001_1.fastq.gz")


# ---------------------------------------------------------------- SRA toolkit


def test_download_run_calls_prefetch_then_fasterq_dump_then_compress(tmp_path: Path) -> None:
    retriever = _sra_retriever(tmp_path)
    fake_run = _fake_run_factory(tmp_path, ("SRR001_1.fastq", "SRR001_2.fastq"))

    with patch("flash_fraq.hcc1395_retrieval.run", side_effect=fake_run) as mock_run:
        record = retriever.download_run("SRR001")

    assert record.source == RETRIEVAL_SOURCE_SRA
    assert record.r1_path == tmp_path / "SRR001" / "SRR001_1.fastq.gz"
    assert record.r2_path == tmp_path / "SRR001" / "SRR001_2.fastq.gz"
    assert record.r1_path.exists()
    assert record.r2_path.exists()
    assert mock_run.call_count == 4  # prefetch, fasterq-dump, compress x2
    # The .sra is deleted once the extraction that needed it succeeds, so one
    # dataset doesn't carry two equally large copies of the same reads while
    # the pair is compressed.
    assert not (tmp_path / "SRR001" / "SRR001.sra").exists()


def test_download_run_compresses_with_pigz_when_available(tmp_path: Path) -> None:
    # Single-threaded gzip over a ~175G run is hours of an otherwise idle
    # multi-core fetch job.
    retriever = HCC1395FastqRetriever(
        HCC1395RetrievalConfig.from_sources(output_dir=tmp_path, threads=12, prefer_ena=False)
    )
    fake_run = _fake_run_factory(tmp_path, ("SRR001_1.fastq", "SRR001_2.fastq"))

    with patch("flash_fraq.hcc1395_retrieval.shutil.which", return_value="/usr/bin/pigz"):
        with patch("flash_fraq.hcc1395_retrieval.run", side_effect=fake_run) as mock_run:
            retriever.download_run("SRR001")

    compress_cmds = [c.args[0] for c in mock_run.call_args_list if c.args[0][0] in ("pigz", "gzip")]
    assert len(compress_cmds) == 2
    assert all(cmd[:4] == ["pigz", "-f", "-p", "12"] for cmd in compress_cmds)


def test_download_run_falls_back_to_gzip_without_pigz(tmp_path: Path) -> None:
    retriever = _sra_retriever(tmp_path)
    fake_run = _fake_run_factory(tmp_path, ("SRR001_1.fastq", "SRR001_2.fastq"))

    with patch("flash_fraq.hcc1395_retrieval.shutil.which", return_value=None):
        with patch("flash_fraq.hcc1395_retrieval.run", side_effect=fake_run) as mock_run:
            retriever.download_run("SRR001")

    compress_cmds = [c.args[0] for c in mock_run.call_args_list if c.args[0][0] in ("pigz", "gzip")]
    assert compress_cmds == [
        ["gzip", "-f", str(tmp_path / "SRR001" / "SRR001_1.fastq")],
        ["gzip", "-f", str(tmp_path / "SRR001" / "SRR001_2.fastq")],
    ]


def test_download_run_passes_threads_to_fasterq_dump(tmp_path: Path) -> None:
    retriever = HCC1395FastqRetriever(
        HCC1395RetrievalConfig.from_sources(output_dir=tmp_path, threads=12, prefer_ena=False)
    )
    fake_run = _fake_run_factory(tmp_path, ("SRR001_1.fastq", "SRR001_2.fastq"))

    with patch("flash_fraq.hcc1395_retrieval.run", side_effect=fake_run) as mock_run:
        retriever.download_run("SRR001")

    fasterq_cmd = next(c.args[0] for c in mock_run.call_args_list if c.args[0][0] == "fasterq-dump")
    assert fasterq_cmd[fasterq_cmd.index("--threads") + 1] == "12"


def test_download_run_normalizes_sra_suffixed_fasterq_dump_output(tmp_path: Path) -> None:
    # fasterq-dump names output after the basename it's handed, so being
    # passed <run>.sra makes it write <run>.sra_1.fastq. Every run used to die
    # at gzip on that mismatch after a full extraction; the pair it actually
    # wrote is resolved and renamed instead.
    retriever = _sra_retriever(tmp_path)
    fake_run = _fake_run_factory(tmp_path, ("SRR001.sra_1.fastq", "SRR001.sra_2.fastq"))

    with patch("flash_fraq.hcc1395_retrieval.run", side_effect=fake_run):
        record = retriever.download_run("SRR001")

    assert record.r1_path == tmp_path / "SRR001" / "SRR001_1.fastq.gz"
    assert record.r2_path == tmp_path / "SRR001" / "SRR001_2.fastq.gz"
    assert record.r1_path.exists()
    assert record.r2_path.exists()


def test_download_run_raises_if_mate_cannot_be_resolved(tmp_path: Path) -> None:
    # A missing or ambiguous mate must name the directory and what was found,
    # rather than surfacing two steps later as a bare gzip failure -- and the
    # .sra must survive it, since a retry would otherwise re-prefetch.
    retriever = _sra_retriever(tmp_path)
    fake_run = _fake_run_factory(tmp_path, ("SRR001.sra_1.fastq",))

    with patch("flash_fraq.hcc1395_retrieval.run", side_effect=fake_run):
        with pytest.raises(RuntimeError, match=r"exactly one \*_2\.fastq"):
            retriever.download_run("SRR001")

    assert (tmp_path / "SRR001" / "SRR001.sra").exists()


def test_download_run_raises_if_fastq_pair_missing(tmp_path: Path) -> None:
    retriever = _sra_retriever(tmp_path)

    with patch("flash_fraq.hcc1395_retrieval.run", return_value=None):
        with pytest.raises(RuntimeError):
            retriever.download_run("SRR001")
