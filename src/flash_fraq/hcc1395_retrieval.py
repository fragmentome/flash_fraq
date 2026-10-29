"""Download HCC1395/HCC1395BL FASTQ pairs, preferring ENA's submitted
`fastq.gz` and falling back to SRA prefetch/fasterq-dump
(docs/phase1-sourcing.md; MAT-23).

Two retrieval paths, chosen per run rather than per cohort:

* **ENA submitted FASTQ** -- when ENA serves the run's originally submitted
  files as a `.fastq.gz` pair, they are downloaded straight to their final
  location. No `.sra` staging copy (which doubles this dataset's isilon
  footprint while it exists) and no decompress/recompress round trip through
  `fasterq-dump` + gzip.
* **SRA toolkit** -- the fallback for every run ENA doesn't serve that way
  (SEQC2 submitted plenty of runs as BAM). Split into two subprocess steps --
  prefetch (downloads the .sra) then fasterq-dump (extracts split _1/_2
  FASTQ) -- mirroring the SRA toolkit's own recommended two-step flow rather
  than fasterq-dump's implicit prefetch, so a failed/interrupted prefetch can
  be retried without re-running the slower extraction step.

Availability is checked per run against ENA's filereport API, not assumed for
the cohort, and which path a run actually took is recorded next to its FASTQ
pair (see `RetrievalRecord`) so a mixed cohort stays auditable.
"""

from __future__ import annotations

import hashlib
import json
import shutil
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from subprocess import CalledProcessError, run

import requests

ENA_FILEREPORT_URL: str = "https://www.ebi.ac.uk/ena/portal/api/filereport"

RETRIEVAL_SOURCE_ENA: str = "ena_submitted_fastq"
RETRIEVAL_SOURCE_SRA: str = "sra_toolkit"


@dataclass(frozen=True)
class EnaSubmittedFastqPair:
    """One run's submitted R1/R2 as served by ENA, with the checksums ENA
    publishes alongside them."""

    r1_url: str
    r2_url: str
    r1_md5: str | None = None
    r2_md5: str | None = None

    @property
    def urls(self) -> tuple[str, str]:
        return self.r1_url, self.r2_url

    @property
    def md5s(self) -> tuple[str | None, str | None]:
        return self.r1_md5, self.r2_md5


@dataclass(frozen=True)
class RetrievalRecord:
    """Which of the two paths a run was actually retrieved by, plus where its
    pair landed. Written to `<run>_retrieval.json` beside the FASTQ pair."""

    run_accession: str
    source: str
    r1_path: Path
    r2_path: Path
    source_urls: tuple[str, ...] = ()

    def to_dict(self) -> dict:
        return {
            "run_accession": self.run_accession,
            "source": self.source,
            "r1_path": str(self.r1_path),
            "r2_path": str(self.r2_path),
            "source_urls": list(self.source_urls),
        }


@dataclass(frozen=True)
class HCC1395RetrievalConfig:
    output_dir: Path
    threads: int = 4
    # Off only for deliberately exercising the SRA fallback on a run ENA does
    # serve; the ENA check is per run, so leaving this on costs one small API
    # call for runs that end up on the toolkit path anyway.
    prefer_ena: bool = True

    @classmethod
    def from_sources(
        cls, output_dir: Path, threads: int = 4, prefer_ena: bool = True
    ) -> "HCC1395RetrievalConfig":
        return cls(output_dir=output_dir, threads=threads, prefer_ena=prefer_ena)


class HCC1395FastqRetriever:
    """Retrieves one SRA run's split FASTQ pair, via ENA's submitted
    `fastq.gz` where available and prefetch + fasterq-dump otherwise."""

    def __init__(self, config: HCC1395RetrievalConfig) -> None:
        self.config = config

    # ---------------------------------------------------------------- ENA

    def resolve_ena_submitted_fastq(self, run_accession: str) -> EnaSubmittedFastqPair | None:
        """Per-run ENA availability check. Returns None -- rather than raising
        -- for every 'not served this way' outcome (submitted as BAM, single
        file, ENA unreachable), because each of those just means this run
        falls back to the SRA toolkit."""
        try:
            response = requests.get(
                ENA_FILEREPORT_URL,
                params={
                    "accession": run_accession,
                    "result": "read_run",
                    "fields": "submitted_ftp,submitted_md5,submitted_format",
                    "format": "tsv",
                },
                timeout=60,
            )
            response.raise_for_status()
        except requests.RequestException:
            return None
        return self.parse_ena_filereport(response.text)

    @staticmethod
    def parse_ena_filereport(filereport_tsv: str) -> EnaSubmittedFastqPair | None:
        lines = [line for line in filereport_tsv.splitlines() if line.strip()]
        if len(lines) < 2:
            return None
        header = lines[0].split("\t")
        if "submitted_ftp" not in header:
            return None
        row = dict(zip(header, lines[1].split("\t")))

        paths = [path for path in row.get("submitted_ftp", "").split(";") if path]
        md5s = [md5 for md5 in row.get("submitted_md5", "").split(";") if md5]
        # Only a submitted pair is usable as-is: a single file is either
        # interleaved or unpaired, and anything else (BAM, CRAM, bz2) would
        # need the same conversion round trip this path exists to avoid.
        if len(paths) != 2 or not all(path.endswith(".fastq.gz") for path in paths):
            return None
        if len(md5s) != 2:
            md5s = [None, None]

        # ENA reports these as bare host/path; it serves the same tree over
        # https, which needs no anonymous-FTP handling in the fetch below.
        urls = [path if "://" in path else f"https://{path}" for path in paths]
        return EnaSubmittedFastqPair(r1_url=urls[0], r2_url=urls[1], r1_md5=md5s[0], r2_md5=md5s[1])

    def _download_ena_pair(
        self, run_accession: str, sample_dir: Path, pair: EnaSubmittedFastqPair
    ) -> tuple[Path, Path]:
        destinations = []
        for mate, (url, expected_md5) in enumerate(zip(pair.urls, pair.md5s), start=1):
            destination = sample_dir / f"{run_accession}_{mate}.fastq.gz"
            # Download to a partial name so an interrupted transfer can never
            # be mistaken for a complete pair by the next run (or by
            # Snakemake, which only sees whether the output exists).
            partial = destination.with_suffix(destination.suffix + ".partial")
            run(
                ["curl", "--location", "--fail", "--silent", "--show-error", "--output", str(partial), url],
                check=True,
            )
            if expected_md5:
                observed_md5 = self._file_md5(partial)
                if observed_md5 != expected_md5:
                    partial.unlink(missing_ok=True)
                    raise RuntimeError(
                        f"{run_accession}: ENA md5 mismatch for mate {mate} ({url}): "
                        f"expected {expected_md5}, got {observed_md5}"
                    )
            partial.rename(destination)
            destinations.append(destination)
        return destinations[0], destinations[1]

    @staticmethod
    def _file_md5(path: Path) -> str:
        digest = hashlib.md5()
        with open(path, "rb") as handle:
            # These files run to tens of GB, so hash in chunks rather than
            # reading one into memory on a shared node.
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        return digest.hexdigest()

    # ---------------------------------------------------------- SRA toolkit

    def _download_via_sra_toolkit(self, run_accession: str, sample_dir: Path) -> tuple[Path, Path]:
        cfg = self.config

        # SEQC2/HCC1395 runs go up to ~175G -- 100G was clipping real runs.
        run(["prefetch", run_accession, "-O", str(cfg.output_dir), "--max-size", "200G"], check=True)

        # prefetch exits 0 even when it skips a run (e.g. over --max-size), so
        # check for the .sra here rather than letting it surface as a confusing
        # gzip failure two steps later.
        sra_path = cfg.output_dir / run_accession / f"{run_accession}.sra"
        if not sra_path.exists():
            raise RuntimeError(f"{run_accession}: prefetch did not produce {sra_path}")

        try:
            run(
                [
                    "fasterq-dump",
                    "--split-files",
                    "--threads", str(cfg.threads),
                    "-O", str(sample_dir),
                    # Without --temp, fasterq-dump scratches to the current working
                    # directory -- the repo checkout under $HOME on this workflow --
                    # and a large run can blow the home quota. Scratch under the
                    # isilon sample dir instead, matching this project's convention
                    # of never writing job output/scratch under REMOTE_REPO_DIR/~.
                    "--temp", str(sample_dir),
                    str(sra_path),
                ],
                check=True,
            )
        except CalledProcessError:
            # A failed/interrupted run leaves its fasterq.tmp.<host>.<pid>
            # scratch dir behind uncleaned -- 58 samples' worth of these
            # (~624G) accumulated from repeated retries and helped exhaust
            # the isilon quota. Clean up before re-raising so a retry
            # doesn't compound the waste.
            for tmp_dir in sample_dir.glob("fasterq.tmp.*"):
                shutil.rmtree(tmp_dir, ignore_errors=True)
            raise

        # fasterq-dump names its output after the basename it was handed, so
        # passing <run>.sra yields <run>.sra_1.fastq -- not the <run>_1.fastq
        # this method returns and gzips. Every run therefore died at `gzip: No
        # such file or directory` after a successful multi-hour extraction
        # (26 files in the .sra_1.fastq form and 0 in the other were found on
        # isilon), which an earlier pass misread as prefetch skipping the run.
        # Resolve whatever pair fasterq-dump actually wrote and normalize the
        # names, rather than depending on how it derives them.
        extracted_mates = []
        for mate in (1, 2):
            extracted = sorted(sample_dir.glob(f"*_{mate}.fastq"))
            if len(extracted) != 1:
                raise RuntimeError(
                    f"{run_accession}: expected exactly one *_{mate}.fastq from fasterq-dump "
                    f"in {sample_dir}, found {[path.name for path in extracted]}"
                )
            extracted_mates.append(extracted[0])

        # Both mates are on disk, so the .sra has served its only purpose and
        # nothing below reads it again. Dropping it here rather than at end of
        # run turns over the sra_fetch disk cap a compression pass sooner,
        # instead of leaving 58 samples' worth of staging copies behind.
        sra_path.unlink()

        r1_path = sample_dir / f"{run_accession}_1.fastq.gz"
        r2_path = sample_dir / f"{run_accession}_2.fastq.gz"
        for mate, extracted in enumerate(extracted_mates, start=1):
            canonical = sample_dir / f"{run_accession}_{mate}.fastq"
            extracted.rename(canonical)
            run(self._compress_command(canonical), check=True)

        if not r1_path.exists() or not r2_path.exists():
            raise RuntimeError(f"{run_accession}: fasterq-dump did not produce split FASTQ pair in {sample_dir}")

        return r1_path, r2_path

    def _compress_command(self, path: Path) -> list[str]:
        """pigz where it exists, gzip otherwise. Compressing a ~175G run's
        FASTQ single-threaded is hours of an otherwise idle multi-core fetch
        job; both write the same gzip container, so the fallback only costs
        speed."""
        if shutil.which("pigz") is not None:
            return ["pigz", "-f", "-p", str(self.config.threads), str(path)]
        return ["gzip", "-f", str(path)]

    # ------------------------------------------------------------- driver

    def download_run(self, run_accession: str) -> RetrievalRecord:
        sample_dir = self.config.output_dir / run_accession
        sample_dir.mkdir(parents=True, exist_ok=True)

        ena_pair = self.resolve_ena_submitted_fastq(run_accession) if self.config.prefer_ena else None
        if ena_pair is not None:
            r1_path, r2_path = self._download_ena_pair(run_accession, sample_dir, ena_pair)
            source, source_urls = RETRIEVAL_SOURCE_ENA, ena_pair.urls
        else:
            r1_path, r2_path = self._download_via_sra_toolkit(run_accession, sample_dir)
            source, source_urls = RETRIEVAL_SOURCE_SRA, ()

        record = RetrievalRecord(
            run_accession=run_accession,
            source=source,
            r1_path=r1_path,
            r2_path=r2_path,
            source_urls=tuple(source_urls),
        )
        self.write_retrieval_record(record, sample_dir)
        return record

    @staticmethod
    def write_retrieval_record(record: RetrievalRecord, sample_dir: Path) -> Path:
        """Per-run sidecar rather than a column appended to the cohort run
        list: fetch jobs run concurrently (sra_fetch=4), so several of them
        rewriting one shared CSV would race, and the run list is a checkpoint
        output Snakemake owns."""
        record_path = sample_dir / f"{record.run_accession}_retrieval.json"
        payload = record.to_dict() | {"retrieved_at": datetime.now(timezone.utc).isoformat()}
        record_path.write_text(json.dumps(payload, indent=2) + "\n")
        return record_path
