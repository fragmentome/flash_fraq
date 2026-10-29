"""Thin coverage for scripts/run_tumor_batch_sampling.py's argument parsing
and plan-loading -- not a re-test of TumorBatchSampler's sampling algorithm,
which tests/test_tumor_batch_sampling.py already covers."""

import csv
import importlib.util
from pathlib import Path

import pytest

SCRIPT_PATH = Path(__file__).resolve().parent.parent / "scripts" / "run_tumor_batch_sampling.py"
_spec = importlib.util.spec_from_file_location("run_tumor_batch_sampling", SCRIPT_PATH)
run_tumor_batch_sampling = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(run_tumor_batch_sampling)


def _write_plan(plan_path: Path, rows: list[dict[str, str]]) -> None:
    fieldnames = list(rows[0].keys())
    with open(plan_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _row(sample_name: str, seed: int, tumor_downsample_fraction: float, tumor_bam_path: str) -> dict[str, str]:
    return {
        "sample_name": sample_name,
        "seed": str(seed),
        "tumor_downsample_fraction": str(tumor_downsample_fraction),
        "tumor_bam_path": tumor_bam_path,
    }


def test_load_targets_builds_one_target_per_row(tmp_path: Path) -> None:
    plan_path = tmp_path / "plan.csv"
    _write_plan(
        plan_path,
        [
            _row("mix_a", 42, 0.25, "/data/tumor.cram"),
            _row("mix_b", 43, 0.5, "/data/tumor.cram"),
        ],
    )

    tumor_bam, targets = run_tumor_batch_sampling._load_targets(plan_path, tmp_path / "out")

    assert tumor_bam == Path("/data/tumor.cram")
    assert len(targets) == 2
    assert targets[0].output_bam == tmp_path / "out" / "mix_a.tumor_sub.bam"
    assert targets[0].seed == 42
    assert targets[0].fraction == pytest.approx(0.25)
    assert targets[1].output_bam == tmp_path / "out" / "mix_b.tumor_sub.bam"


def test_load_targets_skips_pure_background_rows_with_zero_tumor_fraction(tmp_path: Path) -> None:
    """A tumor_downsample_fraction == 0 row is a legitimate pure-background
    negative control (nominal admixture fraction f=0) -- there is no tumor
    read to keep for it, so it must not become a TumorDownsampleTarget at
    all (TumorDownsampleTarget's own `0 < fraction <= 1` validation would
    reject constructing one anyway)."""
    plan_path = tmp_path / "plan.csv"
    _write_plan(
        plan_path,
        [
            _row("mix_zero", 42, 0.0, "/data/tumor.cram"),
            _row("mix_nonzero", 43, 0.5, "/data/tumor.cram"),
        ],
    )

    tumor_bam, targets = run_tumor_batch_sampling._load_targets(plan_path, tmp_path / "out")

    assert tumor_bam == Path("/data/tumor.cram")
    assert len(targets) == 1
    assert targets[0].output_bam == tmp_path / "out" / "mix_nonzero.tumor_sub.bam"


def test_load_targets_rejects_a_plan_with_more_than_one_tumor_bam_path(tmp_path: Path) -> None:
    plan_path = tmp_path / "plan.csv"
    _write_plan(
        plan_path,
        [
            _row("mix_a", 42, 0.25, "/data/tumor_one.cram"),
            _row("mix_b", 43, 0.5, "/data/tumor_two.cram"),
        ],
    )

    with pytest.raises(ValueError, match="expected exactly one tumor_bam_path"):
        run_tumor_batch_sampling._load_targets(plan_path, tmp_path / "out")


def test_load_targets_rejects_an_empty_plan(tmp_path: Path) -> None:
    plan_path = tmp_path / "plan.csv"
    _write_plan(plan_path, [_row("mix_a", 42, 0.25, "/data/tumor.cram")])
    with open(plan_path, "w", newline="") as handle:
        handle.write("sample_name,seed,tumor_downsample_fraction,tumor_bam_path\n")

    with pytest.raises(ValueError, match="no rows"):
        run_tumor_batch_sampling._load_targets(plan_path, tmp_path / "out")


def test_main_runs_the_sampler_and_writes_the_total_reads_file(tmp_path: Path, monkeypatch) -> None:
    plan_path = tmp_path / "plan.csv"
    _write_plan(
        plan_path,
        [
            _row("mix_a", 42, 0.25, "/data/tumor.cram"),
            _row("mix_b", 43, 0.5, "/data/tumor.cram"),
        ],
    )
    output_dir = tmp_path / "out"
    captured_calls = []

    class _FakeSampler:
        def run(self, tumor_bam, reference_fasta, targets, threads=1):
            captured_calls.append((tumor_bam, reference_fasta, list(targets), threads))
            return 12345

    monkeypatch.setattr(run_tumor_batch_sampling, "TumorBatchSampler", _FakeSampler)
    monkeypatch.setattr(run_tumor_batch_sampling, "_raise_open_file_limit", lambda minimum: None)
    monkeypatch.setattr(
        "sys.argv",
        [
            "run_tumor_batch_sampling.py",
            "--plan", str(plan_path),
            "--reference-fasta", str(tmp_path / "ref.fa"),
            "--output-dir", str(output_dir),
            "--threads", "4",
        ],
    )

    run_tumor_batch_sampling.main()

    assert len(captured_calls) == 1
    tumor_bam, reference_fasta, targets, threads = captured_calls[0]
    assert tumor_bam == Path("/data/tumor.cram")
    assert reference_fasta == tmp_path / "ref.fa"
    assert len(targets) == 2
    assert threads == 4

    total_reads_path = output_dir / run_tumor_batch_sampling.TOTAL_READS_FILENAME
    assert total_reads_path.read_text().strip() == "12345"
