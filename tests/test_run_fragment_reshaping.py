"""Thin coverage for scripts/run_fragment_reshaping.py's argument parsing --
not a re-test of FragmentReshaper's pipeline, which
tests/test_fragment_reshaping.py already covers end to end."""

import importlib.util
from pathlib import Path

import pytest

SCRIPT_PATH = Path(__file__).resolve().parent.parent / "scripts" / "run_fragment_reshaping.py"
_spec = importlib.util.spec_from_file_location("run_fragment_reshaping", SCRIPT_PATH)
run_fragment_reshaping = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(run_fragment_reshaping)


def test_main_builds_config_with_target_bam(tmp_path: Path, monkeypatch) -> None:
    captured_configs = []

    class _FakeReshaper:
        def __init__(self, config):
            captured_configs.append(config)

        def run(self):
            class _Verification:
                def to_dict(self):
                    return {"statistic": 0.01, "threshold": 0.15, "passed": True}

            return {
                "output_bam": tmp_path / "out" / "reshaped.reshaped.bam",
                "n_source_fragments": 100,
                "n_reshaped_fragments": 40,
                "verification": _Verification(),
            }

    monkeypatch.setattr(run_fragment_reshaping, "FragmentReshaper", _FakeReshaper)
    monkeypatch.setattr(
        "sys.argv",
        [
            "run_fragment_reshaping.py",
            "--source-bam", str(tmp_path / "hcc1395.cram"),
            "--reference-fasta", str(tmp_path / "ref.fa"),
            "--output-dir", str(tmp_path / "out"),
            "--target-bam", str(tmp_path / "background.cram"),
            "--sample-name", "reshaped",
        ],
    )

    run_fragment_reshaping.main()

    assert len(captured_configs) == 1
    config = captured_configs[0]
    assert config.source_bam == tmp_path / "hcc1395.cram"
    assert config.target_bam == tmp_path / "background.cram"
    assert config.target_histogram_path is None
    assert config.sample_name == "reshaped"
    assert config.min_mapq == 30  # DEFAULT_MIN_MAPQ, not passed explicitly above


def test_min_mapq_flag_is_threaded_into_the_config(tmp_path: Path, monkeypatch) -> None:
    captured_configs = []

    class _FakeReshaper:
        def __init__(self, config):
            captured_configs.append(config)

        def run(self):
            class _Verification:
                def to_dict(self):
                    return {"statistic": 0.01, "threshold": 0.15, "passed": True}

            return {
                "output_bam": tmp_path / "out" / "reshaped.reshaped.bam",
                "n_source_fragments": 100,
                "n_reshaped_fragments": 40,
                "verification": _Verification(),
            }

    monkeypatch.setattr(run_fragment_reshaping, "FragmentReshaper", _FakeReshaper)
    monkeypatch.setattr(
        "sys.argv",
        [
            "run_fragment_reshaping.py",
            "--source-bam", str(tmp_path / "hcc1395.cram"),
            "--reference-fasta", str(tmp_path / "ref.fa"),
            "--output-dir", str(tmp_path / "out"),
            "--target-bam", str(tmp_path / "background.cram"),
            "--min-mapq", "40",
        ],
    )

    run_fragment_reshaping.main()

    assert captured_configs[0].min_mapq == 40


def test_main_builds_config_with_target_histogram_path(tmp_path: Path, monkeypatch) -> None:
    captured_configs = []

    class _FakeReshaper:
        def __init__(self, config):
            captured_configs.append(config)

        def run(self):
            class _Verification:
                def to_dict(self):
                    return {"statistic": 0.01, "threshold": 0.15, "passed": True}

            return {
                "output_bam": tmp_path / "out" / "reshaped.reshaped.bam",
                "n_source_fragments": 100,
                "n_reshaped_fragments": 40,
                "verification": _Verification(),
            }

    monkeypatch.setattr(run_fragment_reshaping, "FragmentReshaper", _FakeReshaper)
    monkeypatch.setattr(
        "sys.argv",
        [
            "run_fragment_reshaping.py",
            "--source-bam", str(tmp_path / "hcc1395.cram"),
            "--reference-fasta", str(tmp_path / "ref.fa"),
            "--output-dir", str(tmp_path / "out"),
            "--target-histogram-path", str(tmp_path / "target_hist.csv"),
        ],
    )

    run_fragment_reshaping.main()

    config = captured_configs[0]
    assert config.target_bam is None
    assert config.target_histogram_path == tmp_path / "target_hist.csv"


def test_target_bam_and_target_histogram_path_are_mutually_exclusive(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(
        "sys.argv",
        [
            "run_fragment_reshaping.py",
            "--source-bam", str(tmp_path / "hcc1395.cram"),
            "--reference-fasta", str(tmp_path / "ref.fa"),
            "--output-dir", str(tmp_path / "out"),
            "--target-bam", str(tmp_path / "background.cram"),
            "--target-histogram-path", str(tmp_path / "target_hist.csv"),
        ],
    )

    with pytest.raises(SystemExit):
        run_fragment_reshaping.main()
