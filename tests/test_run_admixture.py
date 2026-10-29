"""Thin coverage for scripts/run_admixture.py's argument parsing -- not a
re-test of AdmixtureBuilder's pipeline, which tests/test_admixture.py already
covers end to end."""

import importlib.util
from pathlib import Path

import pytest

SCRIPT_PATH = Path(__file__).resolve().parent.parent / "scripts" / "run_admixture.py"
_spec = importlib.util.spec_from_file_location("run_admixture", SCRIPT_PATH)
run_admixture = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(run_admixture)


def test_parse_snv_vaf_path_splits_label_and_path() -> None:
    label, path = run_admixture._parse_snv_vaf_path("sites=/mnt/data/sites.vcf.gz")
    assert label == "sites"
    assert path == Path("/mnt/data/sites.vcf.gz")


def test_parse_snv_vaf_path_rejects_missing_separator() -> None:
    with pytest.raises(ValueError, match="LABEL=PATH"):
        run_admixture._parse_snv_vaf_path("sites")


def test_build_parser_accepts_one_plan_row_worth_of_arguments(tmp_path: Path) -> None:
    parser = run_admixture.ArgumentParser(description=run_admixture.__doc__)
    parser.add_argument("--tumor-bam", required=True, type=Path)
    args = parser.parse_args(["--tumor-bam", str(tmp_path / "t.bam")])
    assert args.tumor_bam == tmp_path / "t.bam"


def test_main_builds_config_from_a_plan_rows_arguments(tmp_path: Path, monkeypatch) -> None:
    captured_configs = []

    class _FakeBuilder:
        def __init__(self, config):
            captured_configs.append(config)

        def run(self, precomputed_tumor_sub=None, precomputed_tumor_total_reads=None):
            return {"sites": 0.1, "clonal": 0.11}

    monkeypatch.setattr(run_admixture, "AdmixtureBuilder", _FakeBuilder)
    monkeypatch.setattr(
        "sys.argv",
        [
            "run_admixture.py",
            "--tumor-bam", str(tmp_path / "tumor.bam"),
            "--healthy-bam", str(tmp_path / "healthy.bam"),
            "--reference-fasta", str(tmp_path / "ref.fa"),
            "--snv-vaf-path", f"sites={tmp_path / 'sites.vcf.gz'}",
            "--snv-vaf-path", f"clonal={tmp_path / 'clonal.vcf.gz'}",
            "--output-dir", str(tmp_path / "out"),
            "--tumor-fraction", "0.1",
            "--target-coverage", "1.0",
            "--tumor-native-coverage", "35",
            "--healthy-native-coverage", "20",
            "--sample-name", "mix_test",
        ],
    )

    run_admixture.main()

    assert len(captured_configs) == 1
    config = captured_configs[0]
    assert config.sample_name == "mix_test"
    assert config.tumor_fraction == pytest.approx(0.1)
    assert set(config.snv_vaf_paths) == {"sites", "clonal"}
    assert config.snv_vaf_paths["sites"] == tmp_path / "sites.vcf.gz"
    assert config.threads == 8  # --threads not passed above -- CLI default


def test_main_threads_the_threads_argument_into_the_config(tmp_path: Path, monkeypatch) -> None:
    captured_configs = []

    class _FakeBuilder:
        def __init__(self, config):
            captured_configs.append(config)

        def run(self, precomputed_tumor_sub=None, precomputed_tumor_total_reads=None):
            return {"sites": 0.1}

    monkeypatch.setattr(run_admixture, "AdmixtureBuilder", _FakeBuilder)
    monkeypatch.setattr(
        "sys.argv",
        [
            "run_admixture.py",
            "--tumor-bam", str(tmp_path / "tumor.bam"),
            "--healthy-bam", str(tmp_path / "healthy.bam"),
            "--reference-fasta", str(tmp_path / "ref.fa"),
            "--snv-vaf-path", f"sites={tmp_path / 'sites.vcf.gz'}",
            "--output-dir", str(tmp_path / "out"),
            "--tumor-fraction", "0.1",
            "--target-coverage", "1.0",
            "--tumor-native-coverage", "35",
            "--healthy-native-coverage", "20",
            "--sample-name", "mix_test",
            "--threads", "16",
        ],
    )

    run_admixture.main()

    assert captured_configs[0].threads == 16


def test_main_passes_precomputed_tumor_sub_and_total_reads_through_to_run(tmp_path: Path, monkeypatch) -> None:
    captured_run_kwargs = []

    class _FakeBuilder:
        def __init__(self, config):
            pass

        def run(self, precomputed_tumor_sub=None, precomputed_tumor_total_reads=None):
            captured_run_kwargs.append(
                {
                    "precomputed_tumor_sub": precomputed_tumor_sub,
                    "precomputed_tumor_total_reads": precomputed_tumor_total_reads,
                }
            )
            return {"sites": 0.1}

    monkeypatch.setattr(run_admixture, "AdmixtureBuilder", _FakeBuilder)
    tumor_sub_bam = tmp_path / "batch" / "mix_test.tumor_sub.bam"
    monkeypatch.setattr(
        "sys.argv",
        [
            "run_admixture.py",
            "--tumor-bam", str(tmp_path / "tumor.bam"),
            "--healthy-bam", str(tmp_path / "healthy.bam"),
            "--reference-fasta", str(tmp_path / "ref.fa"),
            "--snv-vaf-path", f"sites={tmp_path / 'sites.vcf.gz'}",
            "--output-dir", str(tmp_path / "out"),
            "--tumor-fraction", "0.1",
            "--target-coverage", "1.0",
            "--tumor-native-coverage", "35",
            "--healthy-native-coverage", "20",
            "--sample-name", "mix_test",
            "--tumor-sub-bam", str(tumor_sub_bam),
            "--tumor-total-reads", "123456",
        ],
    )

    run_admixture.main()

    assert captured_run_kwargs == [
        {"precomputed_tumor_sub": tumor_sub_bam, "precomputed_tumor_total_reads": 123456}
    ]


def test_main_defaults_precomputed_tumor_arguments_to_none(tmp_path: Path, monkeypatch) -> None:
    captured_run_kwargs = []

    class _FakeBuilder:
        def __init__(self, config):
            pass

        def run(self, precomputed_tumor_sub=None, precomputed_tumor_total_reads=None):
            captured_run_kwargs.append(
                {
                    "precomputed_tumor_sub": precomputed_tumor_sub,
                    "precomputed_tumor_total_reads": precomputed_tumor_total_reads,
                }
            )
            return {"sites": 0.1}

    monkeypatch.setattr(run_admixture, "AdmixtureBuilder", _FakeBuilder)
    monkeypatch.setattr(
        "sys.argv",
        [
            "run_admixture.py",
            "--tumor-bam", str(tmp_path / "tumor.bam"),
            "--healthy-bam", str(tmp_path / "healthy.bam"),
            "--reference-fasta", str(tmp_path / "ref.fa"),
            "--snv-vaf-path", f"sites={tmp_path / 'sites.vcf.gz'}",
            "--output-dir", str(tmp_path / "out"),
            "--tumor-fraction", "0.1",
            "--target-coverage", "1.0",
            "--tumor-native-coverage", "35",
            "--healthy-native-coverage", "20",
            "--sample-name", "mix_test",
        ],
    )

    run_admixture.main()

    assert captured_run_kwargs == [{"precomputed_tumor_sub": None, "precomputed_tumor_total_reads": None}]
