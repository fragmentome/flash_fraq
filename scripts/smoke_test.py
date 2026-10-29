"""End-to-end smoke test on SYNTHETIC data (no cohort data, no cluster, ~1 minute).

Builds a small §9-shaped feature matrix with a known signal, then runs the real
`train_tf_model.py` and `run_tf_eval.py` CLIs on it and checks the held-out
correlation clears a floor. This verifies the install and the train/eval wiring;
the numbers it prints say nothing about real cfDNA performance.

    python scripts/smoke_test.py [--workdir DIR]
"""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
from argparse import ArgumentParser
from pathlib import Path

import numpy
from pandas import DataFrame

SCRIPTS_DIR = Path(__file__).resolve().parent
N_ROWS = 300
N_NOISE_FEATURES = 20
CORRELATION_FLOOR = 0.8


class SyntheticMatrixBuilder:
    """Synthetic stand-in for the real feature matrix: same column names and
    split semantics, a few CNV-like columns that scale with the true fraction,
    and a noisier ichorCNA-like column as the baseline."""

    def __init__(self, n_rows: int = N_ROWS, seed: int = 0) -> None:
        self.n_rows = n_rows
        self.rng = numpy.random.default_rng(seed)

    def build(self) -> DataFrame:
        rng, n = self.rng, self.n_rows
        fraction = numpy.concatenate([numpy.zeros(n // 10), rng.uniform(0.0, 0.5, n - n // 10)])
        rng.shuffle(fraction)
        columns: dict[str, object] = {
            "sample_id": [f"synthetic_{i:04d}" for i in range(n)],
            "regime": ["real_background"] * n,
            "split_group": rng.choice(["train", "validation", "test"], size=n, p=[0.7, 0.15, 0.15]),
            "nominal_tumor_fraction": fraction,
            "empirical_tumor_fraction": fraction,
            "ichorCNA_tumor_fraction_estimate": fraction + rng.normal(0.0, 0.12, n),
        }
        for index in range(1, 6):
            columns[f"cnv_logr_bin_{index}"] = index * fraction + rng.normal(0.0, 0.05, n)
        for index in range(N_NOISE_FEATURES):
            columns[f"gc_corrected_total_coverage_window_{index}"] = rng.normal(10.0, 0.3, n)
        return DataFrame(columns)


def run(command: list[str]) -> None:
    print("$", " ".join(command), flush=True)
    subprocess.run(command, check=True)


def main() -> None:
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--workdir", type=Path, help="Keep artifacts here (default: a temp dir).")
    args = parser.parse_args()

    workdir = args.workdir or Path(tempfile.mkdtemp(prefix="flash_fraq_smoke_"))
    workdir.mkdir(parents=True, exist_ok=True)
    matrix_path = workdir / "feature_matrix.csv"
    model_root = workdir / "models"
    report_path = workdir / "tf_eval_report.json"

    SyntheticMatrixBuilder().build().to_csv(matrix_path, index=False)
    run([sys.executable, str(SCRIPTS_DIR / "train_tf_model.py"),
         "--feature-matrix-path", str(matrix_path), "--output-root", str(model_root),
         "--feature-set", "cnv_only", "--n-estimators", "100"])
    run([sys.executable, str(SCRIPTS_DIR / "run_tf_eval.py"),
         "--feature-matrix-path", str(matrix_path), "--model-root", str(model_root),
         "--feature-set", "cnv_only", "--report-path", str(report_path)])

    result = json.loads(report_path.read_text())["real_background"]
    pearson = result["model_metrics"]["pearson"]
    print(f"\nSYNTHETIC data: n_eval={result['n_eval_rows']} model pearson={pearson:.3f} "
          f"rmse={result['model_metrics']['rmse']:.4f} "
          f"ichorCNA-like rmse={result['ichorcna_metrics']['rmse']:.4f}")
    print(f"artifacts in {workdir}")
    if pearson < CORRELATION_FLOOR:
        sys.exit(f"smoke test FAILED: held-out Pearson {pearson:.3f} < {CORRELATION_FLOOR}")
    print("smoke test OK")


if __name__ == "__main__":
    main()
