"""Build the README figures and the metrics tables behind them.

Inputs are the artifacts the pipeline already writes: the admixture feature
matrix CSV and the trained model joblib from `train_tf_model.py`. Predictions
are recomputed from those two (the eval report does not store per-sample
predictions), restricted to the held-out `--eval-split-groups` rows of one
regime -- the same rows `run_tf_eval.py` scores.

    python scripts/make_figures.py --feature-matrix M.csv --model-root MODELS \\
        --feature-set cnv_only --out-dir docs/figures --results-dir docs/results

Figures whose source data this repo does not persist are skipped with a message
rather than faked:
  * fragment-length distribution -- extraction keeps only median/mode/entropy
    per sample, not the TLEN histogram;
  * training curves -- LightGBM fits are not logged per boosting round.
"""

from __future__ import annotations

import json
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy
from pandas import DataFrame, read_csv
from scipy.stats import pearsonr, spearmanr

from flash_fraq.tf_model import TrainedTumorFractionModel

LABEL = "empirical_tumor_fraction"
BASELINE = "ichorCNA_tumor_fraction_estimate"
DPI = 200
# Okabe-Ito: colorblind-safe.
MODEL_COLOR, BASELINE_COLOR, ACCENT = "#0072B2", "#D55E00", "#009E73"
CATEGORICAL = ["#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00"]
# Bins on the true label; the last edge is open.
BIN_EDGES = (0.0, 0.01, 0.03, 0.10, numpy.inf)
BIN_LABELS = ("<1%", "1-3%", "3-10%", ">10%")
N_BOOTSTRAP = 2000


def apply_style() -> None:
    plt.rcParams.update({
        "figure.dpi": DPI, "savefig.dpi": DPI, "savefig.bbox": "tight",
        "font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "grid.alpha": 0.25, "legend.frameon": False,
    })


def point_metrics(prediction: numpy.ndarray, truth: numpy.ndarray) -> dict[str, float]:
    return {
        "n": int(len(truth)),
        "pearson_r": float(pearsonr(prediction, truth)[0]) if len(truth) > 2 else float("nan"),
        "spearman_rho": float(spearmanr(prediction, truth)[0]) if len(truth) > 2 else float("nan"),
        "mae": float(numpy.mean(numpy.abs(prediction - truth))),
        "rmse": float(numpy.sqrt(numpy.mean((prediction - truth) ** 2))),
    }


def bootstrap_ci(values: numpy.ndarray, rng: numpy.random.Generator) -> tuple[float, float]:
    draws = rng.choice(values, size=(N_BOOTSTRAP, len(values)), replace=True).mean(axis=1)
    return float(numpy.percentile(draws, 2.5)), float(numpy.percentile(draws, 97.5))


class FigureBuilder:
    def __init__(self, rows: DataFrame, prediction: numpy.ndarray, out_dir: Path, color_by: str) -> None:
        self.rows = rows.reset_index(drop=True)
        self.truth = self.rows[LABEL].to_numpy(float)
        self.prediction = prediction
        self.baseline = self.rows[BASELINE].to_numpy(float) if BASELINE in self.rows else None
        self.out_dir = out_dir
        self.color_by = color_by if color_by in self.rows.columns else None
        self.rng = numpy.random.default_rng(0)
        self.written: list[Path] = []

    def _save(self, fig: plt.Figure, name: str) -> None:
        path = self.out_dir / name
        fig.savefig(path)
        plt.close(fig)
        self.written.append(path)

    def _groups(self) -> list[tuple[str, numpy.ndarray]]:
        if self.color_by is None:
            return [("held-out", numpy.ones(len(self.rows), bool))]
        values = self.rows[self.color_by].astype(str)
        return [(name, (values == name).to_numpy()) for name in sorted(values.unique())]

    def predicted_vs_true(self) -> dict[str, float]:
        stats = point_metrics(self.prediction, self.truth)
        fig, ax = plt.subplots(figsize=(4.8, 4.8))
        top = max(self.truth.max(), self.prediction.max()) * 1.05
        ax.plot([0, top], [0, top], color="0.4", lw=1, ls="--", label="identity")
        for (name, mask), color in zip(self._groups(), CATEGORICAL):
            ax.scatter(self.truth[mask], self.prediction[mask], s=16, alpha=0.75, color=color, label=name)
        ax.set(xlim=(0, top), ylim=(0, top), xlabel="True tumor fraction (empirical label)",
               ylabel="Predicted tumor fraction", aspect="equal")
        ax.text(0.04, 0.96, f"Pearson r = {stats['pearson_r']:.3f}\nSpearman ρ = {stats['spearman_rho']:.3f}\n"
                f"MAE = {stats['mae']:.4f}\nn = {stats['n']}",
                transform=ax.transAxes, va="top", fontsize=9)
        ax.legend(loc="lower right")
        self._save(fig, "predicted_vs_true.png")
        return stats

    def bland_altman(self) -> dict[str, float]:
        mean, diff = (self.prediction + self.truth) / 2, self.prediction - self.truth
        bias, spread = float(diff.mean()), float(diff.std(ddof=1))
        fig, ax = plt.subplots(figsize=(5.6, 4))
        for (name, mask), color in zip(self._groups(), CATEGORICAL):
            ax.scatter(mean[mask], diff[mask], s=16, alpha=0.75, color=color, label=name)
        for level, style in ((bias, "-"), (bias + 1.96 * spread, "--"), (bias - 1.96 * spread, "--")):
            ax.axhline(level, color="0.3", lw=1, ls=style)
        ax.set(xlabel="Mean of predicted and true tumor fraction",
               ylabel="Predicted − true tumor fraction")
        ax.set_title(f"bias {bias:+.4f}, 95% limits ±{1.96 * spread:.4f}", fontsize=10)
        ax.legend()
        self._save(fig, "bland_altman.png")
        return {"bias": bias, "loa_lower": bias - 1.96 * spread, "loa_upper": bias + 1.96 * spread}

    def error_by_bin(self) -> list[dict[str, object]]:
        records: list[dict[str, object]] = []
        series = {"model": self.prediction}
        if self.baseline is not None:
            series["ichorCNA"] = self.baseline
        for label, (lo, hi) in zip(BIN_LABELS, zip(BIN_EDGES[:-1], BIN_EDGES[1:])):
            mask = (self.truth >= lo) & (self.truth < hi)
            if mask.sum() < 3:
                continue
            for name, values in series.items():
                abs_err = numpy.abs(values[mask] - self.truth[mask])
                low, high = bootstrap_ci(abs_err, self.rng)
                records.append({"bin": label, "method": name, "n": int(mask.sum()),
                                "mae": float(abs_err.mean()), "mae_ci_low": low, "mae_ci_high": high})
        fig, ax = plt.subplots(figsize=(5.6, 4))
        bins = [label for label in BIN_LABELS if any(r["bin"] == label for r in records)]
        colors = {"model": MODEL_COLOR, "ichorCNA": BASELINE_COLOR}
        width = 0.38
        for offset, name in enumerate(series):
            sel = {r["bin"]: r for r in records if r["method"] == name}
            xs = numpy.arange(len(bins)) + (offset - (len(series) - 1) / 2) * width
            mae = numpy.array([sel[b]["mae"] for b in bins])
            err = numpy.array([[sel[b]["mae"] - sel[b]["mae_ci_low"] for b in bins],
                               [sel[b]["mae_ci_high"] - sel[b]["mae"] for b in bins]])
            ax.bar(xs, mae, width, yerr=err, color=colors[name], label=name, capsize=3)
        counts = {b: next(r["n"] for r in records if r["bin"] == b) for b in bins}
        ax.set_xticks(numpy.arange(len(bins)), [f"{b}\n(n={counts[b]})" for b in bins])
        ax.set(xlabel="True tumor fraction bin", ylabel="Mean absolute error (tumor fraction)")
        ax.legend()
        ax.set_title("95% bootstrap CI over mixtures", fontsize=10)
        self._save(fig, "error_by_tumor_fraction_bin.png")
        return records


def cohort_overview(matrix: DataFrame, out_dir: Path) -> list[Path]:
    split = matrix["split_group"].fillna("none") if "split_group" in matrix else None
    fig, axes = plt.subplots(1, 2, figsize=(8, 3.6))
    counts = matrix["regime"].value_counts()
    axes[0].bar(counts.index, counts.values, color=CATEGORICAL[: len(counts)])
    axes[0].set(xlabel="Mixture regime", ylabel="Mixtures (samples)")
    for x, value in enumerate(counts.values):
        axes[0].text(x, value, str(value), ha="center", va="bottom")
    if split is not None:
        order = ["train", "validation", "test", "none"]
        by_split = split.value_counts().reindex([s for s in order if s in set(split)])
        axes[1].bar(by_split.index, by_split.values, color=CATEGORICAL[: len(by_split)])
        axes[1].set(xlabel="Split group", ylabel="Mixtures (samples)")
        for x, value in enumerate(by_split.values):
            axes[1].text(x, value, str(value), ha="center", va="bottom")
    path = out_dir / "cohort_overview.png"
    fig.savefig(path)
    plt.close(fig)
    return [path]


def main() -> None:
    parser = ArgumentParser(description=__doc__, formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("--feature-matrix", required=True, type=Path)
    parser.add_argument("--model-root", required=True, type=Path)
    parser.add_argument("--regime", default="real_background")
    parser.add_argument("--variant", default="primary")
    parser.add_argument("--feature-set", default="all", choices=["all", "cnv_only", "fragment_only"])
    parser.add_argument("--eval-split-groups", nargs="+", default=["validation", "test"])
    parser.add_argument("--color-by", default="split_group", help="Matrix column used to color points.")
    parser.add_argument("--out-dir", type=Path, default=Path("docs/figures"))
    parser.add_argument("--results-dir", type=Path, default=Path("docs/results"))
    args = parser.parse_args()

    apply_style()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    args.results_dir.mkdir(parents=True, exist_ok=True)

    matrix = read_csv(args.feature_matrix)
    suffix = "" if args.feature_set == "all" else f"_{args.feature_set}"
    model = TrainedTumorFractionModel.load(
        args.model_root / f"tf_model_{args.regime}_{args.variant}{suffix}.joblib")
    rows = matrix[(matrix["regime"] == args.regime) & matrix["split_group"].isin(args.eval_split_groups)]
    if rows.empty:
        raise SystemExit(f"no held-out rows for regime {args.regime!r} in {args.eval_split_groups}")
    prediction = model.predict(rows).to_numpy(float)

    builder = FigureBuilder(rows, prediction, args.out_dir, args.color_by)
    metrics: dict[str, object] = {
        "regime": args.regime, "feature_set": args.feature_set,
        "eval_split_groups": args.eval_split_groups,
        "model": builder.predicted_vs_true(),
        "bland_altman": builder.bland_altman(),
        "error_by_bin": builder.error_by_bin(),
    }
    if BASELINE in rows:
        metrics["ichorCNA"] = point_metrics(rows[BASELINE].to_numpy(float), builder.truth)
    builder.written += cohort_overview(matrix, args.out_dir)

    (args.results_dir / "metrics.json").write_text(json.dumps(metrics, indent=2))
    flat = [{"scope": "overall", "method": "model", **metrics["model"]}]
    if "ichorCNA" in metrics:
        flat.append({"scope": "overall", "method": "ichorCNA", **metrics["ichorCNA"]})
    flat += [{"scope": f"bin {r['bin']}", "method": r["method"], "n": r["n"], "mae": r["mae"],
              "mae_ci_low": r["mae_ci_low"], "mae_ci_high": r["mae_ci_high"]} for r in metrics["error_by_bin"]]
    DataFrame(flat).to_csv(args.results_dir / "metrics.csv", index=False)

    print("wrote:", *[str(p) for p in builder.written], sep="\n  ")
    print("skipped: fragment_length_distribution (no per-sample TLEN histograms are persisted), "
          "training_curves (LightGBM rounds are not logged)")


if __name__ == "__main__":
    main()
