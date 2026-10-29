"""Driver for the phase2 §6 decoy-classifier batch-effect gate.

`--matrix-kind cohort` runs the gate against the pre-admixture cohort matrix
(source label `source_dataset|platform`, fraction from ichorCNA); `admixture`
runs it against the §9 mixed-sample matrix (source label
`platform_tumor|platform_healthy`, fraction from SNV VAFs). The cohort variant is
the weaker of the two -- see DecoyClassifierQcConfig.for_cohort.

Exits non-zero when the gate fails, so a Snakemake run stops rather than leaving
a failing gate to be noticed in the report. `--no-fail-on-flag` records the
verdict without stopping the run.
"""

from __future__ import annotations

import sys
from argparse import ArgumentParser
from pathlib import Path

from flash_fraq.decoy_qc import DecoyClassifierQc, DecoyClassifierQcConfig


def main() -> None:
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--feature-matrix-path", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--matrix-kind", choices=("cohort", "admixture"), default="cohort")
    parser.add_argument(
        "--regime",
        choices=("matched_pair", "real_background"),
        default=None,
        help=(
            "Required for an admixture matrix that carries more than one regime: the gate "
            "refuses to run pooled across matched_pair/real_background (audit finding 4)."
        ),
    )
    parser.add_argument("--correlation-threshold", type=float, default=0.1)
    parser.add_argument("--library-confound-correlation-threshold", type=float, default=0.1)
    parser.add_argument("--fixed-fraction-target", type=float, default=0.0)
    parser.add_argument("--fixed-fraction-tolerance", type=float, default=0.05)
    parser.add_argument("--max-missing-fraction", type=float, default=0.0)
    parser.add_argument(
        "--n-splits",
        type=int,
        default=5,
        help="Folds for check 2/check 3's out-of-fold predictions; skipped if smaller than the data",
    )
    parser.add_argument("--no-fail-on-flag", action="store_true")
    args = parser.parse_args()

    if args.matrix_kind == "admixture" and args.regime is None:
        parser.error(
            "--regime is required with --matrix-kind admixture; the gate refuses to run pooled "
            "across regimes (pass matched_pair or real_background explicitly)"
        )

    overrides = {
        "correlation_threshold": args.correlation_threshold,
        "library_confound_correlation_threshold": args.library_confound_correlation_threshold,
        "fixed_fraction_target": args.fixed_fraction_target,
        "fixed_fraction_tolerance": args.fixed_fraction_tolerance,
        "max_missing_fraction": args.max_missing_fraction,
        "n_splits": args.n_splits,
        "regime": args.regime,
    }
    if args.matrix_kind == "cohort":
        config = DecoyClassifierQcConfig.for_cohort(
            feature_matrix_path=args.feature_matrix_path, output_dir=args.output_dir, **overrides
        )
    else:
        config = DecoyClassifierQcConfig(
            feature_matrix_path=args.feature_matrix_path, output_dir=args.output_dir, **overrides
        )

    qc = DecoyClassifierQc(config)
    report = qc.run()

    print(f"decoy QC report: {qc.report_path()}")
    print(f"fixed-fraction check: {report.fixed_fraction_check}")
    print(f"correlation check: {report.correlation_check}")
    print(f"passed: {report.passed}")

    if not report.passed and not args.no_fail_on_flag:
        sys.exit(1)


if __name__ == "__main__":
    main()
