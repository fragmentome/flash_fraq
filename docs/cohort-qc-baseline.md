# Pre-admixture cohort QC — first run (2026-08-31)

The phase2 §6 decoy-classifier gate, run on the real Snyder + Cristiano samples before any
admixture exists. Target `cohort_qc_complete`; module reference in
[pipeline-modules.md](pipeline-modules.md) §7.

```
snakemake --workflow-profile workflow/profiles/local cohort_qc_complete
```

Snakemake strips `PYTHONPATH` inside `--use-conda` jobs, so run this from a real
`conda activate snakemake` shell (as `scripts/drivers/smoke_local_features_driver.sh` does),
not by putting the snakemake env on `PATH` by hand — the latter leaves the job running the
snakemake env's interpreter, which has no `flash_fraq`.

## Cohort matrix

`/mnt/8tb_hdd/flash_fraq_derived/cohort/cohort_feature_matrix.csv` — 547 samples
(503 Cristiano, 44 Snyder) x 4359 columns, of which 4354 are features and **4346 are
present for every sample**. Whole run: ~34 s.

The 8 ragged columns are all `cnv_logr_bin_*`, each missing for 1-6 samples (bins ichorCNA
filtered out of individual fits). Window and arm blocks come out fully dense at cohort
level.

Two arms never appear as columns: `chr15p` (no 5 Mb windows on it at all) and `chr22p`
(present as a window, but its short:long ratio is undefined in every sample, so the z-score
is NaN everywhere). 44 arm columns, not the §9 placeholder's 39.

## Gate result

| Check | Result |
| --- | --- |
| 1. Source discriminable at f=0 (273 samples, 4351 features) | **accuracy 1.000 vs. 0.890 baseline — yes, perfectly** |
| 2. corr(out-of-fold predicted-source-probability, tumor fraction), 4346 features | max abs 0.0181, threshold 0.10 — **not flagged** |

`passed: true`, but read both lines together before taking that at face value:

- **Check 1 is a maximal fail.** Snyder and Cristiano separate perfectly on fragmentomics
  features at zero ichorCNA tumor fraction. Platform/center signal is unambiguously present
  in the feature matrix. That is what §6 predicts for two cohorts from different centers,
  and it is exactly the axis that must never line up with the label.
- **Check 2 degenerates because check 1 succeeds so completely.** A classifier that
  separates the sources perfectly emits saturated probabilities, so check 2 collapses to
  the plain correlation between source membership and the label: `corr(source indicator,
  ichorCNA fraction) = 0.02197` against check 2's reported 0.01809. The check is measuring
  "do these two cohorts happen to have similar ichorCNA tumor-fraction distributions"
  (they do: mean 0.051 Cristiano vs. 0.046 Snyder) rather than anything about feature
  leakage. The report says so itself:
  `degenerate_with_source_indicator: true`, set when the classifier-derived and
  indicator-derived numbers agree to within `degeneracy_tolerance` (0.005).
- **The label is ichorCNA's own estimate**, not an SNV-VAF empirical fraction, so a pass
  says "source does not predict ichorCNA's CNV call", not "source does not predict tumor
  burden".

So: **necessary, not sufficient.** `fail_on_flag` is false in
`workflow/config/cohort_qc.yaml` for that reason. The real gate is this same code with
`label_column="empirical_tumor_fraction"` (`DecoyClassifierQcConfig.from_sources`) once
admixture exists — and the check-1 result above is advance warning that mixture design has
to keep source membership orthogonal to mixing ratio, because the features will otherwise
carry the source right through.

## Check 2 is out-of-fold (updated 2026-08-31)

The first run of check 2 fitted the classifier on the full matrix and then called
`predict_proba` on those same rows. Over 4354 features and 547 samples a
`LogisticRegression` separates its own training set perfectly, so those probabilities were
saturated and the reported number (0.02199) was the source indicator (0.02198) to four
decimal places — the check carried no information the label's own class means did not.

It now scores every sample out-of-fold: `cross_val_predict(..., method="predict_proba")`
under a shuffled `StratifiedKFold(n_splits=5, random_state=0)`, with the `StandardScaler`
kept inside the pipeline so it refits per fold rather than leaking held-out means into
training. `n_splits` and `fold_random_state` are written into the report, so the number is
reproducible from the report alone; `n_splits` is configurable in
`workflow/config/cohort_qc.yaml` and the check *skips with a reason* (rather than letting
sklearn raise, or silently using fewer folds) if the smallest source has fewer samples than
folds. A skipped check writes `max_abs_correlation: null`, never a 0.0 that could be
misread as a measurement.

| | in-sample (first run) | out-of-fold (5 folds) |
| --- | --- | --- |
| max abs corr | 0.02199 | **0.01809** |
| corr(source indicator, fraction) | 0.02198 | 0.02197 |
| flagged | no | no |

The verdict does not change, and the check is still degenerate — 0.01809 against 0.02197 is
within the 0.005 tolerance, because the two cohorts really are separable enough that even a
model that never saw a row still calls its source almost perfectly. What changed is that
0.01809 now means what the report claims it means. Whole run: ~22 s.

## What this changed downstream (2026-08-31)

Check 1's result drove a design decision, recorded in full in
[MAT-13](todo/todo/admixture-and-empirical-tumor-fraction.md#superseding-decision-2026-08-31--background-demoted-to-the-matched-normal):

**The admixture healthy background moved from Snyder to HCC1395BL, the SEQC2 matched normal.** With
a single tumor source, mixing tumor from SEQC2 into background from another cohort makes tumor
fraction *itself* a batch mixing ratio — raising `f` raises one center's share of the library by
construction. Accuracy 1.000 at `f=0` says the features read that off perfectly, so the model could
hit the regression target without using tumor biology at all. Stratifying backgrounds across Snyder
and Cristiano does not help, because the tumor side never varies. Mixing tumor against its own
matched normal holds batch composition constant across the whole ladder and removes the confound by
construction.

**Snyder stays in this matrix.** It was demoted only as a mixing source; the 44 samples remain in the
cohort feature matrix as the gate's second source, and as a frozen cross-center evaluation set. That
matters for this document specifically: with one source, check 1 is undefined and check 2 skips
outright (`skipped_reason`, `max_abs_correlation: null`), so the gate would go blind.

Also open, and not resolved by the above: matched-pair mixing leaves one biological background
subject, which makes subject-level train/val/test splitting impossible. See MAT-13 for the
two-regime design question that follows from it.
