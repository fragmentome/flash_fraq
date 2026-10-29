# CNV-only interim model — labeled interim result

**Status:** interim, 2026-09-29. This is the "one defensible number if the reshaping route
fails" result named in
[cnv-only-interim-model.md](todo/completed/cnv-only-interim-model.md), not the project's
headline result. It is superseded once fragment-size reshaping lands on real tumor data and the
primary (`all`-feature-set) model is retrained/re-evaluated on reshaped mixtures.

## Caveats (read before the numbers)

- **Unreshaped tumor.** The tumor side of every mixture (HCC1395/HCC1395BL sheared gDNA) was
  never put through fragment-size reshaping. The fragmentomics block is confounded with a
  library-prep signal that scales with tumor fraction (audit finding 1,
  [decoy-qc-library-confound-check.md](todo/completed/decoy-qc-library-confound-check.md)). This
  is exactly why this run is CNV-only (`--feature-set cnv_only`): `cnv_logr_bin_*`,
  `arm_zscore_*`, `gc_corrected_total_coverage_window_*` only — no `short_long_ratio_*` or
  `fragment_size_*` columns entered training.
- **ichorCNA low-TF prior.** ichorCNA ran with the `low_tf` profile (not the previous default),
  chosen because these are admixture mixtures spanning nominal tumor fraction 0-0.5, not real
  high-purity tumors. The `ichorCNA_tumor_fraction_estimate` column used as the baseline/
  comparator throughout reflects this prior.
- **1.0x-equivalent low-pass coverage.** As throughout this project's admixture design.
- **Centromere / chromosome-end mask.** `gc_correction.py` and `feature_matrix.pivot_cnv_bins`
  both drop 1-5 Mb windows/bins within `mask_margin_bp=5,000,000` of a centromere or chromosome
  end before the CNV features ever reach the model (see
  [cnv-only-interim-model.md](todo/completed/cnv-only-interim-model.md)'s 2026-09-28 entry).
  Masked bins are absent from the matrix, not present-as-NaN.
- **Decoy QC's check 3 (library-confound) is expected to flag on these mixtures** — see below.
  That is the pre-reshaping baseline this interim result is measured against, not a surprise.

## 1. Feature matrix

Built with `scripts/run_admixture_mixture_feature_matrix.py` from the 495 unreshaped mixtures
(fragments + gc_correct + ichorCNA `low_tf`, chain finished 2026-09-29 20:07).

| | |
| --- | --- |
| `feature_matrix.csv` | `/mnt/8tb_hdd/flash_fraq_derived/admixture_mixtures/feature_matrix.csv` |
| `feature_matrix_summary.json` | `/mnt/8tb_hdd/flash_fraq_derived/admixture_mixtures/feature_matrix_summary.json` |
| `feature_matrix_real_background.csv` | `/mnt/8tb_hdd/flash_fraq_derived/admixture_mixtures/feature_matrix_real_background.csv` |

Sanity checks:

- **495 rows**: 484 `real_background` + 11 `matched_pair`, matching the phenotype-filtered plan.
- **3,699 columns total**, 3,683 primary feature columns; 2,243 `cnv_logr_bin_*`, 41
  `arm_zscore_*`, 465 `gc_corrected_total_coverage_window_*` columns present (masked bins already
  excluded by construction, not zeroed/NaN'd).
- **NaN counts:** only 4 columns carry any NaN across the whole matrix —
  `split_group` (11 NaN, all and only the `matched_pair` rows, which have no subject-level
  train/validation/test split by design — see the matched_pair descriptive report below), and
  three `cnv_logr_bin_*` columns (13, 2, and 1 NaN rows respectively) from bins not called on a
  handful of individual mixtures. 16 of 495 rows have at least one NaN across the full
  cnv/arm/coverage feature block; the LightGBM trainer handles this natively (native NaN
  support), no imputation was applied.
- **Label (`empirical_tumor_fraction`) range per regime:**

  | Regime | n | mean | std | min | 25% | 50% | 75% | max |
  | --- | --- | --- | --- | --- | --- | --- | --- | --- |
  | `real_background` | 484 | 0.138 | 0.149 | 0.0016 | 0.0259 | 0.0823 | 0.2076 | 0.5277 |
  | `matched_pair` | 11 | 0.130 | 0.150 | 0.0039 | 0.0289 | 0.0748 | 0.1717 | 0.4965 |

  (full split_group breakdown: 341 train / 66 validation / 77 test, all in `real_background`;
  `matched_pair` rows carry no split group, see above.)

## 2. Decoy QC (per regime, `--no-fail-on-flag`) — pre-reshaping baseline

Run with `scripts/run_decoy_qc.py --matrix-kind admixture --regime {real_background,matched_pair}
--no-fail-on-flag` against `feature_matrix.csv`.

Reports: `/mnt/8tb_hdd/flash_fraq_derived/admixture_mixtures/decoy_qc_interim/decoy_qc_report_real_background.json`
and `..._matched_pair.json`.

| Check | `real_background` | `matched_pair` |
| --- | --- | --- |
| Fixed-fraction (check 1) | not discriminable (n=176, accuracy 1.0 == baseline 1.0) | not discriminable (n=4, accuracy 1.0 == baseline 1.0) |
| Fraction-correlation (check 2) | skipped — only one source (`HiSeq4000\|HiSeq2500`) in the matrix, no source signal to leak | skipped — only one source (`HiSeq4000\|HiSeq4000`) |
| **Library-confound (check 3)** | **flagged: True**, Spearman ρ = **0.254** (threshold 0.1), n=484, 933 non-CNV features, 2,746 CNV features, 5-fold residualised | **flagged: True**, Spearman ρ = **0.709** (threshold 0.1), n=11, 933 non-CNV features, 2,749 CNV features, 5-fold residualised |
| **`passed`** | **False** | **False** |

This is exactly the expected pre-reshaping failure mode: the fragment/library-prep block still
carries a tumor-fraction-correlated signal in both regimes (worse in `matched_pair`, where n=11
makes the correlation estimate noisier but the confound itself — cell-line gDNA tumor vs. cfDNA
background, within every mixture — is the same mechanism). It is the reason this run trains
`cnv_only`, and it is the number the reshaping route needs to beat once it lands.

Note: `run_decoy_qc.py`'s output filename previously did not include the `--regime`, so running
both regimes into the same `--output-dir` silently overwrote the first regime's report — fixed in
this session (`src/flash_fraq/decoy_qc.py`, `scripts/run_decoy_qc.py`; see "Code changes" below).

## 3. Model training (`feature_set=cnv_only`)

`scripts/train_tf_model.py --feature-set cnv_only`, hyperparameters from `workflow/config/tf_model.yaml`
(n_estimators=300, num_leaves=31, max_depth=-1, learning_rate=0.05, subsample=0.8,
colsample_bytree=0.8, min_child_samples=10, seed=42), against the full `feature_matrix.csv`
(the trainer fits one booster per regime present, never pooled).

| | |
| --- | --- |
| Model root | `/mnt/8tb_hdd/flash_fraq_derived/models/tf_model_cnv_only_interim/` |
| `real_background` model | `tf_model_real_background_primary_cnv_only.joblib`, trained on 341 rows (the `train` split) |
| `matched_pair` model | **skipped** — no rows labeled `split_group='train'` (by design, no booster fit); its descriptive report comes from `tf_eval` instead (§4) |
| Feature importances | `feature_importances.json` (top features: `cnv_logr_bin_6118`, `cnv_logr_bin_6113`, `cnv_logr_bin_7126`, `cnv_logr_bin_85`, `cnv_logr_bin_13089`, ...) |

## 4. Held-out evaluation (`scripts/run_tf_eval.py --feature-set cnv_only`)

Report: `/mnt/8tb_hdd/flash_fraq_derived/models/tf_model_cnv_only_interim/tf_eval_report.json`.
Eval split groups: `validation` + `test` (n=143, `real_background` only, per `tf_eval.yaml`).

### Headline metrics, `real_background`, held-out (n=143)

| Metric | Model (cnv_only) | ichorCNA (low_tf) |
| --- | --- | --- |
| RMSE | **0.00961** | 0.1273 |
| MAE | **0.00629** | 0.0838 |
| Pearson r | 0.9979 | 0.9817 |
| Spearman ρ | 0.9878 | 0.9108 |
| R² | **0.9958** | 0.2668 |

Model-vs-ichorCNA agreement on the same held-out rows: Pearson 0.981, Spearman 0.913.

### Low-fraction slice (true `empirical_tumor_fraction` ≤ 0.05, n=52)

| Metric | Model (cnv_only) | ichorCNA (low_tf) |
| --- | --- | --- |
| RMSE | **0.00687** | 0.01591 |
| MAE | **0.00541** | 0.01312 |
| Pearson r | 0.847 | -0.159 |
| Spearman ρ | 0.819 | -0.108 |
| R² | **0.682** | -0.702 |

ichorCNA is anti-correlated with truth in the low-fraction slice (uninformative there, as
expected of a CNA-based estimator at low purity); the CNV-only model retains a strong positive
correlation and a materially lower RMSE/MAE.

### Detection floor

`detection_floor = 0.01` (1%) — computed against the nominal-`f=0` null pool
(`detection_floor_null_max=0.0`, one-sided Mann-Whitney U per bin, α=0.05), matching the strict
default and at or below ichorCNA's own published ~3% floor.

### ichorCNA comparison, same held-out mixtures

Reported above alongside the headline metrics (§4.1/§4.2): the model beats ichorCNA's own
tumor-fraction call on the identical 143 held-out rows by RMSE (12.6x smaller), MAE (13.3x
smaller), and — most tellingly — is the only one of the two with positive correlation in the
low-fraction slice.

### `matched_pair` descriptive report (n=11, no held-out split by construction)

`tf_eval` produces this automatically when a regime carries no `split_group` (no booster is
fit for it in §3). All numbers here are descriptive, not a held-out evaluation:

| | RMSE | MAE | Pearson | Spearman |
| --- | --- | --- | --- | --- |
| Label vs. nominal fraction | 0.00425 | 0.00353 | 0.9998 | 1.000 |
| ichorCNA vs. label | 0.0960 | 0.0618 | 0.9657 | 0.9727 |
| Leave-one-rung-out ridge (CNV features, 2,749 cols) | 0.0250 | 0.0164 | 0.9952 | 0.9909 |

## 5. Decision-rule read (informal, interim)

Using the same decision rule stated in
[held-out-evaluation-and-ichorcna-comparison.md](todo/completed/held-out-evaluation-and-ichorcna-comparison.md):
RMSE is materially better than ichorCNA's (0.0096 vs. 0.127) and the detection floor (1%) is at
ichorCNA's own published floor while the low-fraction slice shows the model clearly ahead
(positive vs. negative correlation). On CNV features alone, over unreshaped mixtures, this reads
as a **yes** — but it is not the headline result: the fragmentomics block is disabled here
specifically because it is broken (§2), so this number cannot yet speak to the project's actual
question (whether fragment-ratio features add anything ichorCNA-style CNA features don't). It is
the before/after reference for the reshaping route.

## Code changes made in this session

- `src/flash_fraq/decoy_qc.py`: `DecoyClassifierQc.report_path()` now suffixes the report
  filename with `--regime` when one is set (`decoy_qc_report_<regime>.json`), instead of always
  writing `decoy_qc_report.json`. Bug found while running step 2 exactly as documented (both
  regimes into the same `--output-dir`): the second run silently overwrote the first regime's
  report. `scripts/run_decoy_qc.py` updated to print the actual (possibly regime-scoped) path.
  Cohort-kind runs (`regime=None`) are unaffected — same filename as before.
  New test: `tests/test_decoy_qc.py::test_regime_scoped_runs_do_not_overwrite_each_others_report`.
  Full suite: 258 passed (257 previously + 1 new), no regressions.

## Output paths (summary)

| What | Path |
| --- | --- |
| Feature matrix | `/mnt/8tb_hdd/flash_fraq_derived/admixture_mixtures/feature_matrix.csv` |
| Feature matrix summary | `/mnt/8tb_hdd/flash_fraq_derived/admixture_mixtures/feature_matrix_summary.json` |
| Feature matrix (real_background subset) | `/mnt/8tb_hdd/flash_fraq_derived/admixture_mixtures/feature_matrix_real_background.csv` |
| Decoy QC (real_background) | `/mnt/8tb_hdd/flash_fraq_derived/admixture_mixtures/decoy_qc_interim/decoy_qc_report_real_background.json` |
| Decoy QC (matched_pair) | `/mnt/8tb_hdd/flash_fraq_derived/admixture_mixtures/decoy_qc_interim/decoy_qc_report_matched_pair.json` |
| Models | `/mnt/8tb_hdd/flash_fraq_derived/models/tf_model_cnv_only_interim/` |
| Feature importances | `/mnt/8tb_hdd/flash_fraq_derived/models/tf_model_cnv_only_interim/feature_importances.json` |
| Eval report | `/mnt/8tb_hdd/flash_fraq_derived/models/tf_model_cnv_only_interim/tf_eval_report.json` |
| This write-up | `docs/interim-cnv-only-results.md` |

## Housekeeping note (not deleted, proposed for cleanup)

While diagnosing the decoy-QC overwrite bug (§2), an earlier attempt at working around it wrote
per-regime subdirectories that are now redundant with the fixed flat-file naming:
`/mnt/8tb_hdd/flash_fraq_derived/admixture_mixtures/decoy_qc_interim/real_background/decoy_qc_report.json`
and `.../decoy_qc_interim/matched_pair/decoy_qc_report.json`, plus the original clobbered
top-level `.../decoy_qc_interim/decoy_qc_report.json` (this last one holds only the
`matched_pair` result — the `real_background` run that produced it first was overwritten before
the fix landed). Contents are identical to the canonical `decoy_qc_report_{regime}.json` files
above (or, for the clobbered one, superseded by them); safe to delete, not deleted here.
