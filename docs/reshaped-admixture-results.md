# Reshaped-admixture held-out results

**Status:** 2026-10-05. First held-out evaluation on mixtures whose tumor side was put through
fragment-size reshaping. Supersedes [interim-cnv-only-results.md](interim-cnv-only-results.md) as
the current number. **Decoy-QC gate overridden** for this run — see §2 before quoting any number.
Partial answer to [motivation-and-goal.md](motivation-and-goal.md) §1 in §6: fragment-ratio
features alone do not recover tumor fraction below ~3–5% on this design; no near-diploid test yet (§5).

## 1. What was run

| | |
| --- | --- |
| Tumor source | HCC1395 NV_1+NV_3 pooled, MAPQ≥30, fragment-size reshaped to 0.72x |
| Mixtures | 495: 484 `real_background` (rebuilt) + 11 `matched_pair` (symlinked from the original, *unreshaped* build) |
| Per-mixture features | fragments → GC correction → ichorCNA (`low_tf` profile), 484/484 OK each |
| Feature matrix | `/mnt/8tb_hdd/flash_fraq_derived/admixture_mixtures_reshaped/feature_matrix.csv` (3,681 features after QC) |
| Models / reports | `/mnt/8tb_hdd/flash_fraq_derived/models/tf_model_reshaped/` — `tf_eval_report_{all,cnv_only}.json`, `feature_importances_{all,cnv_only}.json` |
| Driver | `scripts/drivers/reshaped_downstream_chain.sh`, finished 2026-10-05 10:32 with `OVERRIDE_DECOY_GATE=1` |
| Split | `real_background` only: 341 train, 143 held-out (validation + test), split by subject |

## 2. Decoy-QC gate: overridden, and why

The `real_background` gate failed on check 3 (library-confound probe) only: residualised
out-of-fold Spearman ρ = **−0.117** vs. the |0.1| threshold. Checks 1–2 did not flag.

A permutation null shows the check itself is miscalibrated, not the data:

| | ρ |
| --- | --- |
| Observed | −0.117 |
| KFold seed 0–9 (same data) | −0.173 … +0.016 |
| 50 label permutations (no signal) | −0.270 ± 0.044; 50/50 exceed \|0.1\|; 50/50 exceed \|observed\| |

Out-of-fold ridge predictions on an unlearnable target are biased negative, so a fixed |0.1|
threshold flags pure noise. The observed value is *closer to zero* than every permutation.
Evidence and script: `/mnt/8tb_hdd/flash_fraq_derived/admixture_mixtures_reshaped/decoy_qc/library_confound_permutation_null.{txt,py}`.
Follow-up (not done): calibrate check 3 against the permutation null instead of a fixed threshold.

Reshaping did remove the confound it targeted. Raw per-feature Spearman ρ against nominal
tumor fraction:

| Feature | Original build | Reshaped |
| --- | --- | --- |
| `fragment_size_entropy` | 0.976 | −0.072 |
| `fragment_size_median` | 0.722 | −0.016 |

The `matched_pair` gate also fails (ρ = 0.71, `fragment_size_entropy` ρ = −0.99), as expected:
those 11 mixtures still use the unreshaped tumor and n=11 is descriptive only.

## 3. Headline metrics — `real_background`, held-out (n=143)

| Metric | Model `all` | Model `cnv_only` | ichorCNA (`low_tf`) |
| --- | --- | --- | --- |
| RMSE | **0.0104** | 0.0123 | 0.0700 |
| MAE | **0.0042** | 0.0045 | 0.0432 |
| Pearson r | 0.9956 | 0.9940 | 0.9844 |
| Spearman ρ | 0.9876 | 0.9876 | 0.9088 |

### Low-fraction slice (true tumor fraction < 0.05, n=68)

| Metric | Model `all` | Model `cnv_only` | ichorCNA (`low_tf`) |
| --- | --- | --- | --- |
| RMSE | 0.0046 | **0.0045** | 0.0144 |
| MAE | 0.0037 | 0.0037 | 0.0113 |
| Pearson r | 0.938 | 0.938 | 0.478 |
| Spearman ρ | 0.917 | 0.915 | 0.359 |

Disagreement table (19 rows with the largest model–ichorCNA gap): the model is closer to the
true label in 19/19, mostly at nominal f=0.5 where ichorCNA under-calls (~0.18 vs. true ~0.39).

### `all` model: share of total split importance by feature family

| Family | Share |
| --- | --- |
| `cnv_logr_bin_*` | 89.8% |
| `gc_corrected_total_coverage_window_*` | 6.9% |
| `arm_zscore_*` | 2.3% |
| `gc_corrected_short_coverage_window_*` | 0.7% |
| `short_long_ratio_window_*` | 0.3% |
| `fragment_size_*` | 0.0% |

## 4. Interpretation

1. **Fragment features add nothing on top of CNV.** `all` and `cnv_only` are within 0.002 RMSE
   overall and identical on the low-fraction slice; 99% of the `all` model's importance is
   copy-number/coverage. Now that the library-prep confound is gone, the booster has no reason
   to use the fragmentomic block — the CNV block already explains tumor fraction almost fully.
2. **The win over ichorCNA is real but flattering.** ~7x lower RMSE overall, ~3x on the
   low-fraction slice, and ichorCNA's rank correlation collapses below 5% (ρ 0.36) while the
   model holds (ρ 0.92). But every mixture uses the same tumor, HCC1395, so the supervised model
   is learning *that one karyotype's* fixed CNV profile; ichorCNA is unsupervised and fits each
   sample from scratch. This measures "a model trained on this tumor's CNVs" vs. "a generic CNV
   caller", not fragmentomics vs. CNV.
3. **Consistent with the interim CNV-only run.** Unreshaped interim `cnv_only`: RMSE 0.0096 vs.
   ichorCNA 0.127. The reshaped rebuild changes the tumor side's fragment sizes, not its CNVs, so
   the CNV model's performance was expected to carry over, and it does. ichorCNA's RMSE improved
   (0.127 → 0.070); the cause (different tumor pool, MAPQ≥30 filter, relabelled empirical
   fractions) was not investigated.

## 5. What this does and does not answer

§1 asks whether fragment-ratio features recover tumor fraction **where ichorCNA's CNA-based
estimate is unidentifiable**. This run cannot answer that:

- **`all` says nothing about fragments alone**, because CNV features dominate it. A separate
  `fragment_only` feature set (`short_long_ratio_window_*` + `fragment_size_*`; short-fragment
  *coverage* excluded because it tracks copy number) was added and run — see §6.
- **No near-diploid axis.** HCC1395 is highly aneuploid, so CNV signal is always present; the
  "ichorCNA unidentifiable" regime does not exist in this design (also flagged in the eval
  report's own `notes`). Answering §1 as posed needs a design change: e.g. CNV-masked features
  or a near-diploid tumor source.

## 6. Fragment-only results

2026-10-05 11:23, same 341/143 split, same gate override. Features: 465
`short_long_ratio_window_*` + 3 `fragment_size_*` (468 total). Report:
`/mnt/8tb_hdd/flash_fraq_derived/models/tf_model_reshaped/tf_eval_report_fragment_only.json`.

| Metric | `fragment_only` | `cnv_only` | ichorCNA (`low_tf`) |
| --- | --- | --- | --- |
| RMSE, held-out (n=143) | 0.0426 | 0.0123 | 0.0700 |
| MAE | 0.0317 | 0.0045 | 0.0432 |
| Pearson r | 0.928 | 0.994 | 0.984 |
| Spearman ρ | 0.799 | 0.988 | 0.909 |
| RMSE, f < 0.05 (n=68) | 0.0328 | 0.0045 | **0.0144** |
| Spearman ρ, f < 0.05 | 0.322 | 0.915 | 0.359 |
| Detection floor | 3% | — | — |

Top features: `short_long_ratio_window_326` (334), `fragment_size_entropy` (333), then
`short_long_ratio_window_{436,222,153,467}`.

### Reading

1. **Fragment features carry real signal over the full range.** Fragment features alone track tumor fraction
   across 0–0.5 (r 0.93, ρ 0.80) and beat ichorCNA's overall RMSE, which is driven by its
   under-calling at f=0.5. In 10 of the 12 rows where the model and ichorCNA disagree most, the
   fragment-only model is closer to the true label.
2. **They do not recover low tumor fraction.** Below 5% the fragment-only model is *worse*
   than ichorCNA: RMSE 0.033 vs. 0.014, and its rank correlation (ρ 0.32) is no better than
   ichorCNA's (0.36). Its detection floor is 3%, vs. the 1% reported for CNV models. This is where §1's
   question lives, and on this design the answer is **no**: the fragmentomic signal dies
   between ~3% and ~5% tumor fraction, above where CNV-based estimates still work.
3. **Even this is an upper bound.** Two things could inflate it:
   - `short_long_ratio` windows are not fully CNV-blind. Tumor DNA is shorter, and a window's
     tumor share scales with its copy number, so HCC1395's CNVs still modulate the ratio per
     window.
   - `fragment_size_entropy` is the #2 feature even though its raw correlation with f dropped to
     −0.07 after reshaping. Either the relationship is non-monotonic, or reshaping left a residual
     library signature. Not investigated.

   In both cases a truly CNV-free fragment signal would be weaker, not stronger.
4. **Still no near-diploid test.** The §1 regime "ichorCNA unidentifiable" is still absent from
   this design (§5). What this run shows is narrower: *even when CNV signal is available to
   ichorCNA*, fragment-ratio features alone do not reach ichorCNA's low-fraction accuracy.
