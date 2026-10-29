# Results-corruption audit: admixture → feature matrix → decoy QC → TF model → eval

Date: 2026-09-25
Scope: bugs and design defects that change **numbers** in the 495-mixture results
(`/mnt/8tb_hdd/flash_fraq_derived/admixture_mixtures/`, `models/tf_model/`). Code style
and consistency issues are out of scope.

This report deliberately leaves out the mpileup multi-sample-column AD bug
(`compute_empirical_tumor_fraction` read only `fields[9]`). That bug is being fixed
separately, with `relabel_admixture_empirical_tf.py` and its driver.

All numbers below come from the on-disk artifacts as they exist today. Spearman ρ is
measured against `nominal_tumor_fraction`, because the empirical label is currently broken.

## Recompute scope at a glance

| # | Finding | Severity | Recompute needed beyond the relabel |
|---|---------|----------|-------------------------------------|
| 1 | Tumor side is sheared gDNA, never reshaped | **Critical** | All 495 mixture BAMs → fragments → GC → ichorCNA → matrix → everything |
| 2 | `mtdna_fraction_log10` is a tumor-read counter used as a feature | **High** | Matrix column drop + retrain/eval only |
| 3 | Cristiano backgrounds are not verified healthy (the config says `tumor`) | **High** | Row filter if a phenotype table exists; else new backgrounds |
| 4 | Decoy QC gate is structurally unable to detect the real confound | **High** | Redesign the check (cheap to run) |
| 5 | matched_pair model is a constant predictor with no evaluation | **High** | Retrain/eval only (or report it as N/A) |
| 6 | Detection floor null pool `truth <= 0.0` stays empty after the label fix | **Medium** | Eval only |
| 7 | Label estimator: unweighted mean of per-site ratios over low-VAF sites | **Medium** | Relabel only (mpileup already rerun by the fix) |
| 8 | ichorCNA baseline run with the high-TF normal prior | **Medium** | ichorCNA rerun (495×) → eval |
| 9 | Minor: samtools `-s` 6-decimal rounding, chrX in GC/arm features, no MAPQ filter | Low | None needed for v1 |

Checked and found clean (no recompute needed): see the end of this report.

---

## 1. Tumor reads are sheared cell-line gDNA; fragment features encode "fraction of HCC1395 library", not ctDNA biology — CRITICAL

**Where**
- `workflow/config/admixture_design.yaml`: the tumor `bam_path` is the raw
  `HCC1395_WGS_NV_3.sorted.cram`.
- `scripts/drivers/admixture_mixture_fragments_driver.sh` runs fragment extraction
  directly on `dedup.bam`.
- `docs/todo/completed/fragment-size-reshaping.md` says reshaping is "not yet wired into
  admixture_design.yaml's tumor source".

**Effect on results.** Every real_background mixture is sonicated HiSeq4000 gDNA (tumor)
plus nuclease-digested cfDNA (background). As f rises, the size distribution moves toward
the gDNA shear profile, whatever the tumor biology is. Measured across the 484
real_background rows:

| feature | Spearman ρ vs nominal f |
|---|---|
| `fragment_size_entropy` | **0.948** (6.31 at f=0 → 8.02 at f=0.5) |
| `fragment_size_median` | 0.626 (166.95 → 179.16 bp) |
| `short_long_ratio_window_*`, `gc_corrected_short_coverage_window_*` | same mechanism (short/long bins are cfDNA-defined) |

- Every fragmentomics feature (DELFI short/long ratios, short coverage, entropy) measures
  library prep. It does not measure tumor-derived fragmentation.
- Real ctDNA is cfDNA-sized and somewhat shorter. It does **not** broaden the distribution
  toward the gDNA profile. The model's apparent skill therefore will not transfer to real
  plasma, and any "fragmentomics adds signal beyond CNV" claim is confounded.
- In matched_pair, both sides are gDNA, so it has no cfDNA fragmentation signal at all.
  Only CNV/coverage features carry f there.

**Fix / recompute.** Wire `fragment_reshaping` into the tumor source, or restrict the
reported features to CNV-only. Reshaping changes which tumor reads exist, so it forces a
rerun of all mixture BAMs and every stage after them. **The relabel alone cannot rescue
this.** If a full recompute is unaffordable, the minimum honest option is to rerun only
the model/eval on a CNV-only feature subset (`cnv_logr_bin_*`, `arm_zscore_*`,
`gc_corrected_total_coverage_window_*`) and state that fragment features are
uninterpretable in this cohort. Total coverage is also somewhat affected through GC/mappability
differences between the libraries.

## 2. `mtdna_fraction_log10` is a direct f proxy and is the #2 feature — HIGH

**Where:** `admixture.py:compute_mtdna_fraction_log10` writes the value into
`sample_metadata.csv`. It is not excluded in `feature_matrix.py`, `tf_model.py` or
`decoy_qc.py`.

**Evidence**
- HCC1395 gDNA has mtDNA fraction ≈ 10^-2.6 (the matched_pair rows are all −2.61 to −2.70).
  Cristiano cfDNA is ≈ 10^-4.9. Each mixture's mtDNA fraction is therefore just the
  share of HCC1395 reads.
- In real_background, ρ = **0.960** vs nominal f, with a monotone ladder from −4.94 at f=0
  to −3.24 at f=0.5. Within-rung SD shrinks from 0.26 to 0.007.
- It is in both saved models' `feature_columns`. In the real_background model it ranks
  **#2 by split importance** (84, behind one CNV bin at 137).

**Effect.** The model is partly reading "how much cell-line library is in the tube". A
real patient sample has nothing equivalent (tumor cfDNA does not carry 100× mtDNA), so
this inflates every real_background metric. It will keep inflating them after the label fix.

**Fix / recompute.** Add `mtdna_fraction_log10` to the non-feature columns in
`feature_matrix.py`, `tf_model.py` and `decoy_qc.py`, then retrain and re-evaluate. This
is cheap because the column is already in the matrix.

## 3. real_background "healthy" donors are not verified healthy — HIGH (to confirm)

**Where:** `admixture_design.load_background` filters only on
`source_group in ["cristiano_delfi"]`. In `cohort/native_coverage.csv`, all 503 Cristiano
rows carry `tumor_status = tumor`. That is the dataset-level vocabulary from
`admixture_design.yaml`, not a per-sample phenotype. No step filters to non-cancer donors.

**Evidence that it matters.** The Cristiano DELFI release is a mix of cancer patients and
controls. Among the 44 f=0 mixtures, ichorCNA reports TF > 0.1 for 6 of them, with a
maximum of 0.346. Across the native cohort, Cristiano's maximum ichorCNA TF is 0.49 vs 0.25
for Snyder healthy. ichorCNA at 1× is noisy, so this is suggestive but not proof.

**Effect.** For any cancer-patient background, the true ctDNA fraction is
`f + (1−f)·TF_bg`. The empirical label counts only HCC1395 somatic SNVs, so it misses
TF_bg, while the CNV/fragment features see it. The result is label noise at every rung,
worst at f=0: the "negative controls" contaminate the null pool, the low-fraction
metrics and the detection floor.

**Fix / recompute.** Join the DELFI per-sample phenotype (cancer/non-cancer) onto the
background pool.
- If the backgrounds are mixed, the cheapest option is to **drop the cancer-background
  rows** from the matrix and retrain. That leaves ~half the rows and changes no BAMs.
- A redesign is needed only if too few rows remain.

## 4. Decoy QC "passed" measures nothing relevant — HIGH

**Where:** `decoy_qc.py`, run on the pooled `feature_matrix.csv`. Report:
`admixture_mixtures/decoy_qc/decoy_qc_report.json`.

- **Source label = regime.** The source label is `platform_tumor|platform_healthy`, which
  takes exactly two values here: `HiSeq4000|HiSeq2500` (484 real_background) and
  `HiSeq4000|HiSeq4000` (11 matched_pair). The classifier therefore predicts
  matched_pair vs real_background. It does not detect a batch effect inside the regime
  being modeled.
- **Check 1 never ran.** It reports `n_samples=45, accuracy=1.0, baseline=1.0` because the
  matched_pair class had fewer than 2 rows in the f≈0 slice. That is the early-return branch.
- **Check 2 is degenerate.** Its correlation is 0.0082 vs 0.0082 for the plain source
  indicator (`degenerate_with_source_indicator: true`), so it is just
  corr(is_matched_pair, label). With the ladder balanced across regimes, that is ≈0 by
  construction.
- **The real confound is invisible to it.** The actual confound is tumor library (gDNA,
  HiSeq4000) vs background library (cfDNA, HiSeq2500) *within each mixture*, and it scales
  with f (finding 1). That confound is not a per-row "source", so no per-row source
  classifier can see it. The gate would pass even with a fully confounded model.

**Effect.** Nothing downstream should be described as "batch-QC passed". This will not
change after the relabel.

**Fix.** Run it per regime, and add a check that targets the confound. For example, test
how well non-CNV features alone predict f (ρ ≈ 0.95 for entropy and mtDNA today, which
should fail the gate), or run a tumor-library-identity probe.

## 5. matched_pair: no split, constant model, no evaluation — HIGH

**Evidence**
- `split_group` is **NaN** for all 11 matched_pair rows in `feature_matrix.csv`.
  `split_assignment.csv` has no HCC1395BL subject, so `assign_split_groups` maps it to
  NaN. (The design intended `MATCHED_PAIR_SPLIT_GROUP`.)
- `tf_model._train_rows` finds no `train` rows and falls back to all 11. With
  `min_child_samples=10`, LightGBM cannot make a split. In `tf_model_report.json`,
  **every feature importance is 0**, so the model predicts the training mean for any input.
- `tf_eval` produced only `tf_eval_report_real_background.json`. Nothing evaluates matched_pair.

**Effect.** Any matched_pair model number, or a "matched_pair vs real_background"
comparison, is meaningless. The relabel will not change this.

**Fix.** Report matched_pair as a descriptive calibration curve (features/ichorCNA vs
label across 11 rungs), or use leave-one-rung-out with a model sized for n=11. Retraining
alone is cheap.

## 6. Detection floor null pool is defined on the empirical label — MEDIUM

**Where:** `tf_eval.compute_detection_floor`, using `null_mask = truth <= 0.0` with
`truth = eval_rows[empirical_tumor_fraction]`.

- Today no row at any rung has empirical TF ≤ 0, so the pool is empty and
  `detection_floor: null`.
- After the AD-sum fix, f=0 mixtures will still get small **positive** empirical values
  (sequencing errors that match the alt at ~39k sites), so the null pool stays empty and
  the floor stays `null`.
- The bins are also defined on the empirical label, so rows drift across bin edges from
  label noise.

**Fix.** Define the null pool as `nominal_tumor_fraction == 0`, or `truth <= ε` with a
realistic ε. Eval-only recompute.

## 7. Label estimator: unweighted mean of per-site ratios — MEDIUM (survives the AD fix)

**Where:** `admixture.py:compute_empirical_tumor_fraction`. After the column fix, it still
returns `mean(observed_af / vaf_tumor_pure)` over every covered site.

- **Variance.** The sites set has 39,141 SNVs. bwaTVAF has p1 = 0.04, p5 = 0.058 and
  median 0.335, and 4,977 sites have VAF < 0.1. At 1× most covered sites have depth 1–2.
  A single alt read at a VAF-0.04 site adds 25/depth to that site's ratio, and every site
  gets equal weight whatever its depth. This greatly inflates label noise at low f, which
  is the regime that matters.
- **Bias.** bwaTVAF is estimated from the same HCC1395 data and is selected by a calling
  threshold. By Jensen's inequality, E[1/VAF_hat] ≠ 1/VAF, and the error is largest for
  low-VAF sites.
- **Spec mismatch.** The ticket specifies a regression of observed on f·VAF (sites) and
  mean(obs)/mean(VAF) (clonal).
- **Suggested estimator.** A depth-weighted ratio of sums, Σalt / Σ(depth·VAF), is
  unbiased under the same model and has far lower variance.
- `empirical_tumor_fraction` = sites is the regression label, so this noise goes straight
  into training and eval.

**Fix / recompute.** Change the estimator inside the relabel that is already being run.
The same mpileup output serves both, so the extra cost is ≈0 if it is done now. Better
still, make the relabel script keep the per-site (alt, depth, VAF) counts so any
estimator can be recomputed without rerunning mpileup.

## 8. ichorCNA baseline uses the high-TF normal prior — MEDIUM

**Where:** `admixture_mixture_ichor_driver.sh` sets `ICHOR_NORMAL="c(0.5,...,0.9)"`,
`ICHOR_SC_STATES="c(1,3)"` and `ICHOR_PLOIDY="c(2,3)"`.

- For low-TF ULP-WGS, ichorCNA's own recommendation is `normal=c(0.95,0.99,0.995,0.999)`
  with no subclonal states. The high-TF prior plus subclonal states biases estimates at
  f ≤ 0.05.
- Observed means at f = 0, 0.01, 0.02, 0.03 and 0.05 are 0.060, 0.048, 0.050, 0.065 and
  0.054: flat and non-monotone.
- At f=0.5 it reads 0.29 on average.

**Effect.** The "model beats ichorCNA" comparison, especially the low-fraction metrics,
uses a mis-configured baseline. It does not corrupt the model itself, since ichorCNA TF
is excluded from the primary features. It does corrupt the headline comparison.

**Fix / recompute.** Rerun ichorCNA with the low-TF settings. That means 495 runs from the
existing `.wig` files, with no BAM work, followed by eval.

## 9. Low-impact items (no recompute recommended for v1)

- **samtools `-s` precision.** `subsample_bam` formats `seed + fraction` to 6 decimals.
  For the smallest tumor fractions (e.g. matched_pair f=0.01 gives 7.9e-5), that is up to
  ~0.6% relative error in the draw fraction. The empirical label absorbs it once fixed.
- **Downsample-depth check.** `verify_downsample_depth`'s absolute tolerance of 0.25× at
  C=1× would not catch a 25% depth error. Nothing indicates that this happened.
- **chrX.** chrX is included in the GC LOESS fit and the arm z-scores. HCC1395 is female,
  as are some backgrounds, so this adds a small sex-dependent shift. The features are
  within-sample normalized, so the effect is minor.
- **No MAPQ filter** in fragment extraction. It adds low-MAPQ noise equally to all rows
  and does not bias with f.

---

## Checked and clean

Nothing here needs recompute.

- **Label columns not leaked into the saved models.** `feature_columns` in both joblibs
  has 4,355 entries. The only non-bin features are `mtdna_fraction_log10` and
  `fragment_size_{median,mode,entropy}`. `empirical_tumor_fraction_{sites,clonal}` are
  absent.
- **Subject-level split.** 483 unique real_background donors, and no donor spans two
  splits. `EGAF00002727156` appears twice (f=0 and f=0.5), both times in validation.
- **Validation rows are held out.** Validation is not used for early stopping or tuning,
  so eval on validation+test is held-out.
- **Tumor batch sampler.** The crc32 + Weyl-phase rule keeps mates together and is
  independent across rungs.
- **Merge → fixmate → markdup chain.** Code review found no read loss or duplication that
  depends on f. Read counts were not reconciled per mixture.

## What this means for "avoid the recompute"

- **Relabel-only path** (already planned): relabel → matrix → decoy QC → model → eval.
  Fold findings **2, 5, 6 and 7** into it; each is a small code change.
- **Needs a decision before the relabel:** finding **3** (filter cancer backgrounds) and
  finding **4** (redesign decoy QC). Both are cheap to run, but they change which results
  are valid.
- **Can't be fixed without a real recompute:**
  - Finding **1**. The fragmentomics results cannot be interpreted without reshaping the
    tumor side and regenerating all mixtures. The cheapest defensible fallback is a
    CNV-only model on the existing features.
  - Finding **8** needs an ichorCNA-only rerun, which is moderate in cost.
