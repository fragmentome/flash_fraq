# Admixture design — coverage feasibility and the mixture plan (2026-08-31)

> **⚠️ Superseded 2026-08-31 — the background pool changed after this was written.**
>
> Everything below is an accurate record of the plan that was generated and of why each choice was
> made *given a Snyder + Cristiano background pool*. That premise no longer holds: the admixture
> background is now **HCC1395BL**, the SEQC2 matched normal for the tumor arm, and Snyder is demoted
> to a held-out cross-center evaluation set.
>
> The reason is the one this document already half-identifies at "achievable on the background side
> and *not* on the tumor side". Making background *identity* orthogonal to `f` leaves the primary
> confound untouched: with a single tumor source, `f` **is** a batch mixing ratio, because every
> tumor read comes from SEQC2 and every background read does not. Check 1 recovers center at accuracy
> 1.000, so `f` is readable straight off the batch signature. Mixing tumor against its own matched
> normal holds batch composition constant across the ladder and removes that by construction.
>
> Concretely invalidated below: the `C_target` = 0.9x choice and the whole feasibility table (the
> healthy side no longer binds at ~34–47x), the 484-mixture plan over 475 backgrounds, the
> 383/83/81 subject-level split, and `max_source_proportion_deviation` as a meaningful statistic.
>
> Still open: matched-pair mixing leaves **one** background subject, so subject-level splitting is no
> longer possible at all. See
> [MAT-13](todo/todo/admixture-and-empirical-tumor-fraction.md#superseding-decision-2026-08-31--background-demoted-to-the-matched-normal)
> for the two-regime question that has to be settled before the plan is regenerated.

## Reopened decision (2026-09-01): two-regime design, matched-pair regime carries no split

Settled: the plan is built from **two regimes**, evaluated separately and never pooled into
one undifferentiated training table.

- **`matched_pair`** — HCC1395 x HCC1395BL mixtures across the fraction ladder. Tumor and
  background identity are fixed, so this regime gives a clean, unconfounded `f` ladder with
  no batch-composition drift — the property v1's orthogonality constraint was trying to buy
  with a multi-subject pool, achieved here for free because there is only one pool member.
- **`real_background`** — mixtures against the retained multi-subject background pool
  (Snyder demoted to held-out cross-center evaluation per MAT-13; Cristiano retained), kept
  for biological diversity in background composition.

`admixture_plan.csv` carries a `regime` column (`matched_pair` | `real_background`) as a hard
partition key, not incidental metadata — downstream consumers (MAT-17's feature-matrix
assembler) must treat it as a partition and must not merge rows across it.

**Splitting.** `matched_pair` has exactly one background subject (HCC1395BL), so a
subject-level split is impossible by construction, and no non-subject-level substitute
(splitting by rung, by seed, ...) actually prevents leakage — the same tumor/normal pair
appears across every split regardless of how the rows are sliced. The regime is therefore
declared **unsplit**: it carries no `split_group` assignment. Only `real_background` is
split, using the existing source-subject-level scheme described below.

**`C_target` re-chosen: 1.0x, uniformly across both regimes.** The v1 choice of 0.9x was
driven entirely by the healthy-side binding constraint
`(1-f)*C_target/C_healthy <= 1`, capped by Snyder's shallowest libraries (0.75x floor). That
constraint no longer holds anywhere in the new design:

- `matched_pair`'s background is HCC1395BL at a declared ~34-47x — MAT-13 already notes this
  makes the healthy side stop binding outright.
- `real_background`'s pool loses Snyder (demoted to held-out cross-center evaluation, see
  below), leaving **Cristiano alone** as the retention-limiting source. At 1.0x — the top of
  phase0 §1.1's 0.1x-1x ULP-WGS range — Cristiano retains 498/503 (99.0%), well above the
  90% floor, so even this regime's one remaining check does not actually bind at the ceiling.
- The tumor side (HCC1395 placeholder 35x) still never binds: downsample fraction tops out
  at ~0.0286 at 1.0x, same "err toward infeasible, not unbuildable" placeholder logic as v1.

With no candidate below the ceiling failing any retention check, the v1 rule — deepest
passing candidate — trivially selects the ceiling itself. A single shared `C_target` (rather
than one value per regime) keeps the two regimes built to comparable depth, which matters
since they are evaluated separately but need to stay comparable to each other.

Consequence for the config: `background_source_groups` for the `real_background` pool must
narrow to `cristiano_delfi` only — Snyder's demotion has to be enforced there, not just
stated in the docs, since the retention computation above assumes Snyder is out of the pool.

**Regenerated 2026-09-01** with the decisions above, over the real 547-row
`native_coverage.csv` and placeholder HCC1395/HCC1395BL coverage (35x / 34x):

| | |
| --- | --- |
| `C_target` chosen | **1.0x**, as predicted — nothing binds below the ceiling |
| Cristiano retained at `f=0`, `C_target=1.0x` | 498 / 503 (99.0%) |
| `matched_pair` mixtures | 11 (one per rung; `matched_pair_replicates_per_fraction: 1`) |
| `real_background` mixtures | 484 (44 x 11 rungs, 483 distinct Cristiano libraries) |
| `real_background` splits | 341 / 77 / 66 mixtures (train/validation/test), all Cristiano |
| `split_assignment.csv` | 504 subjects: 503 Cristiano + HCC1395 (`SPLIT_GROUP_ANY`); HCC1395BL absent — `matched_pair` is unsplit |
| `max_source_proportion_deviation` | 0.0 (trivially, single-source pool) |

The phase2 §7 step-1 design: which mixtures to build, from which background libraries, at
what depth. Target `admixture_design_complete`; module reference in
[pipeline-modules.md](pipeline-modules.md) §8.

```
snakemake --workflow-profile workflow/profiles/local admixture_design_complete
```

Run it from a real `conda activate snakemake` shell, same caveat as
[cohort-qc-baseline.md](cohort-qc-baseline.md). Whole run: ~6 s.

This is the design, not the build. It ends at `admixture_plan.csv`; MAT-13 consumes that
file and does the actual downsample/merge/dedup. Nothing here needs HCC1395 — the tumor
source enters as one sample id and one coverage number, both declared in
`workflow/config/admixture_design.yaml`.

### 2026-09-17 — re-run with real HCC1395/HCC1395BL coverage, plan structure unchanged

MAT-23's blocker cleared: `hcc1395_merged_coverage_report` measured
`HCC1395_WGS_NV_3` (tumor) at **126.59x** and `HCC1395BL_WGS_NV_3` (matched_pair_background) at
**95.91x**, replacing the 35.0x/34.0x placeholders above. Re-ran `admixture_design_complete`
with no other config change.

Both real numbers landed *above* the placeholders (deliberately chosen low, see the `C_target`
discussion above), so the tumor side still never binds and the plan's structure is identical:
**495 mixtures** (11 `matched_pair` + 484 `real_background`), `C_target` still 1.0x, splits still
341/77/66. Only the `matched_pair` rows' downsample fractions changed (both now tighter, since
the real coverage is 2.5-3x the placeholder). Nothing here needed a design decision — this was
the "re-run, not new design work" the reopened Status line anticipated.

## Outputs

All under `/mnt/8tb_hdd/flash_fraq_derived/cohort/`:

| File | What it is |
| --- | --- |
| `native_coverage.csv` | 547 rows, an `input_manifest` CSV with `native_coverage` populated |
| `admixture_feasibility.csv` | 60,170 rows = 547 samples x 11 rungs x 10 candidate depths |
| `admixture_plan.csv` | 484 planned mixtures, one row per mixed sample |
| `split_assignment.csv` | 548 rows, exactly `subject_id,split_group` |
| `admixture_design_summary.json` | the two constraints below, as numbers |

## Native coverage, harvested

Genome-wide mean depth from the `total` row of each sample's `mosdepth.summary.txt`, gated
on `localfeat_sample_complete.flag`. Quantiles are nearest-rank, so every number below is a
depth some library actually has.

| Cohort | n | min | p25 | median | p75 | max | <1.0x |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Snyder | 44 | 0.75x | 0.99x | 1.18x | 1.36x | 1.58x | 12 |
| Cristiano | 503 | 0.67x | 1.99x | 2.26x | 2.63x | 13.0x | 5 |

## The `C_target` choice: **0.9x**

The healthy-side downsample fraction is `(1-f)*C_target/C_healthy`, largest at `f=0`, so
the lowest rung of the ladder is the binding case and `C_target` is capped by the
lowest-coverage background the design wants to keep. The rule: take the **deepest**
candidate in phase0 §1.1's 0.1x-1x ULP-WGS range that keeps >= 90% of *every* source at
`f=0` (`min_background_retention_per_source`, with a 20-sample floor).

| `C_target` at `f=0` | Snyder retained | Cristiano retained |
| --- | --- | --- |
| 0.9x (**chosen**) | **43 / 44 (97.7%)** | **500 / 503 (99.4%)** |
| 1.0x | 32 / 44 (72.7%) | 498 / 503 (99.0%) |

1.0x fails on Snyder: 12 of its 44 samples sit below 1x, so the deepest rung of the
ichorCNA range would silently discard a quarter of the smaller cohort at exactly the low-`f`
end where the regression label matters most. 0.9x costs one Snyder sample (the 0.75x floor)
and three Cristiano, and buys 10% more reads per mixture than the next step down.

The tumor side never binds: at a declared 35x native coverage its downsample fraction peaks
at 0.0129 across the whole ladder. The placeholder is set to the *low* end of phase2 §7's
35-47x range on purpose — a lower declared coverage inflates the tumor-side fraction, so the
error, when MAT-23 lands, is in the direction of calling a mixture infeasible rather than
planning one that cannot be built.

## The fraction ladder

`[0.0, 0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.3, 0.5]`, declared in
`workflow/config/admixture_design.yaml`, not in code. It includes `f=0` because §6 check 1
is specified there, and is dense below 0.1 because ichorCNA's reported sensitivity floor is
~3% (phase0 §1.1) — that is the region this model exists to beat, so it is where the label
needs resolution. 44 mixtures per rung x 11 rungs = **484 mixtures**, drawing on 475
distinct background libraries.

## Source membership is orthogonal to `f` — by construction

[cohort-qc-baseline.md](cohort-qc-baseline.md)'s check 1 separates Snyder from Cristiano at
accuracy **1.000** against a 0.890 baseline. Source identity is unambiguously legible in the
fragmentomics features, so if the background pool's composition drifted with `f`, a model
could hit the regression target by reading batch identity and never learn tumor signal.

The naive design does exactly that. Per-rung feasible pools *grow* with `f` (the healthy
fraction shrinks), and they grow unevenly — Snyder's 0.75x floor drops out long before
Cristiano's tail does — so pool composition would itself encode `f`. Instead every rung
draws from **one pool: the 543 samples feasible at every rung**, allocated across
(source, split) cells in fixed proportions by largest remainder. Every rung is therefore
90.9% Cristiano / 9.1% Snyder, and `max_source_proportion_deviation` in the summary is
**0.0**, not "small".

Round-robin within each cell, not sampling with replacement: a background library used at
many rungs while another was used at none would put library identity back into correlation
with `f`. As planned, 475 libraries cover 484 mixtures.

## Splits

Assigned at source-subject level *before* mixing (phase0 §5.2), stratified within each
source so both cohorts reach all three splits — a single-source split would turn a
validation score into a cross-cohort generalisation measurement. 383 train / 83 validation
/ 81 test subjects, giving 330 / 77 / 77 mixtures, each split covering all 11 rungs.
Assignment runs over every background subject in the manifest rather than over the eligible
pool, so changing `C_target` cannot move a subject between splits.

**One caveat, and it is structural.** phase0 §5.2 wants no source subject in two splits.
That is achievable on the background side and *not* on the tumor side: v1 has exactly one
tumor source, and it is in every mixture. Binding HCC1395 to one split would collapse the
whole matrix into that split. It is therefore assigned the sentinel
`feature_matrix.SPLIT_GROUP_ANY` (`"all"`), and a mixed sample takes its split from
whichever side is split-bound — in practice the background subject.
`FeatureMatrixAssembler.assign_split_groups` resolves the sentinel and still raises when
two *split-bound* subjects disagree, so the leakage check survives for every subject it can
apply to. What §5.2 actually buys here is "no background subject spans two splits"; the
tumor side's guarantee has to come from §5.2's second bullet — holding out an entire
source/platform — once a second tumor source exists.

## What "healthy background" means here

`background_source_groups: []` in the config admits both cohorts, and most of those
libraries are cancer-patient plasma: 39 of the 44 Snyder samples are the IC0x lung-cancer
series (only BH01/IA0x/IH0x are healthy donors, and just IA01-04 and IH03 are on disk), and
the Cristiano subset is pancreatic cancer. They are backgrounds in the *mixing* sense — the
library the tumor reads are spiked into — not verified tumor-free.

This is a deliberate trade. Restricting to the five genuinely-healthy Snyder libraries would
leave one source, and a single-source background pool makes the orthogonality constraint
above vacuous while handing the model a perfect batch shortcut. Residual ctDNA in a
background instead biases the *nominal* `f` low, which is precisely the error phase2 §7.4
already handles by making the SNV-VAF empirical fraction the regression label and
`nominal_tumor_fraction` mere metadata. Narrowing the pool is one config key if that trade
is ever made the other way.
