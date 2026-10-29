# `src/flash_fraq` module reference and callstack

Reference doc for the preprocessing pipeline modules implemented so far, corresponding to
tickets MAT-13 through MAT-18 (see `docs/todo/`) and `docs/phase2-preprocessing-plan.md`. There is no
Snakefile yet (see CLAUDE.md — pre-implementation status), so the pipeline order below is
inferred from each module's documented inputs/outputs, not from a wired-up DAG.

## Pipeline order (data-dependency order)

```
0. admixture_design.py    (design)  mosdepth summaries -> native_coverage.csv,
                                    admixture_feasibility.csv, admixture_plan.csv,
                                    split_assignment.csv  (decides WHAT step 1 builds)
1. admixture.py           (MAT-13) tumor_bam + healthy_bam -> mixed, deduped BAM
                                    + sample_metadata.csv (empirical_tumor_fraction label)
2. fragment_extraction.py (MAT-14) mixed BAM -> fragment_stats.csv, window_fragment_counts.tsv
3. gc_correction.py       (MAT-15) window_fragment_counts.tsv -> window_ratios.csv, arm_zscores.csv
4. cnv_binning.py         (MAT-16) mixed BAM -> cnv_bins.csv, ichor_estimates.csv (via HMMcopy + ichorCNA)
5. feature_matrix.py      (MAT-17) joins outputs of 1-4 (+ split_assignment.csv) -> one-row-per-sample matrix
6. decoy_qc.py            (MAT-18) feature matrix -> decoy_qc_report.json (batch-effect gate; must pass
                                    before treating any regression result as informative)
```

Steps 1 and 5 need HCC1395 as the tumor source, so neither can run yet. Step 0 does not —
the tumor enters it as a declared sample id and coverage number — so it has already run
over the real 547 samples; see [admixture-design.md](admixture-design.md). The
**pre-admixture branch** below runs today on the real Snyder + Cristiano samples and
short-circuits 1 -> 5, feeding the same step 6 gate:

```
2/3/4 per-sample outputs (already written for 547 samples by workflow/rules/local_features.smk)
   -> cohort_matrix.py  one row per REAL sample -> cohort_feature_matrix.csv
   -> decoy_qc.py       DecoyClassifierQcConfig.for_cohort -> decoy_qc_report.json
```

`__init__.py` is currently empty — no package-level exports yet.

---

## 1. `admixture.py` — in silico tumor/healthy admixture builder

Builds one mixed sample by downsampling a tumor BAM and a healthy BAM to a target tumor
fraction `f`, merging, deduping **once on the merged BAM** (never per-source — see
CLAUDE.md), and deriving the **empirical** SNV-VAF tumor fraction (the regression label,
not the nominal mixing ratio).

**Key types**
- `SomaticSnv` — one known HCC1395 somatic SNV site (`chrom`, `pos`, `ref`, `alt`, `vaf_tumor_pure`).
- `AdmixtureConfig` — frozen config; built via `AdmixtureConfig.from_sources(...)`.
- `AdmixtureBuilder` — orchestrator, holds `self.config`.

**Callstack from `AdmixtureBuilder.run()`:**
```
run()
├── compute_downsample_fractions()        # coverage-normalized fraction for tumor & healthy
├── subsample_bam(tumor_bam, ...)         # samtools view -s
├── subsample_bam(healthy_bam, ...)       # samtools view -s
├── verify_downsample_depth(...)  x2      # samtools flagstat, sanity-check retained fraction
│   └── _flagstat_total_reads()
├── merge_and_dedup([tumor_sub, healthy_sub])
│   └── samtools merge -> sort -n -> fixmate -m -> sort -> markdup -> index
├── load_snvs()                           # reads snv_vaf_path TSV -> list[SomaticSnv]
├── compute_empirical_tumor_fraction(dedup_bam, snvs)
│   └── bcftools mpileup -R sites.bed -a AD  (subprocess), then per-site observed_AF / vaf_tumor_pure, meaned
├── compute_mtdna_fraction_log10(dedup_bam)
│   └── samtools idxstats                 # log10(mtDNA reads / total mapped reads)
└── write_sample_metadata(...)            # -> sample_metadata.csv (one row: label + source/platform/subject IDs)
```

External tools invoked as subprocesses: `samtools` (view, merge, sort, fixmate, markdup, index,
flagstat, idxstats), `bcftools mpileup`.

---

## 2. `fragment_extraction.py` — cfDNA fragment extraction and window intersection

Extracts fragment midpoints from a (mixed, deduped) BAM using `pysam`, intersects them
against DELFI 5 Mb and ichorCNA 1 Mb window BED files via `bedtools`, and summarizes
per-window short/long fragment counts plus genome-wide fragment-size stats.

**Key types**
- `FragmentExtractionConfig` — via `from_sources(...)`.
- `FragmentExtractor` — orchestrator.

Filtering constants: `REQUIRE_FLAGS` (proper pair), `EXCLUDE_FLAGS` (dup/secondary/supplementary),
`READ1_FLAG` (dedupe pair to one row per fragment). DELFI short/long ranges:
`SHORT_FRAGMENT_RANGE = (100, 150)`, `LONG_FRAGMENT_RANGE = (151, 220)`.

**Callstack from `FragmentExtractor.run_full()`** (the entry point that produces both
MAT-14 output files; `run()` is a subset used for just the intersected BEDs):
```
run_full()
├── extract_midpoints()                   # pysam.AlignmentFile.fetch, _passes_filters() per read
│                                          # -> fragment_midpoints.bed (chrom, midpoint, midpoint+1, tlen)
├── intersect_windows(midpoints_bed)
│   ├── _run_bedtools_intersect(midpoints_bed, delfi_windows_bed, ...)  # -> fragments_delfi_5mb.bed
│   └── _run_bedtools_intersect(midpoints_bed, ichor_bins_bed, ...)    # -> fragments_ichor_1mb.bed
├── write_window_fragment_counts(delfi_out)
│   └── aggregate_window_counts(delfi_intersected)
│       ├── _parse_window_index(window_name)
│       └── _window_gc_content(fasta, chrom, start, end)   # pysam.FastaFile.fetch, per window
│                                                            # -> window_fragment_counts.tsv
└── write_fragment_stats(midpoints_bed)
    └── compute_fragment_stats(midpoints_bed)
        └── _shannon_entropy(tlens)                         # -> fragment_stats.csv
                                                              #    (fragment_size_median/mode/entropy)
```

External tools invoked: `bedtools intersect` (subprocess); BAM/FASTA reads via `pysam`.

---

## 3. `gc_correction.py` — DELFI-style LOESS GC-bias correction

Applies LOESS GC correction to the short/long fragment counts produced by
`fragment_extraction.py`, separately for short and long counts (not directly to their
ratio), then recomputes the short:long ratio from the corrected counts. Also derives
chromosome-arm fragmentation z-scores.

**Key types**
- `GcCorrectionConfig` — via `from_sources(...)`; derives `window_ratios_path` /
  `arm_zscore_path` defaults from `output_windows_path`'s parent if not given.
- `GcCorrector` — orchestrator.
- `HG19_CENTROMERE_MB` — approximate hg19 centromere midpoints (Mb) used to bucket 5 Mb
  windows into `p`/`q` arms; flagged in-file as a coarse approximation to validate later.

**Callstack from `GcCorrector.run_full()`:**
```
run_full()
├── run()
│   └── correct(windows)
│       └── _correct_counts(gc_content, counts)  x2   # statsmodels.api.nonparametric.lowess(frac=span)
│                                                       # short_count and long_count independently
│                                                       # -> short_long_ratio, gc_corrected_*_coverage
│                                                       # (also writes output_windows_path, the full diagnostic table)
├── write_window_ratios(corrected)                     # -> window_ratios.csv
└── write_arm_zscores(corrected)
    └── compute_arm_zscores(corrected)
        └── _arm_id(chrom, start)                       # per-window arm bucket via HG19_CENTROMERE_MB
                                                          # -> arm_zscores.csv: sample_id,arm,zscore
                                                          #    (self-normalized z-score against this
                                                          #    sample's own arm distribution; keyed on the
                                                          #    arm label, since the arm set is ragged)
```

No external subprocess calls in this module — pure pandas/statsmodels.

---

## 4. `cnv_binning.py` — HMMcopy + ichorCNA subprocess wrapper

Thin subprocess wrapper around HMMcopy Suite (`readCounter`, `gcCounter`, `mapCounter`)
and ichorCNA's `runIchorCNA.R`, for 1 Mb copy-number binning. Logic is never reimplemented
in Python — these are GPL-licensed external tools invoked as subprocesses (see module
docstring for the licensing rationale).

**Key types**
- `IchorCnaConfig` — via `from_sources(...)`; `__post_init__` derives all intermediate/output
  paths (`read_wig_path`, `gc_wig_path`, `map_wig_path`, `cna_seg_path`, `params_path`,
  `cnv_bins_output_path`, `ichor_estimate_output_path`) from `output_dir` + `sample_id`.
- `IchorCnaRunner` — orchestrator.
- `PARAMS_TUMOR_FRACTION_PREFIX` — the label ichorCNA writes in its `.params.txt` output.

**Callstack from `IchorCnaRunner.run()`:**
```
run()
├── run_read_counter()   -> _run_to_wig -> _run_subprocess(["readCounter", ...])   # -> *.readcounter.wig
├── run_gc_counter()     -> _run_to_wig -> _run_subprocess(["gcCounter", ...])     # -> *.gc.wig
├── run_map_counter()    -> _run_to_wig -> _run_subprocess(["mapCounter", ...])    # -> *.map.wig
├── run_ichor_cna()      -> _run_subprocess(["Rscript", "runIchorCNA.R", ...])     # -> *.cna.seg, *.params.txt
├── parse_cna_bins()                          # read_csv(cna_seg_path, sep="\t")
├── write_cnv_bins(bins)                      # long-format -> cnv_bins.csv (sample_id, bin_index, cnv_logr)
└── write_ichor_estimate(parse_tumor_fraction_estimate())
    └── parse_tumor_fraction_estimate()       # parses "Tumor Fraction:" line from params.txt
                                               # -> ichor_estimates.csv (ichorCNA_tumor_fraction_estimate,
                                               #    ablation-only feature — see feature_matrix.py)
```

`_run_subprocess` centralizes error handling: wraps `CalledProcessError` into a `RuntimeError`
with the tool name and captured stderr.

---

## 5. `feature_matrix.py` — assembles the Phase 3 output-contract feature matrix

Joins the per-sample outputs of modules 1–4 (MAT-13 through MAT-16) into the single
one-row-per-mixed-sample matrix defined by phase2 §9. Also resolves `split_group`
(assigned at the source-subject level, before mixing, per phase0 §5.2) and separates the
ablation-only ichorCNA column from the primary feature set.

**Key types**
- `FeatureMatrixConfig` — via `from_sources(input_dir, output_path)`; all input paths are
  derived from `input_dir` by convention (`sample_metadata.csv`, `fragment_stats.csv`,
  `arm_zscores.csv`, `window_ratios.csv`, `cnv_bins.csv`, `ichor_estimates.csv`,
  `split_assignment.csv`).
- `FeatureMatrixAssembler` — orchestrator.
- `REGRESSION_LABEL_COLUMN = "empirical_tumor_fraction"`, `ABLATION_ONLY_COLUMNS = ("ichorCNA_tumor_fraction_estimate",)`.

**Callstack from `FeatureMatrixAssembler.assemble()` / `.write()`:**
```
assemble()
├── load_sample_metadata()                    # sample_metadata.csv (MAT-13 output)
├── load_split_assignment()                   # split_assignment.csv (subject_id -> split_group)
├── assign_split_groups(sample_metadata, split_assignment)
│                                              # maps tumor_subject_id & healthy_subject_id to split_group,
│                                              # raises ValueError if they disagree (source-subject leak)
└── merge(...) chain over:
    ├── load_fragment_stats()                 # fragment_stats.csv (MAT-14)
    ├── load_arm_zscores_wide()                # arm_zscores.csv (MAT-15) pivoted long -> wide (arm_zscore_chr1p)
    ├── load_window_ratios_wide()              # window_ratios.csv (MAT-15) pivoted -> ratio/coverage columns
    ├── load_cnv_bins_wide()                   # cnv_bins.csv (MAT-16) pivoted -> cnv_logr_bin_N
    └── load_ichor_estimates()                 # ichor_estimates.csv (MAT-16)
                                                # all left-merged on sample_id

write(matrix)                                  # -> output_path CSV
get_primary_feature_columns(matrix)            # excludes IDs/metadata + REGRESSION_LABEL_COLUMN
                                                # + ABLATION_ONLY_COLUMNS -> primary model feature list
```

---

## 6. `decoy_qc.py` — decoy classifier batch-effect QC gate

Implements the phase2 §6 gate: trains a classifier to predict source/platform from the
feature matrix and checks whether that predicted-source signal correlates with the
tumor-fraction label. A positive correlation signals leaked batch effects rather than
real biological signal, and must be resolved before trusting any training result.

**Key types**
- `DecoyClassifierQcConfig` — via `from_sources(feature_matrix_path, output_dir)`;
  `__post_init__` derives `non_feature_columns` (IDs, metadata, label, platform columns)
  to exclude from the classifier's feature set.
- `FixedFractionCheckResult`, `FractionCorrelationCheckResult`, `DecoyClassifierQcReport` —
  result dataclasses.
- `DecoyClassifierQc` — orchestrator.

**Callstack from `DecoyClassifierQc.run()`:**
```
run()
├── load_feature_matrix()                     # read_csv(feature_matrix_path)
├── run_fixed_fraction_check(feature_matrix)   # Check 1 (context only, does not gate `passed`)
│   ├── get_feature_columns(fixed_subset)
│   ├── get_source_label(fixed_subset)         # platform_tumor + "|" + platform_healthy
│   └── sklearn: train_test_split -> LogisticRegression.fit/.score
│                                               # "can classifier discriminate source at fixed tumor fraction?"
├── run_fraction_correlation_check(feature_matrix)   # Check 2 (this is the actual gate)
│   ├── get_feature_columns(feature_matrix)
│   ├── get_source_label(feature_matrix)
│   ├── skipped_correlation_check(reason, ...)  # one source / no features / smallest source < n_splits
│   ├── sklearn: cross_val_predict(pipeline, cv=StratifiedKFold, method="predict_proba")
│   │                                           # OUT-OF-FOLD per-class probability vs. tumor_fraction label,
│   │                                           # correlation_by_source, flagged if max |corr| > threshold
│   └── correlate_with_fraction(...)            # also the plain corr(source indicator, fraction), so a
│                                               #   degenerate check reports itself as degenerate
└── write_report(report)                       # -> decoy_qc_report.json
                                                 #    report.passed = not correlation_check.flagged
```

External libraries: `scikit-learn` (`LogisticRegression`, `train_test_split`,
`StratifiedKFold`, `cross_val_predict`) — the only module in this set that trains a model
rather than just transforming data.

**Check 2 is scored out-of-fold.** It was originally fitted on the full matrix and scored on
the same rows; over 4354 features and 547 samples the classifier separates its training set
perfectly, the probabilities saturate, and the correlation collapses onto plain
`corr(source indicator, fraction)` (0.02199 reported vs. 0.02198 computed directly — see
[cohort-qc-baseline.md](cohort-qc-baseline.md)). Three things follow, and all three are in
the code:

- `cross_val_predict` under a shuffled `StratifiedKFold(n_splits, random_state)`. The
  `StandardScaler` must stay inside the `Pipeline` `build_classifier()` returns, or it
  refits once on the full matrix and leaks held-out means back into training.
- `n_splits` and `fold_random_state` go into the report, so the number is reproducible from
  the report alone. If the smallest source has fewer samples than `n_splits` the check is
  **skipped with a reason** rather than raising or quietly using fewer folds, and a skipped
  check writes `max_abs_correlation: null` — never a 0.0 that reads as measured evidence.
  The one-source and no-usable-features cases skip the same way.
- `source_indicator_correlation_by_source` is reported alongside, and
  `degenerate_with_source_indicator` is true when the two agree within
  `degeneracy_tolerance`. Out-of-fold probabilities can still saturate when the sources are
  genuinely separable (they are, here), so the report says so instead of leaving the next
  reader to rediscover it.

---

## 7. `cohort_matrix.py` — pre-admixture cohort feature matrix

The §9 matrix `feature_matrix.py` builds is one row per *mixed admixture sample*, keyed on
`empirical_tumor_fraction`. Admixture needs HCC1395 as the tumor source, which is not
aligned, so that matrix cannot exist yet. This module builds the closest thing that can:
one row per **real, unmixed sample**, pivoted and joined from the per-sample MAT-14/15/16
tables that `workflow/rules/local_features.smk` already wrote under
`<derived_root>/<dataset>/<sample_id>/`.

It shares `feature_matrix.py`'s pivot helpers (`pivot_window_ratios`, `pivot_cnv_bins`,
`pivot_arm_zscores`, over the generic `pivot_long_to_wide`), so every column that survives
into the real §9 matrix is named identically here.

**Deliberate differences from §9**
- No `empirical_tumor_fraction` — there is no mixture and no SNV-VAF ground truth.
  `ichorCNA_tumor_fraction_estimate` stands in wherever a fraction is needed.
- `source_dataset` + `platform` (one library per row) replace the
  `platform_tumor`/`platform_healthy` pair.
- No `split_group`: splits are assigned at the source-subject level before mixing, which
  has not happened.
- **Arm z-scores are keyed by arm label (`arm_zscore_chr1p`), not by the row's position in
  arm_zscores.csv.** A sample only has windows on the arms its coverage reaches — across
  the real cohort 44, 45 or 46 of the 46 primary-contig arms — so a positional key would
  place different chromosome arms in one column. `arm_zscores.csv` carries the label
  itself; `load_arm_zscores_long` rejects a table written before it did rather than
  falling back to position (back-fill: `scripts/backfill_arm_zscore_labels.py`).

**Key types**
- `CohortDataset(name, platform)` — the platform label is source metadata, not recoverable
  from the derived tables, so it is configured (`workflow/config/cohort_qc.yaml`).
- `CohortSample` — one sample directory and the six tables it must carry.
- `CohortMatrixConfig` — via `from_sources(derived_root, dataset_platforms, output_path)`.
- `CohortFeatureMatrixAssembler` — orchestrator.

**Callstack from `CohortFeatureMatrixAssembler.assemble()`:**
```
assemble()
├── discover_samples()                        # dirs with localfeat_sample_complete.flag AND all six
│                                             #   tables; raises on a sample_id in two datasets
├── load_sample_metadata(samples)             # sample_id, source_dataset, platform, reference_build
├── load_long_table(samples, ...)             # concat per-sample fragment_stats.csv / ichor_estimates.csv
├── load_arm_zscores_long(samples)            # concat arm_zscores.csv; rejects a pre-`arm` table
└── feature_matrix pivots + left merge on sample_id:
    ├── pivot_window_ratios(...)              # short_long_ratio_window_N, gc_corrected_*_window_N
    ├── pivot_cnv_bins(...)                   # cnv_logr_bin_N
    └── pivot_arm_zscores(...)                # arm_zscore_chr1p ...

write(matrix)                                  # -> cohort_feature_matrix.csv
summarise(matrix)                              # n_samples, n_feature_columns,
                                               #   complete_feature_columns (present for EVERY sample),
                                               #   samples_by_source_dataset, tumor-fraction counts
get_feature_columns(matrix)                    # excludes metadata + the ichorCNA estimate (which
                                               #   doubles as the stand-in label -> circular as a feature)
```

The pivot is inherently ragged: per-sample window/arm/bin sets differ (a window whose
GC-corrected ratio is undefined, an ichorCNA-filtered bin), so a missing cell is a genuine
NaN rather than a shifted column. `summarise`'s `complete_feature_columns` is what the
decoy classifier can actually use at `max_missing_fraction: 0.0`.

**Snakemake:** `workflow/rules/cohort_qc.smk`, target `cohort_qc_complete`. Workstation-only,
like `local_features` — ask for it by name:
```
snakemake --workflow-profile workflow/profiles/local cohort_qc_complete
```

---

## Running the §6 gate against the cohort matrix

`DecoyClassifierQcConfig.for_cohort(...)` adapts `decoy_qc.py` to the single-source matrix:
the source label becomes `source_dataset|platform` (`get_source_label` joins any number of
platform columns, not exactly two) and the fraction becomes
`ichorCNA_tumor_fraction_estimate`.

**This is a weaker check than the one Phase 3 needs.** ichorCNA's estimate is derived from
these same samples' coverage rather than from SNV VAFs of a known mixture, so a correlation
here reads as "source predicts ichorCNA's CNV call", not "source predicts true tumor
burden". A pass is necessary, not sufficient; `fail_on_flag` therefore defaults to false in
`workflow/config/cohort_qc.yaml` until the gate runs on `empirical_tumor_fraction`.

Two additions were needed before the gate would run on real (rather than synthetic) data:
- `prepare_features()` — coerces to numeric, drops columns missing in more than
  `max_missing_fraction` of rows (0.0 = keep only cohort-complete columns), median-imputes
  what remains, and drops constant columns. sklearn rejects NaN outright, and the ragged
  pivot guarantees some.
- `build_classifier()` — `StandardScaler` before `LogisticRegression`. Not cosmetic: the §9
  feature block mixes arm z-scores (order 1) with GC-corrected window coverages (order
  1e4), and unscaled L2-penalised logistic regression regularises the small-scale columns
  away and stops converging within `max_iter`.

`run_fraction_correlation_check` also now refuses to run on a constant label column. A
constant label makes every correlation NaN, and `NaN > threshold` is False — the gate would
have reported a clean pass having measured nothing.

---

## 8. `admixture_design.py` — coverage harvest, feasibility table, mixture plan

`admixture.py` (§1) knows how to build one mixed sample once someone tells it which two
libraries to mix, at what `f`, and to what depth. This module is the part that decides
**what to build**, and it needs no tumor reads: HCC1395 enters as one sample id and one
coverage number, declared in `workflow/config/admixture_design.yaml` until MAT-23 lands.

It fills three gaps that nothing else did. `AdmixtureConfig` requires
`tumor_native_coverage` and `healthy_native_coverage` and nothing produced either, even
though `ManifestEntry.native_coverage` was already a field and mosdepth had already run on
all 547 samples. `feature_matrix.assign_split_groups` consumes a `split_assignment`
DataFrame that nothing emitted. And there was no fraction ladder, no target-coverage choice
and no healthy-partner assignment anywhere in the repo.

```
NativeCoverageHarvester(NativeCoverageConfig)
├── harvest()                                  # walk <derived_root>/<dataset>/<sample_id>/,
│   └── parse_mosdepth_summary(path)           #   `total` row col 4; gated on
│                                              #   localfeat_sample_complete.flag
├── write(entries)                             # -> native_coverage.csv, via InputManifestBuilder
└── summarise(entries)                         # nearest-rank per-source coverage distribution

AdmixtureFeasibilityAnalyzer(AdmixtureDesignConfig)
├── load_background()                          # coverage manifest, filtered by
│                                              #   background_source_groups ([] = all)
├── build_table(background)                    # sample x f x C_target, both downsample
│                                              #   fractions + the `<= 1` verdict
├── retention_by_source(table, C)              # how many survive at the BINDING (lowest) rung
└── choose_target_coverage(table)              # deepest C keeping >= N and >= share of every source

MixturePlanner(AdmixtureDesignConfig)
├── assign_split_groups(background)            # -> split_assignment.csv (subject_id,split_group)
├── eligible_background(background, C)         # the ONE pool every rung draws from
├── build_plan(background, C, splits)          # -> admixture_plan.csv
├── validate_plan(plan)                        # asks AdmixtureBuilder itself, not a second copy
└── summarise(plan, splits, C)                 # max_source_proportion_deviation, split sizes

AdmixtureDesigner(AdmixtureDesignConfig).run() # all three stages -> four CSVs + summary JSON
```

**The `<= 1` arithmetic lives in one place.** `admixture.downsample_fractions` was split out
of `AdmixtureBuilder.compute_downsample_fractions` so both callers share one statement: the
builder still validates and raises (right, when it is about to run samtools), while the
feasibility table needs the *offending* value to say by how much a candidate depth misses.
`validate_plan` then re-checks every planned row through `AdmixtureBuilder` itself, so a
plan that MAT-13 would reject is impossible to emit rather than caught at build time.

**Two constraints shape the plan**, both covered in
[admixture-design.md](admixture-design.md) with the numbers from the real run:

- *Coverage.* The healthy downsample fraction `(1-f)*C_target/C_healthy` is largest at
  `f=0`, so the lowest rung binds and `C_target` is capped by the lowest-coverage background
  kept. Chosen: **0.9x** — the deepest candidate in phase0 §1.1's 0.1x-1x range that keeps
  >= 90% of both cohorts at `f=0` (43/44 Snyder, 500/503 Cristiano; 1.0x would drop 12 of
  Snyder's 44). The tumor side never binds.
- *Orthogonality.* §6 check 1 separates the two cohorts at accuracy 1.000
  ([cohort-qc-baseline.md](cohort-qc-baseline.md)), so any drift of pool composition with
  `f` is a batch shortcut to the label. Every rung therefore draws from the single pool of
  samples feasible at *every* rung, allocated across (source, split) cells by largest
  remainder — making `max_source_proportion_deviation` exactly 0.0.

`split_group` is assigned at source-subject level before mixing, stratified within source.
The single tumor source is assigned `feature_matrix.SPLIT_GROUP_ANY` rather than one split,
which it would otherwise collapse the whole matrix into; `assign_split_groups` resolves the
sentinel and still raises when two split-bound subjects disagree. See
[admixture-design.md](admixture-design.md) for what that costs.

**Snakemake:** `workflow/rules/admixture_design.smk`, target `admixture_design_complete`.
Workstation-only, like `local_features` and `cohort_qc`:
```
snakemake --workflow-profile workflow/profiles/local admixture_design_complete
```
