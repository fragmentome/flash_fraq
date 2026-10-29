"""Decoy classifier QC gate (phase2-preprocessing-plan.md §6).

Trains a classifier to predict source/platform of origin from the fragmentomics
feature matrix and checks whether that predicted source signal leaks into the
tumor-fraction label. This must pass before any tumor-fraction regression
result is treated as informative.
"""

import json
from dataclasses import dataclass, field
from pathlib import Path

from numpy import ndarray, unique
from pandas import DataFrame, Series, read_csv, to_numeric
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.model_selection import KFold, StratifiedKFold, cross_val_predict, train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

# The CNV feature block (feature_matrix.pivot_cnv_bins/pivot_arm_zscores and the
# coverage-window pivot) -- the one part of the feature set real ctDNA copy-number
# biology is expected to drive as a function of tumor fraction. Check 3 residualises
# on exactly this block so a genuine CNV-driven correlation between the non-CNV
# (fragmentomics) block and f is not mistaken for a library confound.
CNV_FEATURE_PREFIXES: tuple[str, ...] = ("cnv_logr_bin_", "arm_zscore_", "gc_corrected_total_coverage_window_")


@dataclass(frozen=True)
class DecoyClassifierQcConfig:
    """Configuration for a single decoy classifier QC run."""

    feature_matrix_path: Path
    output_dir: Path
    label_column: str = "empirical_tumor_fraction"
    platform_columns: tuple[str, ...] = ("platform_tumor", "platform_healthy")
    extra_non_feature_columns: tuple[str, ...] = ()
    non_feature_columns: tuple[str, ...] = field(init=False)
    fixed_fraction_target: float = 0.0
    fixed_fraction_tolerance: float = 0.05
    correlation_threshold: float = 0.1
    # Check 3 (library-confound probe) target column and threshold. Separate from
    # `label_column`/`correlation_threshold`: check 2 gates on the fitted label
    # (`empirical_tumor_fraction`), check 3 is specified against the *nominal* design
    # fraction, since it is checking a library-identity leak, not label correctness.
    nominal_fraction_column: str = "nominal_tumor_fraction"
    library_confound_correlation_threshold: float = 0.1
    cnv_feature_prefixes: tuple[str, ...] = CNV_FEATURE_PREFIXES
    # A matrix pooling more than one regime must never be handed to check 3 (or the
    # gate as a whole) implicitly -- see admixture_design's REGIME_MATCHED_PAIR /
    # REGIME_REAL_BACKGROUND. None only tolerated when the matrix carries at most one
    # regime already (e.g. the pre-admixture cohort matrix, which has no regime column).
    regime: str | None = None
    regime_column: str = "regime"
    # Feature columns missing in more than this fraction of rows are dropped rather than
    # imputed. Real per-sample tables are ragged (windows/arms/bins that drop out of one
    # sample's GC-correction or ichorCNA fit), and sklearn rejects NaN outright, so the
    # choice is between dropping a column and inventing values for it. 0.0 keeps only
    # columns every sample has -- raise it to trade completeness for feature count.
    max_missing_fraction: float = 0.0
    test_size: float = 0.3
    # Check 2's out-of-fold folds. Every probability the correlation is computed on comes
    # from a model that never saw that row, so the number is not an artefact of a
    # 4354-feature logistic regression memorising 547 rows. The smallest source must have
    # at least this many samples or the check is skipped rather than run on fewer folds --
    # silently changing the fold count would make two reports incomparable.
    n_splits: int = 5
    # Check 2 is degenerate when the classifier's probabilities carry no more information
    # than the source indicator itself; |classifier corr - indicator corr| under this is
    # treated as "the same number", and the report says so.
    degeneracy_tolerance: float = 0.005
    random_state: int = 0

    def __post_init__(self) -> None:
        non_feature_columns = (
            "sample_id",
            "tumor_source_id",
            "healthy_source_id",
            "tumor_subject_id",
            "healthy_subject_id",
            "reference_build",
            "target_coverage",
            "ichorCNA_tumor_fraction_estimate",
            # See tf_model.NON_FEATURE_COLUMNS: empirical_tumor_fraction_sites *is* the
            # empirical_tumor_fraction label verbatim in this cohort, and
            # empirical_tumor_fraction_clonal is the other cross-estimator on the same
            # SNVs -- neither is a fragmentomics feature, and both leak whichever label
            # this config is checking correlation against.
            "empirical_tumor_fraction_sites",
            "empirical_tumor_fraction_clonal",
            "split_group",
            # audit finding 2: a tumor-read counter (HCC1395 gDNA ~10^-2.6 vs cfDNA
            # ~10^-4.9 mtDNA fraction), never a fragmentomics feature -- see
            # feature_matrix.get_primary_feature_columns / tf_model.NON_FEATURE_COLUMNS.
            "mtdna_fraction_log10",
            self.regime_column,
            self.nominal_fraction_column,
            self.label_column,
            *self.platform_columns,
            *self.extra_non_feature_columns,
        )
        object.__setattr__(self, "non_feature_columns", non_feature_columns)

    @classmethod
    def from_sources(cls, feature_matrix_path: Path, output_dir: Path) -> "DecoyClassifierQcConfig":
        return cls(feature_matrix_path=feature_matrix_path, output_dir=output_dir)

    @classmethod
    def for_cohort(cls, feature_matrix_path: Path, output_dir: Path, **overrides: object) -> "DecoyClassifierQcConfig":
        """Config for `cohort_matrix`'s one-row-per-real-sample table.

        Two adaptations, both forced by there being no admixture yet:

        - The source label is `source_dataset|platform` (one library per row) instead of
          the `platform_tumor|platform_healthy` pair a mixed sample carries.
        - The fraction in check 2 is `ichorCNA_tumor_fraction_estimate`, not
          `empirical_tumor_fraction`. **This weakens the check's meaning**: ichorCNA's
          estimate is itself derived from these samples' coverage, and it is 0.0 for most
          of the cohort, so a correlation here reads as "source predicts ichorCNA's CNV
          call", not "source predicts true tumor burden". Treat a pass as necessary but
          not sufficient; the real gate runs on `empirical_tumor_fraction` once
          admixture exists.
        """
        return cls(
            feature_matrix_path=feature_matrix_path,
            output_dir=output_dir,
            label_column="ichorCNA_tumor_fraction_estimate",
            platform_columns=("source_dataset", "platform"),
            extra_non_feature_columns=("reference_build",),
            **overrides,
        )


@dataclass(frozen=True)
class FixedFractionCheckResult:
    """Check 1: at a fixed tumor fraction, can the decoy classifier still discriminate source?"""

    fraction_tested: float
    n_samples: int
    n_features: int
    accuracy: float
    baseline_accuracy: float
    source_discriminable: bool


@dataclass(frozen=True)
class FractionCorrelationCheckResult:
    """Check 2: does out-of-fold predicted-source-probability correlate with the tumor-fraction label?

    `max_abs_correlation` is None when the check could not run (one source, no usable
    features, or a source smaller than `n_splits`); `skipped_reason` then says why. A
    skipped check is not a measured 0.0 and must not read as one.

    `source_indicator_correlation_by_source` is the plain corr(is-this-source, fraction),
    which needs no classifier at all. When it matches the classifier-derived number within
    `degeneracy_tolerance`, `degenerate_with_source_indicator` is True: the classifier
    separates the sources well enough that its probabilities are saturated, and check 2 is
    reporting the source indicator under another name.
    """

    max_abs_correlation: float | None
    correlation_by_source: dict[str, float]
    source_indicator_correlation_by_source: dict[str, float]
    max_abs_source_indicator_correlation: float | None
    degenerate_with_source_indicator: bool
    n_features: int
    n_splits: int | None
    fold_random_state: int | None
    flagged: bool
    skipped_reason: str | None = None


@dataclass(frozen=True)
class LibraryConfoundCheckResult:
    """Check 3: does the non-CNV (fragmentomics) feature block predict the *nominal*
    tumor fraction beyond what the CNV block already explains?

    The real confound this check exists for (audit finding 4) does not show up as a
    source classifier: it is inside every mixture (tumor library vs background
    library) and scales with f, so no per-row source label carries it. Instead this
    fits an out-of-fold Ridge of the non-CNV block directly against
    `nominal_tumor_fraction`, residualised on an out-of-fold Ridge fit of the CNV
    block against the same target -- a real ctDNA CNV-driven correlation between the
    non-CNV block and f is absorbed by the residual, so what is left is attributable
    to library identity rather than tumor biology. `spearman_rho` is the correlation
    of the non-CNV out-of-fold prediction with that residual.

    Residualisation is skipped (and `residualised` is False) when no CNV feature
    column survives `prepare_features`, e.g. a matrix with no CNV block at all; the
    check then correlates the non-CNV block's out-of-fold prediction directly
    against `nominal_tumor_fraction`.

    `spearman_rho` is None when the check could not run; `skipped_reason` says why,
    mirroring check 2's `max_abs_correlation is None` convention.
    """

    spearman_rho: float | None
    per_feature_spearman_rho: dict[str, float]
    n_samples: int
    n_non_cnv_features: int
    n_cnv_features: int
    n_splits: int | None
    residualised: bool
    flagged: bool
    skipped_reason: str | None = None


@dataclass(frozen=True)
class DecoyClassifierQcReport:
    fixed_fraction_check: FixedFractionCheckResult
    correlation_check: FractionCorrelationCheckResult
    library_confound_check: LibraryConfoundCheckResult
    passed: bool


class DecoyClassifierQc:
    """Batch-effect QC gate: predicts source/platform of origin from fragmentomics
    features and checks whether that signal correlates with the tumor-fraction label.

    Per phase2-preprocessing-plan.md §6, check 2 (predicted-source-probability vs.
    true tumor-fraction correlation) is the concrete red flag that gates `passed`;
    check 1 (fixed-fraction discriminability) is reported for context.
    """

    def __init__(self, config: DecoyClassifierQcConfig) -> None:
        self.config = config

    def load_feature_matrix(self) -> DataFrame:
        return read_csv(self.config.feature_matrix_path)

    def get_feature_columns(self, feature_matrix: DataFrame) -> list[str]:
        return [
            column
            for column in feature_matrix.columns
            if column not in self.config.non_feature_columns
        ]

    def get_source_label(self, feature_matrix: DataFrame) -> Series:
        """The decoy classifier's target: every configured source/platform column, joined.

        Any number of columns, not exactly two: a mixed admixture sample identifies its
        source as the `platform_tumor|platform_healthy` pair, while a real single-source
        cohort sample identifies it as `source_dataset|platform`.
        """
        columns = [feature_matrix[column].astype(str) for column in self.config.platform_columns]
        joined = columns[0]
        for column in columns[1:]:
            joined = joined + "|" + column
        return joined

    def prepare_features(self, feature_matrix: DataFrame) -> DataFrame:
        """Numeric, NaN-free, non-constant feature block for the classifier.

        The synthetic matrices in the tests need none of this, but a matrix assembled
        from real per-sample tables does: window/arm/bin sets differ between samples, so
        the pivot is ragged, and ichorCNA leaves blank `cnv_logr` in bins it filtered
        out. sklearn rejects NaN, and a column that is constant across the cohort (or
        that survives only in a subset) carries no source signal anyway, so dropping it
        makes the reported feature count honest rather than shrinking the check.

        Column selection is done on the *full* matrix by the caller's slice, so the
        fixed-fraction subset and the full-range check can legitimately keep different
        columns -- each reports its own `n_features`.
        """
        features = feature_matrix[self.get_feature_columns(feature_matrix)]
        features = features.apply(to_numeric, errors="coerce")

        missing_fraction = features.isna().mean()
        features = features.loc[:, missing_fraction <= self.config.max_missing_fraction]
        # Only reachable with max_missing_fraction > 0; the median is the cohort's own,
        # so an imputed cell never invents an out-of-range value.
        features = features.fillna(features.median())
        return features.loc[:, features.std(ddof=0) > 0]

    def build_classifier(self) -> Pipeline:
        """Standardise, then fit logistic regression.

        The scaler is not cosmetic here: the §9 feature block mixes arm z-scores (order
        1) with GC-corrected window coverages (order 1e4), and unscaled, L2-penalised
        logistic regression effectively regularises the small-scale columns away and
        stops converging within max_iter.
        """
        return Pipeline(
            [
                ("scale", StandardScaler()),
                ("classify", LogisticRegression(max_iter=1000, random_state=self.config.random_state)),
            ]
        )

    def run_fixed_fraction_check(self, feature_matrix: DataFrame) -> FixedFractionCheckResult:
        label = feature_matrix[self.config.label_column]
        fixed_fraction_mask = (label - self.config.fixed_fraction_target).abs() <= self.config.fixed_fraction_tolerance
        fixed_subset = feature_matrix.loc[fixed_fraction_mask]

        features = self.prepare_features(fixed_subset)
        source_label = self.get_source_label(fixed_subset)

        # A stratified split needs at least one sample per class on each side; below that
        # the check cannot run at all, which is not the same as "source is unpredictable".
        too_few_per_source = source_label.value_counts().min() < 2 if len(source_label) else True
        if source_label.nunique() < 2 or too_few_per_source or features.empty:
            baseline_accuracy = 1.0
            return FixedFractionCheckResult(
                fraction_tested=self.config.fixed_fraction_target,
                n_samples=fixed_subset.shape[0],
                n_features=features.shape[1],
                accuracy=baseline_accuracy,
                baseline_accuracy=baseline_accuracy,
                source_discriminable=False,
            )

        features_train, features_test, source_train, source_test = train_test_split(
            features,
            source_label,
            test_size=self.config.test_size,
            random_state=self.config.random_state,
            stratify=source_label,
        )

        classifier = self.build_classifier()
        classifier.fit(features_train, source_train)
        accuracy = float(classifier.score(features_test, source_test))
        baseline_accuracy = float(source_test.value_counts(normalize=True).max())

        return FixedFractionCheckResult(
            fraction_tested=self.config.fixed_fraction_target,
            n_samples=fixed_subset.shape[0],
            n_features=features.shape[1],
            accuracy=accuracy,
            baseline_accuracy=baseline_accuracy,
            source_discriminable=bool(accuracy > baseline_accuracy),
        )

    def correlate_with_fraction(self, values: ndarray, tumor_fraction: Series) -> float:
        """Pearson correlation, computed the way pandas does it: NaN rows dropped pairwise.

        `tumor_fraction` carries NaN wherever the label failed to parse, and the label is
        the one column this gate cannot afford to silently fill in.
        """
        return float(DataFrame({"values": values, "fraction": tumor_fraction}).corr().iloc[0, 1])

    def skipped_correlation_check(self, reason: str, n_features: int) -> FractionCorrelationCheckResult:
        """A check that could not run. Deliberately not a correlation of 0.0.

        `flagged` is False because nothing was measured, so nothing can be a red flag --
        but `max_abs_correlation` is None rather than a number a reader could mistake for
        evidence, and `skipped_reason` says which precondition failed.
        """
        return FractionCorrelationCheckResult(
            max_abs_correlation=None,
            correlation_by_source={},
            source_indicator_correlation_by_source={},
            max_abs_source_indicator_correlation=None,
            degenerate_with_source_indicator=False,
            n_features=n_features,
            n_splits=None,
            fold_random_state=None,
            flagged=False,
            skipped_reason=reason,
        )

    def run_fraction_correlation_check(self, feature_matrix: DataFrame) -> FractionCorrelationCheckResult:
        """Check 2, out-of-fold.

        Every probability comes from `cross_val_predict` under a shuffled `StratifiedKFold`,
        so no row is scored by a model that saw it. Fitting in-sample instead made this
        check meaningless on the real cohort: a LogisticRegression over 4354 features and
        547 rows separates its own training set perfectly, the probabilities saturate at 0
        and 1, and the correlation collapses onto plain corr(source indicator, fraction) --
        observed as 0.02199 reported against 0.02198 computed directly.

        The scaler stays inside the pipeline `build_classifier()` returns precisely so that
        `cross_val_predict` refits it per fold; a scaler fitted once on the full matrix
        would leak the held-out rows' means back into training.

        The direct source-indicator correlation is reported alongside regardless, because
        out-of-fold probabilities can *still* saturate when the sources are genuinely
        separable. When the two numbers coincide the report says so
        (`degenerate_with_source_indicator`) rather than leaving the next reader to
        rediscover it.
        """
        features = self.prepare_features(feature_matrix)
        source_label = self.get_source_label(feature_matrix)
        tumor_fraction = to_numeric(feature_matrix[self.config.label_column], errors="coerce")

        # A constant label makes every correlation NaN, and `NaN > threshold` is False --
        # the check would report a clean pass while having measured nothing. That is the
        # one failure mode this gate must never have, so refuse to run instead.
        if tumor_fraction.std(ddof=0) == 0 or tumor_fraction.isna().all():
            raise ValueError(
                f"{self.config.label_column!r} has no variance across the matrix; the correlation check "
                f"cannot distinguish leaked source signal from no signal"
            )

        if source_label.nunique() < 2:
            # Nothing to leak: with one source there is no source signal to correlate.
            return self.skipped_correlation_check(
                f"only one source ({source_label.iloc[0]!r}) in the matrix; no source signal to leak",
                features.shape[1],
            )
        if features.empty:
            return self.skipped_correlation_check(
                f"no feature column survives max_missing_fraction={self.config.max_missing_fraction} "
                f"and the constant-column drop",
                features.shape[1],
            )

        # StratifiedKFold needs every class represented in every fold. Letting sklearn
        # raise here would turn a small source into a crashed workflow; quietly lowering
        # n_splits instead would make two reports' numbers incomparable without saying so.
        source_counts = source_label.value_counts()
        smallest_source, smallest_count = str(source_counts.index[-1]), int(source_counts.iloc[-1])
        if smallest_count < self.config.n_splits:
            return self.skipped_correlation_check(
                f"smallest source {smallest_source!r} has {smallest_count} samples, fewer than "
                f"n_splits={self.config.n_splits}; out-of-fold probabilities cannot be computed",
                features.shape[1],
            )

        folds = StratifiedKFold(
            n_splits=self.config.n_splits, shuffle=True, random_state=self.config.random_state
        )
        predicted_probability = cross_val_predict(
            self.build_classifier(), features, source_label, cv=folds, method="predict_proba"
        )

        # cross_val_predict orders its probability columns by np.unique(y), not by any one
        # fold's fitted classes_ -- there is no fitted estimator to ask.
        classes = unique(source_label)
        correlation_by_source = {
            str(source): self.correlate_with_fraction(predicted_probability[:, index], tumor_fraction)
            for index, source in enumerate(classes)
        }
        source_indicator_correlation_by_source = {
            str(source): self.correlate_with_fraction(
                (source_label == source).astype(float).to_numpy(), tumor_fraction
            )
            for source in classes
        }

        max_abs_correlation = max(abs(correlation) for correlation in correlation_by_source.values())
        max_abs_indicator_correlation = max(
            abs(correlation) for correlation in source_indicator_correlation_by_source.values()
        )
        degenerate = all(
            abs(correlation_by_source[source] - source_indicator_correlation_by_source[source])
            <= self.config.degeneracy_tolerance
            for source in correlation_by_source
        )

        return FractionCorrelationCheckResult(
            max_abs_correlation=max_abs_correlation,
            correlation_by_source=correlation_by_source,
            source_indicator_correlation_by_source=source_indicator_correlation_by_source,
            max_abs_source_indicator_correlation=max_abs_indicator_correlation,
            degenerate_with_source_indicator=bool(degenerate),
            n_features=features.shape[1],
            n_splits=self.config.n_splits,
            fold_random_state=self.config.random_state,
            flagged=bool(max_abs_correlation > self.config.correlation_threshold),
        )

    def enforce_single_regime(self, feature_matrix: DataFrame) -> DataFrame:
        """Refuse to run on more than one regime pooled together (audit finding 4;
        `matched_pair`/`real_background` are a hard partition, see
        `admixture_design.REGIME_MATCHED_PAIR`/`REGIME_REAL_BACKGROUND`).

        A matrix with no regime column (the pre-admixture cohort matrix) passes
        through unchanged. A matrix with a regime column and `self.config.regime` set
        is filtered to that regime. A matrix with a regime column, no configured
        regime, and more than one distinct value is refused outright rather than
        silently pooled.
        """
        regime_column = self.config.regime_column
        if regime_column not in feature_matrix.columns:
            return feature_matrix

        if self.config.regime is not None:
            filtered = feature_matrix.loc[feature_matrix[regime_column] == self.config.regime]
            if filtered.empty:
                raise ValueError(f"no rows with {regime_column}={self.config.regime!r} in the matrix")
            return filtered

        distinct_regimes = feature_matrix[regime_column].unique()
        if len(distinct_regimes) > 1:
            raise ValueError(
                f"matrix pools {len(distinct_regimes)} regimes {sorted(map(str, distinct_regimes))}; "
                f"pass an explicit regime rather than running the gate pooled"
            )
        return feature_matrix

    def get_cnv_and_non_cnv_features(self, features: DataFrame) -> tuple[DataFrame, DataFrame]:
        """Split a prepared feature block into its CNV and non-CNV (fragmentomics) halves.

        Runs after `prepare_features` so both halves are already numeric, NaN-free
        and non-constant.
        """
        cnv_columns = [column for column in features.columns if column.startswith(self.config.cnv_feature_prefixes)]
        non_cnv_columns = [column for column in features.columns if column not in cnv_columns]
        return features[cnv_columns], features[non_cnv_columns]

    def build_regressor(self) -> Pipeline:
        """Standardise, then fit ridge regression -- check 3's regressor.

        Mirrors `build_classifier()`'s reasoning: the feature block mixes arm
        z-scores with GC-corrected window coverages on wildly different scales, and
        unscaled ridge regularises the small-scale columns away.
        """
        return Pipeline(
            [
                ("scale", StandardScaler()),
                ("regress", Ridge(random_state=self.config.random_state)),
            ]
        )

    def spearman_with_fraction(self, values: ndarray, fraction: Series) -> float:
        """Spearman rank correlation, NaN rows dropped pairwise -- see
        `correlate_with_fraction`'s docstring for why that matters here too."""
        return float(DataFrame({"values": values, "fraction": fraction}).corr(method="spearman").iloc[0, 1])

    def skipped_library_confound_check(
        self, reason: str, n_non_cnv_features: int, n_cnv_features: int
    ) -> LibraryConfoundCheckResult:
        """A check that could not run -- see `skipped_correlation_check`'s docstring;
        the same "None, not a measured 0.0" discipline applies here."""
        return LibraryConfoundCheckResult(
            spearman_rho=None,
            per_feature_spearman_rho={},
            n_samples=0,
            n_non_cnv_features=n_non_cnv_features,
            n_cnv_features=n_cnv_features,
            n_splits=None,
            residualised=False,
            flagged=False,
            skipped_reason=reason,
        )

    def run_library_confound_check(self, feature_matrix: DataFrame) -> LibraryConfoundCheckResult:
        """Check 3, out-of-fold -- see the class docstring on `LibraryConfoundCheckResult`
        for what this measures and why it residualises on the CNV block.

        Every prediction comes from a fold that never saw that row, for the same
        reason `run_fraction_correlation_check` insists on it: an in-sample ridge
        over a wide feature block would fit the training rows well enough to
        saturate, and the correlation would stop measuring generalisable signal.
        """
        if self.config.nominal_fraction_column not in feature_matrix.columns:
            # The pre-admixture cohort matrix (DecoyClassifierQcConfig.for_cohort) has no
            # design fraction to check a library confound against -- only real admixture
            # mixtures carry one.
            return self.skipped_library_confound_check(
                f"{self.config.nominal_fraction_column!r} not in the matrix; check 3 needs an "
                f"admixture matrix's design fraction",
                0,
                0,
            )

        nominal_fraction = to_numeric(feature_matrix[self.config.nominal_fraction_column], errors="coerce")
        if nominal_fraction.std(ddof=0) == 0 or nominal_fraction.isna().all():
            # Unlike check 2's label, a constant nominal fraction is a legitimate
            # precondition failure rather than a correctness bug (e.g. a fixed-fraction
            # slice of the matrix) -- skip rather than raise, since check 3 is an
            # ancillary probe, not the gate's primary label.
            return self.skipped_library_confound_check(
                f"{self.config.nominal_fraction_column!r} has no variance across the matrix; check 3 "
                f"cannot distinguish a library confound from no signal",
                0,
                0,
            )

        features = self.prepare_features(feature_matrix)
        cnv_features, non_cnv_features = self.get_cnv_and_non_cnv_features(features)

        per_feature_spearman_rho = {
            column: self.spearman_with_fraction(non_cnv_features[column].to_numpy(), nominal_fraction)
            for column in non_cnv_features.columns
        }

        if non_cnv_features.empty:
            return self.skipped_library_confound_check(
                f"no non-CNV feature column survives max_missing_fraction={self.config.max_missing_fraction} "
                f"and the constant-column drop",
                non_cnv_features.shape[1],
                cnv_features.shape[1],
            )
        if len(feature_matrix) < self.config.n_splits:
            return self.skipped_library_confound_check(
                f"{len(feature_matrix)} samples, fewer than n_splits={self.config.n_splits}; "
                f"out-of-fold predictions cannot be computed",
                non_cnv_features.shape[1],
                cnv_features.shape[1],
            )

        folds = KFold(n_splits=self.config.n_splits, shuffle=True, random_state=self.config.random_state)

        residualised = not cnv_features.empty
        target = (
            nominal_fraction - cross_val_predict(self.build_regressor(), cnv_features, nominal_fraction, cv=folds)
            if residualised
            else nominal_fraction
        )
        non_cnv_prediction = cross_val_predict(self.build_regressor(), non_cnv_features, target, cv=folds)
        spearman_rho = self.spearman_with_fraction(non_cnv_prediction, target)

        return LibraryConfoundCheckResult(
            spearman_rho=spearman_rho,
            per_feature_spearman_rho=per_feature_spearman_rho,
            n_samples=len(feature_matrix),
            n_non_cnv_features=non_cnv_features.shape[1],
            n_cnv_features=cnv_features.shape[1],
            n_splits=self.config.n_splits,
            residualised=residualised,
            flagged=bool(abs(spearman_rho) > self.config.library_confound_correlation_threshold),
        )

    def run(self) -> DecoyClassifierQcReport:
        feature_matrix = self.enforce_single_regime(self.load_feature_matrix())
        fixed_fraction_check = self.run_fixed_fraction_check(feature_matrix)
        correlation_check = self.run_fraction_correlation_check(feature_matrix)
        library_confound_check = self.run_library_confound_check(feature_matrix)
        report = DecoyClassifierQcReport(
            fixed_fraction_check=fixed_fraction_check,
            correlation_check=correlation_check,
            library_confound_check=library_confound_check,
            passed=not correlation_check.flagged and not library_confound_check.flagged,
        )
        self.write_report(report)
        return report

    def report_path(self) -> Path:
        # A regime-scoped run (admixture matrix-kind requires --regime, audit finding 4)
        # must not silently overwrite another regime's report when both are written to
        # the same --output-dir -- suffix the filename with the regime in that case.
        filename = f"decoy_qc_report_{self.config.regime}.json" if self.config.regime else "decoy_qc_report.json"
        return self.config.output_dir / filename

    def write_report(self, report: DecoyClassifierQcReport) -> None:
        self.config.output_dir.mkdir(parents=True, exist_ok=True)
        report_path = self.report_path()

        report_path.write_text(
            json.dumps(
                {
                    "fixed_fraction_check": vars(report.fixed_fraction_check),
                    "correlation_check": vars(report.correlation_check),
                    "library_confound_check": vars(report.library_confound_check),
                    "passed": report.passed,
                },
                indent=2,
                default=float,
            )
        )
