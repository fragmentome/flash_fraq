# Phase 0 Technical Spec: Fragmentomics-First ctDNA Tumor Fraction Estimator

**Status:** Draft v1
**Positioning:** Complementary to ichorCNA (CNV-based), not a reimplementation. Fragmentomics (DELFI-style) as the primary signal, CNV as a secondary/confirmatory feature.

---

## 1. Methodology grounding

### 1.1 ichorCNA (Adalsteinsson et al. 2017, *Nat Commun*)

- **Binning:** Genome divided into non-overlapping **1 Mb bins** (`T` total bins); read counts tallied per bin via HMMcopy Suite. Centromeric regions excluded (UCSC gap track).
- **GC/mappability correction:** Two sequential **LOESS regressions** on bin-wise read counts — first against GC fraction, then against mappability (uniqueness) score. A panel-of-normals from ~27 healthy donors is used to further normalize/flag systematically noisy bins.
- **HMM state structure:** 7 total copy-number states — hemizygous deletion (1 copy), neutral (2), gain (3), amplification (4), high-level amplification (5–7) — plus 2 subclonal states (subclonal deletion/gain). Transition matrix is stationary with a very high self-transition probability (~0.99999), i.e. strong prior toward long, contiguous segments. Emission model: log copy ratios follow a **Student's t-distribution** conditioned on state, with the state mean parameterized by tumor fraction and ploidy.
- **Tumor fraction estimation — the actual optimized parameter:** ichorCNA jointly estimates **tumor fraction (normal fraction `n`) and tumor ploidy (`ϕ`)** via **EM**, maximizing the complete-data log-likelihood over the HMM. The observed copy ratio in a bin is modeled as a mixture: `observed CN ∝ 2(1-n) + n·c` (normal diploid background admixed with tumor copy number `c`). EM is run from multiple initializations (`n ∈ {0.35, 0.45, 0.50, 0.65, 0.75, 0.85, 0.95}`, `ϕ ∈ {2,3,4}`) to avoid local optima, and Viterbi decoding gives the final CN segmentation. **Tumor fraction is thus a byproduct of fitting the CN mixture model, not a directly-measured quantity** — it is only identifiable when there's enough CNA signal to constrain the mixture proportion.
- **Their own validation:** Four benchmarking designs — "serial" dilutions of one breast-cancer patient's cfDNA across coverage (0.01×–1×), "merged" mixtures of 44 cancer + 18 healthy donor ULP-WGS samples at 0.1×, "exact tumor fraction" mixtures (0.1–1.0 in 0.1 steps) from 50 patient samples, and CNA spike-in limit-of-detection analysis. Reported sensitivity floor ≈ 3% tumor fraction — this is the number your model should be benchmarked against, not an arbitrary target.

### 1.2 DELFI (Cristiano et al. 2019, *Nature*)

- **Windows:** hg19 autosomes tiled into 100 kb bins, aggregated into **504 windows of 5 Mb** (~2.6 Gb covered), excluding centromeres/telomeres and ENCODE/Duke blacklist regions (~13% of genome excluded).
- **Fragment ratio per window:** Fragments classified **short (100–150 bp)** vs **long (151–220 bp)** by insert size from paired-end alignment; short:long ratio computed per 5 Mb window. GC correction applied separately to short and long counts via **LOESS (span = 3/4)**: fit read count vs. GC%, subtract predicted trend, add back genome-wide median to keep counts on the original scale.
- **Feature set fed to the classifier (no PCA on the primary published DELFI model):** GC-corrected total and short-fragment coverage across all 504 bins (centered/scaled), **z-scores of copy number for 39 chromosome arms**, and log₁₀-transformed mitochondrial DNA representation. Highly correlated bins (r > 0.9) removed during training. *(Note: some derivative/follow-on DELFI-family papers, e.g. lung-cancer fragmentome extensions, do apply PCA — retaining components explaining ~90% of variance — as an alternative dimensionality-reduction step before a logistic-regression-with-LASSO comparator model. The core 2019 Nature paper's primary reported model does not use PCA.)*
- **Model & target:** **Gradient boosted machine** (`gbm`: n.trees=150, interaction.depth=3, shrinkage=0.1). Target is **binary classification** (cancer vs. healthy), not a continuous tumor fraction — this is the key structural difference from ichorCNA's output and from what you're building.
- **Validation:** 10-fold cross-validation repeated 10×, with within-training-fold-only feature/bin selection to avoid leakage.

### 1.3 The structural gap between the two papers (why your project isn't redundant with either)

- ichorCNA outputs a **continuous, mechanistically-interpretable tumor fraction** but needs enough CNA burden to be identifiable — it's blind to tumors with a flat/near-diploid CN profile, and its 1 Mb bins are too coarse to see fragment-size signal.
- DELFI outputs a **binary/probabilistic cancer score**, not tumor fraction, and was never designed or validated as a quantitative fraction estimator.
- **Nobody in the published methodology is doing continuous tumor-fraction regression from fragmentomics features.** That's the gap you're filling — a DELFI-style feature set (short/long ratio) trained toward a continuous tumor-fraction target (ichorCNA's mixture-model estimand), validated against admixture ground truth. This is the one-sentence pitch for the "extension, not clone" framing.

Sources: [Adalsteinsson et al. 2017 (PMC5673918)](https://pmc.ncbi.nlm.nih.gov/articles/PMC5673918/) · [ichorCNA GitHub](https://github.com/broadinstitute/ichorCNA) · [Cristiano et al. 2019 (PMC6774252)](https://pmc.ncbi.nlm.nih.gov/articles/PMC6774252/) · [DELFI dataset, EGA EGAD00001005339](https://ega-archive.org/datasets/EGAD00001005339) · [SEQC2 / HCC1395 reference standard](https://github.com/bioinform/somaticseq/blob/master/docs/seqc2.md)

---

## 2. Feature set

| Feature | Source | Complementary or redundant with CNV? |
|---|---|---|
| Short:long fragment ratio per 5 Mb window (100–150bp / 151–220bp, GC-corrected) | DELFI-style | **Complementary** — captures nucleosome-protection/chromatin signal independent of copy number; present even in copy-neutral tumors |
| GC-corrected total & short-fragment coverage per window | DELFI-style | **Complementary** — depth-of-coverage texture, distinct from binned CN ratio |
| Chromosome-arm-level fragment z-scores (39 arms) | DELFI-style | **Partially redundant** with CNV arm-level gain/loss calls — include as a bridge feature, not core signal |
| Genome-wide fragment size distribution summary stats (median, mode, entropy of the size histogram) | Fragmentomics, novel to this project vs. either paper | **Complementary** — global signal, not spatially binned |
| ichorCNA-style 1 Mb bin log2 copy ratio (GC/mappability-corrected) | CNV | **Redundant with itself only where CNA burden exists** — near-diploid tumors contribute near-zero information here, which is exactly the case fragmentomics should compensate for |
| ichorCNA point estimate of tumor fraction (as an input feature, not the label) | CNV, secondary | Use only as a **sanity-check / stacking feature** in an ablation, not a primary predictor — including it directly muddies the "complementary to ichorCNA" story if it dominates the model |
| Mitochondrial DNA fraction (log10) | DELFI-style | Weak/auxiliary signal, keep for parity with DELFI features |

**Explicitly excluded from v1:** fragment end-motif / 4-mer signal, WPS (windowed protection score) nucleosome positioning, methylation. See §5.

---

## 3. Labeling / ground truth strategy: in silico admixture

**Approach:** Mix reads from tumor WGS at known ratios into real healthy cfDNA background BAMs, downsample to target coverage, and use the known mixing ratio as the regression label.

### 3.1 How this compares to the papers' own validation

- ichorCNA itself validated with **exactly this design** ("merged" and "exact tumor fraction" mixtures — real patient cfDNA mixed with real healthy donor cfDNA at known proportions, §1.1). Your approach is directly precedented, not a methodological stretch — but note ichorCNA mixed **cfDNA into cfDNA** (fragment-size distributions from both sources are already cfDNA-like). You are proposing tumor **WGS** (long, genomic-DNA fragment sizes) into healthy **cfDNA** — see 3.2, this is the single biggest technical risk in the plan.
- DELFI did not use in silico admixture at all — real patient plasma only, binary label. So there's no direct DELFI precedent for continuous-fraction admixture; you're borrowing ichorCNA's validation design and applying it to a fragmentomics feature set.

### 3.2 Technical issues to flag

1. **Fragment-size distribution mismatch (biggest issue).** Tumor tissue/cell-line WGS libraries are sheared mechanically (e.g. Covaris) to ~300–500 bp inserts with a roughly uniform/Gaussian size distribution. Real ctDNA fragments are the product of *in vivo* nuclease digestion and carry the short-fragment enrichment signal that is your entire feature set. Naively admixing genomic WGS reads into cfDNA background will teach the model to detect "sheared library vs. nuclease-digested library," not tumor fraction — a trivial and clinically meaningless classifier. **Mitigation:** either (a) computationally reshape tumor fragment-size distributions to match empirical ctDNA size profiles before mixing (in silico fragmentation), or (b) prioritize tumor **cfDNA-like** sources (e.g. HCC1395 spike-in plasma studies, or patient cfDNA mixtures as ichorCNA did) over raw tumor tissue WGS. This must be resolved before any label is trustworthy.
2. **Coverage matching.** Mixing at the read level changes effective coverage; downstream depth-dependent features (short-fragment coverage counts) need to be coverage-normalized post-mix, and the final merged BAM must be downsampled to a *consistent* target coverage (e.g. 0.1×–1× ULP-WGS matching ichorCNA's operating range) rather than left to float with the sum of two input depths.
3. **Batch/platform effects.** Healthy background and tumor reads will likely come from different sequencing centers/platforms/library prep protocols (e.g. HiSeq2500 DELFI cohort vs. whatever TCGA/SEQC2 platform the tumor reads come from). This introduces GC-bias and fragment-size batch effects that are *not* biological — they will leak into both the fragmentomics and CNV features and can masquerade as tumor-fraction-correlated signal if platform correlates with mixing ratio in your design matrix. **Mitigation:** stratify train/val/test splits so platform is never a proxy for label; run a negative control (mix 0% tumor from two different platforms and confirm the model doesn't score them differently).
4. **Reads-level vs. UMI-level mixing.** If either source uses UMI-tagged libraries, naive read-level mixing before dedup can distort effective molecule counts (duplicate collapse rates differ pre- vs. post-mix). Do the mixing on raw/pre-dedup reads and run the full pipeline (dedup, alignment, filtering) once on the merged file — never mix already-deduplicated BAMs and never dedup each source independently pre-mix.
5. **Purity of the "tumor" reads themselves.** Cell-line or tissue WGS is not itself 100% pure — need to know/estimate the actual tumor cellularity of the source before treating the nominal mixing ratio as ground truth; otherwise your "known ratio" is off by the source's own purity, silently adding label noise.

### 3.3 Candidate data sources

| Source | Type | Access | Notes |
|---|---|---|---|
| DELFI cohort (EGAD00001005339) | Healthy + cancer patient plasma cfDNA, HiSeq2500 | **DAC-controlled** (EGA, formal application + moratorium) | Best-matched healthy cfDNA background; expect weeks of lead time for access |
| HCC1395 / HCC1395BL (SEQC2) | Matched tumor/normal cell line WGS, ~55× | **Open access** (no dbGaP) | Best starting point precisely because it's unrestricted and well-characterized (~40k somatic SNVs); still has the fragment-size-mismatch issue in 3.2 since it's tissue WGS, not cfDNA |
| TCGA WGS BAMs (GDC/dbGaP, phs000178) | Tumor tissue WGS | **Controlled**, requires eRA Commons + dbGaP authorization via a PI | Broad tumor-type coverage but heavier access burden than HCC1395; only pursue if HCC1395 diversity proves insufficient |
| ichorCNA's own panel-of-normals / published healthy ULP-WGS | Healthy cfDNA, low coverage | Check ichorCNA GitHub for released panel | Useful as a second healthy background source for cross-platform validation control (issue 3.3 above) |

**Recommendation for v1:** Start with HCC1395 (open access, fast to obtain) mixed into whichever healthy cfDNA background you can access fastest, explicitly building the fragment-size-reshaping step from day one rather than treating it as a later fix — it's load-bearing for label validity, not a nice-to-have.

---

## 4. Model type: gradient boosting, not deep learning

**Recommendation: feature-based gradient boosting (XGBoost or LightGBM)**, matching DELFI's own choice of `gbm` for its classifier.

Justification given expected sample sizes:
- In silico admixture gives you effectively unlimited *label* combinations, but the number of **independent background samples** (distinct healthy/tumor source BAMs you can mix) is realistically in the tens-to-low-hundreds given access constraints (§3.3) — this is a small-N regime by deep learning standards, even before accounting for the need for held-out subjects (not just held-out mixing ratios) in any valid split.
- The feature set (§2) is already hand-engineered, low-dimensional (hundreds of window-level features, not raw sequence/read-level input), and structurally similar to DELFI's tabular design — there's no case for a representation-learning model (CNN over WPS tracks, transformer over reads) when the input is already a compact, interpretable feature vector.
- Gradient boosting gives you built-in feature importance for the "where's the complementary signal coming from" story that's central to the complementary-to-ichorCNA positioning — directly useful for figures/explanation, unlike a neural net.
- Matches precedent from both papers under review: ichorCNA's own signal is EM/HMM-based (not deep learning), DELFI's classifier is `gbm`. Diverging into deep learning would be a bigger methodological leap than the project needs and would weaken (not strengthen) the "principled extension" narrative.

**Defer:** any sequence-level/deep learning approach (e.g. learned representations of fragment end-motifs or raw coverage tracks) — flag as a stretch/v2 idea only if v1 shows a clear ceiling from hand-engineered features.

---

## 5. Validation strategy

### 5.1 Primary metrics

- **RMSE (and MAE) against known admixture ratio** — your core quantitative metric, directly comparable to ichorCNA's own dilution-series validation design (§1.1, §3.1).
- **Correlation (Pearson/Spearman) vs. ichorCNA's own tumor-fraction calls** on the same public benchmark samples — this is the key "complementary, not competing" evidence: report where your model and ichorCNA agree (validates both), and specifically characterize the **disagreement cases**, especially low-CNA-burden / near-diploid samples where ichorCNA is expected to be uninformative and fragmentomics is expected to add value. That disagreement analysis, not aggregate correlation, is the actual deliverable for the "extension of the state of the art" framing.
- **Detection floor / limit-of-blank:** minimum tumor fraction at which your model's estimate is statistically distinguishable from the 0%-tumor-fraction null distribution (analogous to ichorCNA's reported ~3% floor) — report this explicitly as a comparable number.

### 5.2 Split design (critical, given §3.2 risk)

- Split by **source subject/sample**, not by mixing ratio — no healthy background or tumor source that appears in training may appear (at any ratio) in validation/test. Otherwise the model can memorize source-specific artifacts rather than learning tumor-fraction signal.
- Hold out at least one **entire data source/platform** (e.g. train on HCC1395-based mixtures, test generalization on a DELFI-cohort-based mixture set once/if access is granted) to directly test for the batch-effect failure mode in §3.2.

### 5.3 External benchmark

- If/when public ULP-WGS samples with **both** a published ichorCNA tumor-fraction call and known clinical context are available (e.g. supplementary data from the ichorCNA paper itself, or from later ichorCNA-using studies), run your model on the same raw BAMs and report correlation — this is a real-world check beyond synthetic admixture, and the strongest evidence against a "just reimplementing ichorCNA" critique.

---

## 6. Scope boundaries

**In v1:**
- Fragment-ratio (short:long) and coverage features per 5 Mb window, GC-corrected
- Chromosome-arm z-scores and 1 Mb CNV bin features as secondary/complementary inputs
- Gradient boosting regression → continuous tumor fraction
- In silico admixture labeling (HCC1395-based, reshaped fragment sizes)
- Benchmarking against ichorCNA's own calls on public/benchmark data

**Explicitly deferred (not v1):**
- **Methylation signal** (e.g. WGBS-derived fragmentomics, methylation-informed tissue-of-origin) — genuinely complementary in principle, but a separate assay/data modality with its own access and preprocessing burden; scope out unless you decide otherwise.
- Fragment end-motif (4-mer) features and WPS nucleosome-positioning tracks — real DELFI-family signal, but adds feature engineering complexity beyond the core short:long ratio; candidate for v1.1 if the core feature set underperforms.
- Cancer type classification / tissue-of-origin — this project targets tumor **fraction**, not cancer detection or localization; don't let scope drift toward reproducing DELFI's actual classification task.
- Deep learning / representation learning approaches (§4).
- Any tumor-type-specific model tuning — v1 should be a general tumor-fraction estimator, not per-cancer-type variants.

---

## 7. Open risks / unknowns

1. **Label validity risk (highest severity):** the fragment-size mismatch between tissue-WGS "tumor" reads and true ctDNA (§3.2.1) is not a minor caveat — if unaddressed, the model learns a library-prep artifact, not tumor fraction, and every downstream metric becomes meaningless. This must be solved or explicitly worked around before any training run is treated as informative.
2. **Data access latency:** the best-matched healthy cfDNA source (DELFI EGA cohort) requires DAC approval with unknown turnaround time; TCGA WGS requires dbGaP authorization via a PI. Both can stall the project for weeks; HCC1395's open access mitigates this for the tumor side but not the healthy-background side — identify an open-access healthy cfDNA source early (e.g. published low-pass healthy plasma WGS in SRA/GEO) as a fallback.
3. **Low tumor fraction detection floor:** ichorCNA's own floor is ~3%; fragmentomics signal is generally considered noisier per-window than CNV at very low coverage, so there's a real chance your model's floor is *worse* than ichorCNA's at matched coverage, which complicates (but doesn't invalidate) the complementary-signal story — plan to report this honestly rather than only in favorable regimes.
4. **Small effective sample size for held-out validation** (§4) — with realistically tens of independent subjects, confidence intervals on RMSE/correlation will be wide; don't overclaim precision from a benchmark with few independent sources even if the number of synthetic mixtures is large.
5. **Batch effects being invisible until tested:** platform/protocol differences between tumor and healthy sources could produce a model that looks excellent in-distribution and fails on any new data source — the source-held-out split (§5.2) is the main defense, but this risk can't be fully retired without an external, truly independent benchmark set.
6. **Positioning risk (non-technical but real):** given your role as the organization's ichorCNA domain owner, any internal presentation of this work needs the "complementary, not competing" framing to be backed by the disagreement analysis (§5.1), not just aggregate correlation — a purely aggregate "matches ichorCNA" result reads as a redundant reimplementation rather than an extension.
