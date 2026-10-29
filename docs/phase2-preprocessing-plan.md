# Phase 2: Preprocessing Pipeline Plan

**Status:** Planning only — no alignment runs, no real file processing executed in this pass.
**Consistency check against Phase 0/1:** This plan is scoped to sources Phase 1 rated **usable now** — GSE71378 (Snyder), the Cristiano/DELFI pancreatic BAM subset, HCC1395/HCC1395BL (SEQC2), and 1000 Genomes NYGC (as negative-control only). TCGA (blocked), FinaleDB (unconfirmed), and the ichorCNA PoN-as-fragmentomics-source (invalid per Phase 1 flag #1) are explicitly **not** planned around below. Where this creates a gap vs. Phase 0's original source list, it's called out rather than silently patched.

---

## 0. Reference genome build: what Phase 1 actually handed us

The prompt's premise — "GDC/TCGA is hg38 and older EGA/GEO cfDNA studies are often hg19" — turns out to be **moot for v1**, because TCGA is the one source Phase 1 marked **Blocked**. Once TCGA is out of scope, the actual build landscape is:

| Source | Build | Status |
|---|---|---|
| GSE71378 (Snyder) | **hg19** (path literally says `original_data_hg19`) | Already-aligned BAM, on hand |
| Cristiano DELFI pancreatic subset | **hg19** (DELFI 2019 paper's own windows are defined on hg19; EGA cfDNA studies of this era are hg19 by convention) | Already-aligned BAM, on hand |
| HCC1395/HCC1395BL (SEQC2) | Build-agnostic at the source (raw FASTQ via SRA) — a pre-aligned BAM mirror also exists on NCBI FTP but its build is **not confirmed** by Phase 1 | Not yet downloaded |
| 1000 Genomes NYGC | **GRCh38** (CRAM) | Open, but flagged in Phase 1 as genomic DNA, not cfDNA — negative-control use only |

So the two admixture-relevant, already-in-hand cfDNA backgrounds (Snyder, Cristiano) are **both hg19**, and DELFI's own published feature windows (§1.2 of the spec) are defined on hg19 — the natural common build isn't a compromise, it's already the majority position. The only build tension is the fallback 1000G source (GRCh38), which per Phase 0/1 is not part of the admixture training set anyway.

**Decided:** the whole v1 pipeline standardizes on **hg19** — confirmed, all in-hand/planned sources are hg19. Align HCC1395's raw FASTQ to hg19 (do **not** use the NCBI FTP pre-aligned BAM mirror unless its build is confirmed hg19 — verify via `samtools view -H` and the `@SQ` line lengths before trusting it). Keep 1000G NYGC in its native GRCh38/CRAM for the isolated batch-effect negative-control check only (§6) — never read-level-merge it with anything, so its build never has to reconcile with the rest of the matrix.

The remaining open item is narrower than a build choice: confirm the exact reference FASTA *variant* (decoy sequences, `chr`-prefixed vs. not) behind the existing Snyder/Cristiano hg19 BAMs, so HCC1395's alignment uses the identical contig set — a build-name match ("hg19") isn't sufficient on its own if the decoy/contig list differs, since that's what `samtools merge` checks.

This is a **reference build consistency plan**, not a "realign everything" plan — see §1 for why partial realignment (only the one source that needs it) is safer than a blanket realign.

---

## 1. Pipeline stages, end to end

Two source classes need genuinely different treatment. Being explicit about which stages each skips avoids accidentally re-processing precious already-aligned cfDNA data and introducing realignment artifacts.

### 1a. Already-aligned cfDNA sources (GSE71378) — **skip alignment**

**Correction:** the Cristiano pancreatic subset was moved into this bucket on the assumption it was already-aligned BAM (per phase1-sourcing.md). That's now known to be wrong — only raw FASTQ is on hand for it, at `/home/isilon/HumGenTempData/WGS/Public_data/Cristiano_Velculescu_DELFI/FASTQ/*/<sample_id>_R{1,2}.FASTQ`. Cristiano moves to §1b (raw-FASTQ path) unless a ready BAM turns up elsewhere.

These arrive as BAM, hg19, already on the Institute archive. Stages that still apply:
1. **Verify, don't realign.** Confirm `@SQ` headers match the hg19 reference we'll use for HCC1395 (same contig names/order — `chr1` vs `1` naming mismatches are a common silent breakage point for `samtools merge` later). Confirm proper-pair flags (`0x2`) are set as expected.
2. **Re-index** if `.bai`/`.crai` is missing or stale.
3. Feed directly into dedup (§2), fragment extraction (§3), binning (§4), feature matrix (§5).
4. **Resolved:** the healthy (BH01/IA0x/IH0x) samples are present alongside the cancer subset — the `.../Cancer` directory name in the archive layout is just an example path, not a scope restriction. No separate pull from SRP061633 needed.

### 1b. Raw-FASTQ sources (HCC1395/HCC1395BL, Cristiano pancreatic subset) — **full path**

1. **Retrieval:** `prefetch` + `fasterq-dump` against SRP162370 for HCC1395 (per-sample SRR accessions via SRA Run Selector, not a blind full-project pull — 6 centers × tumor+normal is a lot of data to move indiscriminately). Cristiano's FASTQ is already on hand locally (see §1a correction above) — no retrieval step needed for it, just alignment (steps 2-4 below).
2. **In silico fragment-size reshaping** (Phase 0 §3.2.1) — this happens here, on the raw reads, **before** alignment or mixing. This is the single highest-severity risk item in the whole project (Phase 0 §7.1) and is out of scope for this document's tooling detail beyond flagging where it sits in the pipeline; it needs its own design pass.
3. **Alignment:** BWA-MEM against the hg19 reference (§0), matching whatever exact reference FASTA variant (decoy sequences, contig naming) the Snyder/Cristiano BAMs were aligned against — **open item, see end of doc**: we don't yet know the exact reference build variant used to produce those BAMs, and mismatched decoys/contig sets will silently corrupt `samtools merge` downstream.
4. Then joins the same dedup/extraction/binning/feature path as 1a.

### 2. Deduplication — samtools markdup vs. Picard MarkDuplicates

**Decided: samtools markdup**, chosen for speed — no JVM startup/heap overhead, and it operates directly in the `samtools sort` → `markdup` pipe without an intermediate serialization step, which matters at the TB-scale, multi-source volumes here (§8).

The trade-off worth keeping in mind: samtools markdup's per-library/optical-duplicate handling has a shorter track record than Picard's for heterogeneous multi-read-group BAMs, and HCC1395 spans **6 sequencing centers/platforms** (Phase 1) with different tile geometries (e.g. HiSeq2500 vs. NovaSeq optical-duplicate pixel distances differ). Since Phase 0 §3.2.4 requires dedup to run once on the post-merge BAM (never per-source pre-mix), that merged BAM will carry multiple read groups from different platforms in a single pass — so it's worth spot-checking samtools markdup's duplicate-rate output per read-group against a Picard run on a sample subset early on, rather than assuming parity, given the speed choice.

### 3. Fragment size extraction

Picard's `CollectInsertSizeMetrics` gives aggregate QC summary stats but is **not** sufficient on its own — the feature set needs per-fragment insert sizes binned into genomic windows, which requires per-read `TLEN` extraction, not an aggregate histogram.

**Tooling:** `pysam` (or `samtools view` piped through a small script) over reads passing `-f 0x2 -F 0x400+0x100+0x800` (properly paired, not duplicate, not secondary/supplementary), taking `abs(TLEN)` from read1 only to avoid double-counting each pair. Fragment midpoints intersected against genomic windows via `bedtools intersect` (5 Mb DELFI windows and 1 Mb ichorCNA bins, both derived from the same BAM pass to avoid redundant I/O over multi-TB inputs).

### 4. GC-bias correction and genome binning — what to reuse vs. reimplement

- **1 Mb CN binning + GC/mappability LOESS correction + centromere/gap exclusion:** reuse **HMMcopy Suite** (`readCounter`, `gcCounter`, `mapCounter`) and **ichorCNA**'s own R pipeline (`runIchorCNA.R`) directly, invoked as external tools/subprocesses rather than reimplemented. Both are actively maintained, meant to be run as-is, and this is exactly the CNV/secondary-feature path the spec wants (Phase 0 §2, ichorCNA-style bins).
- **License note:** HMMcopy Suite and ichorCNA are GPL-licensed (Broad Institute). Invoking them as external subprocesses (not statically linking their code into our own package) doesn't trigger GPL's viral distribution requirements for an internal research pipeline. This only becomes a real constraint if the org later wants to distribute this pipeline as closed-source software — flagging so it's a known future consideration, not a blocker now.
- **5 Mb DELFI-style short/long window ratios + their own GC-correction LOESS:** **reimplement**, don't reuse. Unlike ichorCNA, DELFI's 2019 paper code isn't published as a general-use, callable library — there's no equivalent off-the-shelf package to invoke. Implement the LOESS GC-correction (span = 3/4, per Phase 0 §1.2) directly in Python (`statsmodels.api.nonparametric.lowess`) against the same read-level extraction from §3.

### 5. Feature matrix construction

Matches Phase 0 §2 exactly — see §7 (output contract) for the literal schema. No feature is added or dropped from what the spec already specified.

---

## 6. Batch effect risks — concrete combinations

Given the actual sources in hand:

- **GSE71378 (Snyder, single lab/platform) mixed with HCC1395 (6 different centers/platforms).** Fragment-size reshaping (§1b.2) targets the size *distribution* specifically, but library-prep signatures beyond size — GC-bias curve shape, sequencer error/quality profiles — may survive reshaping and correlate with tumor fraction simply because tumor fraction *is* "how much HCC1395 signal is present" in a given mixture. This is the sharpest version of the Phase 0 §3.2.3 risk.
- **Cristiano DELFI pancreatic BAM (different center than Snyder) used as a second healthy/cancer-plasma background alongside GSE71378** — introduces a second, independent center effect if both are pooled into one training set.
- **1000G NYGC (GRCh38, blood gDNA, NovaSeq)** — per Phase 0/1, this is a negative-control-only source and should **never** be mixed into the training set at any nonzero tumor fraction. The concrete risk is subtler: if it's ever used as an additional "0% tumor fraction" healthy control alongside real cfDNA 0% controls, a model could learn "blood-gDNA library prep vs. cfDNA library prep" as a spurious proxy, since that axis would then vary system-wide with class membership. Keep it fully out of the label-bearing training/val/test split; use it only for the standalone QC check below.

**QC check (concrete):** train a decoy classifier (logistic regression or small GBM) to predict **source/platform of origin** from the same feature matrix used for tumor-fraction regression. Two checks, not one:
1. At fixed tumor fraction (especially f=0, pure healthy), can the decoy classifier still discriminate GSE71378 vs. Cristiano vs. HCC1395-derived reads? If yes, platform signal is present in the features.
2. Is the decoy classifier's predicted-source-probability correlated with the true tumor-fraction label across the full training range? A nonzero correlation here is the actual red flag (Phase 0 §3.2.3's failure mode made measurable) — platform-discriminating signal that's *also* label-correlated is what would get learned as spurious tumor-fraction signal.

Run this before any real model training is treated as informative, not as a post-hoc check.

---

## 7. In silico admixture mechanics — concrete steps

Given tumor = HCC1395 (post-reshaping, native coverage ~35–47×) and healthy background = a GSE71378 sample (BH01 pooled or an individual IA0x/IH0x, native coverage varies by sample):

1. **Downsampling ratio calculation.** For a target admixture at tumor fraction `f` and target final coverage `C_target` (ULP-WGS range, matching ichorCNA's 0.1×–1× operating range per Phase 0 §1.1):
   - Tumor-derived target depth = `f × C_target`; downsample fraction for tumor source = `(f × C_target) / C_tumor_native`
   - Healthy-derived target depth = `(1-f) × C_target`; downsample fraction for healthy source = `((1-f) × C_target) / C_healthy_native`
   - Both fractions must be ≤ 1. HCC1395's high native coverage (35–47×) makes the tumor-side fraction small and safe across the full `f` range; the healthy side is more likely to be the binding constraint at low `f` if the chosen GSE71378 sample's native coverage is itself low — check per-sample coverage before picking mixing partners.
2. **Subsampling tool:** `samtools view -s <seed>.<fraction>` per source BAM. Note this is probabilistic, not exact-count, subsampling — after subsampling, verify actual retained depth (via `samtools flagstat` or `mosdepth`) matches the computed target before merging, especially at small fractions where variance is proportionally larger.
3. **Merge:** `samtools merge <mixed>.bam <tumor_subsampled>.bam <healthy_subsampled>.bam`, then sort, then run dedup **once** on the merged file (§2) — per Phase 0 §3.2.4, never merge already-deduplicated BAMs.
4. **Verify empirical tumor fraction against known somatic SNVs.** HCC1395 has a SEQC2 high-confidence somatic call set (~40k SNVs) with known VAF in the pure tumor sample. After mixing, at loci that are genuinely tumor-somatic (absent from HCC1395BL matched normal), the healthy background contributes only reference alleles, so expected VAF in the merged BAM ≈ `f × VAF_tumor_pure`. Use `bcftools mpileup` at those known positions on the merged BAM, compute observed alt-allele frequency, and compare to the expected value.
   - **This addresses Phase 0 §3.2.5 (tumor purity of the source itself)** directly: rather than trusting the nominal downsample-ratio `f` as ground truth, use the **empirical SNV-VAF-derived fraction** as the actual regression label. This absorbs any gap between nominal mixing ratio and true tumor content of the HCC1395 source, so the label reflects genomic reality rather than an assumption of 100%-pure input.

---

## 8. Tooling / environment

- **Alignment/processing:** BWA-MEM, samtools (≥1.15 recommended), Picard, bedtools, SRA toolkit (`prefetch`/`fasterq-dump`)
- **CNV path:** HMMcopy Suite (readCounter/gcCounter/mapCounter), ichorCNA (R)
- **Feature engineering:** Python — `pysam`, `pandas`, `numpy`, `statsmodels` (LOESS)
- **Modeling (Phase 3, not this phase):** XGBoost/LightGBM per Phase 0 §4

**Compute — the workstation runs everything downstream of a BAM.**

Once a BAM exists, every remaining stage for every source — Snyder, Cristiano, HCC1395 — runs locally on the workstation (36 cores, 125 GB RAM, 8 TB HDD). Nothing downstream of alignment runs on the cluster. Cluster scratch cannot hold these BAMs long-term anyway, so this is forced, not merely permitted.

**The one boundary: FASTQ→BAM alignment stays on the HPC** (§1b, HCC1395 + Cristiano). Not a capability limit — the workstation could do it — but the FASTQ lives on isilon, alignment is where the compute genuinely belongs, and moving hundreds of GB of raw reads down a 92 MB/s link to align them locally costs more than it saves. §8c.6 prices the alternative if you want it anyway.

An earlier version of this section said the phase "cannot run on a laptop … not a personal workstation." That was written assuming everything reads raw FASTQ. It is wrong for §1a, which never touches FASTQ, and over-general for the rest.

---

## 8b. Local workstation execution plan — all sources

### 8b.1 Data location

**Superseded by the CRAM-era implementation.** The cluster-side alignment
rules now convert to CRAM and verify it before the workstation ever sees the
sample, so nothing BAM-shaped is mirrored here
any more. The table below reflects what actually ships in
`workflow/config/local_features.yaml`; see `docs/local-features-pipeline.md`
for the full pipeline walkthrough.

| What | Path |
|---|---|
| Mirrored CRAMs (all sources) | `/mnt/8tb_hdd/flash_fraq_crams/<source_dir_name>/CRAM/<sample>/<sample>.sorted.cram(.crai)` |
| Transfer script (re-runnable; skips what is already local) | `scripts/sync_crams_to_hdd.sh` |
| hg19 reference | `/mnt/8tb_hdd/references/genome/hg19/hg19.fa` (+ `.fai`) |
| Derived outputs | `/mnt/8tb_hdd/flash_fraq_derived/<source_dir_name>/<sample>/` |
| SSD staging scratch (per-sample, purged via `temp()` once consumed) | `/mnt/data/flash_fraq_stage/` |

There is no separate BAM sort/scratch path to manage: nothing is decoded to
BAM except the one `temp()` BAM `localfeat_cram_to_bam` produces for
`readCounter`'s benefit, and it never leaves `/mnt/data`.

### 8b.2 Environments already present

| Need | Where |
|---|---|
| samtools, bwa, bedtools, pysam, pandas, numpy, statsmodels | conda env `flash_fraq` |
| HMMcopy Suite — `readCounter`, `gcCounter`, `mapCounter` | conda env `TSO500ichorCNA` |
| ichorCNA R package | conda env `problem2` (`lib/R/library/ichorCNA`) |
| mosdepth | conda env `mosdepth` |

The HMMcopy binaries and the ichorCNA R package are in **different envs** — consolidate into one before §8b.3 step 6. Prior working 1 Mb hg19 ichorCNA driver: `/mnt/8tb_hdd/tso500offtarget_ichorcna/`. Reuse the invocation pattern; its panel-of-normals is TSO500 off-target panel data, not cfDNA WGS, so do not reuse the PoN.

### 8b.3 Stage order

Per source BAM:

1. **Header/reference verification** — `samtools view -H`; assert `@SQ` names/order match `hg19.fa.fai`. Index-only, seconds per file. Closes the reference-variant open item and gates everything else, so run across all samples first.
2. **Index freshness** — rsync preserves source mtimes, so transferred `.bai` files can be **older than their BAM**, and samtools emits `[W::hts_idx_load3] The index file is older than the data file` (already observed on the Snyder mirror). Cosmetic, but `touch` the `.bai` after transfer to keep real staleness visible.
3. **Native coverage** — `mosdepth --by 1000000 --no-per-base`. Closes the BH01-coverage open item and feeds the §7.1 downsampling fractions.
4. **Duplicate-rate characterization** — `samtools markdup` on source BAMs is for QC only. Per Phase 0 §3.2.4 the dedup that counts runs **once on each merged mixture**, which is tiny (§8c.4).
5. **Fragment extraction** — one `pysam` pass emitting per-fragment `abs(TLEN)` + midpoint.
6. **Binning + GC correction** — HMMcopy `readCounter` at 1 Mb; own LOESS for 5 Mb DELFI windows.
7. **Feature matrix** — §9 schema.

**Superseded: steps 3, 5 and 6 no longer share a single pass.** This
constraint held for a direct HDD read, where three separate sequential scans
would triple the dominant cost in §8c. It was retired once CRAMs are staged
onto NVMe SSD first (one HDD read per sample, then arbitrarily many free SSD
reads) — at that point the reason to fuse the stages disappears, and the
implementation runs them as separate, independently-restartable Snakemake
rules instead. See `docs/local-features-pipeline.md` §1 for the full
before/after reasoning.

### 8b.4 Concurrency — the HDD is the only real constraint

Measured on this machine:

| Quantity | Measured |
|---|---|
| BAM decode, 1 thread | 212 MB/s |
| BAM decode, `-@ 4` / `-@ 8` | 933 / 1001 MB/s |
| `pysam` Python-level iteration | 1.41 M reads/s per process |
| Read density | ~25 M reads per GB of BAM |
| HDD sequential, saturated | ~131 MB/s at 98% `%util` |
| HDD O_DIRECT read **while a transfer runs** | 27.5 MB/s (1 reader), 52 MB/s (4 readers) |

Consequences, and they are not the intuitive ones:

- **CPU decode is never the bottleneck.** At 4 threads samtools decodes ~7× faster than the disk delivers.
- **A single `pysam` worker is slower than the disk**: 1.41 M reads/s ÷ 25 M reads/GB ≈ **56 MB/s**. One worker leaves the drive two-thirds idle. Run **3–4 workers** to saturate ~150 MB/s and no more — past that they contend for the same spindle.
- **Never run analysis while a transfer is running.** Read throughput collapsed from ~150 to 27.5 MB/s under concurrent rsync — a 5× penalty. Finish transfers first.
- Give samtools `-@ 4` per process rather than raising sample-level concurrency.

---

## 8c. Estimated runtime

All figures assume **~150 MB/s effective sequential HDD throughput** with no competing transfer, derived from the 131 MB/s write-saturated measurement above. This single constant dominates everything; if the drive benchmarks differently once transfers are idle, rescale linearly.

### 8c.1 Cohort volume

| Source | Samples | BAM volume |
|---|---|---|
| Snyder / GSE71378 | ~60 (59 confirmed) | **~950 GB** (927.5 GB measured for 59) |
| Cristiano DELFI | ~211 | **~2.94 TB** (measured, 423 files) |
| HCC1395 + matched normal | ~12 (6 centers × T/N) | **~1.2 TB** (estimated at 35–47×) |
| **Total** | | **~5.1 TB** |

**One full sequential pass over everything = 5.1 TB ÷ 150 MB/s ≈ 9.5 h.** That is the unit of cost for every stage below. Snyder alone is ~1.8 h per pass.

### 8c.2 Per-stage estimates

| Stage | Passes | Snyder only | Full 5.1 TB |
|---|---|---|---|
| 1. Header verification | 0 (header only) | < 1 min | ~2 min |
| 2. Index freshness / re-index | 0, or 1 if rebuilding | seconds | ~9.5 h *if* indexes must be rebuilt — they should not be |
| 3+5+6. **Combined coverage + fragment + binning pass** | **1** | **~2 h** | **~9.5–11 h** |
| 4. Duplicate-rate QC (subset of samples, not all) | ~3× I/O on those chosen | ~1 h for 5 samples | run on a sample, not the cohort |
| 7. Feature matrix assembly | 0 (operates on extracted tables) | minutes | ~30 min |

Stages 3/5/6 split into three separate passes instead of one: **~29 h instead of ~10 h** for the full cohort. This is the single biggest runtime decision in the phase.

### 8c.3 The admixture trap — read sources once, not once per mixture

The naive implementation of §7 runs `samtools view -s` against the source BAMs for every mixture. **This does not scale, and the failure is severe:** BH01 is 208.9 GB, so 1,000 mixtures drawing from it means 209 TB of reads — **~16 days** on this drive, for the healthy side alone.

Instead: **one pass per source BAM, emitting a reusable read pool.** Mixtures target 0.1×–1× ULP-WGS (Phase 0 §3.2.2), i.e. ~0.3–3 GB per mixed BAM against a ~100 GB source — three orders of magnitude smaller. So extract, once per source, a subsampled pool at the *deepest* fraction any mixture needs, then draw every mixture from that pool. Cost drops from 209 TB to one 9.5 h cohort pass, and every subsequent draw is a few seconds against a file that fits in page cache.

Design this before writing mixture code. It is far more expensive to retrofit.

### 8c.4 Per-mixture cost, once pools exist

| Step | Per mixture |
|---|---|
| Subsample + merge from pools | seconds |
| `samtools sort` + `markdup` on the merged 0.3–3 GB BAM | ~1–2 min |
| `readCounter` 1 Mb | ~1 min |
| `runIchorCNA.R` | ~2–5 min |
| DELFI 5 Mb fragment features | seconds |
| **Total, serial** | **~5–8 min** |

These are small files and CPU-bound, so they parallelize well across 36 cores — unlike the source-BAM passes. At ~12 concurrent, **1,000 mixtures ≈ 7–11 h**.

### 8c.5 End-to-end

| Block | Wall clock |
|---|---|
| Transfer Cristiano (2.94 TB @ 92 MB/s, 4 streams) | ~9 h |
| Transfer HCC1395 BAMs after HPC alignment | ~4 h |
| Header/index verification, all sources | ~15 min |
| Combined coverage + fragment + binning pass, 5.1 TB | ~10 h |
| Source read-pool extraction (§8c.3) | folded into the pass above |
| 1,000 mixtures: merge → dedup → ichorCNA → features | ~7–11 h |
| Feature matrix assembly | ~30 min |
| **Total, transfers excluded** | **~18–22 h of workstation compute** |
| **Total, including transfers, nothing overlapped** | **~31–35 h** |

Call it **two to three days of elapsed time** with realistic restarts and re-runs, of which roughly half is transfer. Transfers and analysis must not overlap (§8b.4).

**Snyder-only, available today:** ~2 h for the combined pass, plus minutes of verification. Everything needed to close both remaining open items is a single afternoon.

### 8c.6 If alignment moves here too

Priced for completeness, not recommended. BWA-MEM on 36 cores runs roughly 10–14 h per 35–47× WGS sample; 12 HCC1395 BAMs ≈ **6–8 days**, plus several hundred GB of FASTQ pulled down first, plus sort/merge scratch on the same HDD already established as the bottleneck. The HPC does this in parallel across nodes next to the data. Keep §1b on the cluster.

---

## 9. Output contract for Phase 3

One row per **mixed admixture sample**:

| Column | Notes |
|---|---|
| `sample_id` | unique per mixed instance |
| `tumor_source_id`, `healthy_source_id` | which raw BAMs contributed |
| `platform_tumor`, `platform_healthy` | sequencing center/platform — required for §6 stratification and decoy-classifier QC |
| `reference_build` | tag, expected `hg19` for all admixture-pair rows |
| `target_coverage` | nominal ULP-WGS depth mixed to |
| `nominal_tumor_fraction` | intended value from the downsample-ratio calculation (§7.1) |
| `empirical_tumor_fraction` | SNV-VAF-derived actual value (§7.4) — **this is the regression label**, not the nominal value |
| `short_long_ratio_window_1..N` | GC-corrected, 5 Mb DELFI windows post-blacklist/correlation filtering |
| `gc_corrected_total_coverage_window_1..N`, `gc_corrected_short_coverage_window_1..N` | per window |
| `arm_zscore_<arm>` | chromosome-arm fragment z-scores, one column per arm label (`arm_zscore_chr1p`). An earlier draft of this table said `arm_zscore_1..39`: the count is not fixed and is not 39 — a sample carries 44, 45 or 46 of the 46 primary-contig arms (chr15p has no 5 Mb windows anywhere; chr22p is NaN in every sample), so only the label is comparable across samples |
| `fragment_size_median`, `fragment_size_mode`, `fragment_size_entropy` | genome-wide summary stats |
| `cnv_logr_bin_1..T` | 1 Mb ichorCNA-style bins, GC/mappability-corrected |
| `ichorCNA_tumor_fraction_estimate` | ablation-only feature per Phase 0 §2 — must be clearly tagged so it's excluded from the primary model by default |
| `mtdna_fraction_log10` | auxiliary |
| `split_group` | source-subject-level train/val/test tag, assigned **before** mixing, per Phase 0 §5.2 |

Window/bin counts (`N`, `T`) are placeholders pending the actual blacklist/correlation-filtering pass — exact counts depend on which windows survive the >0.9-correlation removal step from the DELFI methodology.

---

## Decisions made this round

1. **Reference build:** hg19 across the board, confirmed (§0). Only remaining sub-item: verify the exact FASTA variant (decoy set, contig naming) behind the existing Snyder/Cristiano BAMs before aligning HCC1395, so `samtools merge` doesn't break on contig mismatches — a mechanical verification step, not a decision.
2. **Dedup tool:** samtools markdup, for speed (§2).
3. **GSE71378 healthy-sample availability:** confirmed present; `/Cancer` was just an example path, not a scope limit.
4. **Which GSE71378 sample to use as the healthy background:** recommend **BH01** (the pooled healthy plasma sample) as the primary background over an individual IA0x/IH0x sample. It's the most widely cited reference sample from the Snyder paper, and pooling generally means higher/more stable coverage than a single individual donor — better suited as the fixed healthy partner across many admixture mixtures. Keep one or two individual IA0x/IH0x samples in reserve as a secondary healthy background, specifically for the held-out-source generalization check in Phase 0 §5.2 (train mixtures use BH01, a held-out test slice uses an individual sample, to confirm the model isn't just keying off BH01-specific artifacts).
5. **FinaleDB / full DELFI-cohort DAC path:** skip both for v1. Phase 1's "usable now" sources (GSE71378 + Cristiano pancreatic + HCC1395) already cover the design without any DAC turnaround risk or FinaleDB's unresolved governance question. Revisit only if v1 results indicate the model needs broader cancer-type diversity than the pancreatic subset provides.

## Still open

- **Reference FASTA variant confirmation** for HCC1395 alignment (mechanical check, not a judgment call — do before running BWA-MEM). Resolved by §8b.3 step 1, run locally against the Snyder BAM headers.
- **Native coverage of BH01** (and the reserved IA0x/IH0x samples) — needed to compute the concrete downsampling fractions in §7.1 once actual files are inspected. Unblocked: BH01 is in transfer to the workstation; measure via §8b.3 step 3.
