# The local (workstation) feature-extraction pipeline

A walkthrough of the Snakemake workflow that turns the mirrored hg19 CRAMs into
the per-sample feature tables `feature_matrix.py` joins. This is the CRAM-era
implementation of `phase2-preprocessing-plan.md` §8b.3, and it runs entirely on
the workstation — the cluster's job ends at FASTQ→BAM alignment.

---

## 1. The one idea worth understanding first

The plan (§8b.4) insisted that stages 3, 5 and 6 — coverage, fragment
extraction, CN binning — share **one sequential pass** over each file. That was
not an aesthetic preference. It came from a measurement: the 8 TB HDD sustains
~150 MB/s when nothing else touches it, and collapses to ~27.5 MB/s while a
transfer runs. Three separate passes over a spinning disk projected to ~29 h for
the cohort against ~10 h for a fused pass.

Staging each CRAM onto the NVMe SSD first **inverts that trade-off**. You pay
exactly one HDD-sequential read per sample; every pass after that reads
`/mnt/data`, where repeated reads are effectively free. Once the expensive
resource is no longer the thing being conserved, the reason to fuse the stages
disappears.

So the pipeline deliberately does the opposite of what §8b.3 prescribes: the
rules stay **separate**. That buys back everything the fused design had to give
up — each stage is independently restartable, independently cacheable, and a
failure in CN binning doesn't cost you the fragment pass you already paid for.

That is the single decision the rest of the file follows from. §8b.3's stage
*order* is preserved exactly; only its single-pass *constraint* is retired,
because the premise it rested on no longer holds.

---

## 2. Shape of the DAG

```
localfeat_make_windows ─────────────┐   (once, cohort-wide)
                                    │
localfeat_stage_cram  HDD → SSD  temp(cram, crai)
  ├─ localfeat_verify_cram ── cram_header_verification.json
  │     ├─ localfeat_coverage ────── mosdepth.summary.txt
  │     ├─ localfeat_fragments ◄─────┘ window_fragment_counts.tsv
  │     │      └─ localfeat_gc_correct ─ window_ratios.csv, arm_zscores.csv
  │     └─ localfeat_cram_to_bam  temp(bam, bai)
  │            └─ localfeat_ichor ──── ichor/cnv_bins.csv, ichor_estimates.csv
  └─ localfeat_sample_complete   ← join point; frees the SSD copies
                                   
local_features_complete            ← cohort target
```

Nine rules. Everything except `localfeat_make_windows` and the two completion
rules is per-`(dataset, sample_id)`.

---

## 3. Stage in, analyse, purge

### `localfeat_stage_cram` — the only HDD read

Copies the CRAM and its `.crai` from the mirror to `/mnt/data`. Three details
that aren't obvious:

- **`cp`, not `rsync`.** The source is a local mount, and a partial staged copy
  must never be resumed — Snakemake deletes and re-runs its outputs on failure,
  and a resumable partial would defeat that.
- **`--no-preserve=mode`**, then **`touch` on the `.crai`.** rsync preserved the
  cluster's mtimes on the mirror, so a mirrored index is usually *older* than
  its CRAM and samtools warns "index file is older than data file" on every
  single read. Stamping the staged index keeps that warning meaningful rather
  than constant background noise.
- Both outputs are wrapped in **`temp()`**.

### `temp()` is the purge mechanism

Snakemake deletes a `temp()` output as soon as every job consuming it has
finished. For the staged CRAM that means once `localfeat_coverage`,
`localfeat_fragments` and `localfeat_cram_to_bam` have all run — no extra rule
is needed to trigger it. `localfeat_sample_complete` does *not* take the staged
files as inputs, so it neither holds them open nor releases them; what it
provides is a single per-sample target for `local_features_complete` to depend
on, instead of enumerating seven output paths per sample.

Consequence worth internalising: **peak SSD usage tracks the number of
concurrent samples, not the cohort size.** You never need 118 CRAMs on the SSD.

### `ssd_mb` — the concurrency budget

A custom resource declared by the two rules that write to the SSD:

```python
ssd_mb=lambda wc: local_cram_size_mb(LOCALFEAT_CRAM_ROOT, wc.dataset, wc.sample_id)
# cram_to_bam uses multiplier=3.0 — a decoded BAM runs ~2.5–3× its CRAM
```

with `resources: ssd_mb: 400000` in the local profile. It is **accounting, not
reservation**: Snakemake simply won't start a job whose declared footprint would
push the running total past 400 GB.

The subtle part is *which file gets stat'ed*. `local_cram_size_mb()` measures
the **HDD mirror**, never the staged copy. Sizing off the staged path was my
first attempt and it was silently wrong — the staged file doesn't exist when the
DAG is built, so `input.size_mb` evaluates to `0`, every job declares a
zero-byte footprint, and the budget quietly stops constraining anything. A
disabled limit that still looks present in the config is worse than no limit.

> **Related landmine, now fixed:** the local profile briefly had two separate
> `resources:` blocks. YAML has no key merging — the second silently *replaced*
> the first, dropping `ssd_mb` entirely. Same failure signature: a limit that
> reads as present and enforces nothing. Everything lives in one mapping now.

---

## 4. The per-sample stages

### `localfeat_make_windows` (cohort-wide, runs once)

Generates the DELFI 5 Mb windows and ichorCNA 1 Mb bins from `hg19.fa.fai`.
**Neither BED existed on this machine** — `fragment_extraction.py` has always
taken them as inputs, and only the test suite ever synthesised them.

Two constraints baked into `scripts/make_genome_windows.py`:

- Windows emit in **caller-specified contig order**, not `.fai` order. `hg19.fa.fai`
  is sorted by size (chrX lands between chr7 and chr8), so trusting file order
  would interleave windows and make `window_index` non-monotonic along the
  genome.
- Names are `window<n>` / `bin<n>`, because `FragmentExtractor._parse_window_index()`
  reads the trailing integer off the BED name column.

Primary contigs only — `hg19.fa.fai` also carries `*_hap*`/`*_random`/`chrUn`
entries that have no windows or bins and would otherwise leak in.

Measured output: **621 DELFI windows, 3053 ichor bins.**

### `localfeat_verify_cram` (§8b.3 step 1)

Two checks for two unrelated failure modes:

1. **`samtools quickcheck`** — the truncation gate. `sync_crams_to_hdd.sh` runs
   rsync with `--partial` and no `--partial-dir`, so an interrupted transfer
   leaves partial bytes **at the final filename**. Such a file opens fine and
   decodes for a while before dying mid-stream. Without this gate that surfaces
   as a baffling failure hours into a fragment pass.
2. **@SQ vs. the reference `.fai`** — closes the phase2 "reference FASTA
   variant" open item. A matching build *name* is not enough: HCC1395 must be
   aligned against the identical contig set or the `samtools merge` inside
   `admixture.py` will reject the mix at admixture time. **Order** is compared,
   not just membership, because @SQ order defines the reference-id encoding —
   and as *relative* order over the contigs the CRAM actually carries. A flat
   `cram_names == fai_names` was the first version and it was useless: any CRAM
   listing a legitimate subset of the reference fails it regardless of
   ordering, so the flag was either trivially true or falsely alarming.

Asymmetric severity, deliberately: a contig in the reference but absent from the
CRAM is normal (the header lists only what the aligner emitted). The reverse is
fatal, as is any length mismatch — that means a different assembly patch wearing
the same name — and so is a reordering, which is now enforced rather than merely
recorded in the report.

Result on IA01: **93 contigs, identical order, no length mismatches.**

### `localfeat_coverage` (step 3)

`mosdepth` at 1 Mb, `--no-per-base`. Feeds the §7.1 downsampling fractions and
closes the open BH01-coverage question. Runs in its own `mosdepth` conda env —
hence `--use-conda` is **required**, not optional, on the run command.

### `localfeat_fragments` (step 5, MAT-14)

The pysam TLEN/midpoint pass. The CRAM fix that made this work at all:

```python
pysam.AlignmentFile(path, "rb", reference_filename=str(reference_fasta))
```

htslib needs an explicit reference for CRAM because the @SQ `UR:` tags point at
a **cluster** path that doesn't exist locally.

`fragment_midpoints.bed` and the two intersected BEDs are `temp()` — one line
per fragment is tens of GB on a deep sample, and nothing downstream reads them.

That size is also why `compute_fragment_stats()` streams the file into a TLEN
**histogram** rather than a list. It used to do `midpoints_bed.read_text()` and
build one Python `int` per fragment, which on a deep sample is hundreds of
millions of objects — far past the rule's 16 GB request, and it would have
surfaced as an OOM kill on the first sample deeper than IA01. Median, mode and
entropy are all recoverable from a histogram bounded by the ~10³ distinct
fragment lengths.

Scaling note: on the HDD a single pysam worker managed ~56 MB/s against a
~150 MB/s spindle, so §8b.4 wanted 3–4 workers to saturate it. Reading from SSD
makes each worker **CPU-bound instead of IO-bound** — so scale this with
`jobs`, not with per-job threads.

IA01 result: **median 153 bp, mode 166 bp, entropy 7.796** — a textbook cfDNA
mononucleosome profile, which is the strongest single signal that the CRAM path
is decoding correctly.

### `localfeat_gc_correct` (step 6, DELFI half, MAT-15)

LOESS GC correction (statsmodels, span 0.75) over the small per-window table.
The one stage that never opens an alignment, so it has no staging or IO concern
at all.

### `localfeat_cram_to_bam` + `localfeat_ichor` (step 6, CN half, MAT-16)

`readCounter` is BamTools-based and **cannot open a CRAM** — it dies with
`BgzfStream ERROR: read block failed - invalid block header`. §4 forbids
reimplementing HMMcopy logic (it stays an external subprocess, which is also
what keeps the GPL boundary clean), so the workflow decodes **one temp BAM on
the SSD for readCounter's sole benefit**. Every other consumer reads the CRAM
directly, and the BAM is `temp()`, so it exists only for the lifetime of this
one rule.

`cnv_bins.csv`'s `bin_index` keys on **genomic position**, not row order.
ichorCNA drops bins that fail its mappability/centromere/PoN filters, so the
row count of `.cna.seg` varies between samples — §9 below shows the same sample
yielding 3053 and 2594 bins on two runs. Numbering rows 1..N and letting
`feature_matrix.load_cnv_bins_wide()` pivot that into `cnv_logr_bin_<n>`
columns would align *different loci* across samples with nothing raising an
error. The index is now `chromosome_rank * 1000 + bin ordinal`, so a locus one
sample lost is a NaN in the pivot rather than a shifted neighbour, and
`chrom`/`start` ride along in the CSV.

`gcCounter` and `mapCounter` never run. ichorCNA ships hg19 1 Mb GC and
mappability WIGs; gcCounter would regenerate a byte-identical GC file once per
sample, and **mapCounter cannot run at all** — it needs a mappability bigWig
that does not exist on this machine (confirmed by a full-filesystem search).
Only `readCounter`, which is genuinely per-sample, is required.

---

## 5. Three ichorCNA traps, and why the config looks the way it does

These cost the most time, so they're worth reading before touching
`local_features.yaml`.

**`--libdir` is not an optimisation.** `runIchorCNA.R` on disk is **0.3.2
(2019)**. Without `--libdir` it calls `library(ichorCNA)` and loads the
conda-installed **0.5.1** — a 2019 driver against a 2024 library whose `runEM`
array contract changed. It fails with `incorrect number of subscripts` partway
through the fit. `--libdir` pins the script to its own sources.

I want to flag how I got here, because the first diagnosis was wrong: I
hypothesised `--scStates` defaulting to `"NULL"`, added the README's ULP-WGS
parameters, and **it did not fix it**. Only then did I stop guessing and diff
versions. (The scStates parameters were independently necessary — that pairing
with `estimateScPrevalence=TRUE` produces the *same* error message — so both
fixes are in, but only the version pin resolved the actual failure.)

**`--genomeStyle UCSC` is mandatory.** These CRAMs are chr-prefixed hg19. The
script defaults to NCBI naming, which matches **zero** chromosomes and yields an
empty fit rather than an error — a silent wrong answer, the worst failure mode
available.

**The panel of normals exists, and I initially said it didn't.** I set
`normal_panel: null` after checking only the conda 0.5.1 extdata (which carries
just a 500 kb hg19 and a 1 Mb hg38). The 0.3.2 checkouts ship
`HD_ULP_PoN_1Mb_median_normAutosome_mapScoreFiltered_median.rds` — correct build,
correct bin size, healthy-donor ULP **cfDNA WGS**. It's wired in now. This
matters rather than being cosmetic: the PoN is the healthy baseline ichorCNA
corrects against, and without it low-tumour-fraction estimates are noticeably
noisier — which is exactly this project's operating regime.

To be explicit, since §8b.2 warns about this: `HD_ULP` is **not** the TSO500
off-target PoN in `/mnt/8tb_hdd/tso500offtarget_ichorcna`. That one is built
from hybrid-capture panel off-target reads. This one is cfDNA WGS from healthy
donors, which is what these samples are.

---

## 6. What keeps a still-transferring sample out of the DAG

`discover_local_crams()` only enqueues a sample whose **`.crai` exists
alongside** the CRAM. Because rsync writes partial bytes at the final filename,
globbing `*.sorted.cram` alone would happily enqueue a truncated file.

Two layers, because one isn't enough: rsync orders transfers by size, not by
pair, so an index can occasionally land *before* its CRAM finishes. That
remaining window is what `samtools quickcheck` in `localfeat_verify_cram` is
there to catch.

The sample list is resolved **once at DAG construction**, not per job — the
mirror is append-only while a sync runs, and a mid-run re-glob would grow the
DAG underneath the scheduler.

---

## 7. Conventions you'll trip over otherwise

- **`LOCALFEAT_` prefixes on every global.** Snakemake's `include:` shares one
  Python namespace across all rule files, so an unprefixed name is silently
  reassigned by whichever file is included last. `MAT19_`/`MAT20_` follow the
  same rule.
- **Helper functions live in `common.smk`.** `--lint` complains about function
  definitions inside a rules file.
- **`rule all` deliberately excludes the local target.** It covers cluster-side
  alignment only. Folding the local target in would make `snakemake` on either
  machine try to build targets whose inputs live on the other one. Ask for it by
  name.
- **`keep-going: true`** in the local profile: one truncated CRAM would
  otherwise abort the cohort, revealing bad samples one run at a time.

---

## 8. Running it

Pre-flight checks are mandatory per CLAUDE.md:

```bash
snakemake --lint -s workflow/Snakefile
snakemake -n -p --use-conda -s workflow/Snakefile \
  --workflow-profile workflow/profiles/local local_features_complete
```

Then:

```bash
snakemake --use-conda -s workflow/Snakefile \
  --workflow-profile workflow/profiles/local local_features_complete
```

`--use-conda` is required — `localfeat_coverage` declares `conda: "mosdepth"`.

Note that `snakemake` lives in **its own conda env**, not in `flash_fraq`. Run
it with that env's `PATH`, and export `PYTHONPATH=$PWD/src` so the drivers can
import `flash_fraq`.

**Timing caveat:** hold off while a CRAM sync is still running. §8b.4
measured HDD read throughput collapsing from ~150 to 27.5 MB/s under concurrent
rsync, and `localfeat_stage_cram` hits the same spindle. The staging advantage
largely evaporates if you fight the transfer for it.

---

## 9. Validation status

Smoke-tested end-to-end on real data, sample IA01, before the cohort run:

| Stage | Result |
|---|---|
| Header verification | 93 contigs, identical order, no length mismatches |
| Window generation | 621 DELFI windows, 3053 ichor bins |
| Fragment extraction | median 153, mode 166, entropy 7.796 |
| GC correction | `window_ratios.csv` + `arm_zscores.csv` produced |
| ichorCNA (no PoN) | TF 0, ploidy 1.997 |
| ichorCNA (with PoN) | TF 0.0, ploidy 2.000, 2594 bins |

Tumour fraction 0 is the **correct** answer for a healthy plasma sample; the
ploidy tightening toward exactly 2 is the direction you'd expect from adding a
healthy baseline correction.

Test suite: 47 passing, none failing.

**Cohort run, past the IA01 smoke test:** `local_features_complete` has since
been run against every CRAM mirrored to the workstation at the time —
`/mnt/8tb_hdd/flash_fraq_derived/` holds a `localfeat_sample_complete.flag`
plus the full stage-3/5/6 output set (coverage, fragment_stats,
window_ratios/arm_zscores, cnv_bins/ichor) for:

- **Cristiano_Velculescu_DELFI: 503/503** mirrored samples complete.
- **Snyder_Shendure_WPS: 44/44** mirrored samples complete (63 total are
  aligned on the cluster; the remaining ~19 haven't been CRAM-migrated to the
  workstation yet).

Every present sample directory has its complete flag — `keep-going: true`
didn't silently drop anything from this run.

Caveat on the two rows above: they were measured before `bin_index` was
switched from row order to genomic position. The fit itself is unchanged — only
the labels in `cnv_bins.csv` — but the column names in any feature matrix built
from the old CSVs will not match one built now, so regenerate rather than mix.

One diagnostic worth recording: a chr21-only smoke test made ichorCNA fail with
`span is too small`. That is an artefact of a single-chromosome synthetic file —
every other chromosome empty, too few bins for LOESS — not a pipeline defect. I
escalated to a full-sample test rather than "fixing" it, which is what surfaced
the real `--libdir` bug.

---

## 10. Still open

Three items, all awaiting your call rather than blocked on anything:

1. **`docs/phase2-preprocessing-plan.md` is stale in two places** — §8b.1 still
   lists BAM paths, and §8b.3 still asserts the single-pass requirement that §1
   above retires. I left the file alone because you have uncommitted edits open
   in it.
2. **`ssd_mb` is weaker than §3 implies.** Snakemake holds a resource for the
   duration of the *job*, not the lifetime of the file it produced. The budget
   is released the moment `cp` returns, while the staged CRAM sits on the SSD
   until its last consumer finishes — so what actually bounds peak SSD usage is
   `jobs`, and `ssd_mb` only prevents a burst of simultaneous staging. Real
   enforcement would mean declaring `ssd_mb` on every rule that reads a staged
   file. Worth doing before raising `jobs` much above 4.
3. **Windows with zero fragments are dropped.** `aggregate_window_counts()`
   only emits windows that some fragment intersected, so a sample's table can
   be shorter than 621 rows. `feature_matrix` pivots tolerate this as NaN, but
   it does mean per-sample tables are not a fixed width.
