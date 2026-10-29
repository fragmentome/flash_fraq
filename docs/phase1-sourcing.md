# Phase 1: Data Sourcing Verification

**Status:** Research only — no downloads, no scripts executed. Confirms sourcability, access tier, and links against [phase0-spec.md](phase0-spec.md) §3.3 before any data pull is attempted.

---

## Sourcing table

| Source | Exists / Reachable | Access tier | What's actually in it | Accession / link | Verdict |
|---|---|---|---|---|---|
| **GSE71378** (Snyder et al. 2016) | Yes, live | Fully open, no registration | 60 samples: 1 pooled healthy (BH01) + 8 individual healthy (IA01–08) + 3 more healthy (IH01–03) + 51 lung cancer plasma (IC02–52). Raw FASTQ for all samples via linked SRA project SRP061633. Only 5 samples have processed BigBed tracks — everything else is raw reads only. **Already downloaded locally**: `/home/isilon/humangenetik/ARCHIVE/WGS/Public_data/Snyder_Shendure_WPS/BAM/original_data_hg19/Cancer` — note this path name says "Cancer," so confirm whether the healthy (BH01/IA/IH) samples were also pulled down alongside it or only the IC0x cancer subset. | [GSE71378](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE71378) · [SRP061633](https://www.ncbi.nlm.nih.gov/sra/?term=SRP061633) | **Already on hand** — verify healthy-sample coverage before assuming the full cohort is present |
| **FinaleDB** | Site returned a TLS handshake error on direct fetch (`SSL...TLSV1_ALERT_INTERNAL_ERROR`) — could not independently confirm current liveness of `finaledb.research.cchmc.org` itself. GitHub repos (`epifluidlab/finaledb_portal`, `finaledb_workflow`) exist and are the app's source. | Published as "fully open browse/download" of de-identified reprocessed data; published paper (Zhou et al. 2021, *Bioinformatics*, PMC8388032) confirms downloadable "coverage, fragment size profile, etc." files — i.e. **summary/derived tracks**, not raw BAM/FASTQ. Whether Cristiano/DELFI-derived samples specifically are among the 2,579 reprocessed datasets was not confirmed from available text. | [finaledb.research.cchmc.org](http://finaledb.research.cchmc.org/) · [paper (PMC8388032)](https://pmc.ncbi.nlm.nih.gov/articles/PMC8388032/) · [GitHub](https://github.com/epifluidlab/finaledb_portal) | **Not confirmed — verify manually in a browser first** ⚠️ see governance flag below |
| **EGAS00001003611 / EGAD00001005339** (Cristiano DELFI cohort) | Yes, both confirmed to exist | DAC-controlled | 541 BAM files (4.4 TB), 537 samples, healthy + 7 cancer types, Illumina HiSeq2500. DAC ID `EGAC00001001180`. Use restricted to "Research Purposes" only under Fort Lauderdale Guidelines; no re-identification/commercial use; publication moratorium listed as "0 months" (i.e., no embargo once granted). No approval-turnaround estimate is published on the page itself. **Correction (supersedes the "already downloaded" BAM claim below this table previously made):** the pancreatic cancer subset actually on hand locally is **raw FASTQ, not BAM** — `/home/isilon/HumGenTempData/WGS/Public_data/Cristiano_Velculescu_DELFI/FASTQ/*/<sample_id>_R{1,2}.FASTQ`, one subdirectory per sample. The `/home/isilon/humangenetik/ARCHIVE/WGS/Public_data/Cristiano_Velculescu_DELFI/BAM/pancreatic_cancer/paired_end` path this table previously cited as an already-downloaded BAM subset has not been confirmed to actually contain ready BAMs — whether ready BAMs exist anywhere for this source is still open. | [EGAD00001005339](https://ega-archive.org/datasets/EGAD00001005339) · [EGAS00001003611](https://ega-archive.org/studies/EGAS00001003611) | **Partial, FASTQ only** — pancreatic subset in hand as raw FASTQ; needs alignment (see phase2-preprocessing-plan.md §1b) unless ready BAMs are found elsewhere |
| **TCGA WGS BAMs** (GDC, phs000178) | Yes, program still active | dbGaP-controlled for raw BAM/sequence; MAF/variant-level and clinical/biospecimen data is open | Confirmed unchanged: controlled tier = "individually identifiable... low level genomic sequencing data"; open tier = "high level genomic data... not individually identifiable." Access requires an eRA Commons account, a dbGaP **research project**, and DAC approval — the dbGaP project-creation flow is built around a PI + institutional Signing Official, not an unaffiliated individual. GDC's own docs don't describe a path for hobbyist/non-institutional applicants, and don't publish a turnaround estimate. | [GDC controlled-access docs](https://gdc.cancer.gov/access-data/obtaining-access-controlled-data) | **Blocked** for BAM-level access without institutional/PI sponsorship — not realistically pursuable solo |
| **ICGC PCAWG "open tier"** | Yes, portal exists | Open tier ≠ open BAMs — flag this as a likely misread of the spec's phrasing | ICGC's open tier covers derived calls (simple somatic mutations, expression, etc.), not raw aligned reads. The PCAWG BAMs referenced (e.g. EGAS00001001692, EGAD00001002127) are hosted at EGA and are themselves DACO/DAC-controlled, same access burden as TCGA/DELFI. | [ICGC PCAWG data docs](http://docs.icgc.org/pcawg/data/) | **Not a shortcut — controlled, same as TCGA** |
| **HCC1143 / HCC1954** (legacy CGHub benchmark set) | CGHub itself was retired (~2016); the historical `cghub.ucsc.edu/benchmark_download.html` link is dead infrastructure | Formerly distributed through CGHub; current copies live in GDC's legacy archive under the same **dbGaP-controlled** tier as other TCGA BAMs | Not a live open substitute today, contrary to how it's sometimes cited in older tutorials/slides. | — | **Not worth pursuing** — dead link / same access burden as TCGA |
| **HCC1395 / HCC1395BL** (SEQC2) | Yes, confirmed live | **Fully open**, no dbGaP | Matched tumor/normal breast cancer cell line pair, WGS at ~6 centers/platforms, ~34–47× (spec's ~55× is in the right range depending on center/replicate). Raw reads on SRA `SRP162370`; some BWA-aligned BAMs also mirrored on NCBI FTP. | [SRP162370](https://www.ncbi.nlm.nih.gov/sra/?term=SRP162370) · [SEQC2 docs](https://github.com/bioinform/somaticseq/blob/master/docs/seqc2.md) | **Usable now** — reconfirms spec's own v1 recommendation |
| **1000 Genomes high-coverage WGS (NYGC)** | Yes, live | Fully open, no registration | 3,202 samples (incl. 602 trios) at 30× on GRCh38, CRAM format. Hosted redundantly on EMBL-EBI/NCBI FTP, AWS, AnVIL. This is **blood-derived genomic DNA**, not cfDNA — mechanically sheared library, no *in vivo* nuclease fragmentation signal. | [IGSR data portal](https://www.internationalgenome.org/) · [NYGC announcement](https://www.internationalgenome.org/announcements/3202-samples-at-high-coverage-from-NYGC/) | **Usable now** — but re-flagging the spec's own caveat: this is not a cfDNA proxy, only a fully-open genomic-DNA background if nothing else is available |
| **ichorCNA panel-of-normals** (as a second healthy background) | Repo exists, PoN download/generation documented | Open (bundled with ichorCNA / generatable from own data) | ⚠️ **Conflicts with the spec's own feature requirements**: the PoN is distributed/generated as **WIG files of binned read counts** (1 Mb bins), not raw reads or BAM. It carries no fragment-length/insert-size information, so it cannot supply the short:long ratio or fragment-size-distribution features in §2 — it's only usable as a CNV-side normalization reference, not a fragmentomics healthy-background source. | [ichorCNA GitHub](https://github.com/broadinstitute/ichorCNA) · [PoN wiki](https://github.com/broadinstitute/ichorCNA/wiki/Create-Panel-of-Normals) | **Usable for CNV normalization only — not a fragmentomics data source as implied in spec §3.3** |

---

## Flags relative to phase0-spec.md

1. **§3.3 table entry "ichorCNA's own panel-of-normals... as a second healthy background source"** — this doesn't hold up. The PoN is pre-binned copy-number read-count data (WIG), with no per-fragment insert-size information, so it can't feed the DELFI-style short:long ratio or fragment-size features that are the project's primary signal (§2). It remains valid only for the secondary CNV feature's GC/mappability normalization, not as a second fragmentomics-comparable healthy cohort. Worth correcting in the spec if this is expected to double as a fragmentomics source.
2. **FinaleDB's governance basis is unclear and worth resolving before relying on it**, not just as an access-speed shortcut. If it is indeed re-hosting reprocessed derivatives of DAC-controlled studies (Cristiano/DELFI, per its own paper) without each downstream user going through the original EGA DAC, that's a provenance question worth settling (what data-use agreement covers the reprocessed copies?) before building a pipeline around it — independent of whether the site is even currently reachable, which this pass couldn't confirm.
3. **"ICGC open tier" is not open BAMs.** The spec doesn't explicitly claim this, but it's a common confusion worth flagging pre-emptively: ICGC's open tier is derived-call data; the actual aligned reads are DACO-controlled at EGA, same burden as TCGA. Don't treat it as a lower-friction TCGA alternative.
4. **HCC1143/HCC1954 via the "CGHub benchmark" pathway is dead** (CGHub retired). If these are ever considered as HCC1395 supplements, budget for the same dbGaP path as TCGA, not a quick open download.

No conflicts found with the spec's core methodology sections (§1, §4, §5) — this pass only touched §3.3 sourcing.

---

## How to get SEQC2 (HCC1395/HCC1395BL)

Two open, no-registration paths — commands shown for reference only, not run in this pass:

**Aligned BAM (BWA-MEM), direct from NCBI FTP:**
```bash
wget -r -np -nH --cut-dirs=4 \
  https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/data/
```
(Browse that index first to pick the specific center/platform BAMs needed — it's a multi-center dataset, not a single file pair.)

**Raw FASTQ via SRA toolkit (accession `SRP162370`):**
```bash
prefetch SRP162370
fasterq-dump --split-files <SRR_accession>
```
(`prefetch SRP162370` pulls every run under the project; to fetch a specific sample/run only, look up its individual `SRR...` accession via the SRA Run Selector first: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP162370)

High-confidence somatic call sets / regions (for validating any mixing pipeline against ground truth) are under `https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/release/latest/`.

---

## Recommendation

**Revised — most of v1's data is already local, not just "accessible":**

- **Tumor side:** HCC1395/HCC1395BL (SEQC2) — open, no registration required; not yet downloaded, but the commands above get it in one step. This still requires the in silico fragment-size reshaping step (§3.2.1) before mixing, regardless of source.
- **Healthy cfDNA background:** GSE71378 (Snyder et al. 2016), **already sitting at** `/home/isilon/humangenetik/ARCHIVE/WGS/Public_data/Snyder_Shendure_WPS/BAM/original_data_hg19/Cancer` — confirm whether the healthy (BH01/IA0x/IH0x) samples are actually present there (the directory name suggests it may be the cancer subset only, in which case the healthy samples still need pulling from SRP061633).
- **Bonus, not originally scoped:** the Cristiano/DELFI pancreatic-cancer subset is **already on hand as raw FASTQ** at `/home/isilon/HumGenTempData/WGS/Public_data/Cristiano_Velculescu_DELFI/FASTQ/*/<sample_id>_R{1,2}.FASTQ` — this is real cfDNA cancer-patient data, no DAC round-trip needed for this subset specifically. No ready BAM has been confirmed for this source; it needs alignment (phase2-preprocessing-plan.md §1b) unless ready BAMs turn up elsewhere. Only file the EGA DAC application if broader cancer-type coverage from the full 537-sample cohort turns out to be needed.
- **Fallback/cross-platform validation background:** 1000 Genomes NYGC high-coverage WGS — open, but flagged again as genomic DNA, not cfDNA; use only for a negative-control/batch-effect check (§3.2.3), not as a primary healthy background.
- **Investigate before use, don't build around it yet:** FinaleDB — confirm the site is actually reachable and clarify its data-use basis for reprocessed DAC-derived samples before depending on it as a way to avoid the EGA DAC process.

With Snyder and the Cristiano pancreatic subset already local, and HCC1395 a single open download away, v1 has effectively **zero remaining access-approval latency** — the DAC path for the full DELFI cohort is now optional scope-expansion, not a blocker.
