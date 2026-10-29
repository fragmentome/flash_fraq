# data/

Nothing in this directory is committed except this file. Raw reads, alignments, feature matrices and
model files are large and (for the Cristiano/DELFI cohort) access-controlled; see `.gitignore`.

| Cohort | Role | Accession | Access |
| --- | --- | --- | --- |
| Snyder et al. 2016 (Shendure lab) | healthy / cancer cfDNA backgrounds | [GSE71378](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE71378), reads at SRA [SRP061633](https://www.ncbi.nlm.nih.gov/sra/?term=SRP061633) | open |
| Cristiano et al. 2019 (DELFI) | cfDNA backgrounds, fragment-ratio reference | EGA [EGAS00001003611](https://ega-archive.org/studies/EGAS00001003611) / [EGAD00001005339](https://ega-archive.org/datasets/EGAD00001005339) | **controlled (DAC approval required)** |
| SEQC2 HCC1395 / HCC1395BL | tumor and matched-normal cell lines for admixture | SRA [SRP162370](https://www.ncbi.nlm.nih.gov/sra/?term=SRP162370) | open |

Run-level HCC1395 manifests are tracked in `manifests/hcc1395/`; fetch with
`scripts/fetch_hcc1395_fastq.py`. The reference is hg19 (`scripts/fetch_hg19_reference.sbatch`).
Expected disk use for the full pipeline is multiple TB of FASTQ/BAM; only CRAMs and derived
feature tables are kept long-term. Set the data roots in `workflow/config/*.yaml` for your site.

You do not need any of this for the smoke test (`python scripts/smoke_test.py`), which uses
synthetic data.
