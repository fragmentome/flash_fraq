# Workstation SSD scratch for CRAM analysis

The workstation has one NVMe SSD and one spinning HDD:

- **SSD** (`nvme0n1`, ~954 GB usable): `/` (250 G, usually near-full) and `/mnt/data` (660 G). Free space on `/mnt/data` fluctuates but is typically only a few hundred GB.
- **HDD** (`sda`, 7.3 TB, single spinning disk, not RAID): mounted at `/mnt/8tb_hdd`. This is where CRAM sync pulls completed `*.sorted.cram` (+`.crai`, manifests) from the cluster archive via parallel rsync streams, landing under `/mnt/8tb_hdd/flash_fraq_crams/<dataset>/CRAM/`.

## Why stage on SSD before analysis

The HDD is a single spinning disk, so concurrent/random-access reads (parallel SLURM-style jobs, indexed region queries, GATK/samtools random access) hit head-seek contention and are much slower than on the SSD. CRAM decompression also benefits from lower storage latency on non-sequential reads. A single sequential pass over one CRAM is fine on the HDD; anything with parallelism or random access benefits from the SSD.

## Constraint: SSD is small relative to the HDD archive

`/mnt/data` has on the order of a few hundred GB free — enough for roughly 8-12 CRAMs (~25-40 GB each) at once, fewer once analysis outputs are written alongside them. This is nowhere near enough to hold the full HDD archive, so SSD use must be a **scratch pattern**, not persistent storage:

1. **Stage in**: copy the batch of CRAMs needed for the current analysis run from `/mnt/8tb_hdd/flash_fraq_crams/...` to a scratch dir on `/mnt/data`.
2. **Analyze**: run the workstation-side analysis against the SSD copies.
3. **Stage out / purge**: copy any results back to the HDD (or isilon), then delete the SSD scratch copies to free space for the next batch before staging the next one in.

## How this is implemented

The HDD-to-SSD-and-back leg is **not** a separate script — it is built into the Snakemake workflow, in `workflow/rules/local_features.smk` (config: `workflow/config/local_features.yaml`). Snakemake is a better fit than a batching script here, because the three steps above map onto machinery it already has:

| Step | Mechanism |
|---|---|
| Stage in | `rule localfeat_stage_cram`, HDD mirror → `/mnt/data/flash_fraq_stage/<dataset>/<sample>/` |
| Analyze | `localfeat_verify_cram`, `localfeat_coverage`, `localfeat_fragments`, `localfeat_gc_correct`, `localfeat_ichor` |
| Stage out / purge | staged CRAM and BAM are `temp()`; results are written straight to the HDD at `/mnt/8tb_hdd/flash_fraq_derived/` |

Two consequences worth knowing:

- **Purging is automatic and per-sample, not per-batch.** `temp()` deletes a sample's staged copies as soon as its last consumer finishes, so the SSD holds only the samples actively being worked on rather than a whole batch. There is no manual clean-up step to forget.
- **"Batch by available SSD space" is the `ssd_mb` resource.** Rather than sizing batches by hand, each staging job declares its footprint (measured off the HDD-side CRAM; the BAM is budgeted at 3×) and `--resources ssd_mb=N` caps the total in flight. The local profile sets 400 GB. Raise or lower it there in one place instead of re-deriving a batch size.

Run it with:

```
snakemake -s workflow/Snakefile --workflow-profile workflow/profiles/local local_features_complete
```

### Sizing note

The "roughly 8–12 CRAMs (~25–40 GB each)" figure above is an over-estimate for the data actually on hand: measured Snyder CRAMs are ~1.0–1.3 GB and Cristiano's are of the same order, so the real constraint is far looser than a dozen samples. The `ssd_mb` budget measures each file rather than assuming a per-sample size, so it stays correct either way.
