"""In silico tumor/healthy admixture builder (docs/phase2-preprocessing-plan.md §7).

Downsamples a tumor source BAM and a healthy background BAM to a target
admixture fraction `f`, merges them, dedups the merged BAM once (never
per-source, see CLAUDE.md), and derives the empirical SNV-VAF tumor fraction
by comparing observed alt-allele frequency at known HCC1395 somatic SNV sites
against the expected `f * VAF_tumor_pure`.

`AdmixtureConfig.snv_vaf_paths` maps a label to an SNV site set (TSV for the
synthetic test fixtures, VCF -- read via `bcftools query`, `INFO/bwaTVAF` as
`vaf_tumor_pure` -- for the real SEQC2 call set); `AdmixtureBuilder.run()`
computes one empirical tumor fraction per label and records all of them,
per MAT-13's "compute both estimators, agreement is a correctness check".
"""

from __future__ import annotations

import csv
import logging
from dataclasses import dataclass
from math import log10
from pathlib import Path
from subprocess import run

from pandas import DataFrame

_LOGGER = logging.getLogger(__name__)

MTDNA_CONTIG: str = "chrM"
MTDNA_FRACTION_FLOOR: float = 1e-9

# `samtools view -s` can only ever *remove* reads, so a source library cannot be
# downsampled up to a target depth it does not already have.
MAX_DOWNSAMPLE_FRACTION: float = 1.0


def downsample_fractions(
    tumor_fraction: float,
    target_coverage: float,
    tumor_native_coverage: float,
    healthy_native_coverage: float,
) -> tuple[float, float]:
    """The two phase2 §7 step-1 downsample fractions, unvalidated.

    Split out of `AdmixtureBuilder.compute_downsample_fractions` so that
    `admixture_design` can tabulate the *infeasible* combinations too: the builder
    raises as soon as a fraction exceeds 1, which is right when it is about to run
    samtools, but a feasibility table needs the offending value in order to say by
    how much a candidate target coverage misses. Both callers share this one
    arithmetic statement -- see admixture_design.AdmixtureFeasibilityAnalyzer.
    """
    tumor = (tumor_fraction * target_coverage) / tumor_native_coverage
    healthy = ((1 - tumor_fraction) * target_coverage) / healthy_native_coverage
    return tumor, healthy


@dataclass(frozen=True)
class SomaticSnv:
    chrom: str
    pos: int
    ref: str
    alt: str
    vaf_tumor_pure: float


@dataclass(frozen=True)
class AdmixtureConfig:
    tumor_bam: Path
    healthy_bam: Path
    reference_fasta: Path
    snv_vaf_paths: dict[str, Path]
    output_dir: Path
    tumor_fraction: float
    target_coverage: float
    tumor_native_coverage: float
    healthy_native_coverage: float
    seed: int = 42
    sample_name: str = "mixed"
    tumor_source_id: str = ""
    healthy_source_id: str = ""
    tumor_subject_id: str = ""
    healthy_subject_id: str = ""
    platform_tumor: str = "unknown"
    platform_healthy: str = "unknown"
    reference_build: str = "hg19"
    threads: int = 8

    @classmethod
    def from_sources(
        cls,
        tumor_bam: Path,
        healthy_bam: Path,
        reference_fasta: Path,
        snv_vaf_paths: dict[str, Path],
        output_dir: Path,
        tumor_fraction: float,
        target_coverage: float,
        tumor_native_coverage: float,
        healthy_native_coverage: float,
        seed: int = 42,
        sample_name: str = "mixed",
        tumor_source_id: str | None = None,
        healthy_source_id: str | None = None,
        tumor_subject_id: str | None = None,
        healthy_subject_id: str | None = None,
        platform_tumor: str = "unknown",
        platform_healthy: str = "unknown",
        reference_build: str = "hg19",
        threads: int = 8,
    ) -> "AdmixtureConfig":
        tumor_source_id = tumor_source_id if tumor_source_id is not None else Path(tumor_bam).stem
        healthy_source_id = healthy_source_id if healthy_source_id is not None else Path(healthy_bam).stem
        return cls(
            tumor_bam=tumor_bam,
            healthy_bam=healthy_bam,
            reference_fasta=reference_fasta,
            snv_vaf_paths={label: Path(path) for label, path in snv_vaf_paths.items()},
            output_dir=output_dir,
            tumor_fraction=tumor_fraction,
            target_coverage=target_coverage,
            tumor_native_coverage=tumor_native_coverage,
            healthy_native_coverage=healthy_native_coverage,
            seed=seed,
            sample_name=sample_name,
            tumor_source_id=tumor_source_id,
            healthy_source_id=healthy_source_id,
            tumor_subject_id=tumor_subject_id if tumor_subject_id is not None else tumor_source_id,
            healthy_subject_id=healthy_subject_id if healthy_subject_id is not None else healthy_source_id,
            platform_tumor=platform_tumor,
            platform_healthy=platform_healthy,
            reference_build=reference_build,
            threads=threads,
        )


class AdmixtureBuilder:
    """Drives the downsample -> merge -> dedup -> empirical-VAF pipeline for one
    mixed admixture sample, per phase2-preprocessing-plan.md §7."""

    def __init__(self, config: AdmixtureConfig) -> None:
        self.config = config

    def compute_downsample_fractions(self) -> tuple[float, float]:
        cfg = self.config
        tumor_fraction, healthy_fraction = downsample_fractions(
            tumor_fraction=cfg.tumor_fraction,
            target_coverage=cfg.target_coverage,
            tumor_native_coverage=cfg.tumor_native_coverage,
            healthy_native_coverage=cfg.healthy_native_coverage,
        )
        if tumor_fraction > MAX_DOWNSAMPLE_FRACTION:
            raise ValueError(f"tumor downsample fraction {tumor_fraction} exceeds 1 (native coverage too low)")
        if healthy_fraction > MAX_DOWNSAMPLE_FRACTION:
            raise ValueError(f"healthy downsample fraction {healthy_fraction} exceeds 1 (native coverage too low)")
        return tumor_fraction, healthy_fraction

    def subsample_bam(self, source_bam: Path, fraction: float, output_bam: Path) -> Path:
        if not 0 < fraction <= 1:
            raise ValueError(f"subsample fraction must be in (0, 1], got {fraction}")
        subsample_value = f"{self.config.seed + fraction:.6f}"
        run(
            [
                "samtools", "view", "-b", "-s", subsample_value,
                "-@", str(self.config.threads),
                "--reference", str(self.config.reference_fasta),
                "-o", str(output_bam), str(source_bam),
            ],
            check=True,
        )
        return output_bam

    def verify_downsample_depth(
        self,
        source_bam: Path,
        subsampled_bam: Path,
        expected_fraction: float,
        tolerance: float = 0.25,
        reads_before: int | None = None,
    ) -> float:
        # `reads_before` lets a caller pass a precomputed total (see
        # `run`'s `precomputed_tumor_total_reads`) instead of a fresh
        # `samtools flagstat` over the full source -- for the tumor side that
        # source is the ~88GB CRAM every admixture row shares, and flagstat
        # walks every record just like `view -s` does, so recomputing it per
        # row would reintroduce the same redundant-decode cost the batch
        # sampler (tumor_batch_sampling.py) exists to eliminate.
        reads_before = reads_before if reads_before is not None else self._flagstat_total_reads(source_bam)
        reads_after = self._flagstat_total_reads(subsampled_bam)
        actual_fraction = reads_after / reads_before
        if abs(actual_fraction - expected_fraction) > tolerance:
            raise ValueError(
                f"retained fraction {actual_fraction:.4f} for {subsampled_bam} deviates from "
                f"target {expected_fraction:.4f} by more than tolerance {tolerance}"
            )
        return actual_fraction

    def merge_and_dedup(self, subsampled_bams: list[Path]) -> Path:
        cfg = self.config
        merged = cfg.output_dir / f"{cfg.sample_name}.merged.bam"
        run(
            [
                "samtools", "merge", "-f",
                "-@", str(cfg.threads),
                "--reference", str(cfg.reference_fasta),
                str(merged), *[str(bam) for bam in subsampled_bams],
            ],
            check=True,
        )

        name_sorted = cfg.output_dir / f"{cfg.sample_name}.namesort.bam"
        run(
            ["samtools", "sort", "-n", "-@", str(cfg.threads), "-o", str(name_sorted), str(merged)],
            check=True,
        )

        fixmated = cfg.output_dir / f"{cfg.sample_name}.fixmate.bam"
        run(
            ["samtools", "fixmate", "-m", "-@", str(cfg.threads), str(name_sorted), str(fixmated)],
            check=True,
        )

        coord_sorted = cfg.output_dir / f"{cfg.sample_name}.sorted.bam"
        run(
            ["samtools", "sort", "-@", str(cfg.threads), "-o", str(coord_sorted), str(fixmated)],
            check=True,
        )

        dedup_bam = cfg.output_dir / f"{cfg.sample_name}.dedup.bam"
        run(
            ["samtools", "markdup", "-@", str(cfg.threads), str(coord_sorted), str(dedup_bam)],
            check=True,
        )
        run(["samtools", "index", "-@", str(cfg.threads), str(dedup_bam)], check=True)
        return dedup_bam

    def load_snvs(self, snv_vaf_path: Path) -> list[SomaticSnv]:
        """Somatic SNV sites + pure-tumor VAF, from a TSV (synthetic test fixtures)
        or a VCF (the real SEQC2 call set, `INFO/bwaTVAF`) -- picked by extension so
        both `AdmixtureConfig.snv_vaf_paths` entries can point at either shape."""
        name = Path(snv_vaf_path).name
        if name.endswith(".vcf") or name.endswith(".vcf.gz"):
            return self._load_snvs_vcf(snv_vaf_path)
        return self._load_snvs_tsv(snv_vaf_path)

    def _load_snvs_tsv(self, snv_vaf_path: Path) -> list[SomaticSnv]:
        with open(snv_vaf_path, newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            return [
                SomaticSnv(
                    chrom=row["chrom"],
                    pos=int(row["pos"]),
                    ref=row["ref"],
                    alt=row["alt"],
                    vaf_tumor_pure=float(row["vaf_tumor_pure"]),
                )
                for row in reader
            ]

    def _load_snvs_vcf(self, snv_vaf_path: Path) -> list[SomaticSnv]:
        # bcftools query has no --threads option (it's not I/O-bound the way
        # mpileup/samtools view/sort are -- this file is a small SNV site set,
        # not a multi-GB alignment stream), so there is nothing to wire here.
        result = run(
            [
                "bcftools", "query",
                "-f", "%CHROM\t%POS\t%REF\t%ALT\t%INFO/bwaTVAF\n",
                str(snv_vaf_path),
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        snvs = []
        skipped = 0
        for line in result.stdout.splitlines():
            chrom, pos, ref, alt, vaf_tumor_pure_raw = line.split("\t")
            try:
                vaf_tumor_pure = float(vaf_tumor_pure_raw)
            except ValueError:
                # bcftools query renders a missing INFO field as "." -- confirmed
                # against the real SEQC2 sites.vcf.gz that this is a real, if rare,
                # occurrence (1/39141 sites there) rather than a hypothetical. One
                # missing bwaTVAF annotation should not invalidate the whole row's
                # empirical-tumor-fraction estimate, which already aggregates
                # across every other site -- skip just this site, but count and
                # log it rather than silently dropping it.
                skipped += 1
                continue
            snvs.append(
                SomaticSnv(chrom=chrom, pos=int(pos), ref=ref, alt=alt, vaf_tumor_pure=vaf_tumor_pure)
            )
        if skipped:
            _LOGGER.warning(
                "%s: skipped %d/%d site(s) with a missing/non-numeric INFO/bwaTVAF",
                snv_vaf_path, skipped, skipped + len(snvs),
            )
        return snvs

    def _sort_snvs_by_reference_order(self, snvs: list[SomaticSnv]) -> list[SomaticSnv]:
        """Sorts by the reference's own contig order (.fai first column), then position.

        Plain lexicographic sort would put chr10 before chr2; the .fai/@SQ order
        (chr1, chr2, ..., chr22, chrX, chrY, chrM in this hg19 setup) is what the
        BAM's BAI index is built against, which is what `bcftools mpileup -R`
        needs to seek forward-only.
        """
        fai_path = Path(f"{self.config.reference_fasta}.fai")
        with open(fai_path) as handle:
            contig_order = {line.split("\t", 1)[0]: index for index, line in enumerate(handle)}
        return sorted(snvs, key=lambda snv: (contig_order[snv.chrom], snv.pos))

    def compute_empirical_tumor_fraction(self, merged_dedup_bam: Path, snvs: list[SomaticSnv], label: str) -> float:
        # bcftools mpileup -R needs its regions sorted to match the BAM's index
        # (.fai/@SQ contig order, then position) for efficient indexed seeking;
        # an unsorted region list forces repeated backward seeks and turns a
        # sub-minute scan into many minutes per mixture on spinning disk.
        sorted_snvs = self._sort_snvs_by_reference_order(snvs)
        sites_bed = self.config.output_dir / f"{self.config.sample_name}.sites.bed"
        with open(sites_bed, "w") as handle:
            for snv in sorted_snvs:
                handle.write(f"{snv.chrom}\t{snv.pos - 1}\t{snv.pos}\n")

        result = run(
            [
                "bcftools", "mpileup",
                "-f", str(self.config.reference_fasta),
                "-R", str(sites_bed),
                "-a", "AD",
                "-Ov",
                "--threads", str(self.config.threads),
                str(merged_dedup_bam),
            ],
            check=True,
            capture_output=True,
            text=True,
        )

        snvs_by_site = {(snv.chrom, snv.pos): snv for snv in snvs}
        # Depth-weighted ratio of sums, not a mean of per-site ratios (finding 7 of
        # docs/results-corruption-audit.md): at 1x most covered sites have depth 1-2,
        # so a single alt read at a VAF-0.04 site used to add 25/depth to an
        # unweighted average, swamping every well-covered site and inflating label
        # noise exactly in the low-f regime that matters. Summing alt counts and
        # depth*vaf_tumor_pure separately before dividing gives every read equal
        # weight instead of every site equal weight, which is unbiased under the
        # same binomial-per-read model and has far lower variance.
        sum_alt = 0
        sum_depth_vaf = 0.0
        site_rows: list[tuple[str, int, int, int, float]] = []
        for line in result.stdout.splitlines():
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            chrom, pos, alt_field = fields[0], int(fields[1]), fields[4]
            snv = snvs_by_site.get((chrom, pos))
            if snv is None:
                continue

            # The merged/dedup BAM keeps the tumor and healthy sources' distinct
            # @RG SM tags, so bcftools mpileup groups by sample and emits one
            # column per SM -- not one column for the whole mixture. Reading a
            # single fixed column (previously fields[9]) silently measured only
            # one source's reads at each site, which happened to reproduce
            # ~vaf_tumor_pure at any site a tumor read touched, regardless of
            # the true admixture fraction. Sum AD across every sample column so
            # depth/allele counts reflect the whole merged mixture.
            format_keys = fields[8].split(":")
            ad_index = format_keys.index("AD")
            allele_depths: list[int] | None = None
            for sample_field in fields[9:]:
                sample_ad = [int(x) for x in sample_field.split(":")[ad_index].split(",")]
                if allele_depths is None:
                    allele_depths = sample_ad
                else:
                    allele_depths = [a + b for a, b in zip(allele_depths, sample_ad)]
            depth = sum(allele_depths)
            if depth == 0:
                continue

            alts = alt_field.split(",") if alt_field != "." else []
            alt_count = allele_depths[alts.index(snv.alt) + 1] if snv.alt in alts else 0
            sum_alt += alt_count
            sum_depth_vaf += depth * snv.vaf_tumor_pure
            site_rows.append((chrom, pos, alt_count, depth, snv.vaf_tumor_pure))

        self._write_empirical_tumor_fraction_sites(label, site_rows)
        return sum_alt / sum_depth_vaf

    def _write_empirical_tumor_fraction_sites(
        self, label: str, site_rows: list[tuple[str, int, int, int, float]]
    ) -> Path:
        """Per-site (chrom, pos, alt_count, depth, vaf_tumor_pure), gzipped next to the
        mixture's other per-label outputs -- so a future estimator change (this one
        replaced the unweighted per-site-ratio mean, see the comment above) can be
        recomputed from this file alone, with no mpileup rerun."""
        cfg = self.config
        cfg.output_dir.mkdir(parents=True, exist_ok=True)
        output_path = cfg.output_dir / f"{cfg.sample_name}.empirical_tumor_fraction.{label}.sites.tsv.gz"
        DataFrame(
            site_rows, columns=["chrom", "pos", "alt_count", "depth", "vaf_tumor_pure"]
        ).to_csv(output_path, sep="\t", index=False, compression="gzip")
        return output_path

    def compute_mtdna_fraction_log10(self, dedup_bam: Path) -> float:
        result = run(
            ["samtools", "idxstats", "-@", str(self.config.threads), str(dedup_bam)],
            check=True,
            capture_output=True,
            text=True,
        )
        mtdna_mapped = 0
        total_mapped = 0
        for line in result.stdout.splitlines():
            contig, _length, mapped, _unmapped = line.split("\t")
            mapped = int(mapped)
            total_mapped += mapped
            if contig == MTDNA_CONTIG:
                mtdna_mapped = mapped

        fraction = mtdna_mapped / total_mapped if total_mapped > 0 else 0.0
        return log10(max(fraction, MTDNA_FRACTION_FLOOR))

    def write_sample_metadata(
        self, empirical_tumor_fractions: dict[str, float], mtdna_fraction_log10: float
    ) -> Path:
        """`empirical_tumor_fraction` is `feature_matrix.REGRESSION_LABEL_COLUMN` -- MAT-17's
        assembler and everything downstream of it key off that exact name, so it must keep
        meaning one number. It is set to the *first* `snv_vaf_paths` entry's estimate (the
        config declares the canonical one first, e.g. the full-set regression sites file
        ahead of the smaller clonal-ratio file); every label, including that one, also gets
        its own `empirical_tumor_fraction_<label>` column so MAT-13's cross-estimator
        agreement check has both numbers to compare."""
        cfg = self.config
        cfg.output_dir.mkdir(parents=True, exist_ok=True)
        output_path = cfg.output_dir / "sample_metadata.csv"
        canonical_label = next(iter(empirical_tumor_fractions))
        row = {
            "sample_id": cfg.sample_name,
            "tumor_source_id": cfg.tumor_source_id,
            "healthy_source_id": cfg.healthy_source_id,
            "platform_tumor": cfg.platform_tumor,
            "platform_healthy": cfg.platform_healthy,
            "reference_build": cfg.reference_build,
            "target_coverage": cfg.target_coverage,
            "tumor_subject_id": cfg.tumor_subject_id,
            "healthy_subject_id": cfg.healthy_subject_id,
            "nominal_tumor_fraction": cfg.tumor_fraction,
            "empirical_tumor_fraction": empirical_tumor_fractions[canonical_label],
        }
        for label, fraction in empirical_tumor_fractions.items():
            row[f"empirical_tumor_fraction_{label}"] = fraction
        row["mtdna_fraction_log10"] = mtdna_fraction_log10
        DataFrame([row]).to_csv(output_path, index=False)
        return output_path

    def run(
        self,
        precomputed_tumor_sub: Path | None = None,
        precomputed_tumor_total_reads: int | None = None,
    ) -> dict[str, float]:
        """`precomputed_tumor_sub` lets the tumor-side downsample be produced
        upstream, out-of-band, by `tumor_batch_sampling.TumorBatchSampler` --
        one process decoding the shared tumor CRAM once for every plan row --
        instead of this method's own `subsample_bam` re-decoding the full
        source here. `precomputed_tumor_total_reads` is that same batch pass's
        free byproduct (see its docstring); passing it skips the redundant
        `samtools flagstat` this method would otherwise run against the full
        tumor source in `verify_downsample_depth`. Both are `None` by default,
        so a caller building a single mixture standalone (e.g. the matched_pair
        smoke tests, or a run outside the Snakemake DAG) is unaffected."""
        cfg = self.config
        cfg.output_dir.mkdir(parents=True, exist_ok=True)

        tumor_fraction, healthy_fraction = self.compute_downsample_fractions()

        # tumor_fraction == 0 is a legitimate plan row: a pure-background
        # negative control (nominal admixture fraction f=0) that
        # tf_eval.compute_detection_floor's null pool needs. There is no
        # tumor read to keep at f=0 -- subsample_bam's `0 < fraction <= 1`
        # check would reject it outright, and there is nothing to merge on
        # the tumor side -- so this builds the mixture from the healthy
        # downsample alone. compute_empirical_tumor_fraction then runs
        # unmodified against that healthy-only dedup BAM and comes out near
        # zero on its own (no tumor-specific alt alleles present in a
        # healthy source), which is exactly the "sane" value this row's
        # metadata should carry -- no separate zero-fraction path is needed
        # downstream of merge_and_dedup.
        if tumor_fraction == 0:
            tumor_sub = None
        elif precomputed_tumor_sub is not None:
            tumor_sub = precomputed_tumor_sub
        else:
            tumor_sub = self.subsample_bam(
                cfg.tumor_bam, tumor_fraction, cfg.output_dir / f"{cfg.sample_name}.tumor_sub.bam"
            )

        healthy_sub = self.subsample_bam(
            cfg.healthy_bam, healthy_fraction, cfg.output_dir / f"{cfg.sample_name}.healthy_sub.bam"
        )
        if tumor_sub is not None:
            self.verify_downsample_depth(
                cfg.tumor_bam, tumor_sub, tumor_fraction, reads_before=precomputed_tumor_total_reads
            )
        self.verify_downsample_depth(cfg.healthy_bam, healthy_sub, healthy_fraction)

        subsampled_bams = [healthy_sub] if tumor_sub is None else [tumor_sub, healthy_sub]
        dedup_bam = self.merge_and_dedup(subsampled_bams)

        # MAT-13: compute both SNV-VAF estimators (full-set regression sites,
        # clonal-ratio sites) per mixture -- agreement between them is a free
        # correctness check, not a gate (see the wiring ticket). Insertion
        # order, not sorted: write_sample_metadata takes the first entry as
        # canonical, so the config's ordering decides which one that is.
        empirical_tumor_fractions: dict[str, float] = {}
        for label, snv_vaf_path in cfg.snv_vaf_paths.items():
            snvs = self.load_snvs(snv_vaf_path)
            fraction = self.compute_empirical_tumor_fraction(dedup_bam, snvs, label)
            empirical_tumor_fractions[label] = fraction
            result_path = cfg.output_dir / f"{cfg.sample_name}.empirical_tumor_fraction.{label}.txt"
            result_path.write_text(f"{fraction}\n")

        mtdna_fraction_log10 = self.compute_mtdna_fraction_log10(dedup_bam)
        self.write_sample_metadata(empirical_tumor_fractions, mtdna_fraction_log10)

        return empirical_tumor_fractions

    def _flagstat_total_reads(self, bam: Path) -> int:
        # flagstat has no top-level --reference flag (unlike view/merge) --
        # CRAM decoding is instead selected via --input-fmt-option.
        result = run(
            [
                "samtools", "flagstat",
                "-@", str(self.config.threads),
                "--input-fmt-option", f"reference={self.config.reference_fasta}",
                str(bam),
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        first_line = result.stdout.splitlines()[0]
        return int(first_line.split()[0])
