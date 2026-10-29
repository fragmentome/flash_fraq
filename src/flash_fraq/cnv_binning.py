"""Subprocess wrapper around HMMcopy Suite + ichorCNA for 1 Mb CN binning (phase2-preprocessing-plan.md §4).

`readCounter`, `gcCounter`, `mapCounter` (HMMcopy Suite) and `runIchorCNA.R` (ichorCNA) are
invoked here as external subprocesses only — their logic is never reimplemented or statically
linked into this package.

License note: HMMcopy Suite and ichorCNA are GPL-licensed (Broad Institute). Invoking them as
external subprocesses, rather than linking their code into this package, avoids triggering
GPL's viral distribution requirements for this internal research pipeline. This only becomes a
real constraint if the pipeline is later distributed as closed-source software — not a blocker
today.
"""

from dataclasses import dataclass, field
from pathlib import Path
from subprocess import CalledProcessError, CompletedProcess, run

from pandas import DataFrame, read_csv

# ichorCNA's <id>.params.txt reports the tumor fraction estimate as a labelled line,
# e.g. "Tumor Fraction: 0.15".
PARAMS_TUMOR_FRACTION_PREFIX: str = "Tumor Fraction:"

# bin_index = chromosome_rank * STRIDE + bin ordinal within the chromosome.
# 1000 comfortably clears the ~250 one-megabase bins on the largest human
# chromosome, so the two components never collide.
BIN_INDEX_CHROMOSOME_STRIDE: int = 1000


@dataclass(frozen=True)
class IchorCnaConfig:
    bam_path: Path
    reference_fasta: Path
    mappability_bigwig: Path
    output_dir: Path
    sample_id: str
    bin_size_bp: int = 1_000_000
    quality_threshold: int = 20
    chromosomes: tuple[str, ...] = tuple(f"chr{i}" for i in range(1, 23)) + ("chrX",)
    readcounter_bin: Path = Path("readCounter")
    gccounter_bin: Path = Path("gcCounter")
    mapcounter_bin: Path = Path("mapCounter")
    rscript_bin: Path = Path("Rscript")
    ichorcna_script: Path = Path("runIchorCNA.R")
    # Prebuilt GC/mappability WIGs. When set, run() uses these instead of
    # invoking gcCounter/mapCounter. The ichorCNA R package ships hg19 and
    # hg38 tracks at 1 Mb, so on the standard bin size gcCounter would only
    # regenerate an identical file per sample, and mapCounter needs a
    # mappability bigWig that is not always available at all.
    gc_wig_source: Path | None = None
    map_wig_source: Path | None = None
    # runIchorCNA.R arguments left at its defaults by the original wrapper.
    # genome_style in particular is not optional for chr-prefixed hg19: the
    # script defaults to NCBI naming, which matches zero chromosomes here and
    # yields an empty fit rather than an error.
    # Directory holding the ichorCNA R sources that match runIchorCNA.R's own
    # version, passed as --libdir. Not an optimization -- without it the script
    # falls back to `library(ichorCNA)` and picks up whatever version conda
    # installed. The checkout on this workstation is 0.3.2 (2019) while the
    # conda package is 0.5.1, and that pairing dies inside runEM with
    # "incorrect number of subscripts" because the array contract changed
    # between them. --libdir pins the script to its own sources.
    libdir: Path | None = None
    centromere: Path | None = None
    normal_panel: Path | None = None
    genome_build: str = "hg19"
    genome_style: str = "UCSC"
    # ULP-WGS HMM parameters, as R-vector literals, following the recipe in
    # ichorCNA's own scripts/README.md rather than runIchorCNA.R's bare
    # defaults. This is not tuning for its own sake: the script defaults leave
    # --scStates at "NULL" while --estimateScPrevalence stays TRUE, and that
    # pair makes runEM index a 2-D lambdas array with three subscripts, so it
    # dies with "incorrect number of subscripts" partway through the fit. The
    # multi-value ploidy/normal initializations and the softer transition
    # parameters come from the same recipe, which targets the 0.1-1x ULP-WGS
    # depth range this project operates in (phase0-spec.md §1.1).
    ploidy: str = "c(2,3)"
    normal: str = "c(0.5,0.6,0.7,0.8,0.9)"
    sc_states: str = "c(1,3)"
    # Finding 8 of docs/results-corruption-audit.md: at low TF (scStates left empty,
    # "c()") there is no subclonal state to estimate a prevalence for, so this should
    # be FALSE alongside an empty sc_states. True stays the default here so every
    # existing caller that does not pass it keeps today's behavior (the runIchorCNA.R
    # default this wrapper always relied on, see the sc_states comment above).
    estimate_sc_prevalence: bool = True
    max_cn: int = 5
    txn_e: float = 0.9999
    txn_strength: int = 10000
    read_wig_path: Path = field(init=False)
    gc_wig_path: Path = field(init=False)
    map_wig_path: Path = field(init=False)
    generated_gc_wig_path: Path = field(init=False)
    generated_map_wig_path: Path = field(init=False)
    cna_seg_path: Path = field(init=False)
    params_path: Path = field(init=False)
    cnv_bins_output_path: Path = field(init=False)
    ichor_estimate_output_path: Path = field(init=False)

    def __post_init__(self) -> None:
        object.__setattr__(self, "read_wig_path", self.output_dir / f"{self.sample_id}.readcounter.wig")
        # Two distinct paths on purpose. generated_* is where gcCounter /
        # mapCounter write; gc/map_wig_path is what run_ichor_cna hands to
        # runIchorCNA.R and prefers the prebuilt track when one was supplied.
        # Collapsing them would point _run_to_wig at gc_wig_source, so calling
        # run_gc_counter() with a prebuilt track configured would overwrite the
        # ichorCNA package's shipped extdata WIG -- a shared, read-only asset.
        generated_gc = self.output_dir / f"{self.sample_id}.gc.wig"
        generated_map = self.output_dir / f"{self.sample_id}.map.wig"
        object.__setattr__(self, "generated_gc_wig_path", generated_gc)
        object.__setattr__(self, "generated_map_wig_path", generated_map)
        object.__setattr__(
            self, "gc_wig_path", self.gc_wig_source if self.gc_wig_source is not None else generated_gc
        )
        object.__setattr__(
            self, "map_wig_path", self.map_wig_source if self.map_wig_source is not None else generated_map
        )
        object.__setattr__(self, "cna_seg_path", self.output_dir / f"{self.sample_id}.cna.seg")
        object.__setattr__(self, "params_path", self.output_dir / f"{self.sample_id}.params.txt")
        object.__setattr__(self, "cnv_bins_output_path", self.output_dir / "cnv_bins.csv")
        object.__setattr__(self, "ichor_estimate_output_path", self.output_dir / "ichor_estimates.csv")

    @classmethod
    def from_sources(
        cls,
        bam_path: Path,
        reference_dir: Path,
        output_dir: Path,
        sample_id: str | None = None,
    ) -> "IchorCnaConfig":
        return cls(
            bam_path=bam_path,
            reference_fasta=reference_dir / "hg19.fa",
            mappability_bigwig=reference_dir / "mappability_hg19_1000kb.bw",
            output_dir=output_dir,
            sample_id=sample_id if sample_id is not None else bam_path.stem,
        )


class IchorCnaRunner:
    """Drives readCounter -> gcCounter -> mapCounter -> runIchorCNA.R as external subprocesses."""

    def __init__(self, config: IchorCnaConfig) -> None:
        self.config = config

    def run_read_counter(self) -> Path:
        config = self.config
        cmd = [
            str(config.readcounter_bin),
            "--window", str(config.bin_size_bp),
            "--quality", str(config.quality_threshold),
            "--chromosome", ",".join(config.chromosomes),
            str(config.bam_path),
        ]
        self._run_to_wig(cmd, "readCounter", config.read_wig_path)
        return config.read_wig_path

    def run_gc_counter(self) -> Path:
        config = self.config
        cmd = [
            str(config.gccounter_bin),
            "--window", str(config.bin_size_bp),
            str(config.reference_fasta),
        ]
        self._run_to_wig(cmd, "gcCounter", config.generated_gc_wig_path)
        return config.generated_gc_wig_path

    def run_map_counter(self) -> Path:
        config = self.config
        cmd = [
            str(config.mapcounter_bin),
            "--window", str(config.bin_size_bp),
            str(config.mappability_bigwig),
        ]
        self._run_to_wig(cmd, "mapCounter", config.generated_map_wig_path)
        return config.generated_map_wig_path

    def run_ichor_cna(self) -> Path:
        config = self.config
        cmd = [
            str(config.rscript_bin),
            str(config.ichorcna_script),
            "--id", config.sample_id,
            "--WIG", str(config.read_wig_path),
            "--gcWig", str(config.gc_wig_path),
            "--mapWig", str(config.map_wig_path),
            "--outDir", str(config.output_dir),
            "--genomeBuild", config.genome_build,
            "--genomeStyle", config.genome_style,
            # R-vector literals, matching runIchorCNA.R's own option format.
            "--chrs", f"c({','.join(repr(chrom) for chrom in config.chromosomes)})",
            "--chrTrain", f"c({','.join(repr(c) for c in config.chromosomes if c != 'chrX')})",
            # 1 Mb bins are too coarse to call homozygous deletions; the
            # script's own help flags includeHOMD as the setting to disable
            # at this bin size.
            "--includeHOMD", "FALSE",
            "--ploidy", config.ploidy,
            "--normal", config.normal,
            "--scStates", config.sc_states,
            "--estimateScPrevalence", "TRUE" if config.estimate_sc_prevalence else "FALSE",
            "--maxCN", str(config.max_cn),
            "--txnE", str(config.txn_e),
            "--txnStrength", str(config.txn_strength),
        ]
        if config.libdir is not None:
            cmd += ["--libdir", str(config.libdir)]
        if config.centromere is not None:
            cmd += ["--centromere", str(config.centromere)]
        if config.normal_panel is not None:
            cmd += ["--normalPanel", str(config.normal_panel)]
        self._run_subprocess(cmd, "runIchorCNA.R")
        return config.cna_seg_path

    def resolve_gc_wig(self) -> Path:
        """Prebuilt GC WIG if one was configured, else generate with gcCounter."""
        if self.config.gc_wig_source is not None:
            return self.config.gc_wig_source
        return self.run_gc_counter()

    def resolve_map_wig(self) -> Path:
        """Prebuilt mappability WIG if one was configured, else mapCounter."""
        if self.config.map_wig_source is not None:
            return self.config.map_wig_source
        return self.run_map_counter()

    def parse_cna_bins(self) -> DataFrame:
        return read_csv(self.config.cna_seg_path, sep="\t")

    def parse_tumor_fraction_estimate(self) -> float:
        for line in self.config.params_path.read_text().splitlines():
            if line.startswith(PARAMS_TUMOR_FRACTION_PREFIX):
                return float(line.removeprefix(PARAMS_TUMOR_FRACTION_PREFIX).strip())
        raise ValueError(f"no '{PARAMS_TUMOR_FRACTION_PREFIX}' line found in {self.config.params_path}")

    def _logr_column(self, bins: DataFrame) -> str:
        """ichorCNA prefixes its per-sample columns with the --id it was given,
        so the real header is "<sample_id>.logR", not "logR". The bare name is
        still accepted as a fallback for hand-built or fixture input."""
        prefixed = f"{self.config.sample_id}.logR"
        if prefixed in bins.columns:
            return prefixed
        if "logR" in bins.columns:
            return "logR"
        raise ValueError(
            f"no logR column in {self.config.cna_seg_path}: expected {prefixed!r} or 'logR', "
            f"found {list(bins.columns)}"
        )

    def _bin_indices(self, bins: DataFrame) -> list[int]:
        """A genome-position-derived bin index, stable across samples.

        Numbering the rows of .cna.seg 1..N is wrong: ichorCNA drops bins that
        fail its mappability/centromere/panel-of-normals filters, so the row
        count varies between samples (and even between runs of one sample --
        3053 generated bins vs. 2594 retained with the PoN on IA01). Since
        feature_matrix.load_cnv_bins_wide() pivots bin_index straight into
        cnv_logr_bin_<n> columns, a positional index silently aligns different
        loci across samples: a wrong answer with no error anywhere.

        Encoding the chromosome's rank in `chromosomes` and the bin's ordinal
        within it makes the same locus produce the same index in every sample,
        and leaves a bin that one sample lacks as a NaN in the pivot rather
        than a shifted neighbour.
        """
        missing = [column for column in ("chr", "start") if column not in bins.columns]
        if missing:
            raise ValueError(
                f"{self.config.cna_seg_path} lacks {missing}; cannot derive a position-stable "
                f"bin_index (found {list(bins.columns)})"
            )
        ranks = {chromosome: rank for rank, chromosome in enumerate(self.config.chromosomes)}
        unknown = sorted(set(bins["chr"].astype(str)) - set(ranks))
        if unknown:
            raise ValueError(
                f"{self.config.cna_seg_path} contains chromosomes outside the configured set: {unknown}"
            )
        return [
            ranks[str(chromosome)] * BIN_INDEX_CHROMOSOME_STRIDE + self._bin_ordinal(int(start))
            for chromosome, start in zip(bins["chr"], bins["start"])
        ]

    def _bin_ordinal(self, start: int) -> int:
        """Bin number within its chromosome, agnostic to the start convention.

        readCounter emits 1-based WIG starts (1, 1000001, ...) and ichorCNA
        carries those through to .cna.seg, but assuming that outright is a
        one-character way to be silently wrong: a plain `start // bin` maps a
        1-based 1 and 1000001 to bins 0 and 1 correctly, while a plain
        `(start - 1) // bin` would collapse a 0-based 0 and 1000000 both to 0.
        A start that lands exactly on a bin boundary is 0-based; one offset by
        a single base is 1-based.
        """
        bin_size = self.config.bin_size_bp
        if start % bin_size == 0:
            return start // bin_size
        return (start - 1) // bin_size

    def write_cnv_bins(self, bins: DataFrame) -> Path:
        long_bins = DataFrame(
            {
                "sample_id": self.config.sample_id,
                "bin_index": self._bin_indices(bins),
                "chrom": bins["chr"].astype(str).to_numpy(),
                "start": bins["start"].to_numpy(),
                "cnv_logr": bins[self._logr_column(bins)].to_numpy(),
            }
        )
        long_bins.to_csv(self.config.cnv_bins_output_path, index=False)
        return self.config.cnv_bins_output_path

    def write_ichor_estimate(self, tumor_fraction_estimate: float) -> Path:
        DataFrame(
            [{"sample_id": self.config.sample_id, "ichorCNA_tumor_fraction_estimate": tumor_fraction_estimate}]
        ).to_csv(self.config.ichor_estimate_output_path, index=False)
        return self.config.ichor_estimate_output_path

    def run(self) -> DataFrame:
        self.config.output_dir.mkdir(parents=True, exist_ok=True)
        self.run_read_counter()
        self.resolve_gc_wig()
        self.resolve_map_wig()
        self.run_ichor_cna()
        bins = self.parse_cna_bins()
        self.write_cnv_bins(bins)
        self.write_ichor_estimate(self.parse_tumor_fraction_estimate())
        return bins

    def _run_to_wig(self, cmd: list[str], tool_name: str, wig_path: Path) -> None:
        result = self._run_subprocess(cmd, tool_name)
        wig_path.write_bytes(result.stdout)

    @staticmethod
    def _run_subprocess(cmd: list[str], tool_name: str) -> CompletedProcess[bytes]:
        try:
            return run(cmd, check=True, capture_output=True)
        except CalledProcessError as error:
            stderr = error.stderr.decode() if error.stderr else ""
            raise RuntimeError(f"{tool_name} failed (exit {error.returncode}): {stderr}") from error
