"""hg19 centromere/chromosome-end masking for coverage-sensitive per-window features.

`docs/results-corruption-audit.md` and the reshaping-pilot coverage check
(`/mnt/8tb_hdd/flash_fraq_derived/reshaping_pilot/coverage_check/coverage_bias_report.json`)
found reshaped-tumor coverage is biased near centromeres and chromosome ends. This module
is the single source of truth for "is this bin too close to one of those to trust" --
5 Mb windows (DELFI short/long features, MAT-15) and 1 Mb bins (ichorCNA CN, MAT-16) both
call it before their per-window features reach `gc_correction.py` / `feature_matrix.py`.

Data sources (vendored under `data/`, not read from `/mnt` -- tests must not depend on
cluster/workstation paths):
- `data/hg19_chrom_sizes.tsv`: contig lengths for chr1-22 + chrX, from
  `/mnt/data/reference_genome/hg19/hg19.fa.fai` (the same hg19 build and "chr" contig
  naming `scripts/make_genome_windows.py` builds the DELFI/ichor window BEDs from).
  chrY is intentionally excluded -- neither `local_features.yaml`'s `primary_contigs`
  nor `cnv_binning.IchorCnaConfig.chromosomes` carry it, so no feature ever has a chrY bin.
- `data/hg19_centromeres.tsv`: centromere spans (start, end), from the same UCSC gap
  table `scripts/drivers/admixture_mixture_ichor_driver.sh` already points ichorCNA's
  `--centromere` at (`GRCh37.p13_centromere_UCSC-gapTable.txt`). ichorCNA itself only
  drops a ~1e5 bp flank around this span; this module adds the wider margin the coverage
  bias report calls for.

Deliberately NOT used for chromosome-arm bucketing: `gc_correction.HG19_CENTROMERE_MB`
(a coarser, hand-entered approximate midpoint table) stays as the single source of truth
for `arm_id`/`arm_zscore_*` labels. Diffing every real 5 Mb window's arm assignment under
the two tables turned up 5 differences (chr4:50e6, chr6:60e6, chr7:60e6, chr9:45e6,
chr10:40e6 flip p/q) out of ~600 windows -- not "trivially replaceable" without changing
which arm those windows land in, so this module and `HG19_CENTROMERE_MB` intentionally
disagree by a few Mb and serve different purposes (masking vs. arm bucketing).
"""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path

from pandas import DataFrame, Series, read_csv

DATA_DIR: Path = Path(__file__).parent / "data"
CHROM_SIZES_PATH: Path = DATA_DIR / "hg19_chrom_sizes.tsv"
CENTROMERES_PATH: Path = DATA_DIR / "hg19_centromeres.tsv"

DEFAULT_MASK_MARGIN_BP: int = 5_000_000


@lru_cache(maxsize=None)
def load_chrom_sizes(path: Path = CHROM_SIZES_PATH) -> dict[str, int]:
    """chrom -> contig length (bp), hg19, chr1-22 + chrX."""
    sizes = read_csv(path, sep="\t")
    return dict(zip(sizes["chrom"], sizes["length"].astype(int)))


@lru_cache(maxsize=None)
def load_centromeres(path: Path = CENTROMERES_PATH) -> dict[str, tuple[int, int]]:
    """chrom -> (centromere_start, centromere_end) (bp), UCSC gap table, hg19."""
    spans = read_csv(path, sep="\t")
    return {
        chrom: (int(start), int(end))
        for chrom, start, end in zip(spans["chrom"], spans["start"], spans["end"])
    }


def masked(
    chrom: str,
    start: int,
    end: int,
    margin_bp: int = DEFAULT_MASK_MARGIN_BP,
    chrom_sizes: dict[str, int] | None = None,
    centromeres: dict[str, tuple[int, int]] | None = None,
) -> bool:
    """True if a bin `[start, end)` on `chrom` falls within `margin_bp` of a centromere
    or either chromosome end.

    `margin_bp <= 0` disables masking entirely (returns False for every bin), which is
    how callers implement a "0 disables" config option without a separate code path.
    """
    if margin_bp <= 0:
        return False

    chrom_sizes = chrom_sizes if chrom_sizes is not None else load_chrom_sizes()
    centromeres = centromeres if centromeres is not None else load_centromeres()

    if chrom not in chrom_sizes:
        raise KeyError(f"{chrom!r} not in vendored hg19 chrom sizes ({CHROM_SIZES_PATH})")

    chrom_len = chrom_sizes[chrom]
    if start < margin_bp or end > chrom_len - margin_bp:
        return True

    centromere = centromeres.get(chrom)
    if centromere is not None:
        centromere_start, centromere_end = centromere
        if end > centromere_start - margin_bp and start < centromere_end + margin_bp:
            return True

    return False


def mask_frame(
    df: DataFrame,
    chrom_col: str = "chrom",
    start_col: str = "start",
    end_col: str | None = "end",
    bin_size_bp: int | None = None,
    margin_bp: int = DEFAULT_MASK_MARGIN_BP,
    chrom_sizes: dict[str, int] | None = None,
    centromeres: dict[str, tuple[int, int]] | None = None,
) -> Series:
    """Vectorized `masked()`: a boolean Series aligned to `df`'s index.

    `end_col` is used when present (e.g. `window_fragment_counts.tsv`'s `end`, which
    already clips a trailing partial window to its contig's real length). When a table
    carries only a bin's `start` (e.g. `cnv_bins.csv`'s fixed-size ichorCNA bins),
    pass `end_col=None` and `bin_size_bp` instead -- the bin's end is computed as
    `start + bin_size_bp`, clipped to the contig length the same way
    `scripts/make_genome_windows.py` clips a trailing window.
    """
    if margin_bp <= 0:
        return Series(False, index=df.index)

    chrom_sizes = chrom_sizes if chrom_sizes is not None else load_chrom_sizes()
    centromeres = centromeres if centromeres is not None else load_centromeres()

    unknown = sorted(set(df[chrom_col]) - set(chrom_sizes))
    if unknown:
        raise KeyError(f"{unknown} not in vendored hg19 chrom sizes ({CHROM_SIZES_PATH})")

    chrom_len = df[chrom_col].map(chrom_sizes)
    start = df[start_col]
    if end_col is not None:
        end = df[end_col]
    else:
        if bin_size_bp is None:
            raise ValueError("mask_frame needs either end_col or bin_size_bp")
        end = (start + bin_size_bp).clip(upper=chrom_len)

    near_chrom_edge = (start < margin_bp) | (end > chrom_len - margin_bp)

    centromere_start = df[chrom_col].map(lambda chrom: centromeres.get(chrom, (None, None))[0])
    centromere_end = df[chrom_col].map(lambda chrom: centromeres.get(chrom, (None, None))[1])
    has_centromere = centromere_start.notna()
    near_centromere = has_centromere & (end > centromere_start - margin_bp) & (start < centromere_end + margin_bp)

    return near_chrom_edge | near_centromere.fillna(False)
