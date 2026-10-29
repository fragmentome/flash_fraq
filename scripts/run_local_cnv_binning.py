"""Driver for MAT-16 1 Mb CN binning over one staged sample
(docs/phase2-preprocessing-plan.md §8b.3 step 6, CN half).

Two deviations from IchorCnaRunner.run()'s default four-tool sequence, both
forced by what is actually on this workstation:

- **gcCounter and mapCounter do not run.** ichorCNA ships hg19 1 Mb GC and
  mappability WIGs in its R package extdata, so gcCounter would regenerate a
  byte-identical GC file once per sample, and mapCounter cannot run at all --
  it needs a mappability bigWig that does not exist locally. The shipped WIGs
  are passed straight through instead.
- **The input is a BAM, not the CRAM.** readCounter is BamTools-based and
  cannot open a CRAM; workflow/rules/local_features.smk decodes a temp BAM on
  the SSD for this rule alone.

Invoked per-sample by the localfeat_ichor rule.
"""

from __future__ import annotations

from argparse import ArgumentParser
from pathlib import Path

from flash_fraq.cnv_binning import IchorCnaConfig, IchorCnaRunner


def main() -> None:
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--alignment-path", required=True, type=Path, help="staged .bam (readCounter cannot read CRAM)")
    parser.add_argument("--reference-fasta", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--gc-wig", required=True, type=Path)
    parser.add_argument("--map-wig", required=True, type=Path)
    parser.add_argument("--centromere", required=True, type=Path)
    parser.add_argument("--ichorcna-script", required=True, type=Path)
    parser.add_argument("--ichorcna-libdir", type=Path, default=None)
    parser.add_argument("--readcounter-bin", required=True, type=Path)
    parser.add_argument("--rscript-bin", required=True, type=Path)
    parser.add_argument("--bin-size-bp", type=int, default=1_000_000)
    parser.add_argument("--genome-build", default="hg19")
    parser.add_argument("--genome-style", default="UCSC")
    parser.add_argument("--chromosomes", required=True, help="comma-separated, chr-prefixed")
    parser.add_argument("--normal-panel", type=Path, default=None)
    # R-vector literals, passed through to runIchorCNA.R verbatim.
    parser.add_argument("--ploidy", default="c(2,3)")
    parser.add_argument("--normal", default="c(0.5,0.6,0.7,0.8,0.9)")
    parser.add_argument("--sc-states", default="c(1,3)")
    parser.add_argument(
        "--estimate-sc-prevalence", default="TRUE", choices=["TRUE", "FALSE"],
        help="R logical literal, passed through to runIchorCNA.R's --estimateScPrevalence",
    )
    parser.add_argument("--max-cn", type=int, default=5)
    parser.add_argument("--txn-e", type=float, default=0.9999)
    parser.add_argument("--txn-strength", type=int, default=10000)
    args = parser.parse_args()

    runner = IchorCnaRunner(
        IchorCnaConfig(
            bam_path=args.alignment_path,
            reference_fasta=args.reference_fasta,
            # Unused on this path: mapCounter never runs, the shipped
            # mappability WIG is supplied directly instead.
            mappability_bigwig=args.map_wig,
            output_dir=args.output_dir,
            sample_id=args.sample_id,
            bin_size_bp=args.bin_size_bp,
            chromosomes=tuple(args.chromosomes.split(",")),
            readcounter_bin=args.readcounter_bin,
            rscript_bin=args.rscript_bin,
            ichorcna_script=args.ichorcna_script,
            libdir=args.ichorcna_libdir,
            gc_wig_source=args.gc_wig,
            map_wig_source=args.map_wig,
            centromere=args.centromere,
            genome_build=args.genome_build,
            genome_style=args.genome_style,
            normal_panel=args.normal_panel,
            ploidy=args.ploidy,
            normal=args.normal,
            sc_states=args.sc_states,
            estimate_sc_prevalence=args.estimate_sc_prevalence == "TRUE",
            max_cn=args.max_cn,
            txn_e=args.txn_e,
            txn_strength=args.txn_strength,
        )
    )
    runner.run()

    print(f"cnv bins: {runner.config.cnv_bins_output_path}")
    print(f"ichor estimate: {runner.config.ichor_estimate_output_path}")


if __name__ == "__main__":
    main()
