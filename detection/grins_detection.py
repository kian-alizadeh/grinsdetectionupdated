#!/usr/bin/env python3
# Command-line entry point for GRINS detection helper functions.
# Licensed under GPLv3 or later, consistent with the original GRINS scripts.

from __future__ import annotations

import argparse

from grins_utils import (
    STEP_SIZE,
    WINDOW_SIZE,
    antismash_json_to_gff3,
    bam_to_duplicate_gff3,
    detect_grins_from_bowtie,
    split_sequence_into_windows,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="GRINS detection utilities used by the Nextflow workflow.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    split_parser = subparsers.add_parser(
        "split-windows",
        help="Split FASTA/GenBank records into overlapping FASTA windows.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    split_parser.add_argument("--input", required=True)
    split_parser.add_argument("--format", default="fasta")
    split_parser.add_argument("--w_size", type=int, required=True)
    split_parser.add_argument("--s_size", type=int, required=True)
    split_parser.add_argument("--output", default=None)
    split_parser.set_defaults(func=run_split_windows)

    asmash_parser = subparsers.add_parser(
        "antismash-to-gff3",
        help="Convert antiSMASH JSON region features to GFF3.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    asmash_parser.add_argument("--input", required=True, help="antiSMASH JSON file")
    asmash_parser.add_argument("--output", default="", help="Output GFF3 path")
    asmash_parser.set_defaults(func=run_antismash_to_gff3)

    bam_parser = subparsers.add_parser(
        "bam-to-duplicates",
        help="Merge repeated-window BAM alignments into duplicated-region GFF3.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    bam_parser.add_argument("--input", required=True, help="BAM file from Bowtie2/Samtools")
    bam_parser.add_argument("--output", default="grins.gff3", help="Output GFF3 path")
    bam_parser.add_argument("--w_size", type=int, default=WINDOW_SIZE, help="Window size used during splitting")
    bam_parser.add_argument("--min_size", type=int, default=0, help="Minimum duplicated-region size to write")
    bam_parser.set_defaults(func=run_bam_to_duplicates)

    grins_parser = subparsers.add_parser(
        "detect-grins",
        help="Detect GRINS from antiSMASH GenBank files and duplicated-region GFF3 files.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    grins_parser.add_argument("--seq_input", required=True, help="Folder with GenBank sequences")
    grins_parser.add_argument("--dupl_input", required=True, help="Folder with duplicated-region GFF3 files")
    grins_parser.add_argument("--seq_output", default="./output/genomes_GRINS", help="Annotated GenBank output folder")
    grins_parser.add_argument("--GRINS_output", "--grins_output", dest="grins_output", default="./output/GRINS.gff3", help="GRINS output folder")
    grins_parser.add_argument("--GRINS_BGC_output", "--grins_bgc_output", dest="grins_bgc_output", default="./output/GRINS_BGC.gff3", help="GRINS-in-BGC output folder")
    grins_parser.add_argument("--with_plots", default="no", help="Write skew plots: yes/no")
    grins_parser.add_argument("--plot_output", default="./output/plots", help="Plot output folder")
    grins_parser.add_argument("--summary_output", default="GRINS_detected_in_genomes_and_BGCs.txt", help="Summary table output file")
    grins_parser.add_argument("--skew_w_size", type=int, default=WINDOW_SIZE, help="Window size for GC/TA skew calculations")
    grins_parser.add_argument("--skew_s_size", type=int, default=STEP_SIZE, help="Step size for GC/TA skew calculations")
    grins_parser.set_defaults(func=run_detect_grins)

    return parser


def run_split_windows(args: argparse.Namespace) -> None:
    split_sequence_into_windows(
        input_path=args.input,
        input_format=args.format,
        window_size=args.w_size,
        step_size=args.s_size,
        output_path=args.output,
    )


def run_antismash_to_gff3(args: argparse.Namespace) -> None:
    antismash_json_to_gff3(input_json=args.input, output_gff3=args.output)


def run_bam_to_duplicates(args: argparse.Namespace) -> None:
    bam_to_duplicate_gff3(
        input_bam=args.input,
        output_gff3=args.output,
        window_size=args.w_size,
        min_size=args.min_size,
    )


def run_detect_grins(args: argparse.Namespace) -> None:
    detect_grins_from_bowtie(
        seq_input=args.seq_input,
        dupl_input=args.dupl_input,
        seq_output=args.seq_output,
        grins_output=args.grins_output,
        grins_bgc_output=args.grins_bgc_output,
        plot_output=args.plot_output,
        with_plots=args.with_plots,
        summary_output=args.summary_output,
        skew_window_size=args.skew_w_size,
        skew_step_size=args.skew_s_size,
    )


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
