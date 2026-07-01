#!/usr/bin/env python3

import argparse
from pathlib import Path
from Bio import SeqIO


def parse_args():
    parser = argparse.ArgumentParser(
        description="Split FASTA/GenBank sequences into overlapping windows."
    )
    parser.add_argument("--input", required=True)
    parser.add_argument("--format", default="fasta")
    parser.add_argument("--w_size", type=int, required=True)
    parser.add_argument("--s_size", type=int, required=True)
    parser.add_argument("--output", default=None)
    return parser.parse_args()


def main():
    args = parse_args()

    input_path = Path(args.input)

    if args.output:
        output_path = Path(args.output)
    else:
        output_path = input_path.with_name(input_path.stem + "_windows.fasta")

    records_out = []

    for record in SeqIO.parse(str(input_path), args.format):
        seq_len = len(record.seq)

        for pos in range(0, seq_len, args.s_size):
            end = min(pos + args.w_size, seq_len)

            if (end - pos) < 50:
                break

            window = record[pos:end]

            # IMPORTANT:
            # produce_windows_from_bam.py expects exactly:
            # query_name = reference|window_start|window_end
            window.id = f"{record.id}|{pos}|{end}"
            window.name = window.id
            window.description = ""

            records_out.append(window)

    SeqIO.write(records_out, str(output_path), "fasta")


if __name__ == "__main__":
    main()
