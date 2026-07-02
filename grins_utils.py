#!/usr/bin/env python3
# Copyright (C) 2019-2021 Aleksandra Nivina, Sur Herrera Paredes
# Refactored utility functions for the GRINS detection workflow.
# Licensed under GPLv3 or later, consistent with the original GRINS scripts.

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation, SeqFeature


GRINS_SUMMARY_HEADER = [
    "Genome",
    "Genome length",
    "Number of contigs in the genome",
    "Number of duplications in the genome",
    "Number of BGCSs in the genome",
    "Total length of BGCs in the genome",
    "Total length of T1 PKS BGCs in the genome",
    "Total length of NRPS BGCs in the genome",
    "Total length of terpene BGCs in the genome",
    "Total length of lanthipeptide BGCs in the genome",
    "Total length of ladderane BGCs in the genome",
    "Total length of nucleoside BGCs in the genome",
    "Number of GRINS detected in the genome",
    "Number of GRINS detected in CDSs",
    "Number of GRINS detected in BGCs",
    "Number of GRINS detected in T1 PKS",
    "Number of GRINS detected in NRPS",
    "Number of GRINS detected in terpene clusters",
    "Number of GRINS detected in lanthipeptide clusters",
    "Number of GRINS detected in ladderane clusters",
    "Number of GRINS detected in nucleoside clusters",
]


def ensure_dir(path: str | Path) -> Path:
    """Create a directory if needed and return it as a Path."""
    out = Path(path)
    out.mkdir(parents=True, exist_ok=True)
    return out


def strip_version(accession: str) -> str:
    """Match the original scripts by comparing accessions before the first dot."""
    return accession.split(".", 1)[0]


def first_qualifier(feature: SeqFeature, key: str, default: str = "") -> str:
    values = feature.qualifiers.get(key)
    if not values:
        return default
    return str(values[0])


def use_plots(value: str | bool) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"yes", "y", "true", "1"}


# ---------------------------------------------------------------------------
# split-windows
# ---------------------------------------------------------------------------


def split_sequence_into_windows(
    input_path: str | Path,
    input_format: str = "fasta",
    window_size: int = 150,
    step_size: int = 30,
    output_path: Optional[str | Path] = None,
    min_window_size: int = 50,
) -> Path:
    """Split sequence records into overlapping windows for Bowtie2 alignment.

    The query IDs must remain `record|start|end`, because the BAM parsing step
    uses those fields to reconstruct duplicated regions.
    """
    input_path = Path(input_path)
    if output_path is None:
        output_path = input_path.with_name(input_path.stem + "_windows.fasta")
    else:
        output_path = Path(output_path)

    n_windows = 0
    with output_path.open("w") as out_handle:
        for record in SeqIO.parse(str(input_path), input_format):
            seq_len = len(record.seq)
            for start in range(0, seq_len, step_size):
                end = min(start + window_size, seq_len)
                if (end - start) < min_window_size:
                    break

                window = record[start:end]
                window.id = f"{record.id}|{start}|{end}"
                window.name = window.id
                window.description = ""
                SeqIO.write(window, out_handle, "fasta")
                n_windows += 1

    print(f"Wrote {n_windows} sequence windows to {output_path}")
    return output_path


# ---------------------------------------------------------------------------
# antismash-to-gff3
# ---------------------------------------------------------------------------


def _parse_antismash_location(location: str) -> Tuple[int, int]:
    """Return 1-based inclusive GFF coordinates from an antiSMASH location.

    antiSMASH region locations are usually simple strings like `[0:12345]`.
    This parser accepts that form and other strings containing at least two
    coordinates, while preserving the original start + 1 conversion.
    """
    numbers = [int(x) for x in re.findall(r"\d+", location)]
    if len(numbers) < 2:
        raise ValueError(f"Could not parse antiSMASH location: {location!r}")
    start0, end = numbers[0], numbers[1]
    return min(start0, end) + 1, max(start0, end)


def antismash_json_to_gff3(input_json: str | Path, output_gff3: str | Path = "") -> Path:
    """Convert antiSMASH JSON region features into a small GFF3 file."""
    input_json = Path(input_json)
    if output_gff3 == "" or output_gff3 is None:
        output_gff3 = input_json.with_suffix(".gff3").name
    output_gff3 = Path(output_gff3)

    with input_json.open("r") as in_handle:
        print("Reading antiSMASH JSON file")
        asmash = json.load(in_handle)

    n_regions = 0
    with output_gff3.open("w") as out_handle:
        out_handle.write("##gff-version 3\n")
        for record in asmash.get("records", []):
            record_id = record.get("id", "unknown_record")
            print(f"Processing record {record_id}...")

            for feature in record.get("features", []):
                if feature.get("type") != "region":
                    continue

                n_regions += 1
                start, end = _parse_antismash_location(feature.get("location", ""))
                qualifiers = feature.get("qualifiers", {})
                region_number = qualifiers.get("region_number", [str(n_regions)])[0]
                products = qualifiers.get("product", ["unknown"])

                attributes = ";".join([
                    f"ID=region_{region_number}",
                    "Product=" + ",".join(products),
                ])
                row = [
                    record_id,
                    "antiSMASH",
                    "region",
                    str(start),
                    str(end),
                    ".",
                    "+",
                    ".",
                    attributes,
                ]
                out_handle.write("\t".join(row) + "\n")

    print(f"Found {n_regions} regions.")
    return output_gff3


# ---------------------------------------------------------------------------
# bam-to-duplicates
# ---------------------------------------------------------------------------


class Window:
    """A genomic interval represented with 0-based start and 1-based-style end."""

    def __init__(self, start: int, end: int):
        if not isinstance(start, int) or not isinstance(end, int):
            raise ValueError("start and end must be integers.")
        if start < 0 or end < 0:
            raise ValueError("start and end must be non-negative.")
        if start > end:
            raise ValueError("start cannot be greater than end.")
        self.start = start
        self.end = end


class WindowCollection:
    """Merge overlapping windows for one sequence record."""

    def __init__(self):
        self.windows: List[Window] = []

    def find_overlaps(self, window: Window) -> List[int]:
        overlaps = []
        for idx, existing in enumerate(self.windows):
            if window.start < existing.end and existing.start < window.end:
                overlaps.append(idx)
        return overlaps

    def remove_windows(self, indexes: Sequence[int]) -> None:
        for idx in sorted(indexes, reverse=True):
            del self.windows[idx]

    def add_window(self, window: Window) -> None:
        overlaps = self.find_overlaps(window)
        if not overlaps:
            self.windows.append(window)
            return

        start = window.start
        end = window.end
        for idx in overlaps:
            start = min(start, self.windows[idx].start)
            end = max(end, self.windows[idx].end)

        self.remove_windows(overlaps)
        self.add_window(Window(start, end))


BamWindow = Tuple[str, int, int, str, int, int, int]


def find_bam_windows(input_bam: str | Path, min_alignment_size: int = 150) -> List[BamWindow]:
    """Read a BAM file and collect non-self duplicated window alignments."""
    try:
        import pysam
    except ImportError as exc:
        raise ImportError("pysam is required for the bam-to-duplicates command.") from exc

    bam_windows: List[BamWindow] = []
    with pysam.AlignmentFile(str(input_bam), "rb") as samfile:
        for read in samfile.fetch(until_eof=True):
            if read.is_unmapped or read.reference_name is None:
                continue

            parts = read.query_name.split("|")
            if len(parts) != 3:
                raise ValueError(
                    f"Malformed query name {read.query_name!r}; expected reference|start|end."
                )

            query_ref, query_start, query_end = parts
            query_start = int(query_start)
            query_end = int(query_end)
            alignment_length = query_end - query_start
            reference_start = int(read.reference_start)

            # Discard exact self-maps but keep matches to other locations/records.
            if alignment_length >= min_alignment_size and reference_start > 0:
                if query_ref != read.reference_name or reference_start != query_start:
                    bam_windows.append((
                        query_ref,
                        query_start,
                        query_end,
                        read.reference_name,
                        reference_start,
                        reference_start + min_alignment_size,
                        int(read.mapping_quality),
                    ))

    return bam_windows


def merge_bam_windows(bam_windows: Iterable[BamWindow]) -> Dict[str, WindowCollection]:
    """Merge query and target intervals by sequence record."""
    merged: Dict[str, WindowCollection] = {}
    for query_ref, query_start, query_end, target_ref, target_start, target_end, _mapq in bam_windows:
        merged.setdefault(query_ref, WindowCollection()).add_window(Window(query_start, query_end))
        merged.setdefault(target_ref, WindowCollection()).add_window(Window(target_start, target_end))
    return merged


def write_duplicate_gff3(
    windows_by_record: Dict[str, WindowCollection],
    output_gff3: str | Path,
    min_size: int = 0,
) -> Path:
    """Write merged duplicated regions as GFF3."""
    output_gff3 = Path(output_gff3)
    n_written = 0
    with output_gff3.open("w") as out_handle:
        out_handle.write("##gff-version 3\n")
        for record_id in sorted(windows_by_record):
            for window in windows_by_record[record_id].windows:
                if (window.end - window.start) < min_size:
                    continue
                n_written += 1
                row = [
                    record_id,
                    "bowtie2",
                    "duplicate",
                    str(window.start + 1),
                    str(window.end),
                    ".",
                    "+",
                    ".",
                    f"ID=dup_{n_written}",
                ]
                out_handle.write("\t".join(row) + "\n")

    print(f"Wrote {n_written} duplicated regions.")
    return output_gff3


def bam_to_duplicate_gff3(
    input_bam: str | Path,
    output_gff3: str | Path,
    window_size: int = 150,
    min_size: int = 0,
) -> Path:
    bam_windows = find_bam_windows(input_bam, min_alignment_size=window_size)
    merged_windows = merge_bam_windows(bam_windows)
    return write_duplicate_gff3(merged_windows, output_gff3, min_size=min_size)


# ---------------------------------------------------------------------------
# detect-grins
# ---------------------------------------------------------------------------


def gc_skew(input_dna: str, window: int = 150, step: int = 30) -> List[float]:
    sequence = str(input_dna).upper()
    values = []
    for start in range(0, len(sequence) - window, step):
        fragment = sequence[start:start + window]
        g_count = fragment.count("G")
        c_count = fragment.count("C")
        values.append(float(g_count - c_count) / float(g_count + c_count) if g_count + c_count else 0.0)
    return values


def ta_skew(input_dna: str, window: int = 150, step: int = 30) -> List[float]:
    sequence = str(input_dna).upper()
    values = []
    for start in range(0, len(sequence) - window, step):
        fragment = sequence[start:start + window]
        t_count = fragment.count("T")
        a_count = fragment.count("A")
        values.append(float(t_count - a_count) / float(t_count + a_count) if t_count + a_count else 0.0)
    return values


def abs_skew_means(record, start: int, end: int, window: int = 150, step: int = 30) -> Tuple[float, float]:
    sequence = str(record.seq[start:end]).upper()
    gc_values = []
    ta_values = []
    for pos in range(0, len(sequence) - window, step):
        fragment = sequence[pos:pos + window]
        g_count = fragment.count("G")
        c_count = fragment.count("C")
        t_count = fragment.count("T")
        a_count = fragment.count("A")
        gc_values.append(abs(float(g_count - c_count) / float(g_count + c_count)) if g_count + c_count else 0.0)
        ta_values.append(abs(float(t_count - a_count) / float(t_count + a_count)) if t_count + a_count else 0.0)

    if not gc_values or not ta_values:
        return 0.0, 0.0
    return float(np.mean(gc_values)), float(np.mean(ta_values))


def plot_grins_region(
    plot_folder: str | Path,
    assembly: str,
    accession: str,
    start: int,
    end: int,
    gc_values: Sequence[float],
    ta_values: Sequence[float],
    duplicated_starts: Sequence[int],
    duplicated_ends: Sequence[int],
    grins_starts: Sequence[int],
    grins_ends: Sequence[int],
) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.patches as patches
    import matplotlib.pyplot as plt

    x_values = np.arange(start, start + len(gc_values) * 30, 30)
    _fig, ax1 = plt.subplots(figsize=(10, 4))
    ax1.set_ylim(-1, 1)
    ax1.set_xlim(start, end)
    ax1.grid(False)
    ax1.plot(x_values, gc_values, color="blue", linewidth=0.5)
    ax1.plot(x_values, ta_values, color="red", linewidth=0.5)
    ax1.legend(["GC skew", "TA skew"], loc=4)

    for region_start, region_end in zip(duplicated_starts, duplicated_ends):
        rect = patches.Rectangle(
            (region_start, 1), region_end - region_start, -2,
            linewidth=1, edgecolor="grey", facecolor="grey", alpha=0.3,
        )
        ax1.add_patch(rect)

    for region_start, region_end in zip(grins_starts, grins_ends):
        rect = patches.Rectangle(
            (region_start, 1), region_end - region_start, -2,
            linewidth=1, edgecolor="teal", facecolor="teal", alpha=0.3,
        )
        ax1.add_patch(rect)

    plot_dir = ensure_dir(Path(plot_folder) / assembly)
    plt.title("GRINS in %s, %d-%d kb" % (accession, start / 1000, end / 1000))
    plt.savefig(plot_dir / f"{accession}_GRINS_{int(start / 1000)}-{int(end / 1000)}kb.png")
    plt.close("all")


def read_duplication_locations(dup_gff3: str | Path) -> Tuple[Dict[str, List[Tuple[int, int]]], int]:
    locations: Dict[str, List[Tuple[int, int]]] = {}
    n_dups = 0
    with Path(dup_gff3).open("r") as in_handle:
        for line in in_handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            accession = strip_version(fields[0])
            start = int(fields[3])
            end = int(fields[4])
            locations.setdefault(accession, []).append((start, end))
            n_dups += 1
    return locations, n_dups


def find_sequence_file(seq_dir: str | Path, assembly: str) -> Path:
    seq_dir = Path(seq_dir)
    files = sorted(p for p in seq_dir.iterdir() if p.is_file())

    exact_names = [
        f"{assembly}.gbk",
        f"{assembly}.gb",
        f"{assembly}.genbank",
    ]
    for name in exact_names:
        candidate = seq_dir / name
        if candidate.exists():
            return candidate

    matches = [p for p in files if assembly in p.name]
    if not matches:
        raise FileNotFoundError(f"Could not find a GenBank file for assembly {assembly!r} in {seq_dir}")
    if len(matches) > 1:
        print(f"WARNING: multiple GenBank files matched {assembly!r}; using {matches[0].name}")
    return matches[0]


def _zero_stats() -> Dict[str, int]:
    return {
        "length_genome": 0,
        "n_contigs_genome": 0,
        "n_BGC": 0,
        "length_BGC": 0,
        "n_dups": 0,
        "length_PKS": 0,
        "length_NRPS": 0,
        "length_terpene": 0,
        "length_lanthi": 0,
        "length_ladderane": 0,
        "length_nucleoside": 0,
        "n_GRINS_total": 0,
        "n_GRINS_CDS": 0,
        "n_GRINS_BGC": 0,
        "n_GRINS_PKS": 0,
        "n_GRINS_NRPS": 0,
        "n_GRINS_terpenes": 0,
        "n_GRINS_lanthi": 0,
        "n_GRINS_ladderane": 0,
        "n_GRINS_nucleoside": 0,
    }


def _update_bgc_stats(stats: Dict[str, int], bgc_type: str, bgc_len: int) -> None:
    stats["n_BGC"] += 1
    stats["length_BGC"] += bgc_len
    if bgc_type in {"T1PKS", "transAT-PKS"}:
        stats["length_PKS"] += bgc_len
    elif bgc_type == "NRPS":
        stats["length_NRPS"] += bgc_len
    elif bgc_type == "terpene":
        stats["length_terpene"] += bgc_len
    elif bgc_type == "lanthipeptide":
        stats["length_lanthi"] += bgc_len
    elif bgc_type == "ladderane":
        stats["length_ladderane"] += bgc_len
    elif bgc_type == "nucleoside":
        stats["length_nucleoside"] += bgc_len


def _update_grins_bgc_stats(stats: Dict[str, int], bgc_type: str) -> None:
    stats["n_GRINS_BGC"] += 1
    if bgc_type in {"T1PKS", "transAT-PKS"}:
        stats["n_GRINS_PKS"] += 1
    elif bgc_type == "NRPS":
        stats["n_GRINS_NRPS"] += 1
    elif bgc_type == "terpene":
        stats["n_GRINS_terpenes"] += 1
    elif bgc_type == "lanthipeptide":
        stats["n_GRINS_lanthi"] += 1
    elif bgc_type == "ladderane":
        stats["n_GRINS_ladderane"] += 1
    elif bgc_type == "nucleoside":
        stats["n_GRINS_nucleoside"] += 1


def _contained_with_flank(start: int, end: int, feature_start: int, feature_end: int, flank: int) -> bool:
    return (
        start >= feature_start - flank
        and start < feature_end + flank
        and end > feature_start - flank
        and end <= feature_end + flank
    )


def _plot_record_windows(
    record,
    assembly: str,
    plot_output: str | Path,
    duplicate_locations: Sequence[Tuple[int, int]],
    grins_locations: Sequence[Tuple[int, int]],
) -> None:
    for start in range(0, len(record.seq), 100000):
        end = start + 100000

        dup_starts: List[int] = []
        dup_ends: List[int] = []
        for dup_start, dup_end in duplicate_locations:
            if (start < dup_start < end) or (start < dup_end < end):
                dup_starts.append(max(start, dup_start))
                dup_ends.append(min(end, dup_end))

        grins_starts: List[int] = []
        grins_ends: List[int] = []
        for grins_start, grins_end in grins_locations:
            if (start < grins_start < end) or (start < grins_end < end):
                grins_starts.append(max(start, grins_start))
                grins_ends.append(min(end, grins_end))

        plot_grins_region(
            plot_output,
            assembly,
            record.id,
            start,
            end,
            gc_skew(str(record.seq[start:end])),
            ta_skew(str(record.seq[start:end])),
            dup_starts,
            dup_ends,
            grins_starts,
            grins_ends,
        )


def detect_grins_from_bowtie(
    seq_input: str | Path,
    dupl_input: str | Path,
    seq_output: str | Path = "./output/genomes_GRINS",
    grins_output: str | Path = "./output/GRINS.gff3",
    grins_bgc_output: str | Path = "./output/GRINS_BGC.gff3",
    plot_output: str | Path = "./output/plots",
    with_plots: str | bool = "no",
    summary_output: str | Path = "GRINS_detected_in_genomes_and_BGCs.txt",
    flank_size: int = 300,
    min_grins_size: int = 500,
    gc_threshold: float = 0.15,
    ta_threshold: float = 0.15,
) -> Path:
    """Detect GRINS from antiSMASH GenBank files and duplicated-region GFF3 files."""
    seq_input = Path(seq_input)
    dupl_input = Path(dupl_input)
    seq_output = ensure_dir(seq_output)
    grins_output = ensure_dir(grins_output)
    grins_bgc_output = ensure_dir(grins_bgc_output)
    plot_output = Path(plot_output)
    make_plots = use_plots(with_plots)
    if make_plots:
        ensure_dir(plot_output)

    duplication_files = sorted(
        p for p in dupl_input.iterdir()
        if p.is_file() and p.name.endswith(".duplicated.gff3")
    )
    if not duplication_files:
        raise FileNotFoundError(f"No .duplicated.gff3 files found in {dupl_input}")

    summary_output = Path(summary_output)
    with summary_output.open("w") as summary_handle:
        summary_handle.write("\t".join(GRINS_SUMMARY_HEADER) + "\n")

        for dup_file in duplication_files:
            assembly = dup_file.name.replace(".duplicated.gff3", "")
            stats = _zero_stats()
            duplicate_locations_by_record, stats["n_dups"] = read_duplication_locations(dup_file)
            seq_file = find_sequence_file(seq_input, assembly)
            records = SeqIO.parse(str(seq_file), "gb")

            with (seq_output / f"{assembly}.GRINS.gbk").open("w") as annotated_handle, \
                 (grins_bgc_output / f"{assembly}.GRINS_BGC.gff3").open("w") as grins_bgc_handle, \
                 (grins_output / f"{assembly}.GRINS_CDS.gff3").open("w") as grins_cds_handle, \
                 (grins_output / f"{assembly}.GRINS.gff3").open("w") as grins_handle:

                grins_bgc_handle.write("\t".join(["Record", "GRINS start", "GRINS end", "BGC"]) + "\n")
                grins_cds_handle.write("\t".join(["Record", "GRINS start", "GRINS end", "Locus", "CDS name"]) + "\n")
                grins_handle.write("\t".join(["Record", "GRINS start", "GRINS end"]) + "\n")

                for record in records:
                    stats["length_genome"] += len(record.seq)
                    stats["n_contigs_genome"] += 1
                    record_features = record.features
                    accession = strip_version(record.id)

                    for feature in record_features:
                        if feature.type != "region":
                            continue
                        bgc_start = min(int(feature.location.start), int(feature.location.end))
                        bgc_end = max(int(feature.location.start), int(feature.location.end))
                        bgc_type = first_qualifier(feature, "product", "unknown")
                        _update_bgc_stats(stats, bgc_type, bgc_end - bgc_start)

                    duplicate_locations = duplicate_locations_by_record.get(accession, [])
                    grins_locations: List[Tuple[int, int]] = []

                    for dup_start, dup_end in duplicate_locations:
                        record.features.append(
                            SeqFeature(FeatureLocation(start=dup_start, end=dup_end), type="Duplication")
                        )

                        if dup_end - dup_start < min_grins_size:
                            continue

                        mean_gc_skew, mean_ta_skew = abs_skew_means(record, dup_start, dup_end)
                        if mean_gc_skew < gc_threshold or mean_ta_skew < ta_threshold:
                            continue

                        stats["n_GRINS_total"] += 1
                        record.features.append(
                            SeqFeature(FeatureLocation(start=dup_start, end=dup_end), type="GRINS")
                        )
                        grins_start = int(dup_start)
                        grins_end = int(dup_end)
                        grins_locations.append((grins_start, grins_end))

                        in_bgc = False
                        current_bgc_type = ""
                        for feature in record_features:
                            if feature.type != "region":
                                continue
                            bgc_start = min(int(feature.location.start), int(feature.location.end))
                            bgc_end = max(int(feature.location.start), int(feature.location.end))
                            if _contained_with_flank(grins_start, grins_end, bgc_start, bgc_end, flank_size):
                                in_bgc = True
                                current_bgc_type = first_qualifier(feature, "product", "unknown")
                                _update_grins_bgc_stats(stats, current_bgc_type)
                                break

                        in_cds = False
                        current_cds = None
                        for feature in record_features:
                            if feature.type != "CDS":
                                continue
                            cds_start = min(int(feature.location.start), int(feature.location.end))
                            cds_end = max(int(feature.location.start), int(feature.location.end))
                            if (
                                _contained_with_flank(grins_start, grins_end, cds_start, cds_end, flank_size)
                                and (cds_end - cds_start) < 500000
                            ):
                                in_cds = True
                                current_cds = feature
                                stats["n_GRINS_CDS"] += 1
                                break

                        grins_handle.write("\t".join([
                            record.id, "GRINSdetect", "GRINS", str(grins_start), str(grins_end), ".", "+", ".", ""
                        ]) + "\n")

                        if in_bgc:
                            grins_bgc_handle.write("\t".join([
                                record.id, "GRINSdetect", "GRINS", str(grins_start), str(grins_end), ".", "+", ".", current_bgc_type
                            ]) + "\n")

                        if in_cds and current_cds is not None:
                            locus_tag = first_qualifier(current_cds, "locus_tag", "N/a")
                            gene_function = first_qualifier(current_cds, "gene_functions", "N/a")
                            grins_cds_handle.write("\t".join([
                                record.id,
                                "GRINSdetect",
                                "GRINS",
                                str(grins_start),
                                str(grins_end),
                                ".",
                                "+",
                                ".",
                                f"locus_tag={locus_tag},gene_function={gene_function}",
                            ]) + "\n")

                    SeqIO.write(record, annotated_handle, "genbank")

                    if make_plots:
                        _plot_record_windows(record, assembly, plot_output, duplicate_locations, grins_locations)

            summary_values = [
                assembly,
                str(stats["length_genome"]),
                str(stats["n_contigs_genome"]),
                str(stats["n_dups"]),
                str(stats["n_BGC"]),
                str(stats["length_BGC"]),
                str(stats["length_PKS"]),
                str(stats["length_NRPS"]),
                str(stats["length_terpene"]),
                str(stats["length_lanthi"]),
                str(stats["length_ladderane"]),
                str(stats["length_nucleoside"]),
                str(stats["n_GRINS_total"]),
                str(stats["n_GRINS_CDS"]),
                str(stats["n_GRINS_BGC"]),
                str(stats["n_GRINS_PKS"]),
                str(stats["n_GRINS_NRPS"]),
                str(stats["n_GRINS_terpenes"]),
                str(stats["n_GRINS_lanthi"]),
                str(stats["n_GRINS_ladderane"]),
                str(stats["n_GRINS_nucleoside"]),
            ]
            summary_handle.write("\t".join(summary_values) + "\n")

    return summary_output
