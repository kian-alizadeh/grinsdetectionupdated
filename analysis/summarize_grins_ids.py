#!/usr/bin/env python3
"""
summarize_grins_ids.py

Post-process GRINS workflow outputs after the GRINS ID modification.

Expected workflow layout for each run folder:

    RUN_FOLDER/
        outputs/
            ASSEMBLY/
                genomes_GRINS/ASSEMBLY.GRINS.gbk
                GRINS.gff3/ASSEMBLY.GRINS.gff3
                GRINS.gff3/ASSEMBLY.GRINS_CDS.gff3
                GRINS_BGC.gff3/ASSEMBLY.GRINS_BGC.gff3
                plots/
        res/
            ASSEMBLY.txt

The script can also accept RUN_FOLDER/outputs directly.

Outputs, per assembly/output folder:
    ASSEMBLY_grins_copies.tsv                  One row per GRINS feature/copy.
    ASSEMBLY_grins_groups.tsv                  One row per GRINS ID/family.
    ASSEMBLY_compound_grins_blocks.tsv         One row per recurring adjacent GRINS block.
    ASSEMBLY_compound_grins_block_copies.tsv   One row per observed block copy.
    ASSEMBLY_grins_copies.fasta                Optional FASTA of GRINS copy sequences.

Requires: Biopython
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple
from urllib.parse import unquote

try:
    from Bio import SeqIO
except ImportError as exc:
    raise SystemExit(
        "Biopython is required. Install it in the environment you use for GRINS, e.g. `pip install biopython`."
    ) from exc


DOMAIN_FEATURE_TYPES = {
    "aSDomain",
    "PFAM_domain",
    "CDS_motif",
    "domain",
    "protein_domain",
    "motif",
}

CDS_TYPES = {"CDS"}
BGC_TYPES = {"region"}
GRINS_TYPES = {"GRINS"}

COPY_COLUMNS = [
    "run",
    "assembly",
    "record_id",
    "grins_id",
    "grins_feature_id",
    "copy_number_within_grins_id",
    "start_1based",
    "end_1based_inclusive",
    "length_bp",
    "strand",
    "dup_id",
    "dup_group_id",
    "bgc_region_id",
    "bgc_product",
    "bgc_start_1based",
    "bgc_end_1based_inclusive",
    "distance_to_left_bgc_edge_bp",
    "distance_to_right_bgc_edge_bp",
    "nearest_bgc_edge_distance_bp",
    "overlapping_cds_count",
    "overlapping_cds_locus_tags",
    "overlapping_cds_products",
    "overlapping_domain_count",
    "overlapping_domains",
    "nearby_domains",
    "domain_context",
    "nearest_domain_distance_bp",
    "nearest_cds_distance_bp",
    "sequence_preview_30bp",
]

GROUP_COLUMNS = [
    "run",
    "assembly",
    "grins_id",
    "n_copies",
    "record_ids",
    "copy_coordinates",
    "copy_lengths_bp",
    "min_copy_length_bp",
    "max_copy_length_bp",
    "mean_copy_length_bp",
    "dup_group_ids",
    "n_bgc_copies",
    "bgc_products",
    "n_cds_overlapping_copies",
    "cds_locus_tags",
    "n_domain_overlapping_copies",
    "domain_names",
    "domain_contexts",
    "min_intercopy_distance_bp_same_record",
    "max_span_bp_same_record",
    "arrangement",
    "priority_flags",
]


BLOCK_COLUMNS = [
    "run",
    "assembly",
    "block_id",
    "member_grins_ids",
    "n_member_regions",
    "n_block_copies",
    "records",
    "block_copy_coordinates",
    "member_copy_coordinates",
    "mean_block_length_bp",
    "min_block_length_bp",
    "max_block_length_bp",
    "mean_internal_gap_bp",
    "internal_gaps_by_copy",
    "bgc_products",
    "domain_contexts",
    "domain_names",
    "cds_locus_tags",
    "arrangement",
    "interpretation_hint",
]


BLOCK_COPY_COLUMNS = [
    "run",
    "assembly",
    "block_id",
    "block_copy_number",
    "record_id",
    "block_start_1based",
    "block_end_1based_inclusive",
    "block_length_bp",
    "member_grins_ids",
    "member_grins_feature_ids",
    "member_coordinates",
    "internal_gaps_bp",
    "bgc_products",
    "domain_contexts",
    "domain_names",
    "cds_locus_tags",
]


def strip_version(accession: str) -> str:
    return accession.split(".", 1)[0]


def clean(value: object) -> str:
    if value is None:
        return ""
    text = str(value).replace("\n", " ").replace("\r", " ").replace("\t", " ")
    return re.sub(r"\s+", " ", text).strip()


def clean_join(values: Iterable[object], sep: str = "|") -> str:
    cleaned: List[str] = []
    seen = set()
    for value in values:
        text = clean(value)
        if not text or text in seen:
            continue
        cleaned.append(text)
        seen.add(text)
    return sep.join(cleaned)


def first_qual(feature, keys: Sequence[str], default: str = "") -> str:
    for key in keys:
        values = feature.qualifiers.get(key)
        if values:
            return clean(values[0])
    return default


def all_quals(feature, keys: Sequence[str]) -> List[str]:
    out: List[str] = []
    for key in keys:
        values = feature.qualifiers.get(key)
        if values:
            out.extend(clean(v) for v in values if clean(v))
    return out


def feature_bounds(feature) -> Tuple[int, int]:
    """Return 0-based, half-open [start, end) bounds for simple or compound features."""
    return int(feature.location.start), int(feature.location.end)


def strand_string(feature) -> str:
    strand = getattr(feature.location, "strand", None)
    if strand == 1:
        return "+"
    if strand == -1:
        return "-"
    return "."


def intervals_overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> bool:
    return a_start < b_end and b_start < a_end


def overlap_bp(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    return max(0, min(a_end, b_end) - max(a_start, b_start))


def distance_between_intervals(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    if intervals_overlap(a_start, a_end, b_start, b_end):
        return 0
    if a_end <= b_start:
        return b_start - a_end
    return a_start - b_end


def parse_gff3_attributes(attr_text: str) -> Dict[str, str]:
    attrs: Dict[str, str] = {}
    for part in attr_text.strip().split(";"):
        if not part or "=" not in part:
            continue
        key, value = part.split("=", 1)
        attrs[unquote(key)] = unquote(value)
    return attrs


def read_grins_gff_attrs(gff_path: Optional[Path]) -> Dict[Tuple[str, int, int], Dict[str, str]]:
    """Read GRINS GFF3 attributes keyed by GFF coordinates, 1-based inclusive."""
    attrs_by_coord: Dict[Tuple[str, int, int], Dict[str, str]] = {}
    if not gff_path or not gff_path.exists():
        return attrs_by_coord

    with gff_path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            record_id = fields[0]
            feature_type = fields[2]
            if feature_type not in GRINS_TYPES:
                continue
            try:
                start = int(fields[3])
                end = int(fields[4])
            except ValueError:
                continue
            attrs_by_coord[(record_id, start, end)] = parse_gff3_attributes(fields[8])

    return attrs_by_coord


def gff_attrs_for_feature(
    attrs_by_coord: Dict[Tuple[str, int, int], Dict[str, str]],
    record_id: str,
    start0: int,
    end0: int,
) -> Dict[str, str]:
    """
    Try to match GenBank feature coords to GFF coords.

    GFF3 should be 1-based inclusive, while Biopython feature coordinates are
    0-based half-open. Some GRINS workflow versions accidentally passed the
    GFF start straight into FeatureLocation, so we try a few nearby forms.
    """
    candidate_keys = [
        (record_id, start0 + 1, end0),
        (record_id, start0, end0),
        (record_id, start0 + 1, end0 + 1),
        (strip_version(record_id), start0 + 1, end0),
        (strip_version(record_id), start0, end0),
    ]
    for key in candidate_keys:
        if key in attrs_by_coord:
            return attrs_by_coord[key]
    return {}


def feature_display_name(feature) -> str:
    keys = [
        "label",
        "Name",
        "name",
        "ID",
        "domain",
        "aSDomain",
        "product",
        "gene",
        "locus_tag",
        "protein_id",
        "note",
        "description",
    ]
    vals = all_quals(feature, keys)
    return clean_join(vals[:3], "/") or feature.type


def domain_name(feature) -> str:
    keys = [
        "aSDomain",
        "domain",
        "label",
        "Name",
        "name",
        "ID",
        "product",
        "note",
        "description",
    ]
    vals = all_quals(feature, keys)
    if vals:
        return clean_join(vals[:3], "/")
    return feature.type


def cds_locus_tag(feature) -> str:
    return first_qual(feature, ["locus_tag", "gene", "protein_id", "ID"], default="")


def cds_product(feature) -> str:
    return first_qual(feature, ["product", "gene_functions", "gene_function", "function", "note"], default="")


def bgc_product(feature) -> str:
    vals = all_quals(feature, ["product", "products", "kind", "category"])
    return clean_join(vals, ",") or "unknown"


def bgc_id(feature) -> str:
    return first_qual(feature, ["region_number", "ID", "label", "Name"], default="region")


@dataclass
class AnnotFeature:
    feature_type: str
    start0: int
    end0: int
    strand: str
    name: str
    raw_feature: object

    @property
    def start1(self) -> int:
        return self.start0 + 1

    @property
    def end1(self) -> int:
        return self.end0


@dataclass
class GrinsCopy:
    run: str
    assembly: str
    record_id: str
    grins_id: str
    grins_feature_id: str
    start0: int
    end0: int
    strand: str
    dup_id: str
    dup_group_id: str
    bgc_region_id: str = ""
    bgc_product: str = ""
    bgc_start0: Optional[int] = None
    bgc_end0: Optional[int] = None
    overlapping_cds: List[AnnotFeature] = field(default_factory=list)
    overlapping_domains: List[AnnotFeature] = field(default_factory=list)
    nearby_domains: List[Tuple[AnnotFeature, int]] = field(default_factory=list)
    nearest_domain_distance: Optional[int] = None
    nearest_cds_distance: Optional[int] = None
    domain_context: str = ""
    sequence_preview: str = ""
    sequence: str = ""

    @property
    def start1(self) -> int:
        return self.start0 + 1

    @property
    def end1(self) -> int:
        return self.end0

    @property
    def length(self) -> int:
        return self.end0 - self.start0


def locate_run_outputs(input_path: Path) -> Tuple[str, Path]:
    """
    Return (run_name, outputs_dir).

    Accepts either a run folder containing outputs/ or the outputs/ folder itself.
    """
    input_path = input_path.resolve()

    if (input_path / "outputs").is_dir():
        return input_path.name, input_path / "outputs"

    if input_path.name == "outputs" and input_path.is_dir():
        return input_path.parent.name, input_path

    assembly_like = [
        p for p in input_path.iterdir()
        if p.is_dir() and (p / "genomes_GRINS").is_dir()
    ]
    if assembly_like:
        return input_path.name, input_path

    raise FileNotFoundError(
        f"Could not find workflow outputs in {input_path}. Expected RUN_FOLDER/outputs or an outputs folder."
    )


def find_assembly_dirs(outputs_dir: Path, selected_assemblies: Optional[Sequence[str]] = None) -> List[Path]:
    selected = set(selected_assemblies or [])
    dirs = []
    for path in sorted(outputs_dir.iterdir()):
        if not path.is_dir():
            continue
        if selected and path.name not in selected:
            continue
        if (path / "genomes_GRINS").is_dir():
            dirs.append(path)
    return dirs


def find_first_existing(paths: Sequence[Path]) -> Optional[Path]:
    for path in paths:
        if path.exists():
            return path
    return None


def find_gbk(assembly_dir: Path) -> Path:
    assembly = assembly_dir.name
    candidates = [
        assembly_dir / "genomes_GRINS" / f"{assembly}.GRINS.gbk",
        assembly_dir / "genomes_GRINS" / f"{assembly}.gbk",
    ]

    found = find_first_existing(candidates)
    if found:
        return found

    matches = sorted((assembly_dir / "genomes_GRINS").glob("*.gbk"))
    if not matches:
        raise FileNotFoundError(f"No GenBank file found in {assembly_dir / 'genomes_GRINS'}")
    return matches[0]


def find_grins_gff(assembly_dir: Path) -> Optional[Path]:
    assembly = assembly_dir.name
    candidates = [
        assembly_dir / "GRINS.gff3" / f"{assembly}.GRINS.gff3",
        assembly_dir / "GRINS.gff3" / f"{assembly}.gff3",
    ]

    found = find_first_existing(candidates)
    if found:
        return found

    folder = assembly_dir / "GRINS.gff3"
    if folder.exists():
        matches = sorted(p for p in folder.glob("*.gff3") if not p.name.endswith("_CDS.gff3"))
        return matches[0] if matches else None

    return None


def collect_features(record) -> Tuple[List[AnnotFeature], List[AnnotFeature], List[AnnotFeature]]:
    bgcs: List[AnnotFeature] = []
    cdss: List[AnnotFeature] = []
    domains: List[AnnotFeature] = []

    for feature in record.features:
        start0, end0 = feature_bounds(feature)
        if end0 <= start0:
            continue

        annot = AnnotFeature(
            feature_type=feature.type,
            start0=start0,
            end0=end0,
            strand=strand_string(feature),
            name=feature_display_name(feature),
            raw_feature=feature,
        )

        if feature.type in BGC_TYPES:
            annot.name = bgc_product(feature)
            bgcs.append(annot)
        elif feature.type in CDS_TYPES:
            name_parts = [cds_locus_tag(feature), cds_product(feature)]
            annot.name = clean_join(name_parts, "/") or "CDS"
            cdss.append(annot)
        elif feature.type in DOMAIN_FEATURE_TYPES or "domain" in feature.type.lower():
            annot.name = domain_name(feature)
            domains.append(annot)

    return bgcs, cdss, domains


def choose_best_bgc(start0: int, end0: int, bgcs: Sequence[AnnotFeature]) -> Optional[AnnotFeature]:
    best: Optional[AnnotFeature] = None
    best_overlap = 0
    for bgc in bgcs:
        ov = overlap_bp(start0, end0, bgc.start0, bgc.end0)
        if ov > best_overlap:
            best = bgc
            best_overlap = ov
    return best


def nearby_features(
    start0: int,
    end0: int,
    features: Sequence[AnnotFeature],
    nearby_bp: int,
) -> List[Tuple[AnnotFeature, int]]:
    nearby: List[Tuple[AnnotFeature, int]] = []
    for feature in features:
        dist = distance_between_intervals(start0, end0, feature.start0, feature.end0)
        if dist <= nearby_bp:
            nearby.append((feature, dist))

    nearby.sort(key=lambda item: (item[1], item[0].start0, item[0].end0))
    return nearby


def nearest_distance(start0: int, end0: int, features: Sequence[AnnotFeature]) -> Optional[int]:
    if not features:
        return None
    return min(distance_between_intervals(start0, end0, f.start0, f.end0) for f in features)


def classify_domain_context(
    start0: int,
    end0: int,
    overlapping_domains: Sequence[AnnotFeature],
    domains: Sequence[AnnotFeature],
    nearby_domains_with_dist: Sequence[Tuple[AnnotFeature, int]],
    boundary_bp: int,
) -> str:
    if overlapping_domains:
        if len(overlapping_domains) > 1:
            return "spans_multiple_domains"

        domain = overlapping_domains[0]
        fully_inside = start0 >= domain.start0 and end0 <= domain.end0
        near_edge = min(
            abs(start0 - domain.start0),
            abs(start0 - domain.end0),
            abs(end0 - domain.start0),
            abs(end0 - domain.end0),
        ) <= boundary_bp

        if fully_inside and not near_edge:
            return "inside_single_domain"
        if near_edge:
            return "overlaps_domain_edge"
        return "overlaps_single_domain"

    left = [d for d in domains if d.end0 <= start0]
    right = [d for d in domains if d.start0 >= end0]

    if left and right:
        left_dist = start0 - max(d.end0 for d in left)
        right_dist = min(d.start0 for d in right) - end0
        if left_dist <= boundary_bp or right_dist <= boundary_bp:
            return "near_domain_boundary_or_linker"
        if nearby_domains_with_dist:
            return "interdomain_near_domains"

    if nearby_domains_with_dist:
        return "near_domain"

    return "no_nearby_domain"


def grins_id_from_feature(feature, gff_attrs: Dict[str, str], fallback_number: int) -> Tuple[str, str, str, str]:
    qualifiers = feature.qualifiers

    def q(key: str) -> str:
        vals = qualifiers.get(key)
        if vals:
            return clean(vals[0])
        return ""

    grins_id = (
        q("grins_id")
        or q("Name")
        or q("label")
        or gff_attrs.get("grins_id", "")
        or gff_attrs.get("Name", "")
    )

    grins_feature_id = (
        q("grins_feature_id")
        or q("ID")
        or gff_attrs.get("ID", "")
        or gff_attrs.get("grins_feature_id", "")
    )

    dup_id = q("dup_id") or gff_attrs.get("dup_id", "")
    dup_group_id = q("dup_group_id") or gff_attrs.get("dup_group_id", "")

    if not grins_id:
        grins_id = f"unlabeled_GRINS_{fallback_number}"

    if not grins_feature_id:
        grins_feature_id = f"{grins_id}_part_unknown_{fallback_number}"

    return grins_id, grins_feature_id, dup_id, dup_group_id


def extract_grins_copies_from_gbk(
    gbk_path: Path,
    run_name: str,
    assembly: str,
    nearby_bp: int,
    boundary_bp: int,
    gff_attrs_by_coord: Dict[Tuple[str, int, int], Dict[str, str]],
) -> List[GrinsCopy]:
    copies: List[GrinsCopy] = []
    fallback_number = 0

    for record in SeqIO.parse(str(gbk_path), "genbank"):
        bgcs, cdss, domains = collect_features(record)
        record_seq = str(record.seq)

        for feature in record.features:
            if feature.type not in GRINS_TYPES:
                continue

            fallback_number += 1
            start0, end0 = feature_bounds(feature)
            if end0 <= start0:
                continue

            gff_attrs = gff_attrs_for_feature(gff_attrs_by_coord, record.id, start0, end0)
            grins_id, grins_feature_id, dup_id, dup_group_id = grins_id_from_feature(
                feature, gff_attrs, fallback_number
            )

            overlapping_cdss = [
                cds for cds in cdss
                if intervals_overlap(start0, end0, cds.start0, cds.end0)
            ]

            overlapping_domains = [
                domain for domain in domains
                if intervals_overlap(start0, end0, domain.start0, domain.end0)
            ]

            nearby_domains_with_dist = nearby_features(start0, end0, domains, nearby_bp=nearby_bp)
            best_bgc = choose_best_bgc(start0, end0, bgcs)

            grins_seq = record_seq[start0:end0]

            copy = GrinsCopy(
                run=run_name,
                assembly=assembly,
                record_id=record.id,
                grins_id=grins_id,
                grins_feature_id=grins_feature_id,
                start0=start0,
                end0=end0,
                strand=strand_string(feature),
                dup_id=dup_id,
                dup_group_id=dup_group_id,
                overlapping_cds=overlapping_cdss,
                overlapping_domains=overlapping_domains,
                nearby_domains=nearby_domains_with_dist,
                nearest_domain_distance=nearest_distance(start0, end0, domains),
                nearest_cds_distance=nearest_distance(start0, end0, cdss),
                sequence_preview=grins_seq[:30],
                sequence=grins_seq,
            )

            if best_bgc is not None:
                copy.bgc_region_id = bgc_id(best_bgc.raw_feature)
                copy.bgc_product = bgc_product(best_bgc.raw_feature)
                copy.bgc_start0 = best_bgc.start0
                copy.bgc_end0 = best_bgc.end0

            copy.domain_context = classify_domain_context(
                start0,
                end0,
                overlapping_domains,
                domains,
                nearby_domains_with_dist,
                boundary_bp=boundary_bp,
            )

            copies.append(copy)

    return copies


def copy_to_row(copy: GrinsCopy, copy_number: int) -> Dict[str, object]:
    if copy.bgc_start0 is not None and copy.bgc_end0 is not None:
        bgc_start1 = copy.bgc_start0 + 1
        bgc_end1 = copy.bgc_end0
        dist_left = max(0, copy.start0 - copy.bgc_start0)
        dist_right = max(0, copy.bgc_end0 - copy.end0)
        nearest_edge = min(dist_left, dist_right)
    else:
        bgc_start1 = ""
        bgc_end1 = ""
        dist_left = ""
        dist_right = ""
        nearest_edge = ""

    nearby_domain_text = []
    for domain, dist in copy.nearby_domains[:10]:
        nearby_domain_text.append(f"{domain.name}@{domain.start1}-{domain.end1}:dist={dist}")

    return {
        "run": copy.run,
        "assembly": copy.assembly,
        "record_id": copy.record_id,
        "grins_id": copy.grins_id,
        "grins_feature_id": copy.grins_feature_id,
        "copy_number_within_grins_id": copy_number,
        "start_1based": copy.start1,
        "end_1based_inclusive": copy.end1,
        "length_bp": copy.length,
        "strand": copy.strand,
        "dup_id": copy.dup_id,
        "dup_group_id": copy.dup_group_id,
        "bgc_region_id": copy.bgc_region_id,
        "bgc_product": copy.bgc_product,
        "bgc_start_1based": bgc_start1,
        "bgc_end_1based_inclusive": bgc_end1,
        "distance_to_left_bgc_edge_bp": dist_left,
        "distance_to_right_bgc_edge_bp": dist_right,
        "nearest_bgc_edge_distance_bp": nearest_edge,
        "overlapping_cds_count": len(copy.overlapping_cds),
        "overlapping_cds_locus_tags": clean_join(cds_locus_tag(c.raw_feature) for c in copy.overlapping_cds),
        "overlapping_cds_products": clean_join(cds_product(c.raw_feature) for c in copy.overlapping_cds),
        "overlapping_domain_count": len(copy.overlapping_domains),
        "overlapping_domains": clean_join(
            f"{d.name}@{d.start1}-{d.end1}" for d in copy.overlapping_domains
        ),
        "nearby_domains": clean_join(nearby_domain_text),
        "domain_context": copy.domain_context,
        "nearest_domain_distance_bp": "" if copy.nearest_domain_distance is None else copy.nearest_domain_distance,
        "nearest_cds_distance_bp": "" if copy.nearest_cds_distance is None else copy.nearest_cds_distance,
        "sequence_preview_30bp": copy.sequence_preview,
    }


def group_arrangement(copies: Sequence[GrinsCopy], local_distance_bp: int) -> Tuple[str, str, str]:
    """
    Return arrangement, min_intercopy_distance, max_span.
    Distances/spans are calculated only within the same record/contig.
    """
    if len(copies) <= 1:
        return "single_copy", "", ""

    by_record: Dict[str, List[GrinsCopy]] = defaultdict(list)
    for copy in copies:
        by_record[copy.record_id].append(copy)

    if len(by_record) > 1:
        arrangement = "multi_record_or_multi_contig"
    else:
        arrangement = ""

    distances: List[int] = []
    spans: List[int] = []

    for record_copies in by_record.values():
        ordered = sorted(record_copies, key=lambda c: (c.start0, c.end0))
        if len(ordered) > 1:
            spans.append(max(c.end0 for c in ordered) - min(c.start0 for c in ordered))
            for left, right in zip(ordered, ordered[1:]):
                distances.append(distance_between_intervals(left.start0, left.end0, right.start0, right.end0))

    min_dist = min(distances) if distances else ""
    max_span = max(spans) if spans else ""

    if not arrangement:
        if distances and max(distances) <= local_distance_bp:
            arrangement = "local_or_tandem_repeats"
        else:
            arrangement = "scattered_within_locus"

    return arrangement, str(min_dist), str(max_span)


def priority_flags(copies: Sequence[GrinsCopy]) -> str:
    flags: List[str] = []
    n = len(copies)
    contexts = {c.domain_context for c in copies}

    if n == 1:
        flags.append("single_copy_low_priority")
    elif 2 <= n <= 5:
        flags.append("moderate_copy_number")
    elif n > 5:
        flags.append("high_copy_repeat_candidate")

    if any(c.bgc_product for c in copies):
        flags.append("BGC_associated")

    if sum(1 for c in copies if c.overlapping_cds) >= 2:
        flags.append("repeated_CDS_overlap")

    if sum(1 for c in copies if c.overlapping_domains) >= 2:
        flags.append("repeated_domain_overlap")

    if "spans_multiple_domains" in contexts:
        flags.append("spans_multiple_domains")

    if "overlaps_domain_edge" in contexts or "near_domain_boundary_or_linker" in contexts:
        flags.append("domain_boundary_or_linker_signal")

    mobile_keywords = ("transpos", "integrase", "recombinase")
    for copy in copies:
        products = clean_join(cds_product(c.raw_feature) for c in copy.overlapping_cds).lower()
        if any(keyword in products for keyword in mobile_keywords):
            flags.append("mobile_element_nearby_or_overlap")
            break

    return clean_join(flags)


def summarize_group(copies: Sequence[GrinsCopy], local_distance_bp: int) -> Dict[str, object]:
    first = copies[0]
    lengths = [c.length for c in copies]
    arrangement, min_dist, max_span = group_arrangement(copies, local_distance_bp=local_distance_bp)

    copy_coords = [
        f"{c.record_id}:{c.start1}-{c.end1}"
        for c in sorted(copies, key=lambda c: (c.record_id, c.start0, c.end0))
    ]

    return {
        "run": first.run,
        "assembly": first.assembly,
        "grins_id": first.grins_id,
        "n_copies": len(copies),
        "record_ids": clean_join(c.record_id for c in copies),
        "copy_coordinates": clean_join(copy_coords),
        "copy_lengths_bp": clean_join(lengths),
        "min_copy_length_bp": min(lengths),
        "max_copy_length_bp": max(lengths),
        "mean_copy_length_bp": round(sum(lengths) / len(lengths), 2),
        "dup_group_ids": clean_join(c.dup_group_id for c in copies),
        "n_bgc_copies": sum(1 for c in copies if c.bgc_product),
        "bgc_products": clean_join(c.bgc_product for c in copies),
        "n_cds_overlapping_copies": sum(1 for c in copies if c.overlapping_cds),
        "cds_locus_tags": clean_join(
            cds_locus_tag(cds.raw_feature)
            for copy in copies
            for cds in copy.overlapping_cds
        ),
        "n_domain_overlapping_copies": sum(1 for c in copies if c.overlapping_domains),
        "domain_names": clean_join(
            domain.name
            for copy in copies
            for domain in copy.overlapping_domains
        ),
        "domain_contexts": clean_join(c.domain_context for c in copies),
        "min_intercopy_distance_bp_same_record": min_dist,
        "max_span_bp_same_record": max_span,
        "arrangement": arrangement,
        "priority_flags": priority_flags(copies),
    }


def write_tsv(path: Path, rows: Sequence[Dict[str, object]], columns: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(columns), delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_fasta(path: Path, copies: Sequence[GrinsCopy]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for copy in copies:
            header = (
                f">{copy.run}|{copy.assembly}|{copy.grins_id}|{copy.grins_feature_id}|"
                f"{copy.record_id}:{copy.start1}-{copy.end1}|dup_group={copy.dup_group_id or 'NA'}"
            )
            handle.write(header + "\n")
            seq = copy.sequence.upper()
            for i in range(0, len(seq), 80):
                handle.write(seq[i:i + 80] + "\n")


@dataclass
class BlockOccurrence:
    """One observed copy of a recurring adjacent/nearby GRINS-ID block."""

    pattern: Tuple[str, ...]
    members: List[GrinsCopy]

    @property
    def run(self) -> str:
        return self.members[0].run

    @property
    def assembly(self) -> str:
        return self.members[0].assembly

    @property
    def record_id(self) -> str:
        return self.members[0].record_id

    @property
    def start1(self) -> int:
        return min(m.start1 for m in self.members)

    @property
    def end1(self) -> int:
        return max(m.end1 for m in self.members)

    @property
    def length(self) -> int:
        return self.end1 - self.start1 + 1

    @property
    def coord(self) -> str:
        return f"{self.record_id}:{self.start1}-{self.end1}"

    @property
    def member_coords(self) -> str:
        return ";".join(f"{m.grins_id}@{m.record_id}:{m.start1}-{m.end1}" for m in self.members)

    @property
    def gaps(self) -> List[int]:
        gaps: List[int] = []
        ordered = sorted(self.members, key=lambda m: (m.start0, m.end0, m.grins_id))
        for left, right in zip(ordered, ordered[1:]):
            gaps.append(max(0, right.start1 - left.end1 - 1))
        return gaps

    @property
    def gaps_text(self) -> str:
        return ",".join(str(g) for g in self.gaps)


def split_multi(value: object) -> List[str]:
    """Split pipe/comma/semicolon-delimited annotation strings into unique tokens."""
    text = clean(value)
    if not text:
        return []
    return [clean(part) for part in re.split(r"[|,;]+", text) if clean(part)]


def uniq_join(values: Iterable[object], sep: str = "|") -> str:
    out: List[str] = []
    seen = set()
    for value in values:
        for part in split_multi(value):
            if part and part not in seen:
                out.append(part)
                seen.add(part)
    return sep.join(out)


def gap_between_copies(left: GrinsCopy, right: GrinsCopy) -> int:
    return max(0, right.start1 - left.end1 - 1)


def copy_domain_names(copy: GrinsCopy) -> str:
    overlapping = clean_join(f"{d.name}@{d.start1}-{d.end1}" for d in copy.overlapping_domains)
    if overlapping:
        return overlapping
    nearby = []
    for domain, dist in copy.nearby_domains[:10]:
        nearby.append(f"{domain.name}@{domain.start1}-{domain.end1}:dist={dist}")
    return clean_join(nearby)


def copy_cds_locus_tags(copy: GrinsCopy) -> str:
    return clean_join(cds_locus_tag(c.raw_feature) for c in copy.overlapping_cds)


def build_candidate_block_occurrences(
    copies: Sequence[GrinsCopy],
    max_gap: int,
    min_members: int,
    max_members: int,
    min_distinct_ids: int,
) -> List[BlockOccurrence]:
    """Generate contiguous nearby GRINS-ID patterns within each record."""
    by_record: Dict[str, List[GrinsCopy]] = defaultdict(list)
    for copy in copies:
        by_record[copy.record_id].append(copy)

    occurrences: List[BlockOccurrence] = []

    for _record_id, record_copies in by_record.items():
        ordered = sorted(record_copies, key=lambda c: (c.start0, c.end0, c.grins_id))
        n = len(ordered)

        for i in range(n):
            members = [ordered[i]]
            for j in range(i + 1, min(n, i + max_members)):
                prev = ordered[j - 1]
                current = ordered[j]

                if gap_between_copies(prev, current) > max_gap:
                    break

                members.append(current)

                if len(members) < min_members:
                    continue

                pattern = tuple(m.grins_id for m in members)
                if len(set(pattern)) < min_distinct_ids:
                    continue

                occurrences.append(BlockOccurrence(pattern=pattern, members=list(members)))

    return occurrences


def remove_exact_duplicate_block_occurrences(occurrences: Sequence[BlockOccurrence]) -> List[BlockOccurrence]:
    seen = set()
    out: List[BlockOccurrence] = []
    for occ in occurrences:
        key = (
            occ.pattern,
            occ.record_id,
            tuple((m.grins_id, m.start1, m.end1) for m in occ.members),
        )
        if key in seen:
            continue
        seen.add(key)
        out.append(occ)
    return out


def group_block_occurrences(
    occurrences: Sequence[BlockOccurrence],
    min_instances: int,
) -> Dict[Tuple[str, ...], List[BlockOccurrence]]:
    grouped: Dict[Tuple[str, ...], List[BlockOccurrence]] = defaultdict(list)
    for occ in occurrences:
        grouped[occ.pattern].append(occ)
    return {pattern: occs for pattern, occs in grouped.items() if len(occs) >= min_instances}


def is_subsequence(short: Tuple[str, ...], long: Tuple[str, ...]) -> bool:
    if len(short) >= len(long):
        return False
    for i in range(0, len(long) - len(short) + 1):
        if long[i:i + len(short)] == short:
            return True
    return False


def apply_maximal_block_filter(
    grouped: Dict[Tuple[str, ...], List[BlockOccurrence]]
) -> Dict[Tuple[str, ...], List[BlockOccurrence]]:
    """Remove shorter patterns if they are contained in longer repeated patterns."""
    keep: Dict[Tuple[str, ...], List[BlockOccurrence]] = {}
    patterns = list(grouped)
    for pattern in patterns:
        if any(is_subsequence(pattern, other) for other in patterns):
            continue
        keep[pattern] = grouped[pattern]
    return keep


def infer_block_arrangement(occs: Sequence[BlockOccurrence], local_distance_bp: int) -> str:
    if len(occs) <= 1:
        return "single_block_copy"

    by_record: Dict[str, List[BlockOccurrence]] = defaultdict(list)
    for occ in occs:
        by_record[occ.record_id].append(occ)

    if len(by_record) > 1:
        return "repeated_across_records_or_contigs"

    ordered = sorted(next(iter(by_record.values())), key=lambda o: (o.start1, o.end1))
    adjacent_gaps = [max(0, right.start1 - left.end1 - 1) for left, right in zip(ordered, ordered[1:])]
    if adjacent_gaps and max(adjacent_gaps) <= local_distance_bp:
        return "local_or_tandem_repeated_block"
    return "scattered_repeated_block"


def block_interpretation_hint(occs: Sequence[BlockOccurrence]) -> str:
    all_domains = uniq_join(copy_domain_names(m) for occ in occs for m in occ.members).lower()
    all_contexts = uniq_join(m.domain_context for occ in occs for m in occ.members).lower()
    pattern = occs[0].pattern if occs else ()

    hints: List[str] = []
    if len(pattern) >= 2:
        hints.append("recurring_adjacent_grins_ids")
    if "ks" in all_domains:
        hints.append("KS_associated")
    if "at" in all_domains:
        hints.append("AT_associated")
    if "pks" in all_domains:
        hints.append("PKS_domain_associated")
    if "nrps" in all_domains:
        hints.append("NRPS_domain_associated")
    if "boundary" in all_contexts or "linker" in all_contexts or "edge" in all_contexts:
        hints.append("possible_domain_boundary_or_linker_signal")
    if len(occs) >= 2:
        hints.append("candidate_duplicated_block")

    return clean_join(hints)


def summarize_compound_block(
    block_id: str,
    pattern: Tuple[str, ...],
    occs: Sequence[BlockOccurrence],
    local_distance_bp: int,
) -> Dict[str, object]:
    first = occs[0]
    lengths = [occ.length for occ in occs]
    all_gaps = [gap for occ in occs for gap in occ.gaps]

    return {
        "run": first.run,
        "assembly": first.assembly,
        "block_id": block_id,
        "member_grins_ids": "+".join(pattern),
        "n_member_regions": len(pattern),
        "n_block_copies": len(occs),
        "records": clean_join(occ.record_id for occ in occs),
        "block_copy_coordinates": clean_join(occ.coord for occ in occs),
        "member_copy_coordinates": clean_join(occ.member_coords for occ in occs),
        "mean_block_length_bp": round(sum(lengths) / len(lengths), 2) if lengths else "",
        "min_block_length_bp": min(lengths) if lengths else "",
        "max_block_length_bp": max(lengths) if lengths else "",
        "mean_internal_gap_bp": round(sum(all_gaps) / len(all_gaps), 2) if all_gaps else "",
        "internal_gaps_by_copy": clean_join(occ.gaps_text for occ in occs),
        "bgc_products": clean_join(m.bgc_product for occ in occs for m in occ.members),
        "domain_contexts": clean_join(m.domain_context for occ in occs for m in occ.members),
        "domain_names": clean_join(copy_domain_names(m) for occ in occs for m in occ.members),
        "cds_locus_tags": clean_join(copy_cds_locus_tags(m) for occ in occs for m in occ.members),
        "arrangement": infer_block_arrangement(occs, local_distance_bp=local_distance_bp),
        "interpretation_hint": block_interpretation_hint(occs),
    }


def compound_block_copy_row(
    block_id: str,
    copy_number: int,
    occ: BlockOccurrence,
) -> Dict[str, object]:
    return {
        "run": occ.run,
        "assembly": occ.assembly,
        "block_id": block_id,
        "block_copy_number": copy_number,
        "record_id": occ.record_id,
        "block_start_1based": occ.start1,
        "block_end_1based_inclusive": occ.end1,
        "block_length_bp": occ.length,
        "member_grins_ids": "+".join(occ.pattern),
        "member_grins_feature_ids": clean_join(m.grins_feature_id for m in occ.members),
        "member_coordinates": occ.member_coords,
        "internal_gaps_bp": occ.gaps_text,
        "bgc_products": clean_join(m.bgc_product for m in occ.members),
        "domain_contexts": clean_join(m.domain_context for m in occ.members),
        "domain_names": clean_join(copy_domain_names(m) for m in occ.members),
        "cds_locus_tags": clean_join(copy_cds_locus_tags(m) for m in occ.members),
    }


def summarize_compound_blocks_from_copies(
    copies: Sequence[GrinsCopy],
    max_gap: int,
    min_instances: int,
    min_members: int,
    max_members: int,
    min_distinct_ids: int,
    local_distance_bp: int,
    maximal_only: bool,
) -> Tuple[List[Dict[str, object]], List[Dict[str, object]]]:
    """Find recurring adjacent/nearby GRINS-ID patterns from extracted GRINS copies."""
    candidate_occs = build_candidate_block_occurrences(
        copies=copies,
        max_gap=max_gap,
        min_members=min_members,
        max_members=max_members,
        min_distinct_ids=min_distinct_ids,
    )
    candidate_occs = remove_exact_duplicate_block_occurrences(candidate_occs)
    grouped = group_block_occurrences(candidate_occs, min_instances=min_instances)

    if maximal_only:
        grouped = apply_maximal_block_filter(grouped)

    sorted_items = sorted(
        grouped.items(),
        key=lambda item: (
            -len(item[0]),
            -len(item[1]),
            item[1][0].record_id if item[1] else "",
            item[1][0].start1 if item[1] else 0,
            item[0],
        ),
    )

    block_rows: List[Dict[str, object]] = []
    block_copy_rows: List[Dict[str, object]] = []

    for idx, (pattern, occs) in enumerate(sorted_items, start=1):
        block_id = f"BLOCK_{idx}"
        occs_sorted = sorted(occs, key=lambda occ: (occ.record_id, occ.start1, occ.end1))
        block_rows.append(
            summarize_compound_block(
                block_id=block_id,
                pattern=pattern,
                occs=occs_sorted,
                local_distance_bp=local_distance_bp,
            )
        )
        for copy_number, occ in enumerate(occs_sorted, start=1):
            block_copy_rows.append(compound_block_copy_row(block_id, copy_number, occ))

    return block_rows, block_copy_rows


def process_run(
    input_path: Path,
    output_dir: Path,
    nearby_bp: int,
    boundary_bp: int,
    local_distance_bp: int,
    selected_assemblies: Optional[Sequence[str]],
    write_sequences: bool,
    compound_max_gap: int,
    compound_min_instances: int,
    compound_min_members: int,
    compound_max_members: int,
    compound_min_distinct_ids: int,
    compound_maximal_only: bool,
) -> Tuple[
    List[Dict[str, object]],
    List[Dict[str, object]],
    List[Dict[str, object]],
    List[Dict[str, object]],
]:
    """
    Process one workflow run folder.

    New behavior:
        - Input is the whole workflow output folder.
        - Each assembly/output folder gets its own files:
              ASSEMBLY_grins_copies.tsv
              ASSEMBLY_grins_groups.tsv
              ASSEMBLY_compound_grins_blocks.tsv
              ASSEMBLY_compound_grins_block_copies.tsv
              ASSEMBLY_grins_copies.fasta, if --write_fasta is used
    """
    run_name, outputs_dir = locate_run_outputs(input_path)
    assembly_dirs = find_assembly_dirs(outputs_dir, selected_assemblies=selected_assemblies)

    if not assembly_dirs:
        raise FileNotFoundError(f"No assembly output folders found in {outputs_dir}")

    all_copy_rows: List[Dict[str, object]] = []
    all_group_rows: List[Dict[str, object]] = []
    all_block_rows: List[Dict[str, object]] = []
    all_block_copy_rows: List[Dict[str, object]] = []

    for assembly_dir in assembly_dirs:
        assembly = assembly_dir.name
        gbk_path = find_gbk(assembly_dir)
        grins_gff = find_grins_gff(assembly_dir)
        gff_attrs = read_grins_gff_attrs(grins_gff)

        copies = extract_grins_copies_from_gbk(
            gbk_path=gbk_path,
            run_name=run_name,
            assembly=assembly,
            nearby_bp=nearby_bp,
            boundary_bp=boundary_bp,
            gff_attrs_by_coord=gff_attrs,
        )

        if not copies:
            print(f"WARNING: no GRINS features found in {gbk_path}", file=sys.stderr)

        # One row per individual GRINS copy.
        copy_counts: Dict[Tuple[str, str, str], int] = defaultdict(int)
        copy_rows: List[Dict[str, object]] = []

        for copy in sorted(copies, key=lambda c: (c.grins_id, c.record_id, c.start0, c.end0)):
            key = (copy.run, copy.assembly, copy.grins_id)
            copy_counts[key] += 1
            copy_rows.append(copy_to_row(copy, copy_number=copy_counts[key]))

        # One row per GRINS ID/group.
        groups: Dict[Tuple[str, str, str], List[GrinsCopy]] = defaultdict(list)
        for copy in copies:
            groups[(copy.run, copy.assembly, copy.grins_id)].append(copy)

        group_rows = [
            summarize_group(group_copies, local_distance_bp=local_distance_bp)
            for _key, group_copies in sorted(groups.items(), key=lambda item: item[0])
        ]

        # One row per recurring adjacent/nearby compound GRINS block, plus one row
        # per observed copy of each block. This is a higher-level interpretation
        # layer on top of direct GRINS sequence-family IDs.
        block_rows, block_copy_rows = summarize_compound_blocks_from_copies(
            copies=copies,
            max_gap=compound_max_gap,
            min_instances=compound_min_instances,
            min_members=compound_min_members,
            max_members=compound_max_members,
            min_distinct_ids=compound_min_distinct_ids,
            local_distance_bp=local_distance_bp,
            maximal_only=compound_maximal_only,
        )

        safe_assembly = re.sub(r"[^A-Za-z0-9_.-]+", "_", assembly)

        copies_out = output_dir / f"{safe_assembly}_grins_copies.tsv"
        groups_out = output_dir / f"{safe_assembly}_grins_groups.tsv"
        blocks_out = output_dir / f"{safe_assembly}_compound_grins_blocks.tsv"
        block_copies_out = output_dir / f"{safe_assembly}_compound_grins_block_copies.tsv"

        write_tsv(copies_out, copy_rows, COPY_COLUMNS)
        write_tsv(groups_out, group_rows, GROUP_COLUMNS)
        write_tsv(blocks_out, block_rows, BLOCK_COLUMNS)
        write_tsv(block_copies_out, block_copy_rows, BLOCK_COPY_COLUMNS)

        if write_sequences:
            fasta_out = output_dir / f"{safe_assembly}_grins_copies.fasta"
            write_fasta(fasta_out, copies)

        all_copy_rows.extend(copy_rows)
        all_group_rows.extend(group_rows)
        all_block_rows.extend(block_rows)
        all_block_copy_rows.extend(block_copy_rows)

        print(f"Processed assembly: {assembly}")
        print(f"  GRINS copies: {len(copy_rows)}")
        print(f"  GRINS IDs/groups: {len(group_rows)}")
        print(f"  compound GRINS blocks: {len(block_rows)}")
        print(f"  compound block copies: {len(block_copy_rows)}")
        print(f"  wrote: {copies_out}")
        print(f"  wrote: {groups_out}")
        print(f"  wrote: {blocks_out}")
        print(f"  wrote: {block_copies_out}")

        if write_sequences:
            print(f"  wrote: {output_dir / f'{safe_assembly}_grins_copies.fasta'}")

    print(f"\nFinished run: {run_name}")
    print(f"  assemblies processed: {len(assembly_dirs)}")
    print(f"  total GRINS copies: {len(all_copy_rows)}")
    print(f"  total GRINS IDs/groups: {len(all_group_rows)}")
    print(f"  total compound GRINS blocks: {len(all_block_rows)}")
    print(f"  total compound block copies: {len(all_block_copy_rows)}")

    return all_copy_rows, all_group_rows, all_block_rows, all_block_copy_rows


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Summarize GRINS IDs from GRINS workflow outputs.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "inputs",
        nargs="+",
        help="One or more run folders containing outputs/, or outputs/ folders directly.",
    )

    parser.add_argument(
        "--output_dir",
        default="grins_id_summaries",
        help="Folder where summary TSV files will be written.",
    )

    parser.add_argument(
        "--assemblies",
        nargs="*",
        default=None,
        help="Optional assembly folder names to include, e.g. AB088224.2 CP035319.1.",
    )

    parser.add_argument(
        "--nearby_bp",
        type=int,
        default=3000,
        help="Distance used to report nearby domains around each GRINS copy.",
    )

    parser.add_argument(
        "--boundary_bp",
        type=int,
        default=300,
        help="Distance used to call a GRINS near a domain boundary/linker.",
    )

    parser.add_argument(
        "--local_distance_bp",
        type=int,
        default=10000,
        help="Max adjacent-copy distance for calling a GRINS group local/tandem.",
    )

    parser.add_argument(
        "--write_fasta",
        action="store_true",
        help="Write FASTA sequences for all GRINS copies.",
    )

    parser.add_argument(
        "--compound_max_gap",
        type=int,
        default=1000,
        help="Maximum gap between adjacent GRINS copies to consider them part of one compound block.",
    )

    parser.add_argument(
        "--compound_min_instances",
        type=int,
        default=2,
        help="Minimum number of repeated block copies needed to report a compound block.",
    )

    parser.add_argument(
        "--compound_min_members",
        type=int,
        default=2,
        help="Minimum number of adjacent GRINS regions in a compound block.",
    )

    parser.add_argument(
        "--compound_max_members",
        type=int,
        default=5,
        help="Maximum number of adjacent GRINS regions to include in candidate block patterns.",
    )

    parser.add_argument(
        "--compound_min_distinct_ids",
        type=int,
        default=2,
        help="Minimum number of distinct GRINS IDs required in a compound block.",
    )

    parser.add_argument(
        "--compound_maximal_only",
        action="store_true",
        help="Suppress shorter repeated patterns if they are contained in longer repeated patterns.",
    )

    parser.add_argument(
        "--combined",
        action="store_true",
        help="Also write combined TSV files across all input runs.",
    )

    return parser


def main() -> None:
    args = build_parser().parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    combined_copy_rows: List[Dict[str, object]] = []
    combined_group_rows: List[Dict[str, object]] = []
    combined_block_rows: List[Dict[str, object]] = []
    combined_block_copy_rows: List[Dict[str, object]] = []

    for input_name in args.inputs:
        copy_rows, group_rows, block_rows, block_copy_rows = process_run(
            input_path=Path(input_name),
            output_dir=output_dir,
            nearby_bp=args.nearby_bp,
            boundary_bp=args.boundary_bp,
            local_distance_bp=args.local_distance_bp,
            selected_assemblies=args.assemblies,
            write_sequences=args.write_fasta,
            compound_max_gap=args.compound_max_gap,
            compound_min_instances=args.compound_min_instances,
            compound_min_members=args.compound_min_members,
            compound_max_members=args.compound_max_members,
            compound_min_distinct_ids=args.compound_min_distinct_ids,
            compound_maximal_only=args.compound_maximal_only,
        )
        combined_copy_rows.extend(copy_rows)
        combined_group_rows.extend(group_rows)
        combined_block_rows.extend(block_rows)
        combined_block_copy_rows.extend(block_copy_rows)

    if args.combined and len(args.inputs) > 1:
        write_tsv(output_dir / "combined_grins_copies.tsv", combined_copy_rows, COPY_COLUMNS)
        write_tsv(output_dir / "combined_grins_groups.tsv", combined_group_rows, GROUP_COLUMNS)
        write_tsv(output_dir / "combined_compound_grins_blocks.tsv", combined_block_rows, BLOCK_COLUMNS)
        write_tsv(output_dir / "combined_compound_grins_block_copies.tsv", combined_block_copy_rows, BLOCK_COPY_COLUMNS)
        print(f"Wrote combined summaries to {output_dir}")


if __name__ == "__main__":
    main()
