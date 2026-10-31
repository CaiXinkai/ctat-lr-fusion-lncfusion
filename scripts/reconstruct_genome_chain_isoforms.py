#!/usr/bin/env python3
"""Reconstruct fusion-isoform core sequences directly from hg38.

This is intentionally separate from ``build_fusion_isoform_sequences.py``.
It does not choose transcript ends or cut annotated transcripts at minimap2
alignment endpoints.  Instead it:

1. reads each global step-01 HC isoform entity exactly once (preferred input);
2. reads the representative block chain embedded in that global table;
3. tests adjacent block transitions and uniquely selects the transition whose
   two projected sides agree best with the named GTF partner genes;
4. projects every fusion-contig block back to its partner locus using that
   local transition and the canonical genomic breakpoints;
5. extracts those intervals from the reference genome in transcript order;
6. never opens per-sample folders or NR sequences in global mode.  The legacy
   sample-call mode remains available only for backward compatibility.

The result is an observed, reference-corrected *core* sequence.  Internal
junctions and the fusion join are represented, but unobserved 5'/3' transcript
ends are not invented.  Therefore these records must not automatically be
described as uniquely proven full-length transcripts.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path


# CTAT preliminary tables can contain extremely long read/annotation fields.
# Python's csv module otherwise rejects fields larger than 128 KiB.
_csv_limit = sys.maxsize
while True:
    try:
        csv.field_size_limit(_csv_limit)
        break
    except OverflowError:
        _csv_limit //= 10


def clean(value) -> str:
    return "" if value is None else str(value).strip()


def open_text(path: Path, mode: str = "rt"):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode, encoding="utf-8", newline="")
    return path.open(mode.replace("t", ""), encoding="utf-8", newline="")


def integer(value, default=0):
    try:
        return int(float(clean(value).replace(",", "")))
    except (TypeError, ValueError):
        return default


def normalize_chain(value: str) -> str:
    text = clean(value).upper().replace("CHR", "")
    text = re.sub(r"\s+", "", text)
    text = re.sub(r"(?<=\d),(?=\d)", "", text)
    text = text.replace("=>", ">").replace("->", ">")
    return re.sub(r"[;|]+", ";", text).strip(";")


def parse_breakpoint(value: str):
    chrom, pos, strand = clean(value).rsplit(":", 2)
    return chrom.removeprefix("chr"), int(pos.replace(",", "")), strand


def close_breakpoint(a: str, b: str, tolerance: int) -> bool:
    try:
        ca, pa, sa = parse_breakpoint(a)
        cb, pb, sb = parse_breakpoint(b)
        return ca == cb and sa == sb and abs(pa - pb) <= tolerance
    except Exception:
        return False


def find_one(directory: Path, patterns: list[str]):
    for pattern in patterns:
        matches = sorted(directory.glob(pattern))
        if matches:
            return matches[0]
    return None


def read_fasta(path: Path) -> dict[str, str]:
    result: dict[str, list[str]] = {}
    name = ""
    with open_text(path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:].split()[0].split("|", 1)[0]
                result.setdefault(name, [])
            elif name:
                result[name].append(line.upper().replace("U", "T"))
    return {name: "".join(parts) for name, parts in result.items()}


def read_rows(path: Path) -> list[dict]:
    with open_text(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def parse_gtf_attributes(text: str) -> dict[str, str]:
    result = {}
    for key, value in re.findall(r'(\S+)\s+"([^"]+)"', text):
        result[key] = value
    return result


def load_target_transcripts(gtf: Path, target_genes: set[str]):
    """Load complete transcript exon chains for fusion partners in this run."""
    aliases = {clean(g) for g in target_genes if clean(g)}
    aliases |= {g.split(".", 1)[0] for g in list(aliases)}
    transcripts = {}
    gene_to_transcripts = defaultdict(set)
    with open_text(gtf) as handle:
        for raw in handle:
            if not raw or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "exon":
                continue
            attrs = parse_gtf_attributes(fields[8])
            keys = {clean(attrs.get("gene_name")), clean(attrs.get("gene_id"))}
            keys |= {key.split(".", 1)[0] for key in list(keys) if key}
            matched = keys & aliases
            if not matched:
                continue
            transcript_id = clean(attrs.get("transcript_id"))
            if not transcript_id:
                continue
            tags = set(re.findall(r'tag\s+"([^"]+)"', fields[8]))
            record = (fields[0].removeprefix("chr"), int(fields[3]), int(fields[4]), fields[6])
            model = transcripts.setdefault(transcript_id, {
                "transcript_id": transcript_id, "chrom": record[0], "strand": record[3],
                "exons": set(), "tags": set(), "genes": set(),
            })
            model["exons"].add((record[1], record[2]))
            model["tags"].update(tags)
            for key in matched:
                model["genes"].add(key)
                gene_to_transcripts[key].add(transcript_id)
    for model in transcripts.values():
        ordered = sorted(model["exons"], reverse=model["strand"] == "-")
        model["exons"] = ordered
    return {gene: [transcripts[tid] for tid in sorted(ids)]
            for gene, ids in gene_to_transcripts.items()}


def sample_ids(path: Path | None) -> set[str]:
    if path is None:
        return set()
    with open_text(path) as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return {clean(row.get("sample_id") or row.get("Sample") or row.get("sample")) for row in rows} - {""}


def load_hc(path: Path, selected: set[str]):
    result: dict[str, list[dict]] = defaultdict(list)
    with open_text(path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            sample = clean(row.get("sample_id"))
            if selected and sample not in selected:
                continue
            if clean(row.get("feature_level")) not in {"", "isoform"}:
                continue
            if clean(row.get("main_analysis_eligible_after_step01")).upper() == "FALSE":
                continue
            parts = clean(row.get("breakpoint_id") or row.get("event_id")).split("|")
            if len(parts) >= 5:
                row["left_breakpoint"] = clean(row.get("left_breakpoint")) or parts[3]
                row["right_breakpoint"] = clean(row.get("right_breakpoint")) or parts[4]
            result[sample].append(row)
    return result


def load_global_isoforms(path: Path, selected: set[str]):
    """Load one-row-per-feature entities without expanding them back to samples."""
    entities = {}
    with open_text(path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            feature_id = clean(row.get("feature_id"))
            if not feature_id:
                raise ValueError("Global isoform table contains a row without feature_id")
            if feature_id in entities:
                raise ValueError(f"Duplicate feature_id in global isoform table: {feature_id}")
            entities[feature_id] = dict(row)
    return entities


def match_nr(hc: dict, rows: list[dict], tolerance: int):
    core = ("fusion_name", "left_breakpoint", "right_breakpoint", "splice_type")
    exact = [row for row in rows if all(clean(row.get(key)) == clean(hc.get(key)) for key in core)
             and clean(row.get("internal_junction_chain")) == clean(hc.get("internal_junction_chain"))]
    if len(exact) == 1:
        return exact[0], "EXACT"
    candidates = [row for row in rows if all(clean(row.get(key)) == clean(hc.get(key)) for key in core)]
    chain = normalize_chain(hc.get("internal_junction_chain", ""))
    normalized = [row for row in candidates if chain and normalize_chain(row.get("internal_junction_chain", "")) == chain]
    if len(normalized) == 1:
        return normalized[0], "NORMALIZED_CHAIN"
    fuzzy = [row for row in rows
             if clean(row.get("fusion_name")) == clean(hc.get("fusion_name"))
             and clean(row.get("splice_type")) == clean(hc.get("splice_type"))
             and close_breakpoint(row.get("left_breakpoint", ""), hc.get("left_breakpoint", ""), tolerance)
             and close_breakpoint(row.get("right_breakpoint", ""), hc.get("right_breakpoint", ""), tolerance)]
    normalized = [row for row in fuzzy if chain and normalize_chain(row.get("internal_junction_chain", "")) == chain]
    if len(normalized) == 1:
        return normalized[0], "FUZZY_BREAKPOINT_NORMALIZED_CHAIN"
    if len(candidates) == 1:
        return candidates[0], "UNIQUE_EVENT"
    if len(fuzzy) == 1:
        return fuzzy[0], "FUZZY_BREAKPOINT"
    return None, "NO_UNIQUE_NR_MATCH" if candidates or fuzzy else "NO_NR_MATCH"


def parse_blocks(value: str):
    blocks = []
    for token in clean(value).split(";"):
        token = token.strip()
        match = re.search(r"(?:^|:)(\d+)-(\d+):([+-])$", token)
        if not match:
            continue
        start, end, strand = int(match.group(1)), int(match.group(2)), match.group(3)
        blocks.append((min(start, end), max(start, end), strand))
    return blocks


def preliminary_index(rows: list[dict]):
    index: dict[tuple[str, str, str], list[dict]] = defaultdict(list)
    for row in rows:
        fusion = clean(row.get("#FusionName") or row.get("FusionName") or row.get("fusion_name"))
        left = clean(row.get("LeftBreakpoint") or row.get("left_breakpoint"))
        right = clean(row.get("RightBreakpoint") or row.get("right_breakpoint"))
        index[(fusion, left, right)].append(row)
    return index


def local_breakpoint_signature(row: dict):
    """Return the only preliminary fields needed for coordinate projection."""
    left = integer(row.get("LeftLocalBreakpoint") or row.get("left_local_breakpoint"), -1)
    right = integer(row.get("RightLocalBreakpoint") or row.get("right_local_breakpoint"), -1)
    return left, right


def resolve_equivalent_preliminary(rows: list[dict], method: str):
    """Collapse duplicate preliminary rows that encode one local mapping.

    CTAT may report the same event on several preliminary rows because other
    evidence columns differ.  Those rows are not ambiguous for sequence
    extraction if their two local breakpoint anchors are identical.
    """
    usable = [row for row in rows if min(local_breakpoint_signature(row)) >= 0]
    signatures = {local_breakpoint_signature(row) for row in usable}
    if len(signatures) == 1:
        return usable[0], method if len(usable) == 1 else method + "_DEDUP_LOCAL_BREAKPOINTS"
    return None, "NO_UNIQUE_LOCAL_BREAKPOINT_MAPPING"


def preliminary_candidates(info: dict, index: dict, tolerance: int):
    key = (clean(info.get("fusion_name")), clean(info.get("left_breakpoint")), clean(info.get("right_breakpoint")))
    rows = index.get(key, [])
    if rows:
        return rows, "EXACT"
    matches = []
    for (fusion, left, right), candidates in index.items():
        if fusion == key[0] and close_breakpoint(left, key[1], tolerance) and close_breakpoint(right, key[2], tolerance):
            matches.extend(candidates)
    return (matches, "FUZZY_BREAKPOINT") if matches else ([], "NO_PRELIMINARY_EVENT")


def match_preliminary(info: dict, index: dict, tolerance: int, reference=None,
                      transcripts=None, terminal_tolerance: int = 5):
    """Resolve preliminary rows without using NR nucleotide sequence.

    If several local breakpoint signatures remain, every signature is projected
    independently.  Rescue is allowed only when exactly one projection is
    structurally valid, or when all valid candidates yield the identical
    sequence and projected block chain.
    """
    rows, source = preliminary_candidates(info, index, tolerance)
    if not rows:
        return None, source
    event, method = resolve_equivalent_preliminary(rows, source)
    if event is not None:
        return event, method
    if reference is None:
        return None, "NO_UNIQUE_LOCAL_BREAKPOINT_MAPPING"
    unique = {}
    for row in rows:
        unique.setdefault(local_breakpoint_signature(row), row)
    valid = []
    for signature, row in unique.items():
        try:
            sequence, reason, metrics = reconstruct(
                info, row, reference, transcripts, terminal_tolerance)
        except Exception:
            continue
        if sequence is not None:
            valid.append((row, sequence, metrics.get("projected_block_chain", ""), signature))
    if len(valid) == 1:
        return valid[0][0], source + "_UNIQUE_VALID_PROJECTION"
    products = {(sequence, chain) for _, sequence, chain, _ in valid}
    if valid and len(products) == 1:
        return valid[0][0], source + "_IDENTICAL_PROJECTION_RESCUE"
    return None, "NO_UNIQUE_LOCAL_BREAKPOINT_MAPPING"


class IndexedFasta:
    def __init__(self, fasta: Path):
        self.path = fasta
        fai = Path(str(fasta) + ".fai")
        if not fai.exists():
            raise FileNotFoundError(f"Missing FASTA index: {fai}; run: samtools faidx {fasta}")
        self.index = {}
        with fai.open(encoding="utf-8") as handle:
            for line in handle:
                chrom, length, offset, bases, width = line.rstrip("\n").split("\t")[:5]
                record = (int(length), int(offset), int(bases), int(width))
                self.index[chrom] = record
                self.index.setdefault(chrom.removeprefix("chr"), record)
        self.handle = fasta.open("rb")

    def fetch(self, chrom: str, start: int, end: int) -> str:
        chrom_key = chrom if chrom in self.index else chrom.removeprefix("chr")
        if chrom_key not in self.index:
            raise KeyError(f"Chromosome absent from FASTA index: {chrom}")
        length, offset, line_bases, line_width = self.index[chrom_key]
        if start < 1 or end > length or start > end:
            raise ValueError(f"Reference interval out of range: {chrom}:{start}-{end}")
        zero = start - 1
        byte_start = offset + (zero // line_bases) * line_width + zero % line_bases
        needed = end - start + 1
        extra = (needed // line_bases + 2) * (line_width - line_bases)
        self.handle.seek(byte_start)
        sequence = re.sub(rb"\s+", b"", self.handle.read(needed + extra))[:needed]
        if len(sequence) != needed:
            raise ValueError(f"Short reference fetch: {chrom}:{start}-{end}")
        return sequence.decode("ascii").upper()

    def close(self):
        self.handle.close()


COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


def map_local_interval(block, local_anchor: int, breakpoint: str, transcript_direction: int, reference: IndexedFasta):
    local_start, local_end, _ = block
    chrom, genomic_anchor, genomic_strand = parse_breakpoint(breakpoint)
    slope = 1 if genomic_strand == "+" else -1
    g1 = genomic_anchor + slope * (local_start - local_anchor)
    g2 = genomic_anchor + slope * (local_end - local_anchor)
    start, end = min(g1, g2), max(g1, g2)
    sequence = reference.fetch(chrom, start, end)
    genomic_traversal = transcript_direction * slope
    if genomic_traversal < 0:
        sequence = revcomp(sequence)
    traversal = "+" if genomic_traversal > 0 else "-"
    return sequence, f"chr{chrom}:{start}-{end}:{traversal}", {
        "chrom": chrom, "start": start, "end": end, "traversal": traversal,
        "genomic_strand": genomic_strand,
    }


def gene_transcripts(transcripts: dict, gene: str):
    gene = clean(gene)
    return transcripts.get(gene, transcripts.get(gene.split(".", 1)[0], []))


def transcript_priority(model: dict) -> int:
    tags = {tag.lower() for tag in model.get("tags", set())}
    if any("mane_select" in tag for tag in tags):
        return 3
    if any("ensembl_canonical" in tag for tag in tags):
        return 2
    if "basic" in tags:
        return 1
    return 0


def overlap_size(block: dict, exon) -> int:
    return max(0, min(block["end"], exon[1]) - max(block["start"], exon[0]) + 1)


def compatible_terminal_candidates(side_blocks: list[dict], gene: str, transcripts: dict,
                                   terminal_index: int, tolerance: int):
    """Match observed partner blocks to consecutive exons of each transcript."""
    candidates = []
    for model in gene_transcripts(transcripts, gene):
        if not side_blocks or model["chrom"] != side_blocks[0]["chrom"]:
            continue
        if model["strand"] != side_blocks[0]["genomic_strand"]:
            continue
        exons = model["exons"]
        width = len(side_blocks)
        for offset in range(0, len(exons) - width + 1):
            chosen = exons[offset:offset + width]
            overlaps = [overlap_size(block, exon) for block, exon in zip(side_blocks, chosen)]
            if any(value <= 0 for value in overlaps):
                continue
            # Every boundary between observed blocks must agree with the same
            # transcript's consecutive exon boundaries (small coordinate
            # normalization differences are tolerated).
            boundary_error = 0
            valid = True
            for left_block, right_block, left_exon, right_exon in zip(
                    side_blocks, side_blocks[1:], chosen, chosen[1:]):
                if model["strand"] == "+":
                    error = abs(left_block["end"] - left_exon[1]) + abs(right_block["start"] - right_exon[0])
                else:
                    error = abs(left_block["start"] - left_exon[0]) + abs(right_block["end"] - right_exon[1])
                if error > 2 * tolerance:
                    valid = False
                    break
                boundary_error += error
            if not valid:
                continue
            terminal_exon = chosen[terminal_index]
            candidates.append({
                "model": model, "terminal_exon": terminal_exon,
                "priority": transcript_priority(model),
                "overlap": sum(overlaps), "boundary_error": boundary_error,
                "match_mode": "CHAIN_COMPATIBLE_TRANSCRIPT",
            })
    return candidates


def terminal_block_fallback_candidates(side_blocks: list[dict], gene: str, transcripts: dict,
                                       terminal_index: int):
    """Find an annotated terminal exon without requiring an annotated internal chain.

    This is needed for genuine non-reference internal junction chains. It only
    supplies the outer boundary of the already observed terminal block and
    never changes or inserts an internal exon.
    """
    if not side_blocks:
        return []
    block = side_blocks[terminal_index]
    candidates = []
    seen = set()
    for model in gene_transcripts(transcripts, gene):
        if model["chrom"] != block["chrom"] or model["strand"] != block["genomic_strand"]:
            continue
        for exon in model["exons"]:
            overlap = overlap_size(block, exon)
            if overlap <= 0:
                continue
            key = (model["transcript_id"], exon)
            if key in seen:
                continue
            seen.add(key)
            candidates.append({
                "model": model, "terminal_exon": exon,
                "priority": transcript_priority(model), "overlap": overlap,
                "boundary_error": 0, "match_mode": "TERMINAL_BLOCK_FALLBACK",
            })
    return candidates


def complete_outer_terminal(block: dict, side_blocks: list[dict], gene: str, which: str,
                            transcripts: dict, reference: IndexedFasta, tolerance: int):
    """Complete an outer boundary using a junction-chain-compatible transcript."""
    terminal_index = 0 if which == "first" else -1
    candidates = compatible_terminal_candidates(side_blocks, gene, transcripts,
                                                 terminal_index, tolerance)
    if not candidates:
        candidates = terminal_block_fallback_candidates(
            side_blocks, gene, transcripts, terminal_index)
    if not candidates:
        return block, "NOVEL_TERMINAL_BOUNDARY", "", "NO_ANNOTATED_TERMINAL_EXON"
    ranked = sorted(candidates, key=lambda item: (
        item["priority"], -item["boundary_error"], item["overlap"]), reverse=True)
    best_rank = (ranked[0]["priority"], ranked[0]["boundary_error"], ranked[0]["overlap"])
    best = [item for item in ranked if
            (item["priority"], item["boundary_error"], item["overlap"]) == best_rank]
    boundaries = set()
    for item in best:
        exon = item["terminal_exon"]
        boundary = (exon[0] if (which == "first" and block["traversal"] == "+")
                    or (which == "last" and block["traversal"] == "-") else exon[1])
        boundaries.add(boundary)
    if len(boundaries) != 1:
        return block, "AMBIGUOUS_TERMINAL_TRANSCRIPT", "", ranked[0]["match_mode"]
    boundary = next(iter(boundaries))
    completed = dict(block)
    if (which == "first" and block["traversal"] == "+") or (which == "last" and block["traversal"] == "-"):
        completed["start"] = boundary
    else:
        completed["end"] = boundary
    sequence = reference.fetch(completed["chrom"], completed["start"], completed["end"])
    if completed["traversal"] == "-":
        sequence = revcomp(sequence)
    completed["sequence"] = sequence
    status = "EXACT_TERMINAL_EXON" if completed == block else "ANNOTATION_COMPLETED_TERMINAL_EXON"
    transcript_ids = ",".join(sorted({item["model"]["transcript_id"] for item in best}))
    return completed, status, transcript_ids, ranked[0]["match_mode"]


def reconstruct(info: dict, event: dict, reference: IndexedFasta, transcripts=None,
                terminal_tolerance: int = 5):
    blocks = parse_blocks(info.get("representative_block_chain", ""))
    if not blocks:
        return None, "NO_REPRESENTATIVE_BLOCK_CHAIN", {}
    left_local = integer(event.get("LeftLocalBreakpoint") or event.get("left_local_breakpoint"), -1)
    right_local = integer(event.get("RightLocalBreakpoint") or event.get("right_local_breakpoint"), -1)
    if left_local < 0 or right_local < 0 or left_local == right_local:
        return None, "INVALID_LOCAL_BREAKPOINTS", {}
    direction = -1 if clean(info.get("strand")) == "-" else 1
    ordered = sorted(blocks, key=lambda block: block[0], reverse=direction < 0)
    midpoint = (left_local + right_local) / 2
    left_is_low = left_local < right_local
    labels = ["left" if ((block[0] + block[1]) / 2 < midpoint) == left_is_low else "right" for block in ordered]
    transition_positions = [index for index, (a, b) in enumerate(zip(labels, labels[1:]), 1) if a != b]
    if len(transition_positions) != 1:
        return None, "MULTIPLE_OR_MISSING_PARTNER_TRANSITIONS", {"block_labels": ",".join(labels)}
    transition_index = transition_positions[0]
    parts, projected, projected_blocks = [], [], []
    for block, side in zip(ordered, labels):
        if side == "left":
            sequence, locus, details = map_local_interval(block, left_local, info["left_breakpoint"], direction, reference)
        else:
            sequence, locus, details = map_local_interval(block, right_local, info["right_breakpoint"], direction, reference)
        details.update({"sequence": sequence, "side": side})
        parts.append(sequence)
        projected.append(f"{side}:{locus}")
        projected_blocks.append(details)
    terminal_status = "NOT_REQUESTED"
    terminal_transcripts = ""
    if transcripts is not None and projected_blocks:
        first_gene = info.get("left_gene") if projected_blocks[0]["side"] == "left" else info.get("right_gene")
        last_gene = info.get("left_gene") if projected_blocks[-1]["side"] == "left" else info.get("right_gene")
        first_side_blocks = [block for block in projected_blocks if block["side"] == projected_blocks[0]["side"]]
        last_side_blocks = [block for block in projected_blocks if block["side"] == projected_blocks[-1]["side"]]
        projected_blocks[0], first_status, first_tx, first_mode = complete_outer_terminal(
            projected_blocks[0], first_side_blocks, first_gene, "first", transcripts,
            reference, terminal_tolerance)
        projected_blocks[-1], last_status, last_tx, last_mode = complete_outer_terminal(
            projected_blocks[-1], last_side_blocks, last_gene, "last", transcripts,
            reference, terminal_tolerance)
        parts = [block["sequence"] for block in projected_blocks]
        projected = [f"{block['side']}:chr{block['chrom']}:{block['start']}-{block['end']}:{block['traversal']}"
                     for block in projected_blocks]
        terminal_status = f"{first_status}|{last_status}"
        terminal_transcripts = f"{first_tx}|{last_tx}"
        terminal_match_modes = f"{first_mode}|{last_mode}"

    # Canonical database orientation is always the step-01 partner order:
    # left_gene -> right_gene.  A '-' CTAT mini-contig can encode the identical
    # molecule as right(-) -> left(-); reverse-complement the complete product
    # and every block-level representation instead of releasing both forms.
    original_partner_order = f"{labels[0]}_TO_{labels[-1]}"
    orientation_normalized = labels[0] == "right" and labels[-1] == "left"
    if orientation_normalized:
        ordered = list(reversed(ordered))
        labels = list(reversed(labels))
        direction *= -1
        transition_index = len(projected_blocks) - transition_index
        normalized_blocks = []
        for block in reversed(projected_blocks):
            normalized = dict(block)
            normalized["sequence"] = revcomp(block["sequence"])
            normalized["traversal"] = "+" if block["traversal"] == "-" else "-"
            normalized_blocks.append(normalized)
        projected_blocks = normalized_blocks
        parts = [block["sequence"] for block in projected_blocks]
        if transcripts is not None:
            terminal_status = "|".join(reversed(terminal_status.split("|")))
            terminal_transcripts = "|".join(reversed(terminal_transcripts.split("|")))
            terminal_match_modes = "|".join(reversed(terminal_match_modes.split("|")))
    projected = [f"{block['side']}:chr{block['chrom']}:{block['start']}-{block['end']}:{block['traversal']}"
                 for block in projected_blocks]

    # These are two different biological coordinates.  CTAT local anchors
    # describe the structural breakpoint in the mini-reference, whereas the
    # fusion RNA joins the fusion-facing boundaries of the last left-partner
    # block and first right-partner block.  Alternative splicing may therefore
    # place the RNA junction away from one or both structural anchors.
    left_facing_local = (ordered[transition_index - 1][1] if direction > 0
                         else ordered[transition_index - 1][0])
    right_facing_local = (ordered[transition_index][0] if direction > 0
                          else ordered[transition_index][1])
    junction_offset = sum(len(part) for part in parts[:transition_index])
    left_join_block = projected_blocks[transition_index - 1]
    right_join_block = projected_blocks[transition_index]
    left_rna_coordinate = (left_join_block["end"] if left_join_block["traversal"] == "+"
                           else left_join_block["start"])
    right_rna_coordinate = (right_join_block["start"] if right_join_block["traversal"] == "+"
                            else right_join_block["end"])
    left_direct = left_facing_local == left_local
    right_direct = right_facing_local == right_local
    junction_relationship = (
        "DIRECT_AT_BOTH_STRUCTURAL_BREAKPOINTS" if left_direct and right_direct else
        "LEFT_DIRECT_RIGHT_SPLICED" if left_direct else
        "LEFT_SPLICED_RIGHT_DIRECT" if right_direct else
        "SPLICED_AROUND_BOTH_STRUCTURAL_BREAKPOINTS"
    )
    return "".join(parts), "PASS", {
        "junction_offset": junction_offset,
        "block_count": len(parts),
        "partner_order": f"{labels[0]}_TO_{labels[-1]}",
        "original_partner_order": original_partner_order,
        "orientation_normalized_to_left_to_right": str(orientation_normalized).upper(),
        "five_prime_partner": labels[0],
        "three_prime_partner": labels[-1],
        "projected_block_chain": ";".join(projected),
        "local_left_breakpoint": left_local,
        "local_right_breakpoint": right_local,
        "local_rna_junction_left": left_facing_local,
        "local_rna_junction_right": right_facing_local,
        "left_rna_junction_distance_from_structural_breakpoint": abs(left_facing_local - left_local),
        "right_rna_junction_distance_from_structural_breakpoint": abs(right_facing_local - right_local),
        "rna_junction_relationship_to_structural_breakpoint": junction_relationship,
        "rna_fusion_junction": (
            f"chr{left_join_block['chrom']}:{left_rna_coordinate}:{left_join_block['traversal']}"
            f"->chr{right_join_block['chrom']}:{right_rna_coordinate}:{right_join_block['traversal']}"
        ),
        "terminal_exon_status": terminal_status,
        "terminal_transcript_ids": terminal_transcripts,
        "terminal_match_modes": terminal_match_modes if transcripts is not None else "",
    }


def write_fasta(path: Path, records):
    with gzip.open(path, "wt", encoding="utf-8") as out:
        for name, sequence in records:
            out.write(f">{name}\n")
            for start in range(0, len(sequence), 80):
                out.write(sequence[start:start + 80] + "\n")


def projected_chain_blocks(value: str):
    blocks = []
    pattern = re.compile(r"^(left|right):chr([^:]+):(\d+)-(\d+):([+-])$")
    for token in clean(value).split(";"):
        match = pattern.match(token)
        if match:
            blocks.append({"side": match.group(1), "chrom": match.group(2),
                           "start": int(match.group(3)), "end": int(match.group(4)),
                           "traversal": match.group(5)})
    return blocks


def projection_annotation_score(metrics: dict, info: dict, transcripts: dict):
    """Score a local-anchor hypothesis against the named partner annotations."""
    matched_blocks = total_overlap = strand_matches = 0
    projected = projected_chain_blocks(metrics.get("projected_block_chain", ""))
    for block in projected:
        gene = info.get("left_gene") if block["side"] == "left" else info.get("right_gene")
        best_overlap = 0
        best_strand = False
        for model in gene_transcripts(transcripts, gene):
            if clean(model.get("chrom")).removeprefix("chr") != clean(block["chrom"]).removeprefix("chr"):
                continue
            for exon in model.get("exons", []):
                # GTF exon tuples use genomic start/end in positions 0/1.
                overlap = max(0, min(block["end"], exon[1]) - max(block["start"], exon[0]) + 1)
                if overlap > best_overlap:
                    best_overlap = overlap
                    model_strand = model.get("strand", "")
                    best_strand = (not model_strand or model_strand == block["traversal"])
        if best_overlap:
            matched_blocks += 1
            total_overlap += best_overlap
            strand_matches += int(best_strand)
    return matched_blocks, total_overlap, strand_matches


def resolve_global_local_breakpoints(info: dict, reference: IndexedFasta, transcripts: dict):
    """Infer the fusion-contig anchors from the retained representative chain.

    Every adjacent block boundary is tested.  The winning boundary must give a
    unique best projection onto the two named GTF partner genes.  No sample,
    preliminary table, NR nucleotide sequence, or sequence alignment is used.
    """
    blocks = parse_blocks(info.get("representative_block_chain", ""))
    if len(blocks) < 2:
        return None, "NO_REPRESENTATIVE_BLOCK_CHAIN", {}
    direction = -1 if clean(info.get("strand")) == "-" else 1
    ordered = sorted(blocks, key=lambda block: block[0], reverse=direction < 0)
    candidates = []
    for index in range(1, len(ordered)):
        previous, following = ordered[index - 1], ordered[index]
        left_local, right_local = ((previous[1], following[0]) if direction > 0
                                   else (previous[0], following[1]))
        event = {"LeftLocalBreakpoint": left_local, "RightLocalBreakpoint": right_local}
        try:
            sequence, reason, metrics = reconstruct(info, event, reference, None)
        except Exception:
            continue
        if sequence is None:
            continue
        score = projection_annotation_score(metrics, info, transcripts)
        candidates.append((score, left_local, right_local, sequence,
                           metrics.get("projected_block_chain", "")))
    if not candidates:
        return None, "NO_VALID_LOCAL_BREAKPOINT_PROJECTION", {}
    candidates.sort(reverse=True)
    best_score = candidates[0][0]
    best = [candidate for candidate in candidates if candidate[0] == best_score]
    products = {(candidate[3], candidate[4]) for candidate in best}
    if best_score[0] == 0:
        return None, "NO_ANNOTATION_SUPPORTED_LOCAL_BREAKPOINT", {
            "local_breakpoint_candidate_count": len(candidates)}
    if len(products) != 1:
        return None, "AMBIGUOUS_LOCAL_BREAKPOINT_PROJECTION", {
            "local_breakpoint_candidate_count": len(candidates),
            "best_local_breakpoint_candidate_count": len(best),
            "best_annotation_matched_blocks": best_score[0],
        }
    chosen = best[0]
    return {"LeftLocalBreakpoint": chosen[1], "RightLocalBreakpoint": chosen[2]}, \
        "UNIQUE_GTF_SUPPORTED_BLOCK_TRANSITION", {
            "local_breakpoint_candidate_count": len(candidates),
            "best_local_breakpoint_candidate_count": len(best),
            "best_annotation_matched_blocks": best_score[0],
            "best_annotation_overlap_bases": best_score[1],
        }


def run_global_feature_mode(args):
    """Reconstruct exactly one sequence attempt per existing global HC entity."""
    entities = load_global_isoforms(args.global_isoforms, set())
    required = {"feature_id", "left_gene", "right_gene", "left_breakpoint",
                "right_breakpoint", "representative_block_chain", "representative_strand",
                "left_local_breakpoint", "right_local_breakpoint", "local_breakpoint_source"}
    missing = required - set(next(iter(entities.values()), {}))
    if missing:
        raise SystemExit(
            "Global isoform table lacks direct reconstruction fields: " +
            ", ".join(sorted(missing)) +
            ". Re-run the updated 01.build_standardized_fusion_keys.py first."
        )
    target_genes = {clean(row.get(key)) for row in entities.values()
                    for key in ("left_gene", "right_gene") if clean(row.get(key))}
    print(f"Loading GTF exons for {len(target_genes):,} fusion-partner identifiers...", flush=True)
    transcripts = load_target_transcripts(args.annotation_gtf, target_genes)
    reference = IndexedFasta(args.reference_fasta)
    records, audit = [], []
    for index, (feature_id, source) in enumerate(sorted(entities.items()), 1):
        info = dict(source)
        info["strand"] = clean(source.get("representative_strand")) or "."
        base = {"feature_id": feature_id, "fusion_name": clean(info.get("fusion_name")),
                "entity_level": "global_isoform"}
        if not clean(info.get("representative_block_chain")):
            audit.append({**base, "status": "FAILED", "reason": "NO_REPRESENTATIVE_BLOCK_CHAIN"})
            continue
        direct_left = integer(info.get("left_local_breakpoint"), -1)
        direct_right = integer(info.get("right_local_breakpoint"), -1)
        if direct_left >= 0 and direct_right >= 0 and direct_left != direct_right:
            event = {"LeftLocalBreakpoint": direct_left,
                     "RightLocalBreakpoint": direct_right}
            anchor_method = ("DIRECT_STEP01_LOCAL_BREAKPOINT:"
                             + (clean(info.get("local_breakpoint_source")) or "UNKNOWN_SOURCE"))
            anchor_metrics = {"local_breakpoint_candidate_count": 1,
                              "best_local_breakpoint_candidate_count": 1}
        else:
            event, anchor_method, anchor_metrics = resolve_global_local_breakpoints(
                info, reference, transcripts)
        if event is None:
            audit.append({**base, **anchor_metrics, "status": "FAILED",
                          "reason": anchor_method})
            continue
        try:
            sequence, reason, metrics = reconstruct(
                info, event, reference, transcripts, args.terminal_boundary_tolerance)
        except Exception as exc:
            audit.append({**base, **anchor_metrics, "status": "FAILED",
                          "reason": "REFERENCE_EXTRACTION_ERROR", "detail": str(exc)})
            continue
        if sequence is None:
            audit.append({**base, **anchor_metrics, **metrics, "status": "FAILED", "reason": reason})
            continue
        header = (f"{feature_id}|entity=GLOBAL_HC_ISOFORM|class=GENOME_CHAIN_REFERENCE_SEQUENCE"
                  f"|junction_offset={metrics['junction_offset']}"
                  f"|rna_fusion_junction={metrics['rna_fusion_junction']}"
                  f"|junction_relationship={metrics['rna_junction_relationship_to_structural_breakpoint']}"
                  f"|orientation_normalized={metrics['orientation_normalized_to_left_to_right']}"
                  f"|sample_count={integer(info.get('sample_count'))}"
                  f"|total_support_reads={integer(info.get('total_support_reads'))}")
        records.append((header, sequence))
        audit.append({**base, **anchor_metrics, **metrics, "status": "PASS", "reason": "",
                      "local_breakpoint_method": anchor_method,
                      "sequence_length": len(sequence)})
        if index % 500 == 0 or index == len(entities):
            print(f"[{index:,}/{len(entities):,} global isoforms] completed", flush=True)
    reference.close()
    write_fasta(args.output_dir / "global_hc_isoform_reference_sequences.fa.gz", records)
    fields = sorted({key for row in audit for key in row})
    with gzip.open(args.output_dir / "global_hc_isoform_reference_sequences.audit.tsv.gz",
                   "wt", encoding="utf-8", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(audit)
    summary = {
        "input_mode": "GLOBAL_ISOFORM_ENTITIES",
        "global_isoform_input_entities": len(entities),
        "sequence_attempts": len(audit),
        "sequences": len(records),
        "unique_reconstructed_sequences": len({sequence for _, sequence in records}),
        "status": dict(Counter(row.get("status", "") for row in audit)),
        "rna_junction_relationship_to_structural_breakpoint": dict(Counter(
            row.get("rna_junction_relationship_to_structural_breakpoint", "")
            for row in audit if row.get("rna_junction_relationship_to_structural_breakpoint"))),
        "orientation_normalized_to_left_to_right": dict(Counter(
            row.get("orientation_normalized_to_left_to_right", "")
            for row in audit if row.get("orientation_normalized_to_left_to_right"))),
        "failure_reasons": dict(Counter(row.get("reason", "") for row in audit
                                         if row.get("status") != "PASS")),
        "important_note": ("One reconstruction attempt is made directly per existing step-01 global "
                           "feature_id; sample directories and NR sequences are not consulted."),
    }
    (args.output_dir / "global_hc_isoform_reference_sequences.summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False))
    return summary


def validate_against_nr(output_dir: Path, reconstructed: dict[str, str], observed: dict[str, str],
                        minimap2: str, threads: int, max_secondary: int):
    """Validate NR observations against deduplicated reconstructed cores.

    Deduplication prevents minimap2 from assigning a query to an equivalent
    isoform target carrying another sample/feature ID.
    """
    target = output_dir / ".genome_chain_targets.fa"
    query = output_dir / ".nr_observations.fa"
    paf = output_dir / "genome_chain_core_to_nr.paf"
    sequence_to_cluster = {}
    cluster_sequences = {}
    expected_cluster = {}
    for name, sequence in reconstructed.items():
        cluster = sequence_to_cluster.get(sequence)
        if cluster is None:
            cluster = "GC_" + hashlib.sha256(sequence.encode()).hexdigest()[:24]
            sequence_to_cluster[sequence] = cluster
            cluster_sequences[cluster] = sequence
        expected_cluster[name] = cluster
    with target.open("w", encoding="utf-8") as out:
        for cluster, sequence in cluster_sequences.items():
            out.write(f">{cluster}\n{sequence}\n")
    with query.open("w", encoding="utf-8") as out:
        for name, sequence in observed.items():
            if sequence:
                out.write(f">{name}\n{sequence}\n")
    command = [minimap2, "-x", "map-ont", "-c", "--secondary=yes", "-N", str(max_secondary),
               "-t", str(threads), str(target), str(query)]
    with paf.open("w", encoding="utf-8") as out:
        subprocess.run(command, stdout=out, check=True)
    best = {}
    with paf.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12 or fields[5] != expected_cluster.get(fields[0]):
                continue
            qlen, qstart, qend = int(fields[1]), int(fields[2]), int(fields[3])
            tlen, tstart, tend = int(fields[6]), int(fields[7]), int(fields[8])
            nmatch, block = int(fields[9]), int(fields[10])
            metrics = {"validation_identity": nmatch / block if block else 0,
                       "validation_nr_coverage": (qend - qstart) / qlen if qlen else 0,
                       "validation_core_coverage": (tend - tstart) / tlen if tlen else 0,
                       "validation_mapq": int(fields[11])}
            rank = (nmatch, metrics["validation_identity"], metrics["validation_core_coverage"])
            if fields[0] not in best or rank > best[fields[0]][0]:
                best[fields[0]] = (rank, metrics)
    target.unlink(missing_ok=True)
    query.unlink(missing_ok=True)
    return ({name: value[1] for name, value in best.items()}, expected_cluster,
            len(cluster_sequences))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--metadata-root", type=Path,
                        help="Required only for legacy --hc-isoform-calls mode")
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument("--global-isoforms", type=Path,
                        help="Preferred one-row-per-feature step-01 global isoform table")
    inputs.add_argument("--hc-isoform-calls", type=Path,
                        help="Legacy sample-level evidence table")
    parser.add_argument("--sample-list", type=Path)
    parser.add_argument("--reference-fasta", type=Path, required=True)
    parser.add_argument("--annotation-gtf", type=Path, required=True,
                        help="CTAT-matched GTF used to complete only the two outer terminal exons")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--breakpoint-tolerance", type=int, default=5)
    parser.add_argument("--minimap2", default="minimap2")
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--validation-secondary", type=int, default=100)
    parser.add_argument("--min-validation-identity", type=float, default=0.85)
    parser.add_argument("--min-core-coverage", type=float, default=0.90)
    parser.add_argument("--audit-against-nr", action="store_true",
                        help="Optionally compare reconstructed cores with NR FASTA; never changes reconstruction PASS/FAIL")
    parser.add_argument("--representative-min-fraction", type=float, default=0.75,
                        help="Minimum fraction for a unique modal sequence to represent a cross-sample feature")
    parser.add_argument("--terminal-boundary-tolerance", type=int, default=5,
                        help="Allowed bp difference when matching projected splice boundaries to GTF")
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    if args.global_isoforms:
        run_global_feature_mode(args)
        return
    if args.metadata_root is None:
        parser.error("--metadata-root is required with --hc-isoform-calls")
    selected = sample_ids(args.sample_list)
    hc_by_sample = load_hc(args.hc_isoform_calls, selected)
    global_entities = {}
    input_mode = "LEGACY_SAMPLE_CALLS"
    target_genes = set()
    for rows in hc_by_sample.values():
        for row in rows:
            target_genes.update(filter(None, (clean(row.get("left_gene")), clean(row.get("right_gene")))))
            fusion = clean(row.get("fusion_name"))
            if "--" in fusion:
                target_genes.update(fusion.split("--", 1))
    print(f"Loading GTF exons for {len(target_genes):,} fusion-partner identifiers...", flush=True)
    transcripts = load_target_transcripts(args.annotation_gtf, target_genes)
    reference = IndexedFasta(args.reference_fasta)
    audit, all_records = [], []
    reconstructed_by_id, observed_by_id = {}, {}

    for sample in sorted(hc_by_sample):
        core = args.metadata_root / sample / "core"
        nr_meta_path = find_one(core, ["*.fusion_isoforms.nr_metadata.tsv", "*nr_metadata.tsv"])
        nr_fa_path = find_one(core, ["*.fusion_isoforms.nr.fa", "*nr.fa"])
        preliminary_path = find_one(core, ["*.ctat-LR-fusion.fusion_predictions.preliminary.tsv", "*preliminary.tsv"])
        if not nr_meta_path or not preliminary_path:
            for hc in hc_by_sample[sample]:
                audit.append({"sample_id": sample, "feature_id": clean(hc.get("feature_id")),
                              "status": "FAILED", "reason": "MISSING_CTAT_INPUT"})
            continue
        nr_rows = read_rows(nr_meta_path)
        nr_sequences = read_fasta(nr_fa_path) if args.audit_against_nr and nr_fa_path else {}
        prelim = preliminary_index(read_rows(preliminary_path))
        for hc in hc_by_sample[sample]:
            nr, association = match_nr(hc, nr_rows, args.breakpoint_tolerance)
            base = {"sample_id": sample, "feature_id": clean(hc.get("feature_id")),
                    "fusion_name": clean(hc.get("fusion_name")), "association_method": association,
                    "sample_support_reads": integer(hc.get("sample_support_reads"))}
            if nr is None:
                audit.append({**base, "status": "FAILED", "reason": association})
                continue
            info = dict(nr)
            for key in ("feature_id", "fusion_name", "left_gene", "right_gene",
                        "left_breakpoint", "right_breakpoint"):
                info[key] = clean(hc.get(key)) or clean(info.get(key))
            if (not clean(info.get("left_gene")) or not clean(info.get("right_gene"))) and "--" in info["fusion_name"]:
                fusion_left, fusion_right = info["fusion_name"].split("--", 1)
                info["left_gene"] = clean(info.get("left_gene")) or fusion_left
                info["right_gene"] = clean(info.get("right_gene")) or fusion_right
            nr_id = clean(info.get("nr_isoform_id"))
            event, event_method = match_preliminary(
                info, prelim, args.breakpoint_tolerance, reference, transcripts,
                args.terminal_boundary_tolerance)
            base.update({"nr_isoform_id": nr_id, "preliminary_match_method": event_method})
            if event is None:
                audit.append({**base, "status": "FAILED", "reason": event_method})
                continue
            try:
                sequence, reason, metrics = reconstruct(
                    info, event, reference, transcripts, args.terminal_boundary_tolerance)
            except Exception as exc:
                audit.append({**base, "status": "FAILED", "reason": "REFERENCE_EXTRACTION_ERROR",
                              "detail": str(exc)})
                continue
            if sequence is None:
                audit.append({**base, **metrics, "status": "FAILED", "reason": reason})
                continue
            nr_sequence = nr_sequences.get(nr_id, "")
            sequence_id = f"{sample}|{nr_id}|{clean(hc.get('feature_id'))}"
            header = (f"{sequence_id}|fusion={info['fusion_name']}|class=GENOME_CHAIN_CORE"
                      f"|junction_offset={metrics['junction_offset']}"
                      f"|rna_fusion_junction={metrics['rna_fusion_junction']}"
                      f"|junction_relationship={metrics['rna_junction_relationship_to_structural_breakpoint']}"
                      f"|orientation_normalized={metrics['orientation_normalized_to_left_to_right']}"
                      f"|support={base['sample_support_reads']}")
            all_records.append((header, sequence))
            reconstructed_by_id[sequence_id] = sequence
            if args.audit_against_nr:
                observed_by_id[sequence_id] = nr_sequence
            audit.append({**base, **metrics, "status": "PASS", "reason": "",
                          "sequence_id": sequence_id, "sequence_class": "GENOME_CHAIN_CORE",
                          "sequence_length": len(sequence),
                          "nr_sequence_length": len(nr_sequence) if args.audit_against_nr else "",
                          "length_difference_vs_nr": (len(sequence) - len(nr_sequence)
                                                      if args.audit_against_nr and nr_sequence else "")})
        print(f"[{len(audit):,} calls] {sample}: completed", flush=True)

    reference.close()
    write_fasta(args.output_dir / "genome_chain_fusion_isoform_core.fa.gz", all_records)
    unique_core_count = len(set(reconstructed_by_id.values()))

    # Sample files above are coordinate-translation evidence only. Collapse
    # them here to exactly one global entity per integrated step-01 feature_id.
    # Only the highest structural-quality projections are allowed to define
    # the feature sequence; lower-quality sample observations remain audit
    # evidence and cannot veto a better resolved feature-level reconstruction.
    rows_by_feature = defaultdict(list)
    row_by_sequence_id = {}
    for row in audit:
        if row.get("status") == "PASS" and row.get("sequence_id"):
            rows_by_feature[row.get("feature_id", "")].append(row)
            row_by_sequence_id[row["sequence_id"]] = row
    record_by_id = {header.split("|fusion=", 1)[0]: (header, sequence)
                    for header, sequence in all_records}
    representative_records, conflict_rows = [], []
    representative_methods = Counter()

    def terminal_quality(row):
        states = clean(row.get("terminal_exon_status")).split("|")
        state_score = {
            "EXACT_TERMINAL_EXON": 3,
            "ANNOTATION_COMPLETED_TERMINAL_EXON": 3,
            "NOVEL_TERMINAL_BOUNDARY": 2,
            "AMBIGUOUS_TERMINAL_TRANSCRIPT": 0,
        }
        modes = clean(row.get("terminal_match_modes")).split("|")
        mode_score = {
            "CHAIN_COMPATIBLE_TRANSCRIPT": 2,
            "TERMINAL_BLOCK_FALLBACK": 1,
            "NO_ANNOTATED_TERMINAL_EXON": 0,
        }
        return (sum(state_score.get(value, 0) for value in states),
                sum(mode_score.get(value, 0) for value in modes))

    for feature_id, rows in rows_by_feature.items():
        best_quality = max(terminal_quality(row) for row in rows)
        defining_rows = [row for row in rows if terminal_quality(row) == best_quality]
        variants = defaultdict(list)
        for row in defining_rows:
            variants[reconstructed_by_id[row["sequence_id"]]].append(row)
        ranked = sorted(variants.items(), key=lambda item: (
            len(item[1]), sum(integer(row.get("sample_support_reads")) for row in item[1]), len(item[0])),
            reverse=True)
        top_sequence, top_rows = ranked[0]
        fraction = len(top_rows) / len(defining_rows)
        tied = len(ranked) > 1 and len(ranked[1][1]) == len(top_rows)
        all_best_ambiguous = all("AMBIGUOUS_TERMINAL_TRANSCRIPT" in
                                 clean(row.get("terminal_exon_status"))
                                 for row in defining_rows)
        if all_best_ambiguous:
            method = "AMBIGUOUS_TERMINAL_TRANSCRIPT"
        elif len(ranked) == 1:
            method = "EXACT_CROSS_SAMPLE_CONSENSUS"
        elif not tied and fraction >= args.representative_min_fraction:
            method = "UNIQUE_MODAL_SEQUENCE"
        else:
            method = "UNRESOLVED_SEQUENCE_CONFLICT"
        representative_methods[method] += 1
        chosen_id = (top_rows[0]["sequence_id"]
                     if method not in {"UNRESOLVED_SEQUENCE_CONFLICT", "AMBIGUOUS_TERMINAL_TRANSCRIPT"} else "")
        if chosen_id:
            _, sequence = record_by_id[chosen_id]
            entity = global_entities.get(feature_id, {})
            total_support = (integer(entity.get("total_support_reads")) if entity else
                             sum(integer(row.get("sample_support_reads")) for row in rows))
            entity_sample_count = integer(entity.get("sample_count")) if entity else len(rows)
            representative_records.append((
                f"{feature_id}|entity=GLOBAL_HC_ISOFORM|representative_method={method}"
                f"|defining_evidence_fraction={fraction:.6f}|sample_count={entity_sample_count}"
                f"|total_support_reads={total_support}", sequence))
        for rank, (sequence, variant_rows) in enumerate(ranked, 1):
            conflict_rows.append({
                "feature_id": feature_id,
                "feature_status": method,
                "variant_rank": rank,
                "variant_sequence_length": len(sequence),
                "variant_sample_count": len(variant_rows),
                "defining_sample_count": len(defining_rows),
                "feature_sample_count": len(rows),
                "variant_sample_fraction": f"{len(variant_rows) / len(defining_rows):.6f}",
                "terminal_quality": f"{best_quality[0]}|{best_quality[1]}",
                "variant_support_reads": sum(integer(row.get("sample_support_reads")) for row in variant_rows),
                "representative": chosen_id in {row["sequence_id"] for row in variant_rows},
                "example_sequence_id": variant_rows[0]["sequence_id"],
                "sample_ids": ",".join(sorted({row.get("sample_id", "") for row in variant_rows})),
            })
        for row in rows:
            row["feature_sequence_status"] = method
            row["feature_sequence_variant_count"] = len(ranked)
            row["feature_representative_sample_fraction"] = f"{fraction:.6f}"
            row["is_feature_representative"] = row.get("sequence_id") == chosen_id
            row["feature_defining_evidence"] = row in defining_rows
            row["feature_terminal_quality"] = f"{terminal_quality(row)[0]}|{terminal_quality(row)[1]}"
    write_fasta(args.output_dir / "genome_chain_fusion_isoform_core.representative.fa.gz",
                representative_records)
    write_fasta(args.output_dir / "global_hc_isoform_reference_sequences.fa.gz",
                representative_records)
    with gzip.open(args.output_dir / "genome_chain_fusion_isoform_core.feature_sequence_variants.tsv.gz",
                   "wt", encoding="utf-8", newline="") as out:
        fields = list(conflict_rows[0]) if conflict_rows else ["feature_id", "feature_status"]
        writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(conflict_rows)
    validation = {}
    if args.audit_against_nr:
        validation, validation_clusters, _ = validate_against_nr(
            args.output_dir, reconstructed_by_id, observed_by_id,
            args.minimap2, args.threads, args.validation_secondary)
        for row in audit:
            sequence_id = row.get("sequence_id", "")
            if not sequence_id:
                continue
            metrics = validation.get(sequence_id)
            row["nr_audit_cluster_id"] = validation_clusters.get(sequence_id, "")
            if metrics is None:
                row["nr_audit_status"] = ("MISSING_NR_SEQUENCE" if not observed_by_id.get(sequence_id)
                                          else "NO_EXPECTED_CLUSTER_ALIGNMENT")
                continue
            row.update({key.replace("validation_", "nr_audit_"): f"{value:.6f}" if isinstance(value, float) else value
                        for key, value in metrics.items()})
            identity_pass = metrics["validation_identity"] >= args.min_validation_identity
            coverage_pass = metrics["validation_core_coverage"] >= args.min_core_coverage
            row["nr_audit_status"] = ("CONCORDANT" if identity_pass and coverage_pass else
                                      "LOW_IDENTITY_AND_CORE_COVERAGE" if not identity_pass and not coverage_pass else
                                      "LOW_IDENTITY" if not identity_pass else "LOW_CORE_COVERAGE")
    fields = sorted({key for row in audit for key in row})
    with gzip.open(args.output_dir / "genome_chain_fusion_isoform_core.audit.tsv.gz", "wt", encoding="utf-8", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader(); writer.writerows(audit)
    summary = {"input_mode": input_mode,
               "global_isoform_input_entities": len(global_entities) if global_entities else "",
               "coordinate_evidence_attempts": len(audit), "sequences": len(all_records),
               "unique_reconstructed_core_sequences": unique_core_count,
               "features_with_reconstructed_sequence": len(rows_by_feature),
               "feature_representative_sequences": len(representative_records),
               "feature_sequence_status": dict(representative_methods),
               "nr_audit_enabled": args.audit_against_nr,
               "status": dict(Counter(row.get("status", "") for row in audit)),
               "nr_audit_status": dict(Counter(row.get("nr_audit_status", "") for row in audit
                                                if row.get("nr_audit_status"))),
               "terminal_exon_status": dict(Counter(row.get("terminal_exon_status", "") for row in audit
                                                     if row.get("terminal_exon_status"))),
               "terminal_match_modes": dict(Counter(row.get("terminal_match_modes", "") for row in audit
                                                     if row.get("terminal_match_modes"))),
               "rna_junction_relationship_to_structural_breakpoint": dict(Counter(
                   row.get("rna_junction_relationship_to_structural_breakpoint", "")
                   for row in audit if row.get("rna_junction_relationship_to_structural_breakpoint"))),
               "orientation_normalized_to_left_to_right": dict(Counter(
                   row.get("orientation_normalized_to_left_to_right", "")
                   for row in audit if row.get("orientation_normalized_to_left_to_right"))),
               "failure_reasons": dict(Counter(row.get("reason", "") for row in audit if row.get("status") != "PASS")),
               "important_note": "PASS is determined only by deterministic chain projection and reference extraction. NR comparison is optional audit evidence and never gates sequence release. Unobserved transcript ends are not inferred."}
    (args.output_dir / "genome_chain_fusion_isoform_core.summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False))


if __name__ == "__main__":
    main()


