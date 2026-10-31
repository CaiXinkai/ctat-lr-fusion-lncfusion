#!/usr/bin/env python3
"""Convert CTAT final calls and fusion-contig alignments into sequence entities.

Only events present in ctat-LR-fusion.fusion_predictions.tsv are admitted.
The preliminary table supplies local fusion-contig anchors for coordinate
projection; it is never used as a source of additional fusion calls.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from collections import defaultdict
from pathlib import Path


def clean(value):
    return "" if value is None else str(value).strip()


def first(row, *names):
    for name in names:
        if clean(row.get(name)):
            return clean(row[name])
    return ""


def read_tsv(path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def fusion_name(row):
    value = first(row, "#FusionName", "FusionName", "fusion_name")
    if value:
        return value
    left = first(row, "LeftGene", "left_gene")
    right = first(row, "RightGene", "right_gene")
    return f"{left}--{right}" if left and right else ""


def local_signature(row):
    left = first(row, "LeftLocalBreakpoint", "left_local_breakpoint")
    right = first(row, "RightLocalBreakpoint", "right_local_breakpoint")
    try:
        return str(int(left)), str(int(right))
    except ValueError:
        return "", ""


def preliminary_anchors(rows):
    index = defaultdict(set)
    for row in rows:
        key = (
            fusion_name(row),
            first(row, "LeftBreakpoint", "left_breakpoint"),
            first(row, "RightBreakpoint", "right_breakpoint"),
        )
        signature = local_signature(row)
        if all(key) and all(signature):
            index[key].add(signature)
    return index


def parse_attrs(text):
    attrs = {}
    for item in text.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            attrs[key] = value
    return attrs


def read_blocks(path):
    by_read = defaultdict(list)
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "cDNA_match":
                continue
            attrs = parse_attrs(fields[8])
            target = clean(attrs.get("Target")).split()
            if not target:
                continue
            by_read[target[0]].append(
                (fields[0], int(fields[3]), int(fields[4]), fields[6])
            )
    return by_read


def read_ids(row):
    value = first(row, "LR_accessions", "LongReads", "long_reads")
    return [item.strip() for item in value.replace(";", ",").split(",") if item.strip()]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--final-tsv", type=Path, required=True)
    parser.add_argument("--preliminary-tsv", type=Path, required=True)
    parser.add_argument("--fusion-gff3", type=Path, required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    csv.field_size_limit(sys.maxsize)
    final_rows = read_tsv(args.final_tsv)
    anchors = preliminary_anchors(read_tsv(args.preliminary_tsv))
    blocks_by_read = read_blocks(args.fusion_gff3)
    entities = {}
    missing_reads = []

    for event_number, row in enumerate(final_rows, 1):
        fusion = fusion_name(row)
        left_gene = first(row, "LeftGene", "left_gene")
        right_gene = first(row, "RightGene", "right_gene")
        left_bp = first(row, "LeftBreakpoint", "left_breakpoint")
        right_bp = first(row, "RightBreakpoint", "right_breakpoint")
        splice_type = first(row, "SpliceType", "splice_type") or "."
        if not all((fusion, left_gene, right_gene, left_bp, right_bp)):
            raise SystemExit(f"Final CTAT row {event_number} lacks fusion/gene/breakpoint fields")
        ids = read_ids(row)
        if not ids:
            raise SystemExit(f"Final CTAT row {event_number} lacks LR_accessions")
        anchor_set = anchors.get((fusion, left_bp, right_bp), set())
        local_left, local_right = next(iter(anchor_set)) if len(anchor_set) == 1 else ("", "")
        anchor_source = "CTAT_PRELIMINARY_EXACT" if len(anchor_set) == 1 else "GTF_INFERRED"

        for read_id in ids:
            blocks = blocks_by_read.get(read_id, [])
            expected = [block for block in blocks if block[0] == fusion]
            if expected:
                blocks = expected
            if not blocks:
                missing_reads.append(read_id)
                continue
            strands = {block[3] for block in blocks if block[3] in {"+", "-"}}
            strand = next(iter(strands)) if len(strands) == 1 else "."
            ordered = sorted(blocks, key=lambda item: (item[1], item[2], item[0], item[3]))
            chain = ";".join(f"{seq}:{start}-{end}:{block_strand}" for seq, start, end, block_strand in ordered)
            key = (fusion, left_bp, right_bp, splice_type, strand, chain, local_left, local_right)
            entity = entities.setdefault(key, {
                "fusion_name": fusion,
                "left_gene": left_gene,
                "right_gene": right_gene,
                "left_breakpoint": left_bp,
                "right_breakpoint": right_bp,
                "splice_type": splice_type,
                "representative_block_chain": chain,
                "representative_strand": strand,
                "left_local_breakpoint": local_left,
                "right_local_breakpoint": local_right,
                "local_breakpoint_source": anchor_source,
                "support_read_ids": set(),
                "ctat_final_row_numbers": set(),
            })
            entity["support_read_ids"].add(read_id)
            entity["ctat_final_row_numbers"].add(str(event_number))

    fields = [
        "feature_id", "sample_id", "fusion_name", "left_gene", "right_gene",
        "left_breakpoint", "right_breakpoint", "splice_type",
        "representative_block_chain", "representative_strand",
        "left_local_breakpoint", "right_local_breakpoint", "local_breakpoint_source",
        "sample_count", "total_support_reads", "support_read_ids", "ctat_final_row_numbers",
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for key, entity in sorted(entities.items()):
            digest = hashlib.sha256((args.sample + "\t" + "\t".join(key)).encode()).hexdigest()[:20]
            reads = sorted(entity.pop("support_read_ids"))
            rows = sorted(entity.pop("ctat_final_row_numbers"), key=int)
            writer.writerow({
                "feature_id": f"{args.sample}|CTATISO|{digest}",
                "sample_id": args.sample,
                **entity,
                "sample_count": 1,
                "total_support_reads": len(reads),
                "support_read_ids": ",".join(reads),
                "ctat_final_row_numbers": ",".join(rows),
            })

    print(json.dumps({
        "ctat_final_rows": len(final_rows),
        "isoform_entities": len(entities),
        "supporting_reads_without_gff3_blocks": len(set(missing_reads)),
        "policy": "Final CTAT rows only; no support threshold and no fuzzy NR clustering",
    }, ensure_ascii=False))


if __name__ == "__main__":
    main()

