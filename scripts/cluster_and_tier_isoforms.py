#!/usr/bin/env python3
"""Cluster exact CTAT isoform entities and emit progressively stricter tiers.

The exact entity table is never modified.  Clustering follows the earlier
pan-cancer principles: directed fusion pair, splice type, compatible genomic
breakpoints, compatible internal junction chain, terminal subset collapsing,
and no use of minimap2 alignment direction.  Each cluster remains traceable to
all source feature IDs and read IDs.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path


BLOCK_RE = re.compile(r"^(.*):(\d+)-(\d+):([+.-])$")
BP_RE = re.compile(r"^(.*):(\d+):([+.-])$")


def clean(value):
    return "" if value is None else str(value).strip()


def integer(value, default=0):
    try:
        return int(float(clean(value)))
    except ValueError:
        return default


def number(value, default=0.0):
    try:
        return float(clean(value))
    except ValueError:
        return default


def read_tsv(path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def parse_chain(text):
    blocks = []
    for token in filter(None, map(str.strip, text.split(";"))):
        match = BLOCK_RE.match(token)
        if not match:
            raise ValueError(f"Invalid block-chain token: {token}")
        blocks.append((match.group(1), int(match.group(2)), int(match.group(3))))
    if not blocks:
        raise ValueError("Empty block chain")
    return tuple(blocks)


def internal_chain(blocks):
    """Return junctions only; transcript start/end coordinates are excluded."""
    return tuple((left[0], left[2], right[0], right[1])
                 for left, right in zip(blocks, blocks[1:]))


def breakpoint(value):
    match = BP_RE.match(clean(value))
    return (match.group(1), int(match.group(2)), match.group(3)) if match else None


def breakpoint_compatible(left, right, tolerance):
    a, b = breakpoint(left), breakpoint(right)
    return bool(a and b and a[0] == b[0] and a[2] == b[2]
                and abs(a[1] - b[1]) <= tolerance)


def junction_compatible(left, right, tolerance):
    return (left[0] == right[0] and left[2] == right[2]
            and abs(left[1] - right[1]) <= tolerance
            and abs(left[3] - right[3]) <= tolerance)


def chain_relation(representative, candidate, tolerance, collapse_terminal_subset=True):
    """Match equal chains or a contiguous terminally truncated subset."""
    if len(candidate) > len(representative):
        return None
    starts = range(len(representative) - len(candidate) + 1)
    if not collapse_terminal_subset and len(candidate) != len(representative):
        return None
    for offset in starts:
        section = representative[offset:offset + len(candidate)]
        if all(junction_compatible(a, b, tolerance) for a, b in zip(section, candidate)):
            relation = "EQUAL" if len(candidate) == len(representative) else "TERMINAL_SUBSET"
            return relation
    return None


def compatible(seed, candidate, junction_tolerance, breakpoint_tolerance,
               collapse_terminal_subset):
    return (
        seed["fusion_name"] == candidate["fusion_name"]
        and seed["splice_type"] == candidate["splice_type"]
        and breakpoint_compatible(seed["left_breakpoint"], candidate["left_breakpoint"], breakpoint_tolerance)
        and breakpoint_compatible(seed["right_breakpoint"], candidate["right_breakpoint"], breakpoint_tolerance)
        and chain_relation(seed["internal_chain"], candidate["internal_chain"],
                           junction_tolerance, collapse_terminal_subset) is not None
    )


def load_fasta(path):
    records, name, chunks = {}, None, []
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(chunks)
                name = line[1:].split("|entity=", 1)[0].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if name is not None:
        records[name] = "".join(chunks)
    return records


def write_tsv(path, rows, fields):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t",
                                lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_fasta(path, rows, sequences):
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        for row in rows:
            sequence = sequences.get(row["representative_feature_id"])
            if sequence:
                handle.write(f">{row['cluster_id']}|fusion={row['fusion_name']}"
                             f"|support={row['total_support_reads']}\n{sequence}\n")


def event_index(final_rows, max_promiscuity, min_pct_dominant,
                min_ref_reads, min_nonref_reads, min_ffpm):
    partners = defaultdict(set)
    for row in final_rows:
        left, right = clean(row.get("LeftGene")), clean(row.get("RightGene"))
        if left and right:
            partners[left].add(right)
            partners[right].add(left)
    result = {}
    for row_number, row in enumerate(final_rows, 1):
        fusion = clean(row.get("#FusionName") or row.get("FusionName"))
        if not fusion:
            fusion = f"{clean(row.get('LeftGene'))}--{clean(row.get('RightGene'))}"
        splice = clean(row.get("SpliceType")) or "."
        reads = integer(row.get("num_LR"))
        required = min_ref_reads if splice == "ONLY_REF_SPLICE" else min_nonref_reads
        frac = number(row.get("frac_dom_iso"), 1.0) * 100.0
        reasons = []
        if reads < required:
            reasons.append(f"event_reads_lt_{required}")
        if number(row.get("LR_FFPM")) < min_ffpm:
            reasons.append(f"event_ffpm_lt_{min_ffpm:g}")
        if frac < min_pct_dominant:
            reasons.append(f"dominant_fraction_lt_{min_pct_dominant:g}pct")
        left, right = clean(row.get("LeftGene")), clean(row.get("RightGene"))
        if max(len(partners[left]), len(partners[right])) > max_promiscuity:
            reasons.append(f"partner_count_gt_{max_promiscuity}")
        key = (fusion, clean(row.get("LeftBreakpoint")),
               clean(row.get("RightBreakpoint")), splice)
        result[key] = {"pass": not reasons, "reasons": ";".join(reasons),
                       "row": row_number, "num_LR": reads,
                       "LR_FFPM": clean(row.get("LR_FFPM")),
                       "frac_dom_iso": clean(row.get("frac_dom_iso"))}
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--entities", type=Path, required=True)
    parser.add_argument("--sequences", type=Path, required=True)
    parser.add_argument("--final-tsv", type=Path, required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--junction-tolerance", type=int, default=5)
    parser.add_argument("--breakpoint-tolerance", type=int, default=5)
    parser.add_argument("--no-collapse-terminal-subset", action="store_true")
    parser.add_argument("--min-ref-reads", type=int, default=1)
    parser.add_argument("--min-nonref-reads", type=int, default=2)
    parser.add_argument("--min-ffpm", type=float, default=0.0)
    parser.add_argument("--max-promiscuity", type=int, default=10)
    parser.add_argument("--min-pct-dominant", type=float, default=20.0)
    args = parser.parse_args()
    if min(args.junction_tolerance, args.breakpoint_tolerance) < 0:
        parser.error("Tolerances must be non-negative")

    csv.field_size_limit(sys.maxsize)
    exact = read_tsv(args.entities)
    sequences = load_fasta(args.sequences)
    final_rows = read_tsv(args.final_tsv)
    events = event_index(final_rows, args.max_promiscuity, args.min_pct_dominant,
                         args.min_ref_reads, args.min_nonref_reads, args.min_ffpm)
    prepared = []
    for row in exact:
        item = dict(row)
        item["blocks"] = parse_chain(row["representative_block_chain"])
        item["internal_chain"] = internal_chain(item["blocks"])
        item["support"] = integer(row.get("total_support_reads"))
        item["sequence_length"] = len(sequences.get(row["feature_id"], ""))
        prepared.append(item)

    grouped = defaultdict(list)
    for row in prepared:
        grouped[(row["fusion_name"], row["splice_type"])].append(row)
    clusters = []
    for _, group in sorted(grouped.items()):
        remaining = set(range(len(group)))
        ordered = sorted(range(len(group)), key=lambda i: (
            -len(group[i]["internal_chain"]), -group[i]["sequence_length"],
            -group[i]["support"], group[i]["feature_id"]))
        while remaining:
            seed_index = next(i for i in ordered if i in remaining)
            seed = group[seed_index]
            members = [i for i in ordered if i in remaining and compatible(
                seed, group[i], args.junction_tolerance, args.breakpoint_tolerance,
                not args.no_collapse_terminal_subset)]
            for index in members:
                remaining.remove(index)
            clusters.append((seed, [group[i] for i in members]))

    output_rows = []
    for seed, members in clusters:
        member_ids = sorted(x["feature_id"] for x in members)
        read_ids = sorted({rid for x in members
                           for rid in clean(x.get("support_read_ids")).split(",") if rid})
        digest = hashlib.sha256((args.sample + "\t" + "\t".join(member_ids)).encode()).hexdigest()[:20]
        event_key = (seed["fusion_name"], seed["left_breakpoint"],
                     seed["right_breakpoint"], seed["splice_type"])
        event = events.get(event_key, {"pass": False, "reasons": "event_not_found",
                                       "row": "", "num_LR": "", "LR_FFPM": "",
                                       "frac_dom_iso": ""})
        support = len(read_ids)
        output_rows.append({
            "cluster_id": f"{args.sample}|PANCLUSTER|{digest}",
            "sample_id": args.sample, "fusion_name": seed["fusion_name"],
            "left_gene": seed["left_gene"], "right_gene": seed["right_gene"],
            "left_breakpoint": seed["left_breakpoint"],
            "right_breakpoint": seed["right_breakpoint"],
            "splice_type": seed["splice_type"],
            "representative_feature_id": seed["feature_id"],
            "representative_block_chain": seed["representative_block_chain"],
            "internal_junction_chain": ";".join(
                f"{a}:{b}->{c}:{d}" for a, b, c, d in seed["internal_chain"]),
            "member_entity_count": len(members), "total_support_reads": support,
            "member_feature_ids": ",".join(member_ids),
            "support_read_ids": ",".join(read_ids),
            "ctat_final_row_number": event["row"], "ctat_event_num_LR": event["num_LR"],
            "ctat_event_LR_FFPM": event["LR_FFPM"],
            "ctat_event_frac_dom_iso": event["frac_dom_iso"],
            "pass_ctat_event_filter": str(event["pass"]).upper(),
            "ctat_event_filter_reason": event["reasons"],
            "pass_support_ge2": str(support >= 2).upper(),
            "pass_pancancer_hc3": str(event["pass"] and support >= 3).upper(),
            "pass_conservative_ge5": str(event["pass"] and support >= 5).upper(),
            "pass_conservative_ref_ge5": str(event["pass"] and support >= 5
                                                and seed["splice_type"] == "ONLY_REF_SPLICE").upper(),
        })

    output_rows.sort(key=lambda x: (x["fusion_name"], -integer(x["total_support_reads"]), x["cluster_id"]))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    fields = list(output_rows[0]) if output_rows else [
        "cluster_id", "sample_id", "fusion_name", "total_support_reads"]
    tier_specs = [
        ("01_clustered_all", lambda row: True),
        ("02_clustered_ctat_event_pass", lambda row: row["pass_ctat_event_filter"] == "TRUE"),
        ("03_clustered_support_ge2", lambda row: row["pass_ctat_event_filter"] == "TRUE" and integer(row["total_support_reads"]) >= 2),
        ("04_pancancer_hc3", lambda row: row["pass_pancancer_hc3"] == "TRUE"),
        ("05_conservative_ge5", lambda row: row["pass_conservative_ge5"] == "TRUE"),
        ("06_conservative_ref_ge5", lambda row: row["pass_conservative_ref_ge5"] == "TRUE"),
    ]
    summary = [{"tier": "00_raw_exact_entities", "isoform_count": len(exact),
                "fusion_count": len({x["fusion_name"] for x in exact}),
                "support_reads": sum(integer(x.get("total_support_reads")) for x in exact),
                "definition": "Unmodified exact block-chain entities"}]
    definitions = {
        "01_clustered_all": "5-bp breakpoint/junction tolerance; terminal-subset collapse; alignment strand ignored",
        "02_clustered_ctat_event_pass": "Tier 01 + pan-cancer CTAT event prefilter",
        "03_clustered_support_ge2": "Tier 02 + clustered isoform support >=2",
        "04_pancancer_hc3": "Tier 02 + configured pan-cancer single-sample support >=3",
        "05_conservative_ge5": "Tier 02 + support >=5",
        "06_conservative_ref_ge5": "Tier 05 + ONLY_REF_SPLICE",
    }
    for tier, predicate in tier_specs:
        selected = [row for row in output_rows if predicate(row)]
        write_tsv(args.output_dir / f"{args.sample}.{tier}.tsv", selected, fields)
        write_fasta(args.output_dir / f"{args.sample}.{tier}.fa.gz", selected, sequences)
        summary.append({"tier": tier, "isoform_count": len(selected),
                        "fusion_count": len({x["fusion_name"] for x in selected}),
                        "support_reads": sum(integer(x["total_support_reads"]) for x in selected),
                        "definition": definitions[tier]})
    write_tsv(args.output_dir / f"{args.sample}.isoform_tier_summary.tsv", summary,
              ["tier", "isoform_count", "fusion_count", "support_reads", "definition"])

    by_fusion = []
    for fusion in sorted({x["fusion_name"] for x in exact}):
        raw = [x for x in exact if x["fusion_name"] == fusion]
        row = {"fusion_name": fusion, "00_raw_exact_entities": len(raw)}
        for tier, predicate in tier_specs:
            row[tier] = sum(x["fusion_name"] == fusion and predicate(x) for x in output_rows)
        by_fusion.append(row)
    write_tsv(args.output_dir / f"{args.sample}.isoform_tier_counts_by_fusion.tsv", by_fusion,
              ["fusion_name", "00_raw_exact_entities"] + [x[0] for x in tier_specs])
    (args.output_dir / f"{args.sample}.isoform_tier_summary.json").write_text(
        json.dumps({"sample": args.sample, "settings": vars(args) | {
            "entities": str(args.entities), "sequences": str(args.sequences),
            "final_tsv": str(args.final_tsv), "output_dir": str(args.output_dir)},
                    "tiers": summary}, ensure_ascii=False, indent=2, default=str) + "\n",
        encoding="utf-8")
    print(json.dumps({row["tier"]: row["isoform_count"] for row in summary}, ensure_ascii=False))


if __name__ == "__main__":
    main()
