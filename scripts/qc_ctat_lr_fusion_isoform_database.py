#!/usr/bin/env python3
"""
QC summaries for a CTAT-LR-fusion isoform database.

Checks included:
  1. Fusion event/gene-pair support distribution from CTAT-LR-fusion TSV.
  2. Number of NR isoforms per fusion event/breakpoint.
  3. SpliceType composition by NR isoform support count.
  4. Candidate high-confidence subsets by support threshold and optional splice type.
"""

import argparse
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path


def set_csv_field_size_limit():
    limit = sys.maxsize
    while True:
        try:
            csv.field_size_limit(limit)
            return
        except OverflowError:
            limit //= 10


def parse_args():
    parser = argparse.ArgumentParser(description="QC CTAT-LR-fusion isoform database outputs.")
    parser.add_argument("--fusion_tsv", required=True, help="CTAT-LR-fusion fusion_predictions TSV.")
    parser.add_argument("--nr_metadata", required=True, help="fusion_isoforms.nr_metadata.tsv.")
    parser.add_argument("--nr_support", required=True, help="fusion_isoform_read_support.tsv.")
    parser.add_argument("--raw_manifest", required=True, help="fusion_isoforms.raw_manifest.tsv.")
    parser.add_argument("--outdir", required=True, help="QC output directory.")
    parser.add_argument(
        "--high_conf_support",
        type=int,
        default=2,
        help="Support threshold for high-confidence NR isoforms. Default: 2.",
    )
    parser.add_argument(
        "--top_n",
        type=int,
        default=50,
        help="Number of top rows to write for ranked QC tables. Default: 50.",
    )
    return parser.parse_args()


def read_tsv(path):
    with open(path, "r", newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def write_tsv(path, rows, fieldnames):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def get_int(row, name, default=0):
    value = row.get(name, "")
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return default


def event_key(row):
    return (
        row.get("fusion_name", "."),
        row.get("left_breakpoint", "."),
        row.get("right_breakpoint", "."),
    )


def qc_event_support(fusion_tsv, outdir, top_n):
    rows = []
    support_counter = Counter()

    for idx, row in enumerate(read_tsv(fusion_tsv), start=1):
        left_gene = row.get("LeftGene", ".")
        right_gene = row.get("RightGene", ".")
        fusion_name = f"{left_gene}--{right_gene}"
        lr_accessions = [x.strip() for x in row.get("LR_accessions", "").split(",") if x.strip()]
        support = get_int(row, "num_LR", len(lr_accessions))
        if support == 0:
            support = len(lr_accessions)
        support_counter[support] += 1
        rows.append(
            {
                "event_rank": idx,
                "fusion_name": fusion_name,
                "left_gene": left_gene,
                "right_gene": right_gene,
                "left_breakpoint": row.get("LeftBreakpoint", "."),
                "right_breakpoint": row.get("RightBreakpoint", "."),
                "splice_type": row.get("SpliceType", "."),
                "num_LR": support,
                "LR_accession_count": len(lr_accessions),
            }
        )

    dist_rows = [
        {"event_support_num_LR": support, "event_count": count}
        for support, count in sorted(support_counter.items())
    ]
    write_tsv(outdir / "01_event_support_distribution.tsv", dist_rows, ["event_support_num_LR", "event_count"])

    top_rows = sorted(rows, key=lambda r: (int(r["num_LR"]), int(r["LR_accession_count"])), reverse=True)[:top_n]
    write_tsv(
        outdir / "01_top_events_by_num_LR.tsv",
        top_rows,
        [
            "event_rank",
            "fusion_name",
            "left_gene",
            "right_gene",
            "left_breakpoint",
            "right_breakpoint",
            "splice_type",
            "num_LR",
            "LR_accession_count",
        ],
    )

    return len(rows), support_counter


def qc_isoforms_per_event(nr_metadata, outdir, top_n):
    event_to_isoforms = defaultdict(list)
    event_to_support_sum = Counter()

    for row in read_tsv(nr_metadata):
        key = event_key(row)
        support = get_int(row, "support_read_count")
        event_to_isoforms[key].append(row)
        event_to_support_sum[key] += support

    rows = []
    isoform_count_dist = Counter()
    for key, isoforms in event_to_isoforms.items():
        fusion_name, left_bp, right_bp = key
        isoform_count = len(isoforms)
        isoform_count_dist[isoform_count] += 1
        supports = sorted((get_int(row, "support_read_count") for row in isoforms), reverse=True)
        rows.append(
            {
                "fusion_name": fusion_name,
                "left_breakpoint": left_bp,
                "right_breakpoint": right_bp,
                "nr_isoform_count": isoform_count,
                "total_nr_support_reads": event_to_support_sum[key],
                "max_isoform_support": supports[0] if supports else 0,
                "support_distribution_within_event": ",".join(map(str, supports)),
            }
        )

    dist_rows = [
        {"nr_isoforms_per_event": isoform_count, "event_count": count}
        for isoform_count, count in sorted(isoform_count_dist.items())
    ]
    write_tsv(outdir / "02_nr_isoforms_per_event_distribution.tsv", dist_rows, ["nr_isoforms_per_event", "event_count"])

    top_rows = sorted(
        rows,
        key=lambda r: (int(r["nr_isoform_count"]), int(r["total_nr_support_reads"]), int(r["max_isoform_support"])),
        reverse=True,
    )[:top_n]
    write_tsv(
        outdir / "02_top_events_by_nr_isoform_count.tsv",
        top_rows,
        [
            "fusion_name",
            "left_breakpoint",
            "right_breakpoint",
            "nr_isoform_count",
            "total_nr_support_reads",
            "max_isoform_support",
            "support_distribution_within_event",
        ],
    )

    return len(event_to_isoforms), isoform_count_dist


def qc_splice_type_by_support(nr_metadata, raw_manifest, outdir):
    nr_counter = Counter()
    raw_counter = Counter()

    for row in read_tsv(nr_metadata):
        support = get_int(row, "support_read_count")
        splice_type = row.get("splice_type", ".")
        nr_counter[(support, splice_type)] += 1

    for row in read_tsv(raw_manifest):
        splice_type = row.get("splice_type", ".")
        raw_counter[splice_type] += 1

    nr_rows = [
        {"support_read_count": support, "splice_type": splice_type, "nr_isoform_count": count}
        for (support, splice_type), count in sorted(nr_counter.items(), key=lambda x: (x[0][0], x[0][1]))
    ]
    write_tsv(
        outdir / "03_nr_support_by_splice_type.tsv",
        nr_rows,
        ["support_read_count", "splice_type", "nr_isoform_count"],
    )

    raw_rows = [
        {"splice_type": splice_type, "raw_isoform_count": count}
        for splice_type, count in sorted(raw_counter.items())
    ]
    write_tsv(outdir / "03_raw_splice_type_distribution.tsv", raw_rows, ["splice_type", "raw_isoform_count"])

    return nr_counter, raw_counter


def qc_high_confidence(nr_metadata, outdir, high_conf_support):
    rows_all = list(read_tsv(nr_metadata))

    high_conf = [
        row for row in rows_all if get_int(row, "support_read_count") >= high_conf_support
    ]
    high_conf_ref = [
        row
        for row in high_conf
        if row.get("splice_type") == "ONLY_REF_SPLICE"
    ]
    singleton_non_ref = [
        row
        for row in rows_all
        if get_int(row, "support_read_count") == 1 and row.get("splice_type") == "INCL_NON_REF_SPLICE"
    ]

    fields = [
        "nr_isoform_id",
        "fusion_name",
        "left_gene",
        "right_gene",
        "left_breakpoint",
        "right_breakpoint",
        "splice_type",
        "fusion_contig",
        "strand",
        "exon_block_count",
        "representative_read_id",
        "representative_sequence_length",
        "support_read_count",
        "internal_junction_chain",
        "representative_block_chain",
    ]

    write_tsv(outdir / f"04_high_conf_support_ge_{high_conf_support}.tsv", high_conf, fields)
    write_tsv(outdir / f"04_high_conf_support_ge_{high_conf_support}_ONLY_REF_SPLICE.tsv", high_conf_ref, fields)
    write_tsv(outdir / "04_singleton_INCL_NON_REF_SPLICE.tsv", singleton_non_ref, fields)

    summary = [
        {"category": "all_nr_isoforms", "count": len(rows_all)},
        {"category": f"support_ge_{high_conf_support}", "count": len(high_conf)},
        {"category": f"support_ge_{high_conf_support}_ONLY_REF_SPLICE", "count": len(high_conf_ref)},
        {"category": "singleton_INCL_NON_REF_SPLICE", "count": len(singleton_non_ref)},
    ]
    write_tsv(outdir / "04_high_confidence_summary.tsv", summary, ["category", "count"])

    return summary


def main():
    set_csv_field_size_limit()

    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    event_count, event_support_counter = qc_event_support(args.fusion_tsv, outdir, args.top_n)
    nr_event_count, isoform_count_dist = qc_isoforms_per_event(args.nr_metadata, outdir, args.top_n)
    nr_splice_counter, raw_splice_counter = qc_splice_type_by_support(args.nr_metadata, args.raw_manifest, outdir)
    high_conf_summary = qc_high_confidence(args.nr_metadata, outdir, args.high_conf_support)

    print(f"QC output directory: {outdir}")
    print(f"Fusion events in fusion TSV: {event_count}")
    print(f"Fusion events with NR isoforms: {nr_event_count}")
    print(f"Event support distribution: {dict(sorted(event_support_counter.items()))}")
    print(f"NR isoforms per event distribution: {dict(sorted(isoform_count_dist.items()))}")
    print(f"Raw splice type distribution: {dict(sorted(raw_splice_counter.items()))}")
    print(f"High-confidence summary: {high_conf_summary}")


if __name__ == "__main__":
    main()

