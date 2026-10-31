#!/usr/bin/env python3
"""
Build a non-redundant fusion isoform database from CTAT-LR-fusion outputs.

Definition used here:
  raw isoform = one LR_accession/read-supported fusion transcript
  NR isoform  = raw isoforms with compatible fusion breakpoints and internal
                splice-junction chains collapsed together. Terminal 5'/3'
                block boundary differences are tolerated, and small alignment
                jitter at fusion breakpoints / internal junctions can be merged
                with user-defined bp tolerances.

Required CTAT-LR-fusion files:
  1. ctat-LR-fusion.fusion_predictions.tsv or preliminary.tsv
  2. fusion_intermediates_dir/chimeric_read_candidates.transcripts.fa.revised.fasta
  3. fusion_intermediates_dir/LR-FI.mm2.fusion_transcripts.gff3

The GTF coordinates are on LR-FI synthetic fusion contigs, not hg38 chromosomes.
Use the output GTF together with fusion_intermediates_dir/LR-FI_targets.fa.
"""

import argparse
import csv
import sys
from collections import defaultdict
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
    parser = argparse.ArgumentParser(
        description="Generate raw and non-redundant fusion isoform FASTA/GTF files from CTAT-LR-fusion outputs."
    )
    parser.add_argument("--fusion_tsv", required=True, help="CTAT-LR-fusion fusion_predictions TSV.")
    parser.add_argument(
        "--read_fasta",
        required=True,
        help="chimeric_read_candidates.transcripts.fa.revised.fasta from fusion_intermediates_dir.",
    )
    parser.add_argument(
        "--fusion_gff3",
        required=True,
        help="LR-FI.mm2.fusion_transcripts.gff3 from fusion_intermediates_dir.",
    )
    parser.add_argument("--outdir", required=True, help="Output directory.")
    parser.add_argument(
        "--min_nr_support",
        type=int,
        default=1,
        help="Minimum number of raw reads required to keep a non-redundant isoform. Default: 1.",
    )
    parser.add_argument(
        "--id_prefix",
        default="CTATLRF",
        help="Prefix for non-redundant isoform IDs. Default: CTATLRF.",
    )
    parser.add_argument(
        "--junction_tolerance",
        type=int,
        default=5,
        help="Maximum bp difference allowed when comparing internal splice junction coordinates. Default: 5.",
    )
    parser.add_argument(
        "--breakpoint_tolerance",
        type=int,
        default=5,
        help="Maximum bp difference allowed when comparing fusion breakpoint coordinates. Default: 5.",
    )
    parser.add_argument(
        "--no_collapse_terminal_subset",
        action="store_true",
        help=(
            "Disable collapsing reads whose internal junction chain is a contiguous subset of another read. "
            "By default, terminally truncated/ISM-like reads are collapsed into the longer compatible isoform."
        ),
    )
    parser.add_argument(
        "--keep_multicontig_blocks",
        action="store_true",
        help=(
            "Keep all GFF3 cDNA_match blocks for a read, including alignments to other synthetic fusion contigs. "
            "By default, only blocks on the expected contig matching fusion_name are used."
        ),
    )
    parser.add_argument(
        "--ignore_alignment_strand_for_clustering",
        action="store_true",
        help=(
            "Ignore minimap2 alignment strand when clustering isoforms. Useful when plus/minus mappings to the same "
            "synthetic fusion contig represent the same transcript structure."
        ),
    )
    parser.add_argument(
        "--allow_mixed_strand_reads",
        action="store_true",
        help=(
            "Keep raw reads whose cDNA_match blocks contain both + and - strands on the expected fusion contig. "
            "By default these biologically inconsistent reads are excluded from raw GTF and NR clustering."
        ),
    )
    return parser.parse_args()


def fasta_iter(path):
    name = None
    parts = []
    with open(path, "r") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(parts)
                name = line[1:].strip()
                parts = []
            else:
                parts.append(line.strip())
        if name is not None:
            yield name, "".join(parts)


def wrap(seq, width=60):
    for i in range(0, len(seq), width):
        yield seq[i : i + width]


def parse_gff3_attrs(attr_text):
    attrs = {}
    for item in attr_text.split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            attrs[key] = value
    return attrs


def gtf_quote(value):
    return str(value).replace('"', '\\"')


def gtf_attrs(pairs):
    return " ".join(f'{key} "{gtf_quote(value)}";' for key, value in pairs)


def first_present(row, names, default="."):
    for name in names:
        value = row.get(name)
        if value not in (None, ""):
            return value
    return default


def load_fusion_events(fusion_tsv):
    raw_isoforms = {}
    reads_to_raw_ids = defaultdict(list)

    with open(fusion_tsv, "r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"LeftGene", "RightGene", "LR_accessions"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise SystemExit(f"ERROR: fusion TSV is missing required columns: {', '.join(sorted(missing))}")

        for event_num, row in enumerate(reader, start=1):
            left_gene = row["LeftGene"]
            right_gene = row["RightGene"]
            fusion_name = f"{left_gene}--{right_gene}"
            event_id = f"event{event_num:06d}"

            meta = {
                "event_id": event_id,
                "fusion_name": fusion_name,
                "left_gene": left_gene,
                "right_gene": right_gene,
                "left_breakpoint": first_present(row, ["LeftBreakpoint"]),
                "right_breakpoint": first_present(row, ["RightBreakpoint"]),
                "left_local_breakpoint": first_present(row, ["LeftLocalBreakpoint"]),
                "right_local_breakpoint": first_present(row, ["RightLocalBreakpoint"]),
                "splice_type": first_present(row, ["SpliceType"]),
                "num_LR": first_present(row, ["num_LR", "Num_LR"]),
                "LR_FFPM": first_present(row, ["LR_FFPM"]),
                "annots": first_present(row, ["annots"]),
            }

            for read_id in row["LR_accessions"].split(","):
                read_id = read_id.strip()
                if not read_id:
                    continue
                raw_id = f"{fusion_name}|{event_id}|{read_id}"
                raw_isoforms[raw_id] = dict(meta, read_id=read_id, raw_id=raw_id)
                reads_to_raw_ids[read_id].append(raw_id)

    return raw_isoforms, reads_to_raw_ids


def load_read_sequences(read_fasta, reads_to_raw_ids):
    seq_by_raw_id = {}

    for header, seq in fasta_iter(read_fasta):
        read_id = header.split()[0]
        if read_id not in reads_to_raw_ids:
            continue
        for raw_id in reads_to_raw_ids[read_id]:
            seq_by_raw_id[raw_id] = seq.upper()

    return seq_by_raw_id


def load_gff3_blocks(fusion_gff3, reads_to_raw_ids, raw_isoforms, keep_multicontig_blocks):
    blocks_by_raw_id = defaultdict(list)

    with open(fusion_gff3, "r") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            seqname, source, feature, start, end, score, strand, phase, attrs_text = fields
            if feature != "cDNA_match":
                continue

            attrs = parse_gff3_attrs(attrs_text)
            target = attrs.get("Target", "")
            target_parts = target.split()
            if not target_parts:
                continue

            read_id = target_parts[0]
            if read_id not in reads_to_raw_ids:
                continue

            block = {
                "seqname": seqname,
                "source": source,
                "feature": feature,
                "start": int(start),
                "end": int(end),
                "score": score,
                "strand": strand,
                "phase": phase,
                "target_start": target_parts[1] if len(target_parts) > 1 else ".",
                "target_end": target_parts[2] if len(target_parts) > 2 else ".",
            }
            for raw_id in reads_to_raw_ids[read_id]:
                expected_contig = raw_isoforms[raw_id]["fusion_name"]
                if not keep_multicontig_blocks and seqname != expected_contig:
                    continue
                blocks_by_raw_id[raw_id].append(block)

    return blocks_by_raw_id


def filter_mixed_strand_reads(blocks_by_raw_id, allow_mixed_strand_reads):
    if allow_mixed_strand_reads:
        return blocks_by_raw_id, set()

    filtered = defaultdict(list)
    excluded = set()
    for raw_id, blocks in blocks_by_raw_id.items():
        strands = {block["strand"] for block in blocks if block["strand"] in {"+", "-"}}
        if len(strands) > 1:
            excluded.add(raw_id)
            continue
        filtered[raw_id] = blocks
    return filtered, excluded


def sorted_blocks_for_gtf(blocks):
    return sorted(blocks, key=lambda block: (block["seqname"], block["start"], block["end"]))


def block_chain(blocks):
    blocks_sorted = sorted_blocks_for_gtf(blocks)
    return ";".join(
        f'{b["seqname"]}:{b["start"]}-{b["end"]}:{b["strand"]}' for b in blocks_sorted
    )


def contig_and_strand(blocks):
    seqnames = sorted({b["seqname"] for b in blocks})
    strands = sorted({b["strand"] for b in blocks})
    seqname = seqnames[0] if len(seqnames) == 1 else ",".join(seqnames)
    strand = strands[0] if len(strands) == 1 else "."
    return seqname, strand


def internal_junction_chain(blocks):
    """Return the internal splice-junction chain used for NR isoform collapse.

    Junctions are represented by adjacent fusion-contig block boundaries.
    Terminal block starts/ends are intentionally excluded, so reads differing
    only at the 5' or 3' ends collapse together.
    """
    blocks_sorted = sorted_blocks_for_gtf(blocks)
    if len(blocks_sorted) <= 1:
        return "NO_INTERNAL_JUNCTION"

    junctions = []
    for left, right in zip(blocks_sorted, blocks_sorted[1:]):
        junctions.append(
            f'{left["seqname"]}:{left["end"]}->{right["seqname"]}:{right["start"]}:{left["strand"]}/{right["strand"]}'
        )
    return ";".join(junctions)


def parse_breakpoint(value):
    parts = str(value).rsplit(":", 2)
    if len(parts) != 3:
        return (str(value), None, ".")
    chrom, pos, strand = parts
    try:
        return (chrom, int(pos), strand)
    except ValueError:
        return (str(value), None, ".")


def compatible_breakpoint(left, right, tolerance):
    left_chrom, left_pos, left_strand = parse_breakpoint(left)
    right_chrom, right_pos, right_strand = parse_breakpoint(right)
    if left_pos is None or right_pos is None:
        return left == right
    return (
        left_chrom == right_chrom
        and left_strand == right_strand
        and abs(left_pos - right_pos) <= tolerance
    )


def internal_junction_signature(blocks, ignore_alignment_strand=False):
    blocks_sorted = sorted_blocks_for_gtf(blocks)
    signature = []
    for left, right in zip(blocks_sorted, blocks_sorted[1:]):
        left_strand = "." if ignore_alignment_strand else left["strand"]
        right_strand = "." if ignore_alignment_strand else right["strand"]
        signature.append(
            (
                left["seqname"],
                left["end"],
                right["seqname"],
                right["start"],
                left_strand,
                right_strand,
            )
        )
    return signature


def compatible_junction_pair(left, right, tolerance):
    if left[0] != right[0] or left[2] != right[2] or left[4] != right[4] or left[5] != right[5]:
        return False
    return abs(left[1] - right[1]) <= tolerance and abs(left[3] - right[3]) <= tolerance


def compatible_junction_signature(left_sig, right_sig, tolerance, collapse_terminal_subset):
    if len(left_sig) == len(right_sig):
        return all(
            compatible_junction_pair(left, right, tolerance)
            for left, right in zip(left_sig, right_sig)
        )

    if not collapse_terminal_subset:
        return False

    shorter, longer = (left_sig, right_sig) if len(left_sig) < len(right_sig) else (right_sig, left_sig)
    if not shorter:
        return True

    for offset in range(len(longer) - len(shorter) + 1):
        window = longer[offset : offset + len(shorter)]
        if all(
            compatible_junction_pair(short_item, long_item, tolerance)
            for short_item, long_item in zip(shorter, window)
        ):
            return True
    return False


def cluster_signature(raw_id, meta, blocks, ignore_alignment_strand=False):
    fusion_contig, strand = contig_and_strand(blocks)
    if ignore_alignment_strand:
        strand = "."
    return {
        "raw_id": raw_id,
        "fusion_name": meta["fusion_name"],
        "left_breakpoint": meta["left_breakpoint"],
        "right_breakpoint": meta["right_breakpoint"],
        "splice_type": meta["splice_type"],
        "fusion_contig": fusion_contig,
        "strand": strand,
        "junction_count": max(len(blocks) - 1, 0),
        "junction_signature": internal_junction_signature(blocks, ignore_alignment_strand),
    }


def compatible_cluster_signature(left, right, junction_tolerance, breakpoint_tolerance, collapse_terminal_subset):
    if (
        left["fusion_name"] != right["fusion_name"]
        or left["splice_type"] != right["splice_type"]
        or left["fusion_contig"] != right["fusion_contig"]
        or left["strand"] != right["strand"]
    ):
        return False
    if not compatible_breakpoint(left["left_breakpoint"], right["left_breakpoint"], breakpoint_tolerance):
        return False
    if not compatible_breakpoint(left["right_breakpoint"], right["right_breakpoint"], breakpoint_tolerance):
        return False
    return compatible_junction_signature(
        left["junction_signature"],
        right["junction_signature"],
        junction_tolerance,
        collapse_terminal_subset,
    )


def write_raw_fasta(path, raw_isoforms, seq_by_raw_id):
    with open(path, "w") as out:
        for raw_id in sorted(raw_isoforms):
            seq = seq_by_raw_id.get(raw_id)
            if not seq:
                continue
            meta = raw_isoforms[raw_id]
            header = (
                f"{raw_id} read_id={meta['read_id']} fusion={meta['fusion_name']} "
                f"left_breakpoint={meta['left_breakpoint']} right_breakpoint={meta['right_breakpoint']} "
                f"splice_type={meta['splice_type']}"
            )
            out.write(f">{header}\n")
            for line in wrap(seq):
                out.write(line + "\n")


def write_gtf_record(out, seqname, source, feature, start, end, strand, attrs):
    out.write(
        "\t".join(
            [seqname, source, feature, str(start), str(end), ".", strand, ".", gtf_attrs(attrs)]
        )
        + "\n"
    )


def write_raw_gtf(path, raw_isoforms, blocks_by_raw_id):
    with open(path, "w") as out:
        for raw_id in sorted(blocks_by_raw_id):
            blocks = blocks_by_raw_id[raw_id]
            if not blocks:
                continue
            meta = raw_isoforms[raw_id]
            seqnames = sorted({b["seqname"] for b in blocks})
            strands = sorted({b["strand"] for b in blocks})
            seqname = seqnames[0] if len(seqnames) == 1 else ",".join(seqnames)
            strand = strands[0] if len(strands) == 1 else "."
            tx_start = min(b["start"] for b in blocks)
            tx_end = max(b["end"] for b in blocks)

            attrs = [
                ("gene_id", meta["fusion_name"]),
                ("transcript_id", raw_id),
                ("fusion_name", meta["fusion_name"]),
                ("event_id", meta["event_id"]),
                ("read_id", meta["read_id"]),
                ("left_gene", meta["left_gene"]),
                ("right_gene", meta["right_gene"]),
                ("left_breakpoint", meta["left_breakpoint"]),
                ("right_breakpoint", meta["right_breakpoint"]),
                ("splice_type", meta["splice_type"]),
            ]
            write_gtf_record(out, seqname, "CTAT_LR_Fusion", "transcript", tx_start, tx_end, strand, attrs)

            blocks_sorted = sorted_blocks_for_gtf(blocks)
            exon_numbers = list(range(1, len(blocks_sorted) + 1))
            if strand == "-":
                exon_numbers = list(reversed(exon_numbers))

            for block, exon_number in zip(blocks_sorted, exon_numbers):
                exon_attrs = attrs + [
                    ("exon_number", exon_number),
                    ("target_start", block["target_start"]),
                    ("target_end", block["target_end"]),
                ]
                write_gtf_record(
                    out,
                    block["seqname"],
                    "CTAT_LR_Fusion",
                    "exon",
                    block["start"],
                    block["end"],
                    block["strand"],
                    exon_attrs,
                )


def cluster_raw_isoforms(
    raw_isoforms,
    blocks_by_raw_id,
    junction_tolerance,
    breakpoint_tolerance,
    collapse_terminal_subset,
    ignore_alignment_strand,
):
    signatures = {}
    buckets = defaultdict(list)

    for raw_id, meta in raw_isoforms.items():
        blocks = blocks_by_raw_id.get(raw_id, [])
        if not blocks:
            continue
        sig = cluster_signature(raw_id, meta, blocks, ignore_alignment_strand)
        signatures[raw_id] = sig
        bucket_key = (
            sig["fusion_name"],
            sig["splice_type"],
            sig["fusion_contig"],
            sig["strand"],
        )
        buckets[bucket_key].append(raw_id)

    clusters = []
    for raw_ids in buckets.values():
        raw_ids = sorted(raw_ids)
        parent = {raw_id: raw_id for raw_id in raw_ids}

        def find(raw_id):
            while parent[raw_id] != raw_id:
                parent[raw_id] = parent[parent[raw_id]]
                raw_id = parent[raw_id]
            return raw_id

        def union(left, right):
            left_root = find(left)
            right_root = find(right)
            if left_root != right_root:
                parent[right_root] = left_root

        for i, left_id in enumerate(raw_ids):
            left_sig = signatures[left_id]
            for right_id in raw_ids[i + 1 :]:
                right_sig = signatures[right_id]
                if compatible_cluster_signature(
                    left_sig,
                    right_sig,
                    junction_tolerance,
                    breakpoint_tolerance,
                    collapse_terminal_subset,
                ):
                    union(left_id, right_id)

        components = defaultdict(list)
        for raw_id in raw_ids:
            components[find(raw_id)].append(raw_id)
        clusters.extend(sorted(sorted(component) for component in components.values()))

    return sorted(clusters, key=lambda component: component[0])


def make_nr_ids(clusters, raw_isoforms, id_prefix):
    nr_records = []
    fusion_counts = defaultdict(int)

    def cluster_sort_key(raw_ids):
        meta = raw_isoforms[raw_ids[0]]
        return (
            meta["fusion_name"],
            meta["left_breakpoint"],
            meta["right_breakpoint"],
            meta["splice_type"],
            raw_ids[0],
        )

    for raw_ids in sorted(clusters, key=cluster_sort_key):
        raw_ids = sorted(raw_ids)
        fusion_name = raw_isoforms[raw_ids[0]]["fusion_name"]
        fusion_counts[fusion_name] += 1
        safe_fusion = fusion_name.replace("--", "_").replace("/", "_").replace("|", "_")
        nr_id = f"{id_prefix}_{safe_fusion}_isoform{fusion_counts[fusion_name]:04d}"
        nr_records.append((nr_id, raw_ids))

    return nr_records


def choose_representative(raw_ids, seq_by_raw_id):
    found = [raw_id for raw_id in raw_ids if raw_id in seq_by_raw_id]
    if not found:
        return raw_ids[0], ""
    representative = max(found, key=lambda raw_id: (len(seq_by_raw_id[raw_id]), raw_id))
    return representative, seq_by_raw_id[representative]


def write_nr_outputs(
    out_fa,
    out_gtf,
    out_support,
    out_meta,
    nr_records,
    raw_isoforms,
    blocks_by_raw_id,
    seq_by_raw_id,
    min_support,
    junction_tolerance,
    breakpoint_tolerance,
    collapse_terminal_subset,
):
    kept = 0
    with open(out_fa, "w") as fa, open(out_gtf, "w") as gtf, open(out_support, "w") as support, open(out_meta, "w") as meta_out:
        support.write(
            "\t".join(["nr_isoform_id", "support_read_count", "representative_raw_id", "support_raw_ids", "support_read_ids"])
            + "\n"
        )
        meta_out.write(
            "\t".join(
                [
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
                    "representative_raw_id",
                    "representative_read_id",
                    "representative_sequence_length",
                    "support_read_count",
                    "clustering_rule",
                    "junction_tolerance_bp",
                    "breakpoint_tolerance_bp",
                    "collapse_terminal_subset",
                    "internal_junction_chain",
                    "representative_block_chain",
                ]
            )
            + "\n"
        )

        for nr_id, raw_ids in nr_records:
            if len(raw_ids) < min_support:
                continue

            kept += 1
            rep_raw_id, rep_seq = choose_representative(raw_ids, seq_by_raw_id)
            rep_meta = raw_isoforms[rep_raw_id]
            rep_blocks = blocks_by_raw_id[rep_raw_id]
            clustering_chain = internal_junction_chain(rep_blocks)
            blocks_sorted = sorted_blocks_for_gtf(rep_blocks)
            seqname, strand = contig_and_strand(rep_blocks)
            tx_start = min(b["start"] for b in rep_blocks)
            tx_end = max(b["end"] for b in rep_blocks)
            read_ids = [raw_isoforms[raw_id]["read_id"] for raw_id in raw_ids]

            if rep_seq:
                header = (
                    f"{nr_id} fusion={rep_meta['fusion_name']} support={len(raw_ids)} "
                    f"representative_read={rep_meta['read_id']} left_breakpoint={rep_meta['left_breakpoint']} "
                    f"right_breakpoint={rep_meta['right_breakpoint']} splice_type={rep_meta['splice_type']}"
                )
                fa.write(f">{header}\n")
                for line in wrap(rep_seq):
                    fa.write(line + "\n")

            attrs = [
                ("gene_id", rep_meta["fusion_name"]),
                ("transcript_id", nr_id),
                ("fusion_name", rep_meta["fusion_name"]),
                ("left_gene", rep_meta["left_gene"]),
                ("right_gene", rep_meta["right_gene"]),
                ("left_breakpoint", rep_meta["left_breakpoint"]),
                ("right_breakpoint", rep_meta["right_breakpoint"]),
                ("splice_type", rep_meta["splice_type"]),
                ("clustering_rule", "same_fusion_splice_type_strand_contig_fuzzy_breakpoint_and_internal_junction_chain"),
                ("junction_tolerance_bp", junction_tolerance),
                ("breakpoint_tolerance_bp", breakpoint_tolerance),
                ("collapse_terminal_subset", collapse_terminal_subset),
                ("internal_junction_chain", clustering_chain),
                ("support_read_count", len(raw_ids)),
                ("representative_raw_id", rep_raw_id),
                ("support_raw_ids", ",".join(raw_ids)),
            ]
            write_gtf_record(gtf, seqname, "CTAT_LR_Fusion_NR", "transcript", tx_start, tx_end, strand, attrs)

            exon_numbers = list(range(1, len(blocks_sorted) + 1))
            if strand == "-":
                exon_numbers = list(reversed(exon_numbers))
            for block, exon_number in zip(blocks_sorted, exon_numbers):
                exon_attrs = attrs + [("exon_number", exon_number)]
                write_gtf_record(
                    gtf,
                    block["seqname"],
                    "CTAT_LR_Fusion_NR",
                    "exon",
                    block["start"],
                    block["end"],
                    block["strand"],
                    exon_attrs,
                )

            support.write(
                "\t".join([nr_id, str(len(raw_ids)), rep_raw_id, ",".join(raw_ids), ",".join(read_ids)])
                + "\n"
            )
            meta_out.write(
                "\t".join(
                    [
                        nr_id,
                        rep_meta["fusion_name"],
                        rep_meta["left_gene"],
                        rep_meta["right_gene"],
                        rep_meta["left_breakpoint"],
                        rep_meta["right_breakpoint"],
                        rep_meta["splice_type"],
                        seqname,
                        strand,
                        str(len(rep_blocks)),
                        rep_raw_id,
                        rep_meta["read_id"],
                        str(len(rep_seq)),
                        str(len(raw_ids)),
                        "same_fusion_splice_type_strand_contig_fuzzy_breakpoint_and_internal_junction_chain",
                        str(junction_tolerance),
                        str(breakpoint_tolerance),
                        str(collapse_terminal_subset),
                        clustering_chain,
                        block_chain(rep_blocks),
                    ]
                )
                + "\n"
            )

    return kept


def write_raw_manifest(path, raw_isoforms, seq_by_raw_id, blocks_by_raw_id, mixed_strand_excluded):
    with open(path, "w") as out:
        out.write(
            "\t".join(
                [
                    "raw_id",
                    "read_id",
                    "fusion_name",
                    "event_id",
                    "left_gene",
                    "right_gene",
                    "left_breakpoint",
                    "right_breakpoint",
                    "splice_type",
                    "in_read_fasta",
                    "in_fusion_gff3",
                    "mixed_strand_excluded",
                    "sequence_length",
                    "exon_block_count",
                    "internal_junction_chain",
                    "block_chain",
                ]
            )
            + "\n"
        )
        for raw_id in sorted(raw_isoforms):
            meta = raw_isoforms[raw_id]
            seq = seq_by_raw_id.get(raw_id, "")
            blocks = blocks_by_raw_id.get(raw_id, [])
            out.write(
                "\t".join(
                    [
                        raw_id,
                        meta["read_id"],
                        meta["fusion_name"],
                        meta["event_id"],
                        meta["left_gene"],
                        meta["right_gene"],
                        meta["left_breakpoint"],
                        meta["right_breakpoint"],
                        meta["splice_type"],
                        "yes" if seq else "no",
                        "yes" if blocks else "no",
                        "yes" if raw_id in mixed_strand_excluded else "no",
                        str(len(seq)),
                        str(len(blocks)),
                        internal_junction_chain(blocks) if blocks else ".",
                        block_chain(blocks) if blocks else ".",
                    ]
                )
                + "\n"
            )


def main():
    set_csv_field_size_limit()

    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    raw_fa = outdir / "fusion_isoforms.raw_reads.fa"
    raw_gtf = outdir / "fusion_isoforms.raw_blocks.gtf"
    raw_manifest = outdir / "fusion_isoforms.raw_manifest.tsv"
    nr_fa = outdir / "fusion_isoforms.nr.fa"
    nr_gtf = outdir / "fusion_isoforms.nr.gtf"
    nr_support = outdir / "fusion_isoform_read_support.tsv"
    nr_meta = outdir / "fusion_isoforms.nr_metadata.tsv"

    raw_isoforms, reads_to_raw_ids = load_fusion_events(args.fusion_tsv)
    seq_by_raw_id = load_read_sequences(args.read_fasta, reads_to_raw_ids)
    blocks_by_raw_id = load_gff3_blocks(
        args.fusion_gff3,
        reads_to_raw_ids,
        raw_isoforms,
        args.keep_multicontig_blocks,
    )
    blocks_by_raw_id, mixed_strand_excluded = filter_mixed_strand_reads(
        blocks_by_raw_id,
        args.allow_mixed_strand_reads,
    )

    write_raw_fasta(raw_fa, raw_isoforms, seq_by_raw_id)
    write_raw_gtf(raw_gtf, raw_isoforms, blocks_by_raw_id)
    write_raw_manifest(raw_manifest, raw_isoforms, seq_by_raw_id, blocks_by_raw_id, mixed_strand_excluded)

    clusters = cluster_raw_isoforms(
        raw_isoforms,
        blocks_by_raw_id,
        args.junction_tolerance,
        args.breakpoint_tolerance,
        not args.no_collapse_terminal_subset,
        args.ignore_alignment_strand_for_clustering,
    )
    nr_records = make_nr_ids(clusters, raw_isoforms, args.id_prefix)
    kept_nr = write_nr_outputs(
        nr_fa,
        nr_gtf,
        nr_support,
        nr_meta,
        nr_records,
        raw_isoforms,
        blocks_by_raw_id,
        seq_by_raw_id,
        args.min_nr_support,
        args.junction_tolerance,
        args.breakpoint_tolerance,
        not args.no_collapse_terminal_subset,
    )

    print(f"Raw LR_accession isoforms: {len(raw_isoforms)}")
    print(f"Raw isoforms found in read FASTA: {len(seq_by_raw_id)}")
    print(f"Raw isoforms found in fusion GFF3: {len(blocks_by_raw_id)}")
    print(f"Mixed-strand raw reads excluded: {len(mixed_strand_excluded)}")
    print(f"NR isoform clusters before min-support filter: {len(nr_records)}")
    print(f"NR isoform clusters written: {kept_nr}")
    print(f"Junction tolerance: {args.junction_tolerance} bp")
    print(f"Breakpoint tolerance: {args.breakpoint_tolerance} bp")
    print(f"Collapse terminal junction subset: {not args.no_collapse_terminal_subset}")
    print(f"Keep multi-contig GFF3 blocks: {args.keep_multicontig_blocks}")
    print(f"Ignore alignment strand for clustering: {args.ignore_alignment_strand_for_clustering}")
    print(f"Allow mixed-strand reads: {args.allow_mixed_strand_reads}")
    print(f"Raw FASTA: {raw_fa}")
    print(f"Raw GTF: {raw_gtf}")
    print(f"Raw manifest: {raw_manifest}")
    print(f"NR FASTA: {nr_fa}")
    print(f"NR GTF: {nr_gtf}")
    print(f"NR read support: {nr_support}")
    print(f"NR metadata: {nr_meta}")


if __name__ == "__main__":
    main()

