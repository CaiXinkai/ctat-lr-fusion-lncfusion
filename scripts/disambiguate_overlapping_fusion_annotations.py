#!/usr/bin/env python3
"""Resolve overlapping fusion-gene annotations using genomic long-read evidence.

This does not redefine fusion RNAs by gene name. It writes one row per genomic
breakpoint/splice event, keeps every compatible annotation, and selects a gene
only when supporting reads cover exon sequence or splice junctions unique to it.
"""

import argparse
import csv
import re
import subprocess
from collections import defaultdict


def args():
    p = argparse.ArgumentParser()
    p.add_argument("--raw-manifest", required=True)
    p.add_argument("--nr-support", required=True)
    p.add_argument("--bam", required=True)
    p.add_argument("--gtf", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--min-unique-bases", type=int, default=20)
    p.add_argument("--junction-tolerance", type=int, default=5)
    return p.parse_args()


def attrs(text):
    return dict(re.findall(r'(\S+)\s+"([^"]*)"', text))


def merge_intervals(items):
    out = []
    for start, end in sorted(items):
        if not out or start > out[-1][1]:
            out.append([start, end])
        else:
            out[-1][1] = max(out[-1][1], end)
    return [tuple(x) for x in out]


def overlap_len(blocks, intervals):
    total = 0
    for a, b in blocks:
        for c, d in intervals:
            if d <= a:
                continue
            if c >= b:
                break
            total += max(0, min(b, d) - max(a, c))
    return total


def subtract_intervals(source, masks):
    result = []
    masks = merge_intervals(masks)
    for start, end in source:
        pieces = [(start, end)]
        for mstart, mend in masks:
            next_pieces = []
            for a, b in pieces:
                if mend <= a or mstart >= b:
                    next_pieces.append((a, b))
                else:
                    if a < mstart:
                        next_pieces.append((a, mstart))
                    if mend < b:
                        next_pieces.append((mend, b))
            pieces = next_pieces
        result.extend(pieces)
    return merge_intervals(result)


def parse_bp(value):
    chrom, pos, strand = value.rsplit(":", 2)
    return chrom, int(pos), strand


def load_retained_raw_ids(path):
    retained = set()
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            retained.update(x for x in row["support_raw_ids"].split(",") if x)
    return retained


def load_manifest(path, retained_raw_ids):
    groups = defaultdict(list)
    reads = set()
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["raw_id"] not in retained_raw_ids:
                continue
            if row.get("in_fusion_gff3") != "yes" or row.get("mixed_strand_excluded") == "yes":
                continue
            key = (row["left_breakpoint"], row["right_breakpoint"], row["splice_type"])
            groups[key].append(row)
            reads.add(row["read_id"])
    return groups, reads


def load_gtf(path, wanted):
    exons = defaultdict(list)
    transcript_exons = defaultdict(list)
    with open(path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) != 9 or f[2] != "exon":
                continue
            at = attrs(f[8])
            names = {at.get("gene_name"), at.get("gene_id")}
            hit = names & wanted
            if not hit:
                continue
            # Convert 1-based inclusive GTF coordinates to 0-based half-open.
            for name in hit:
                exons[(name, f[0], f[6])].append((int(f[3]) - 1, int(f[4])))
                transcript_exons[(name, f[0], f[6], at.get("transcript_id", "."))].append(
                    (int(f[3]) - 1, int(f[4]))
                )
    junctions = defaultdict(set)
    for (name, chrom, strand, _), tx_exons in transcript_exons.items():
        ordered = sorted(tx_exons)
        junctions[(name, chrom, strand)].update(
            (ordered[i][1], ordered[i + 1][0]) for i in range(len(ordered) - 1)
        )
    return ({key: merge_intervals(value) for key, value in exons.items()}, junctions)


def cigar_blocks(pos0, cigar):
    blocks, junctions = [], []
    ref = pos0
    block_start = ref
    for number, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        n = int(number)
        if op in "M=X":
            ref += n
        elif op == "D":
            ref += n
        elif op == "N":
            if ref > block_start:
                blocks.append((block_start, ref))
            junctions.append((ref, ref + n))
            ref += n
            block_start = ref
    if ref > block_start:
        blocks.append((block_start, ref))
    return blocks, junctions


def load_alignments(bam, wanted_reads):
    alignments = defaultdict(list)
    proc = subprocess.Popen(
        ["samtools", "view", bam], text=True, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE, bufsize=1024 * 1024,
    )
    for line in proc.stdout:
        f = line.rstrip("\n").split("\t")
        if len(f) < 6 or f[0] not in wanted_reads or f[2] == "*" or f[5] == "*":
            continue
        flag = int(f[1])
        if flag & 0x100 or flag & 0x400:
            continue
        blocks, junctions = cigar_blocks(int(f[3]) - 1, f[5])
        alignments[f[0]].append((f[2], "-" if flag & 16 else "+", blocks, junctions))
    stderr = proc.stderr.read()
    if proc.wait() != 0:
        raise SystemExit(f"samtools view failed: {stderr}")
    return alignments


def gene_score(gene, chrom, strand, alternatives, reads, exons, annotated_junctions, alignments, tolerance):
    gene_exons = exons.get((gene, chrom, strand), [])
    other_exons = []
    for other in alternatives:
        if other != gene:
            other_exons.extend(exons.get((other, chrom, strand), []))
    unique_exons = subtract_intervals(gene_exons, other_exons)

    gene_junctions = annotated_junctions.get((gene, chrom, strand), set())
    other_junctions = set()
    for other in alternatives:
        if other == gene:
            continue
        other_junctions.update(annotated_junctions.get((other, chrom, strand), set()))
    unique_junctions = gene_junctions - other_junctions

    unique_bases = total_bases = junction_hits = 0
    evidence_reads = set()
    for read in reads:
        for aln_chrom, aln_strand, blocks, junctions in alignments.get(read, []):
            if aln_chrom != chrom or aln_strand != strand:
                continue
            ub = overlap_len(blocks, unique_exons)
            tb = overlap_len(blocks, gene_exons)
            jh = sum(
                any(abs(a - c) <= tolerance and abs(b - d) <= tolerance for c, d in unique_junctions)
                for a, b in junctions
            )
            unique_bases += ub
            total_bases += tb
            junction_hits += jh
            if ub or jh:
                evidence_reads.add(read)
    return unique_bases, junction_hits, len(evidence_reads), total_bases


def resolve_side(genes, bp, reads, exons, annotated_junctions, alignments, min_bases, tolerance):
    chrom, _, strand = parse_bp(bp)
    genes = sorted(set(genes))
    scores = {
        gene: gene_score(gene, chrom, strand, genes, reads, exons, annotated_junctions, alignments, tolerance)
        for gene in genes
    }
    # A single annotation is not an overlap ambiguity. Keep it as the primary
    # label even when the retained fusion fragment does not cover an exon that
    # is unique relative to another (non-existent) candidate.
    if len(genes) == 1:
        return genes[0], "single_candidate_no_disambiguation_needed", scores
    informative = [g for g, s in scores.items() if s[0] >= min_bases or s[1] > 0]
    if not informative:
        return ".", "ambiguous_no_unique_read_evidence", scores
    ranked = sorted(informative, key=lambda g: (scores[g][1], scores[g][0], scores[g][2], scores[g][3]), reverse=True)
    best = ranked[0]
    if len(ranked) > 1 and scores[ranked[0]][:2] == scores[ranked[1]][:2]:
        return ".", "ambiguous_tied_unique_evidence", scores
    return best, "resolved_by_unique_exon_or_junction", scores


def format_scores(scores):
    return ";".join(
        f"{g}:unique_bp={s[0]},unique_junc={s[1]},evidence_reads={s[2]},total_exonic_bp={s[3]}"
        for g, s in sorted(scores.items())
    )


def main():
    a = args()
    retained_raw_ids = load_retained_raw_ids(a.nr_support)
    groups, wanted_reads = load_manifest(a.raw_manifest, retained_raw_ids)
    wanted_genes = {r[side] for rows in groups.values() for r in rows for side in ("left_gene", "right_gene")}
    exons, annotated_junctions = load_gtf(a.gtf, wanted_genes)
    alignments = load_alignments(a.bam, wanted_reads)
    fields = [
        "fusion_rna_event_id", "left_breakpoint", "right_breakpoint", "splice_type",
        "unique_supporting_reads", "left_candidate_genes", "right_candidate_genes",
        "resolved_left_gene", "resolved_right_gene", "left_resolution", "right_resolution",
        "primary_annotation", "all_compatible_annotations", "left_evidence_scores", "right_evidence_scores",
    ]
    with open(a.output, "w", newline="") as out:
        writer = csv.DictWriter(out, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for number, (key, rows) in enumerate(sorted(groups.items()), 1):
            left = sorted({r["left_gene"] for r in rows}); right = sorted({r["right_gene"] for r in rows})
            read_ids = {r["read_id"] for r in rows}
            lg, ls, lscore = resolve_side(left, key[0], read_ids, exons, annotated_junctions, alignments, a.min_unique_bases, a.junction_tolerance)
            rg, rs, rscore = resolve_side(right, key[1], read_ids, exons, annotated_junctions, alignments, a.min_unique_bases, a.junction_tolerance)
            primary = f"{lg}--{rg}" if lg != "." and rg != "." else "."
            writer.writerow({
                "fusion_rna_event_id": f"FRNA{number:07d}", "left_breakpoint": key[0],
                "right_breakpoint": key[1], "splice_type": key[2],
                "unique_supporting_reads": len(read_ids), "left_candidate_genes": ",".join(left),
                "right_candidate_genes": ",".join(right), "resolved_left_gene": lg,
                "resolved_right_gene": rg, "left_resolution": ls, "right_resolution": rs,
                "primary_annotation": primary,
                "all_compatible_annotations": ";".join(f"{x}--{y}" for x in left for y in right),
                "left_evidence_scores": format_scores(lscore), "right_evidence_scores": format_scores(rscore),
            })
    print(f"Wrote {len(groups)} breakpoint-defined fusion RNA events: {a.output}")


if __name__ == "__main__":
    main()

