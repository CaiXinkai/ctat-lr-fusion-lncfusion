#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'USAGE'
Align multiple long-read cDNA FASTQ files to the LncFusion reference and run
the CTAT-LR-Fusion isoform workflow.

Required:
  --samples FILE          Tab-separated: sample_id, condition, fastq
  --lncfusion-lib DIR     LncFusion Lib root or CTAT genome library directory
  --outdir DIR            Analysis output directory

Optional:
  --threads N             Default: 32
  --scripts-dir DIR       Default: repository scripts directory
  --min-unique-bases N    Default: 20

This default uses `minimap2 -ax splice`, suitable for ordinary long-read cDNA.
For ONT direct RNA or PacBio Iso-Seq, review and modify the minimap2 command.
USAGE
}

die() { echo "ERROR: $*" >&2; exit 1; }

samples=""
lncfusion_lib=""
outdir=""
threads=32
min_unique_bases=20
repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
scripts_dir="$repo_root/scripts"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --samples) samples="${2:-}"; shift 2 ;;
        --lncfusion-lib) lncfusion_lib="${2:-}"; shift 2 ;;
        --outdir) outdir="${2:-}"; shift 2 ;;
        --threads) threads="${2:-}"; shift 2 ;;
        --scripts-dir) scripts_dir="${2:-}"; shift 2 ;;
        --min-unique-bases) min_unique_bases="${2:-}"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) die "Unknown argument: $1" ;;
    esac
done

[[ -s "$samples" ]] || die "Sample table not found: $samples"
[[ -n "$lncfusion_lib" ]] || die "--lncfusion-lib is required"
[[ -n "$outdir" ]] || die "--outdir is required"
command -v minimap2 >/dev/null || die "minimap2 not found"
command -v samtools >/dev/null || die "samtools not found"

if [[ -d "$lncfusion_lib/LncFusion_v1_CTAT/ctat_genome_lib_build_dir" ]]; then
    ctat_lib="$lncfusion_lib/LncFusion_v1_CTAT/ctat_genome_lib_build_dir"
elif [[ -d "$lncfusion_lib/ctat_genome_lib_build_dir" ]]; then
    ctat_lib="$lncfusion_lib/ctat_genome_lib_build_dir"
else
    ctat_lib="$lncfusion_lib"
fi

ref="$ctat_lib/ref_genome.fa"
[[ -s "$ref" ]] || die "Missing reference FASTA: $ref"
mkdir -p "$outdir/bam" "$outdir/logs"

tail -n +2 "$samples" |
while IFS=$'\t' read -r sample_id condition fastq extra; do
    [[ -n "$sample_id" ]] || continue
    [[ "$sample_id" =~ ^[A-Za-z0-9._-]+$ ]] || die "Unsafe sample ID: $sample_id"
    [[ -s "$fastq" ]] || die "FASTQ not found: $fastq"
    bam="$outdir/bam/${sample_id}_aln.bam"
    log="$outdir/logs/${sample_id}.minimap2.log"

    if [[ ! -s "$bam" ]]; then
        echo "Aligning $sample_id ($condition)" >&2
        minimap2 -ax splice -t "$threads" "$ref" "$fastq" 2> "$log" |
            samtools sort -@ "$threads" -o "$bam" -
        samtools index -@ "$threads" "$bam"
    fi

    bash "$scripts_dir/run_ctat_lr_fusion_with_lncfusion_database.sh" \
        --bam "$bam" \
        --outdir "$outdir/results" \
        --lncfusion-lib "$ctat_lib" \
        --scripts-dir "$scripts_dir" \
        --threads "$threads" \
        --min-unique-bases "$min_unique_bases" \
        --keep-ctat-workdir
done

