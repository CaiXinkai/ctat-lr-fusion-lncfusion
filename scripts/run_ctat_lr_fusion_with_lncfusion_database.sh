#!/usr/bin/env bash
set -euo pipefail

# Launch the existing CTAT-LR-fusion isoform pipeline with the CTAT database
# distributed by LncFusion (LncBook v2 + GENCODE v42).

usage() {
    cat <<'USAGE'
Usage:
  bash run_ctat_lr_fusion_with_lncfusion_database.sh \
    --bam /path/sample_aln.bam [--bam /path/another_aln.bam] \
    --outdir /path/output \
    [--lncfusion-lib /data3/LncFusion/Lib] \
    [--min-unique-bases 20] \
    [other options accepted by run_ctat_lr_fusion_isoform_database_cleanup.sh]

The value of --lncfusion-lib can be either:
  1. /data3/LncFusion/Lib
  2. /data3/LncFusion/Lib/LncFusion_v1_CTAT
  3. /data3/LncFusion/Lib/LncFusion_v1_CTAT/ctat_genome_lib_build_dir

Default:
  --lncfusion-lib /data3/LncFusion/Lib

Example:
  bash run_ctat_lr_fusion_with_lncfusion_database.sh \
    --bam '/data3/fu_iso/output/*_aln.bam' \
    --outdir /data3/fu_iso/output/CTAT_LR_LncFusion_batch \
    --threads 32 \
    --vis \
    --keep-ctat-workdir
USAGE
}

die() {
    echo "ERROR: $*" >&2
    exit 1
}

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
base_runner="$script_dir/run_ctat_lr_fusion_isoform_database_cleanup.sh"
lncfusion_lib="/data3/LncFusion/Lib"
outdir=""
min_unique_bases=20
forward_args=()
bam_patterns=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --lncfusion-lib)
            [[ $# -ge 2 ]] || die "--lncfusion-lib requires a path"
            lncfusion_lib="$2"
            shift 2
            ;;
        --bam)
            [[ $# -ge 2 ]] || die "--bam requires a path or quoted glob"
            bam_patterns+=("$2")
            forward_args+=("$1" "$2")
            shift 2
            ;;
        --outdir)
            [[ $# -ge 2 ]] || die "--outdir requires a path"
            outdir="$2"
            forward_args+=("$1" "$2")
            shift 2
            ;;
        --min-unique-bases)
            [[ $# -ge 2 ]] || die "--min-unique-bases requires an integer"
            min_unique_bases="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            forward_args+=("$1")
            shift
            ;;
    esac
done

[[ -s "$base_runner" ]] || die "Base pipeline script not found: $base_runner"
[[ "${#bam_patterns[@]}" -gt 0 ]] || die "At least one --bam is required"
[[ -n "$outdir" ]] || die "--outdir is required"
[[ "$min_unique_bases" =~ ^[0-9]+$ ]] || die "--min-unique-bases must be a non-negative integer"
disambiguator="$script_dir/disambiguate_overlapping_fusion_annotations.py"
[[ -s "$disambiguator" ]] || die "Annotation disambiguation script not found: $disambiguator"

# Accept the LncFusion Lib root, LncFusion_v1_CTAT, or the final CTAT directory.
if [[ -d "$lncfusion_lib/LncFusion_v1_CTAT/ctat_genome_lib_build_dir" ]]; then
    ctat_genome_lib="$lncfusion_lib/LncFusion_v1_CTAT/ctat_genome_lib_build_dir"
elif [[ -d "$lncfusion_lib/ctat_genome_lib_build_dir" ]]; then
    ctat_genome_lib="$lncfusion_lib/ctat_genome_lib_build_dir"
elif [[ -d "$lncfusion_lib" && "$(basename "$lncfusion_lib")" == "ctat_genome_lib_build_dir" ]]; then
    ctat_genome_lib="$lncfusion_lib"
else
    die "Cannot locate LncFusion CTAT genome library below: $lncfusion_lib"
fi

ctat_genome_lib="$(realpath "$ctat_genome_lib")"
ref_fasta="$ctat_genome_lib/ref_genome.fa"
ref_gtf="$ctat_genome_lib/ref_annot.gtf"

[[ -s "$ref_fasta" ]] || die "Missing CTAT reference FASTA: $ref_fasta"
[[ -s "$ref_gtf" ]] || die "Missing CTAT reference GTF: $ref_gtf"

command -v samtools >/dev/null 2>&1 || die "samtools not found in PATH"

# Validate the BAM sequence dictionary against the LncFusion reference. This
# catches chr1/1 naming differences and reference builds with different lengths.
tmp_dir="$(mktemp -d)"
trap 'rm -rf -- "$tmp_dir"' EXIT

if [[ -s "$ref_fasta.fai" ]]; then
    cut -f1,2 "$ref_fasta.fai" | sort -k1,1 > "$tmp_dir/reference.contigs.tsv"

    shopt -s nullglob
    expanded_bams=()
    for pattern in "${bam_patterns[@]}"; do
        matches=( $pattern )
        [[ "${#matches[@]}" -gt 0 ]] || die "No BAM matched: $pattern"
        expanded_bams+=("${matches[@]}")
    done
    shopt -u nullglob

    for bam in "${expanded_bams[@]}"; do
        [[ -s "$bam" ]] || die "BAM not found or empty: $bam"
        samtools view -H "$bam" |
            awk -F '\t' '/^@SQ/ {
                name=""; seq_len="";
                for (i=1; i<=NF; i++) {
                    if ($i ~ /^SN:/) name=substr($i,4);
                    if ($i ~ /^LN:/) seq_len=substr($i,4);
                }
                if (name != "" && seq_len != "") print name "\t" seq_len;
            }' |
            sort -k1,1 > "$tmp_dir/bam.contigs.tsv"

        # Different hg38 distributions can contain different ALT/random contigs.
        # That alone is not fatal. Fail only if identically named contigs have
        # conflicting lengths, or if the dictionaries share no contig names.
        awk 'NR==FNR {ref_len[$1]=$2; next}
             ($1 in ref_len) && ref_len[$1] != $2 {
                 print $1 "\tBAM=" $2 "\tLncFusion=" ref_len[$1]
             }' \
            "$tmp_dir/reference.contigs.tsv" "$tmp_dir/bam.contigs.tsv" \
            > "$tmp_dir/conflicting.contigs.tsv"

        shared_count="$(awk 'NR==FNR {ref[$1]=1; next} ($1 in ref) {n++} END {print n+0}' \
            "$tmp_dir/reference.contigs.tsv" "$tmp_dir/bam.contigs.tsv")"

        if [[ "$shared_count" -eq 0 ]]; then
            die "BAM and LncFusion reference share no contig names (for example, 1 versus chr1): $bam"
        fi
        if [[ -s "$tmp_dir/conflicting.contigs.tsv" ]]; then
            echo "Conflicting lengths for identically named BAM/reference contigs: $bam" >&2
            head -20 "$tmp_dir/conflicting.contigs.tsv" >&2
            die "Re-align this sample to the LncFusion ref_genome.fa before using --LR_bam"
        fi

        if ! diff -q "$tmp_dir/reference.contigs.tsv" "$tmp_dir/bam.contigs.tsv" >/dev/null; then
            ref_count="$(wc -l < "$tmp_dir/reference.contigs.tsv")"
            bam_count="$(wc -l < "$tmp_dir/bam.contigs.tsv")"
            echo "WARNING: BAM and LncFusion contain different contig sets, but all $shared_count shared contigs have matching lengths." >&2
            echo "         BAM contigs: $bam_count; LncFusion contigs: $ref_count. Continuing." >&2
        else
            echo "BAM/reference sequence dictionaries match: $bam" >&2
        fi
    done
else
    echo "WARNING: $ref_fasta.fai is absent; BAM/reference dictionary validation was skipped." >&2
fi

echo "LncFusion CTAT library: $ctat_genome_lib" >&2
echo "Annotation: $ref_gtf" >&2
echo "Reference:  $ref_fasta" >&2

bash "$base_runner" \
    --ctat-genome-lib "$ctat_genome_lib" \
    --scripts-dir "$script_dir" \
    "${forward_args[@]}"

# Create one breakpoint-defined event table per sample. Gene labels remain
# annotations; support is counted from unique read IDs, so overlapping labels
# cannot inflate fusion RNA counts during downstream multi-sample integration.
outdir="$(realpath "$outdir")"
shopt -s nullglob
for pattern in "${bam_patterns[@]}"; do
    matches=( $pattern )
    for bam in "${matches[@]}"; do
        base="$(basename "$bam")"
        if [[ "$base" == *_aln.bam ]]; then
            sample="${base%_aln.bam}"
        else
            sample="${base%.bam}"
        fi
        core="$outdir/$sample/core"
        manifest="$core/${sample}.fusion_isoforms.raw_manifest.tsv"
        nr_support="$core/${sample}.fusion_isoform_read_support.tsv"
        disambiguated="$core/${sample}.fusion_rna_events.disambiguated.tsv"
        if [[ ! -s "$manifest" ]]; then
            echo "WARNING: skipping annotation disambiguation; manifest missing: $manifest" >&2
            continue
        fi
        if [[ ! -s "$nr_support" ]]; then
            echo "WARNING: skipping annotation disambiguation; NR support table missing: $nr_support" >&2
            continue
        fi
        echo "Disambiguating overlapping annotations for $sample" >&2
        python "$disambiguator" \
            --raw-manifest "$manifest" \
            --nr-support "$nr_support" \
            --bam "$bam" \
            --gtf "$ref_gtf" \
            --output "$disambiguated" \
            --min-unique-bases "$min_unique_bases" \
            --junction-tolerance 5
    done
done
shopt -u nullglob

