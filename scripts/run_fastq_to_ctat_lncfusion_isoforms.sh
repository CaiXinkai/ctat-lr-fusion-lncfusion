#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage:
  bash scripts/run_fastq_to_ctat_lncfusion_isoforms.sh \
    --fastq reads.fastq.gz --sample SAMPLE --outdir /path/output \
    --lncfusion-lib /data3/LncFusion/Lib --threads 32 \
    [--platform ont-cdna|ont-directrna|pacbio-isoseq] [--vis] [--resume]

The final deliverables are CTAT's final fusion_predictions.tsv plus reference-
reconstructed isoform sequences. No post-CTAT support, confidence, annotation,
or event filter is applied. The CTAT preliminary table is coordinate evidence
only and can never add a call absent from the final CTAT table.

--resume reuses a completed BAM and/or CTAT core outputs. It is intended for
recovering runs that failed during CTAT's post-result IGV packaging stage.
EOF
}

die() { echo "ERROR: $*" >&2; exit 1; }
need() { command -v "$1" >/dev/null 2>&1 || die "$1 not found in PATH"; }

fastq=""; sample=""; outdir=""; lncfusion_lib="/data3/LncFusion/Lib"
threads=8; platform="ont-cdna"; vis=0; resume=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --fastq) fastq="$2"; shift 2 ;;
        --sample) sample="$2"; shift 2 ;;
        --outdir) outdir="$2"; shift 2 ;;
        --lncfusion-lib) lncfusion_lib="$2"; shift 2 ;;
        --threads) threads="$2"; shift 2 ;;
        --platform) platform="$2"; shift 2 ;;
        --vis) vis=1; shift ;;
        --resume) resume=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) die "Unknown argument: $1" ;;
    esac
done

[[ -s "$fastq" ]] || die "--fastq is missing or empty"
[[ -n "$sample" && "$sample" =~ ^[A-Za-z0-9._-]+$ ]] || die "Invalid --sample"
[[ -n "$outdir" ]] || die "--outdir is required"
[[ "$threads" =~ ^[1-9][0-9]*$ ]] || die "--threads must be a positive integer"
need minimap2; need samtools; need ctat-LR-fusion; need python3; need gzip

if [[ -d "$lncfusion_lib/LncFusion_v1_CTAT/ctat_genome_lib_build_dir" ]]; then
    ctat_lib="$lncfusion_lib/LncFusion_v1_CTAT/ctat_genome_lib_build_dir"
elif [[ -d "$lncfusion_lib/ctat_genome_lib_build_dir" ]]; then
    ctat_lib="$lncfusion_lib/ctat_genome_lib_build_dir"
elif [[ "$(basename "$lncfusion_lib")" == "ctat_genome_lib_build_dir" ]]; then
    ctat_lib="$lncfusion_lib"
else
    die "Cannot resolve LncFusion CTAT library below: $lncfusion_lib"
fi
ctat_lib="$(realpath "$ctat_lib")"
ref="$ctat_lib/ref_genome.fa"; gtf="$ctat_lib/ref_annot.gtf"
[[ -s "$ref" && -s "$gtf" ]] || die "Missing ref_genome.fa or ref_annot.gtf in $ctat_lib"
[[ -s "$ref.fai" ]] || samtools faidx "$ref"

case "$platform" in
    ont-cdna) mm2=(minimap2 -ax splice) ;;
    ont-directrna) mm2=(minimap2 -ax splice -uf -k14) ;;
    pacbio-isoseq) mm2=(minimap2 -ax splice:hq -uf) ;;
    *) die "Unsupported --platform: $platform" ;;
esac
splice_bed="$gtf.mm2.splice.bed"
[[ -s "$splice_bed" ]] && mm2+=(--junc-bed "$splice_bed")

sample_dir="$(realpath -m "$outdir")/$sample"
work="$sample_dir/work"; final="$sample_dir/final"
mkdir -p "$work" "$final"
bam="$work/${sample}.lncfusion_hg38.sorted.bam"
if [[ "$resume" -eq 1 && -s "$bam" && -s "$bam.bai" ]]; then
    echo "[$sample] Reusing existing BAM (--resume): $bam" >&2
else
    echo "[$sample] FASTQ -> coordinate-sorted BAM" >&2
    "${mm2[@]}" -t "$threads" "$ref" "$fastq" 2>"$work/minimap2.log" |
        samtools sort -@ "$threads" -o "$bam" -
    samtools index -@ "$threads" "$bam"
fi

if [[ ! -s "$ctat_lib/ref_genome.fa.mm2" ]]; then
    echo "[$sample] Preparing CTAT minimap2 reference index once" >&2
    ctat-LR-fusion -T "$fastq" --genome_lib_dir "$ctat_lib" \
        --output "$work/reference_prep" --CPU "$threads" --prep_reference_only
fi

ctat_out="$work/ctat_lr_fusion"
final_tsv="$ctat_out/ctat-LR-fusion.fusion_predictions.tsv"
prelim_tsv="$ctat_out/ctat-LR-fusion.fusion_predictions.preliminary.tsv"
fusion_gff3="$ctat_out/fusion_intermediates_dir/LR-FI.mm2.fusion_transcripts.gff3"
ctat_args=(--LR_bam "$bam" --genome_lib_dir "$ctat_lib" --output "$ctat_out" --CPU "$threads" --min_J 1)
[[ "$vis" -eq 1 ]] && ctat_args+=(--vis)

ctat_exit=0
late_igv_failure_rescued=0
if [[ "$resume" -eq 1 && -s "$final_tsv" && -s "$prelim_tsv" && -s "$fusion_gff3" ]]; then
    echo "[$sample] Reusing completed CTAT core outputs (--resume)" >&2
else
    echo "[$sample] Running CTAT-LR-Fusion with the LncFusion CTAT library" >&2
    set +e
    ctat-LR-fusion "${ctat_args[@]}"
    ctat_exit=$?
    set -e
    if [[ "$ctat_exit" -ne 0 ]]; then
        if [[ -s "$final_tsv" && -s "$prelim_tsv" && -s "$fusion_gff3" ]]; then
            late_igv_failure_rescued=1
            echo "WARNING: CTAT exited with status $ctat_exit after writing all core deliverables." >&2
            echo "         Treating this as a post-result IGV packaging failure and continuing." >&2
        else
            die "CTAT failed with status $ctat_exit before all core deliverables were created"
        fi
    fi
fi

[[ -s "$final_tsv" ]] || die "CTAT final table missing: $final_tsv"
[[ -s "$prelim_tsv" ]] || die "CTAT coordinate table missing: $prelim_tsv"
[[ -s "$fusion_gff3" ]] || die "CTAT fusion GFF3 missing: $fusion_gff3"
cp -f "$final_tsv" "$final/${sample}.ctat-LR-fusion.fusion_predictions.tsv"
[[ -s "$ctat_out/ctat-LR-fusion.fusion_predictions.abridged.tsv" ]] && \
    cp -f "$ctat_out/ctat-LR-fusion.fusion_predictions.abridged.tsv" \
        "$final/${sample}.ctat-LR-fusion.fusion_predictions.abridged.tsv"

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
entities="$work/${sample}.ctat_final_isoform_entities.tsv"
python3 "$script_dir/prepare_ctat_isoform_entities.py" \
    --final-tsv "$final_tsv" --preliminary-tsv "$prelim_tsv" \
    --fusion-gff3 "$fusion_gff3" --sample "$sample" --output "$entities" \
    | tee "$work/isoform_entity_summary.json"

if [[ "$(wc -l < "$entities")" -eq 1 ]]; then
    gzip -c /dev/null > "$final/${sample}.fusion_isoforms.reference.fa.gz"
    printf 'feature_id\tstatus\treason\n' | gzip -c > "$final/${sample}.fusion_isoforms.audit.tsv.gz"
    printf '{"ctat_final_rows":0,"isoform_sequences":0}\n' > "$final/${sample}.fusion_isoforms.summary.json"
else
    sequence_dir="$work/sequence_reconstruction"
    python3 "$script_dir/reconstruct_genome_chain_isoforms.py" \
        --global-isoforms "$entities" --reference-fasta "$ref" \
        --annotation-gtf "$gtf" --output-dir "$sequence_dir" --threads "$threads"
    cp -f "$sequence_dir/global_hc_isoform_reference_sequences.fa.gz" \
        "$final/${sample}.fusion_isoforms.reference.fa.gz"
    cp -f "$sequence_dir/global_hc_isoform_reference_sequences.audit.tsv.gz" \
        "$final/${sample}.fusion_isoforms.audit.tsv.gz"
    cp -f "$sequence_dir/global_hc_isoform_reference_sequences.summary.json" \
        "$final/${sample}.fusion_isoforms.summary.json"
fi
cp -f "$entities" "$final/${sample}.fusion_isoforms.entities.tsv"
{
    printf 'ctat_exit_status\t%s\n' "$ctat_exit"
    printf 'late_igv_failure_rescued\t%s\n' "$late_igv_failure_rescued"
    printf 'ctat_final_table\t%s\n' "$final_tsv"
    printf 'ctat_preliminary_coordinate_source\t%s\n' "$prelim_tsv"
    printf 'ctat_fusion_gff3_coordinate_source\t%s\n' "$fusion_gff3"
} > "$final/${sample}.ctat_run_status.tsv"

echo "[$sample] Complete. Final deliverables:" >&2
find "$final" -maxdepth 1 -type f -printf '  %f\n' | sort >&2

