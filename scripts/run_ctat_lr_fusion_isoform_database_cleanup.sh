#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'USAGE'
Batch run CTAT-LR-fusion -> fusion isoform database -> QC from existing BAM files.

Input BAM names should be like:
  SRR32957131_aln.bam

The sample name is inferred by removing the trailing "_aln.bam".
Each sample is written to:
  OUTDIR/SAMPLE/

Core outputs are collected with sample prefixes in:
  OUTDIR/SAMPLE/core/

Required:
  --bam BAM                  Repeat this option for multiple BAMs.
                             A quoted glob is also OK: --bam '/data3/fu_iso/output/*_aln.bam'
  --outdir OUTDIR            Parent output directory.

Common optional:
  --ctat-genome-lib DIR      Default: $CTAT_GENOME_LIB
  --threads N                Default: 32
  --scripts-dir DIR          Default: /data3/fu_iso/output
  --min-nr-support N         Default: 1
  --high-conf-support N      Default: 2
  --id-prefix PREFIX         Default: CTATLRF
  --junction-tolerance N     Default: 5
  --breakpoint-tolerance N   Default: 5
  --min-mapped-reads N       Skip BAMs with fewer mapped reads. Default: 100
  --vis                      Run CTAT-LR-fusion with --vis. Default: off.
                             Note: --vis can fail in CTAT's visualization-only step for some reports.
  --quiet                    Print only pipeline step messages to terminal; keep full command output in logs.
  --rerun-complete           Re-run samples even if core outputs already exist.
  --keep-ctat-workdir        Keep CTAT_LR work directory after collecting core outputs.
  --keep-ctat-intermediates  Keep CTAT fusion_intermediates_dir and checkpoint directories.
  --no-cleanup               Do not delete any generated intermediate files.

Examples:
  bash run_ctat_lr_fusion_isoform_database_cleanup.sh \
    --bam /data3/fu_iso/output/SRR32957131_aln.bam \
    --bam /data3/fu_iso/output/SRR32957132_aln.bam \
    --outdir /data3/fu_iso/output/fusion_isoform_batch \
    --ctat-genome-lib "$CTAT_GENOME_LIB" \
    --threads 32

  bash run_ctat_lr_fusion_isoform_database_cleanup.sh \
    --bam '/data3/fu_iso/output/*_aln.bam' \
    --outdir /data3/fu_iso/output/fusion_isoform_batch \
    --ctat-genome-lib "$CTAT_GENOME_LIB"
USAGE
}

die() {
    echo "ERROR: $*" >&2
    exit 1
}

log() {
    echo "[$(date '+%F %T')] $*" >&2
}

sample_from_bam() {
    local bam_path="$1"
    local base
    base="$(basename "$bam_path")"
    if [[ "$base" == *_aln.bam ]]; then
        echo "${base%_aln.bam}"
    elif [[ "$base" == *.bam ]]; then
        echo "${base%.bam}"
    else
        die "BAM filename must end with .bam: $bam_path"
    fi
}

safe_rm_path() {
    local sample_out="$1"
    local target="$2"
    [[ -e "$target" ]] || return 0

    local resolved_sample_out
    local resolved_target
    resolved_sample_out="$(realpath "$sample_out")"
    resolved_target="$(realpath "$target")"

    case "$resolved_target" in
        "$resolved_sample_out"/*) rm -rf "$resolved_target" ;;
        *) die "Refusing to delete outside sample output directory: $resolved_target" ;;
    esac
}

copy_required() {
    local src="$1"
    local dst="$2"
    [[ -s "$src" ]] || die "Required output missing or empty: $src"
    cp -f "$src" "$dst"
}

copy_optional_glob() {
    local pattern="$1"
    local dest_dir="$2"
    local sample="$3"
    local path
    shopt -s nullglob
    for path in $pattern; do
        cp -f "$path" "$dest_dir/${sample}.$(basename "$path")"
    done
    shopt -u nullglob
}

run_command() {
    local message="$1"
    local logfile="$2"
    shift 2

    log "$message"
    if [[ "$quiet" -eq 1 ]]; then
        "$@" >> "$logfile" 2>&1
    else
        "$@"
    fi
}

record_skip() {
    local table="$1"
    local sample="$2"
    local reason="$3"
    local bam="$4"
    local mapped_reads="$5"
    printf "%s\t%s\t%s\t%s\n" "$sample" "$reason" "$mapped_reads" "$bam" >> "$table"
}

core_outputs_complete() {
    local core_out="$1"
    local sample="$2"

    [[ -s "$core_out/${sample}.fusion_isoforms.raw_reads.fa" ]] || return 1
    [[ -s "$core_out/${sample}.fusion_isoforms.raw_blocks.gtf" ]] || return 1
    [[ -s "$core_out/${sample}.fusion_isoforms.raw_manifest.tsv" ]] || return 1
    [[ -s "$core_out/${sample}.fusion_isoforms.nr.fa" ]] || return 1
    [[ -s "$core_out/${sample}.fusion_isoforms.nr.gtf" ]] || return 1
    [[ -s "$core_out/${sample}.fusion_isoform_read_support.tsv" ]] || return 1
    [[ -s "$core_out/${sample}.fusion_isoforms.nr_metadata.tsv" ]] || return 1
    compgen -G "$core_out/${sample}.ctat-LR-fusion.fusion_predictions*.tsv" > /dev/null || return 1
    return 0
}

bam_patterns=()
outdir=""
ctat_genome_lib="${CTAT_GENOME_LIB:-}"
threads=32
scripts_dir="/data3/fu_iso/output"
min_nr_support=1
high_conf_support=2
id_prefix="CTATLRF"
junction_tolerance=5
breakpoint_tolerance=5
min_mapped_reads=100
run_vis=0
quiet=0
rerun_complete=0
keep_ctat_workdir=0
keep_ctat_intermediates=0
no_cleanup=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --bam) bam_patterns+=("${2:-}"); shift 2 ;;
        --outdir) outdir="${2:-}"; shift 2 ;;
        --ctat-genome-lib) ctat_genome_lib="${2:-}"; shift 2 ;;
        --threads) threads="${2:-}"; shift 2 ;;
        --scripts-dir) scripts_dir="${2:-}"; shift 2 ;;
        --min-nr-support) min_nr_support="${2:-}"; shift 2 ;;
        --high-conf-support) high_conf_support="${2:-}"; shift 2 ;;
        --id-prefix) id_prefix="${2:-}"; shift 2 ;;
        --junction-tolerance) junction_tolerance="${2:-}"; shift 2 ;;
        --breakpoint-tolerance) breakpoint_tolerance="${2:-}"; shift 2 ;;
        --min-mapped-reads) min_mapped_reads="${2:-}"; shift 2 ;;
        --vis) run_vis=1; shift ;;
        --quiet) quiet=1; shift ;;
        --rerun-complete) rerun_complete=1; shift ;;
        --keep-ctat-workdir) keep_ctat_workdir=1; shift ;;
        --keep-ctat-intermediates) keep_ctat_intermediates=1; shift ;;
        --no-cleanup) no_cleanup=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) die "Unknown argument: $1" ;;
    esac
done

[[ "${#bam_patterns[@]}" -gt 0 ]] || die "At least one --bam is required"
[[ -n "$outdir" ]] || die "--outdir is required"
[[ -n "$ctat_genome_lib" ]] || die "--ctat-genome-lib is required, or export CTAT_GENOME_LIB first"
[[ -d "$ctat_genome_lib" ]] || die "CTAT genome library directory not found: $ctat_genome_lib"
[[ -s "$scripts_dir/make_ctat_lr_fusion_isoform_database.py" ]] || die "Missing make script in: $scripts_dir"
[[ -s "$scripts_dir/qc_ctat_lr_fusion_isoform_database.py" ]] || die "Missing QC script in: $scripts_dir"

command -v ctat-LR-fusion >/dev/null 2>&1 || die "ctat-LR-fusion not found in PATH"
command -v python >/dev/null 2>&1 || die "python not found in PATH"
command -v samtools >/dev/null 2>&1 || die "samtools not found in PATH"

mkdir -p "$outdir"
outdir="$(realpath "$outdir")"
cd "$outdir" || die "Cannot cd to output directory: $outdir"
skipped_table="$outdir/skipped_samples.tsv"
printf "sample_id\treason\tmapped_reads\tbam\n" > "$skipped_table"

bams=()
shopt -s nullglob
for pattern in "${bam_patterns[@]}"; do
    matches=( $pattern )
    if [[ "${#matches[@]}" -eq 0 ]]; then
        die "No BAM matched: $pattern"
    fi
    bams+=( "${matches[@]}" )
done
shopt -u nullglob

log "Found ${#bams[@]} BAM file(s)"

for bam in "${bams[@]}"; do
    [[ -s "$bam" ]] || die "BAM not found or empty: $bam"

    input_bam="$(realpath "$bam")"
    sample="$(sample_from_bam "$input_bam")"
    sample_out="$outdir/$sample"
    run_log_dir="$sample_out/logs"
    ctat_out="$sample_out/CTAT_LR"
    ctat_intermediates="$ctat_out/fusion_intermediates_dir"
    db_out="$ctat_out/fusion_isoform_database"
    qc_out="$db_out/qc"
    core_out="$sample_out/core"

    mkdir -p "$run_log_dir" "$core_out"

    sample_log="$run_log_dir/${sample}.run.log"
    {
        log "========== Sample: $sample =========="
        log "BAM: $input_bam"
        log "Sample output: $sample_out"

        if [[ "$rerun_complete" -eq 0 ]] && core_outputs_complete "$core_out" "$sample"; then
            log "Skipping sample because core outputs already exist"
            continue
        fi

        mapped_reads="$(samtools view -c -F 4 "$input_bam")"
        log "Mapped reads: $mapped_reads"
        if [[ "$mapped_reads" -lt "$min_mapped_reads" ]]; then
            log "Skipping sample because mapped reads < $min_mapped_reads"
            record_skip "$skipped_table" "$sample" "mapped_reads_below_${min_mapped_reads}" "$input_bam" "$mapped_reads"
            continue
        fi

        fusion_tsv="$ctat_out/ctat-LR-fusion.fusion_predictions.preliminary.tsv"
        final_fusion_tsv="$ctat_out/ctat-LR-fusion.fusion_predictions.tsv"

        if [[ -s "$final_fusion_tsv" || -s "$fusion_tsv" ]]; then
            log "Skipping CTAT-LR-fusion because fusion prediction TSV already exists"
        else
            ctat_cmd=(
                ctat-LR-fusion
                --LR_bam "$input_bam"
                --genome_lib_dir "$ctat_genome_lib"
                --output "$ctat_out"
                --CPU "$threads"
                --min_J 1
            )
            if [[ "$run_vis" -eq 1 ]]; then
                ctat_cmd+=(--vis)
            fi

            run_command "Running CTAT-LR-fusion" "$sample_log" "${ctat_cmd[@]}"
        fi

        if [[ ! -s "$fusion_tsv" ]]; then
            fusion_tsv="$final_fusion_tsv"
        fi
        if [[ ! -s "$fusion_tsv" ]]; then
            log "Skipping sample because CTAT fusion prediction TSV was not generated"
            record_skip "$skipped_table" "$sample" "missing_ctat_fusion_tsv" "$input_bam" "$mapped_reads"
            continue
        fi

        read_fasta="$ctat_intermediates/chimeric_read_candidates.transcripts.fa.revised.fasta"
        fusion_gff3="$ctat_intermediates/LR-FI.mm2.fusion_transcripts.gff3"
        if [[ ! -s "$read_fasta" ]]; then
            log "Skipping sample because CTAT read FASTA was not generated: $read_fasta"
            record_skip "$skipped_table" "$sample" "missing_ctat_read_fasta" "$input_bam" "$mapped_reads"
            continue
        fi
        if [[ ! -s "$fusion_gff3" ]]; then
            log "Skipping sample because CTAT fusion GFF3 was not generated: $fusion_gff3"
            record_skip "$skipped_table" "$sample" "missing_ctat_fusion_gff3" "$input_bam" "$mapped_reads"
            continue
        fi

        run_command "Building fusion isoform database" "$sample_log" \
            python "$scripts_dir/make_ctat_lr_fusion_isoform_database.py" \
            --fusion_tsv "$fusion_tsv" \
            --read_fasta "$read_fasta" \
            --fusion_gff3 "$fusion_gff3" \
            --outdir "$db_out" \
            --min_nr_support "$min_nr_support" \
            --id_prefix "$id_prefix" \
            --junction_tolerance "$junction_tolerance" \
            --breakpoint_tolerance "$breakpoint_tolerance" \
            --ignore_alignment_strand_for_clustering

        run_command "Running database QC" "$sample_log" \
            python "$scripts_dir/qc_ctat_lr_fusion_isoform_database.py" \
            --fusion_tsv "$fusion_tsv" \
            --nr_metadata "$db_out/fusion_isoforms.nr_metadata.tsv" \
            --nr_support "$db_out/fusion_isoform_read_support.tsv" \
            --raw_manifest "$db_out/fusion_isoforms.raw_manifest.tsv" \
            --outdir "$qc_out" \
            --high_conf_support "$high_conf_support"

        log "Collecting sample-prefixed core outputs"
        copy_required "$fusion_tsv" "$core_out/${sample}.$(basename "$fusion_tsv")"
        copy_required "$db_out/fusion_isoforms.raw_reads.fa" "$core_out/${sample}.fusion_isoforms.raw_reads.fa"
        copy_required "$db_out/fusion_isoforms.raw_blocks.gtf" "$core_out/${sample}.fusion_isoforms.raw_blocks.gtf"
        copy_required "$db_out/fusion_isoforms.raw_manifest.tsv" "$core_out/${sample}.fusion_isoforms.raw_manifest.tsv"
        copy_required "$db_out/fusion_isoforms.nr.fa" "$core_out/${sample}.fusion_isoforms.nr.fa"
        copy_required "$db_out/fusion_isoforms.nr.gtf" "$core_out/${sample}.fusion_isoforms.nr.gtf"
        copy_required "$db_out/fusion_isoform_read_support.tsv" "$core_out/${sample}.fusion_isoform_read_support.tsv"
        copy_required "$db_out/fusion_isoforms.nr_metadata.tsv" "$core_out/${sample}.fusion_isoforms.nr_metadata.tsv"
        copy_optional_glob "$qc_out/*" "$core_out" "$sample"

        cp -f "$sample_log" "$core_out/${sample}.run.log"

        find "$core_out" -maxdepth 1 -type f | sort > "$sample_out/${sample}.core_files.retained.txt"

        if [[ "$no_cleanup" -eq 0 ]]; then
            log "Cleaning non-core generated files inside sample output directory"

            if [[ "$keep_ctat_workdir" -eq 0 ]]; then
                safe_rm_path "$sample_out" "$ctat_out"
            elif [[ "$keep_ctat_intermediates" -eq 0 ]]; then
                safe_rm_path "$sample_out" "$ctat_intermediates"
                safe_rm_path "$sample_out" "$ctat_out/__checkpts"
                safe_rm_path "$sample_out" "$ctat_out/__checkpts_phase2"
                safe_rm_path "$sample_out" "$ctat_out/__mm2_prep_chkpts"
            fi
        else
            log "Cleanup disabled by --no-cleanup"
        fi

        log "Done: $sample"
        log "Core outputs: $core_out"
    } > >(tee -a "$sample_log") 2>&1
done

log "All samples finished. Parent output: $outdir"

