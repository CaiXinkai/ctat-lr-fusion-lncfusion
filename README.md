# CTAT-LR-Fusion with LncFusion annotation

Reproducible wrappers for detecting long-read fusion RNAs with CTAT-LR-Fusion
and the expanded LncFusion annotation library (LncBook v2 + GENCODE v42,
GRCh38/hg38).

This repository does **not** redistribute the multi-gigabyte LncFusion library,
container images, reference genome, or sequencing data. Download those resources
from their original providers and pass their local paths to the scripts.

## What this workflow does

1. Uses the CTAT genome library bundled with LncFusion:
   `LncFusion_v1_CTAT/ctat_genome_lib_build_dir`.
2. Runs CTAT-LR-Fusion on existing long-read genome-aligned BAM files.
3. Builds raw and non-redundant fusion RNA isoform FASTA/GTF/metadata outputs.
4. Performs QC by read support and splice type.
5. Collapses identical genomic breakpoint events and disambiguates overlapping
   gene annotations using read coverage of candidate-specific exons and splice
   junctions.

Fusion RNAs are treated as breakpoint- and splice-structure-defined events.
Gene-pair labels are annotations and are not used as the sole event identity.

## Requirements

- Linux
- Python 3.8 or newer (standard library only for the included Python scripts)
- CTAT-LR-Fusion in `PATH`
- `samtools`
- `minimap2` when starting from FASTQ
- LncFusion `Lib` directory
- sufficient disk space for BAM, CTAT intermediates, and the minimap2 index

Test availability:

```bash
ctat-LR-fusion --help >/dev/null
samtools --version | head -1
minimap2 --version
python --version
```

## LncFusion library layout

The workflow accepts any of these paths through `--lncfusion-lib`:

```text
/path/LncFusion/Lib
/path/LncFusion/Lib/LncFusion_v1_CTAT
/path/LncFusion/Lib/LncFusion_v1_CTAT/ctat_genome_lib_build_dir
```

The resolved CTAT directory must contain:

```text
ref_genome.fa
ref_genome.fa.fai
ref_annot.gtf
```

CTAT-LR-Fusion also needs a minimap2 index and splice BED. Build them once:

```bash
CTAT_LIB=/path/LncFusion/Lib/LncFusion_v1_CTAT/ctat_genome_lib_build_dir
EXAMPLE_BAM=/path/example_aln.bam

ctat-LR-fusion \
  --LR_bam "$EXAMPLE_BAM" \
  --genome_lib_dir "$CTAT_LIB" \
  --output /path/reference_prep \
  --CPU 32 \
  --prep_reference_only
```

This creates resources including:

```text
ref_genome.fa.mm2
ref_annot.gtf.mm2.splice.bed
```

## Run existing BAM files

```bash
bash scripts/run_ctat_lr_fusion_with_lncfusion_database.sh \
  --bam '/data/project/bam/*_aln.bam' \
  --outdir /data/project/ctat_lncfusion \
  --lncfusion-lib /data/reference/LncFusion/Lib \
  --threads 32 \
  --min-unique-bases 20 \
  --keep-ctat-workdir
```

The BAM reference dictionary is checked against the LncFusion FASTA. Different
ALT/random contig sets are allowed when all shared contigs have identical
lengths. Conflicting lengths or incompatible chromosome names stop the run.

## Run multiple FASTQ files

Copy and edit the sample table:

```bash
cp examples/samples.tsv.example samples.tsv
```

Then run:

```bash
bash examples/run_multiple_fastq.sh \
  --samples samples.tsv \
  --lncfusion-lib /data/reference/LncFusion/Lib \
  --outdir /data/project/lncfusion_analysis \
  --threads 32
```

The example aligns each FASTQ to the exact LncFusion `ref_genome.fa`, sorts and
indexes the BAM, then calls the BAM workflow. Its default minimap2 preset is
appropriate for ordinary long-read cDNA. Review `examples/run_multiple_fastq.sh`
before using direct-RNA or PacBio Iso-Seq data.

## Main outputs

For sample `SAMPLE`, retained files are written below:

```text
OUTDIR/SAMPLE/core/
```

Important files:

- `SAMPLE.ctat-LR-fusion.fusion_predictions*.tsv`: CTAT predictions.
- `SAMPLE.fusion_isoforms.nr.fa`: representative non-redundant isoform sequences.
- `SAMPLE.fusion_isoforms.nr.gtf`: isoforms on synthetic fusion contigs.
- `SAMPLE.fusion_isoforms.nr_metadata.tsv`: isoform metadata.
- `SAMPLE.fusion_isoform_read_support.tsv`: raw-read membership per isoform.
- `SAMPLE.fusion_rna_events.disambiguated.tsv`: breakpoint-defined events with
  unique-read support and overlapping-gene resolution.

Resolution states include:

```text
single_candidate_no_disambiguation_needed
resolved_by_unique_exon_or_junction
ambiguous_no_unique_read_evidence
ambiguous_tied_unique_evidence
```

For multi-sample integration, match events by chromosome, strand, and breakpoint
coordinates with an explicitly documented tolerance (for example 5 bp). Do not
sum support across alternative gene labels; use `unique_supporting_reads`.

## Suggested confidence tiers

An initial ranking can use:

- higher confidence: at least two unique supporting reads and `ONLY_REF_SPLICE`;
- intermediate: at least two reads and `INCL_NON_REF_SPLICE`;
- exploratory: singleton events.

Also review same-chromosome neighboring genes, read-through transcription,
paralogs/pseudogenes, repetitive sequence, internal oligo(dT) priming, mapping
quality, and IGV read evidence. These tiers are prioritization rules, not proof
of a genomic rearrangement.

## Reproducibility

Record software versions and reference checksums for every run:

```bash
ctat-LR-fusion --version 2>&1 | head
samtools --version | head -1
minimap2 --version
sha256sum "$CTAT_LIB/ref_genome.fa" "$CTAT_LIB/ref_annot.gtf"
```

Use the same CTAT-LR-Fusion version, LncFusion library, tolerances, and filtering
rules across all conditions.

## Included scripts

- `run_ctat_lr_fusion_with_lncfusion_database.sh`: LncFusion-aware entry point.
- `run_ctat_lr_fusion_isoform_database_cleanup.sh`: BAM batch workflow.
- `make_ctat_lr_fusion_isoform_database.py`: non-redundant isoform construction.
- `qc_ctat_lr_fusion_isoform_database.py`: QC summaries.
- `disambiguate_overlapping_fusion_annotations.py`: event-level gene annotation
  disambiguation using genome-aligned reads.

## Method description

Long-read fusion RNAs are detected using CTAT-LR-Fusion with a GRCh38 CTAT
genome library derived from the LncFusion reference, which integrates LncBook v2
and GENCODE v42 annotations. Candidate fusion transcripts are collapsed into
non-redundant isoforms using compatible fusion breakpoints and internal splice
junction chains. Fusion RNA events are subsequently consolidated by genomic
breakpoint pair and splice class. When a breakpoint overlaps multiple gene
models, supporting genome-aligned reads are evaluated for candidate-specific
exonic sequence and splice-junction evidence; unresolved alternatives are
retained as ambiguous annotations. Read support is counted using unique read
identifiers to prevent inflation by overlapping gene labels.

