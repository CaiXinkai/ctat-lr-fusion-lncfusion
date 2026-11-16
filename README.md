# FASTQ → BAM → CTAT-LR-Fusion (LncFusion) → isoform sequences

A single-entry workflow for long-read fusion RNA detection with
[CTAT-LR-Fusion](https://github.com/TrinityCTAT/CTAT-LR-fusion) and the expanded
LncFusion CTAT library (LncBook v2 + GENCODE v42, GRCh38/hg38).

The workflow keeps a lossless unfiltered layer: no exact entity or reconstructed
sequence is deleted. It also writes a separate `isoform_tiers/` comparison
layer using the earlier pan-cancer clustering principles and progressively
stricter support tiers. These derivative files never replace the raw layer.

## Result policy

The fusion-call deliverable is exclusively:

```text
SAMPLE.ctat-LR-fusion.fusion_predictions.tsv
```

This is CTAT-LR-Fusion's official final table. The preliminary table and fusion
GFF3 remain under `work/` and are used only as coordinate evidence for sequence
reconstruction. They never add an event absent from the final table.

CTAT itself applies its official internal logic when creating the final table.
Removing that logic would require using the preliminary table, which is not a
final CTAT result.

## Requirements

- Linux
- Python 3.8+
- CTAT-LR-Fusion in `PATH`
- minimap2
- samtools
- gzip
- LncFusion `Lib` directory

## One-command run

```bash
bash scripts/run_fastq_to_ctat_lncfusion_isoforms.sh \
  --fastq /data/project/fastq/SAMPLE.fastq.gz \
  --sample SAMPLE \
  --outdir /data/project/ctat_lncfusion \
  --lncfusion-lib /data3/LncFusion/Lib \
  --platform ont-cdna \
  --threads 32 \
  --vis
```

Platforms: `ont-cdna`, `ont-directrna`, or `pacbio-isoseq`.

### Resume a post-result IGV failure

CTAT-LR-Fusion writes its final prediction table before preparing the IGV
bundle. If IGV preparation fails after the final table, resume without
re-aligning the FASTQ or rerunning CTAT:

```bash
bash scripts/run_fastq_to_ctat_lncfusion_isoforms.sh \
  --fastq /data/project/fastq/SAMPLE.fastq.gz \
  --sample SAMPLE \
  --outdir /data/project/ctat_lncfusion \
  --lncfusion-lib /data3/LncFusion/Lib \
  --platform ont-cdna \
  --threads 32 \
  --resume
```

Resume requires the final table, preliminary coordinate table, and fusion GFF3
to be present and nonempty. The recovery decision is recorded in
`SAMPLE.ctat_run_status.tsv`.

For multiple samples, run the command once per FASTQ with a unique `--sample`.
The sample directories are independent and suitable for job arrays.

## What the command does

1. Resolves `LncFusion_v1_CTAT/ctat_genome_lib_build_dir`.
2. Aligns FASTQ to that library's exact `ref_genome.fa` with a
   platform-appropriate splice-aware minimap2 preset.
3. Sorts and indexes the BAM.
4. Runs CTAT-LR-Fusion in BAM kickstart mode.
5. Copies CTAT's official final table into `final/`.
6. Constructs one entity for each distinct exact block chain supported by the
   final CTAT calls. This lossless table is never filtered or overwritten.
7. Projects CTAT fusion-contig blocks back to hg38, extracts reference intervals
   in transcript order, and uses the matched LncFusion GTF only to complete
   compatible outer terminal-exon boundaries.
8. Separately clusters entities with 5-bp breakpoint/junction tolerance,
   terminal-subset collapsing, and alignment-strand-independent matching, then
   writes progressively stricter comparison tiers.

## Outputs

```text
OUTDIR/SAMPLE/
├── final/
│   ├── SAMPLE.ctat-LR-fusion.fusion_predictions.tsv
│   ├── SAMPLE.ctat-LR-fusion.fusion_predictions.abridged.tsv
│   ├── SAMPLE.fusion_isoforms.reference.fa.gz
│   ├── SAMPLE.fusion_isoforms.entities.tsv
│   ├── SAMPLE.fusion_isoforms.audit.tsv.gz
│   ├── SAMPLE.fusion_isoforms.summary.json
│   ├── isoform_tiers/
│   │   ├── SAMPLE.isoform_tier_summary.tsv
│   │   ├── SAMPLE.isoform_tier_counts_by_fusion.tsv
│   │   ├── SAMPLE.01_clustered_all.tsv
│   │   ├── SAMPLE.04_pancancer_hc3.tsv
│   │   └── ...matching representative FASTA files...
│   └── SAMPLE.ctat_run_status.tsv
└── work/
    ├── SAMPLE.lncfusion_hg38.sorted.bam
    ├── SAMPLE.lncfusion_hg38.sorted.bam.bai
    ├── ctat_lr_fusion/
    └── sequence_reconstruction/
```

The abridged table is copied only when CTAT creates it. The BAM and CTAT
directory are reproducibility intermediates; the final CTAT table is the fusion
detection result.

## Post-detection condition integration

`scripts/integrate_compare_fusion_isoforms.py` implements the CTL/E2 integration
rule used in the accompanying MCF7 analysis. It is deliberately separate from
fusion detection and does not alter any CTAT result.

Prepare one TSV by concatenating the per-sample `01_clustered_all.tsv` files and
add a `group` column (for example `CTL` or `E2`). Then run:

```bash
python scripts/integrate_compare_fusion_isoforms.py \
  --input all_samples.clustered.tsv \
  --group-a CTL \
  --group-b E2 \
  --mapped-reads-a 7936508 \
  --mapped-reads-b 9345794 \
  --output-prefix results/MCF7_CTL_vs_E2
```

The integration rule:

1. removes non-monotonic/overlapping junction chains;
2. matches only the same directed fusion and splice type;
3. trims the first and last junction only for core matching when a chain has
   more than two junctions;
4. applies 5-bp breakpoint and junction tolerance;
5. permits a contiguous terminal subset during core matching;
6. retains the representative full chain and all source cluster support;
7. verifies support-read conservation after integration.

Exploratory expression candidates require total isoform reads >=5 and
absolute library-depth-normalized log2 fold change >=1. Within-fusion usage
candidates require at least 10 fusion reads in each group, at least 3 reads for
the isoform, and an absolute usage change >=20 percentage points. A 0.5-read
pseudocount is used only for the log2 fold change.

The outputs are:

```text
PREFIX.core_isoforms.tsv
PREFIX.expression_candidates.tsv
PREFIX.isoform_usage_changes.tsv
PREFIX.dominant_switches.tsv
PREFIX.summary.json
```

These are effect-size and prioritization tables. With one library per condition
they are not statistical differential-expression results: no valid p-value or
FDR can be estimated. Biological replicates must be kept separate for a formal
count model.

## Isoform sequence meaning

The FASTA is reference-guided and exon-chain-defined. Each sequence is rebuilt
from hg38 using the CTAT-observed block chain and normalized to the CTAT
left-gene → right-gene orientation. It is not the longest NR read.

The method preserves the observed internal splice structure and fusion junction
but does not invent unobserved transcript ends. These are therefore
**reference-corrected fusion isoform core sequences**, not automatically
experimentally proven complete 5′–3′ transcripts.

## Reproducibility

```bash
ctat-LR-fusion --version 2>&1 | head
minimap2 --version
samtools --version | head -1
sha256sum \
  /data3/LncFusion/Lib/LncFusion_v1_CTAT/ctat_genome_lib_build_dir/ref_genome.fa \
  /data3/LncFusion/Lib/LncFusion_v1_CTAT/ctat_genome_lib_build_dir/ref_annot.gtf
```

## Sources

- [LncFusion](https://github.com/CZhouLab/LncFusion)
- [CTAT-LR-Fusion](https://github.com/TrinityCTAT/CTAT-LR-fusion)
- [minimap2](https://github.com/lh3/minimap2)

