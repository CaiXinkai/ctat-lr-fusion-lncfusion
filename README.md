# FASTQ → BAM → CTAT-LR-Fusion (LncFusion) → isoform sequences

A single-entry workflow for long-read fusion RNA detection with
[CTAT-LR-Fusion](https://github.com/TrinityCTAT/CTAT-LR-fusion) and the expanded
LncFusion CTAT library (LncBook v2 + GENCODE v42, GRCh38/hg38).

The workflow deliberately has no post-CTAT filtering: no minimum supporting-read
threshold, no confidence tier, no overlapping-gene disambiguation, no NR
clustering, and no pan-cancer integration.

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
   final CTAT calls. No fuzzy/tolerance-based NR clustering is performed.
7. Projects CTAT fusion-contig blocks back to hg38, extracts reference intervals
   in transcript order, and uses the matched LncFusion GTF only to complete
   compatible outer terminal-exon boundaries.

## Outputs

```text
OUTDIR/SAMPLE/
├── final/
│   ├── SAMPLE.ctat-LR-fusion.fusion_predictions.tsv
│   ├── SAMPLE.ctat-LR-fusion.fusion_predictions.abridged.tsv
│   ├── SAMPLE.fusion_isoforms.reference.fa.gz
│   ├── SAMPLE.fusion_isoforms.entities.tsv
│   ├── SAMPLE.fusion_isoforms.audit.tsv.gz
│   └── SAMPLE.fusion_isoforms.summary.json
└── work/
    ├── SAMPLE.lncfusion_hg38.sorted.bam
    ├── SAMPLE.lncfusion_hg38.sorted.bam.bai
    ├── ctat_lr_fusion/
    └── sequence_reconstruction/
```

The abridged table is copied only when CTAT creates it. The BAM and CTAT
directory are reproducibility intermediates; the final CTAT table is the fusion
detection result.

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

