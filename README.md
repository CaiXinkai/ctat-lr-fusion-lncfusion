# FASTQ → BAM → CTAT-LR-Fusion（LncFusion 库）→ fusion isoform 序列

本项目是一个用于长读长 fusion RNA 检测的单入口流程，使用
[CTAT-LR-Fusion](https://github.com/TrinityCTAT/CTAT-LR-fusion) 和扩展后的
LncFusion CTAT 参考库（LncBook v2 + GENCODE v42，GRCh38/hg38）。

流程保留无损的未过滤层：不删除任何精确 isoform entity 或重建序列。
同时在独立的 `isoform_tiers/` 目录中，按先前的泛癌聚类原则生成不同支持
强度的比较层。这些衍生结果不会覆盖原始无损层。

## 结果定义

fusion 检测的最终结果唯一指向：

```text
SAMPLE.ctat-LR-fusion.fusion_predictions.tsv
```

该文件是 CTAT-LR-Fusion 的官方最终表。preliminary 表和 fusion GFF3 保留在
`work/` 中，只用于提供序列重建所需的坐标证据，不会添加任何未出现在
最终表中的 fusion 事件。

CTAT 在生成最终表时会应用自身的官方内部逻辑。如果绕过这些逻辑，就必须使用
preliminary 表，而它不是 CTAT 的最终结果。

## 运行要求

- Linux
- Python 3.8+
- `PATH` 中可调用 CTAT-LR-Fusion
- minimap2
- samtools
- gzip
- LncFusion `Lib` 目录

## 单命令运行

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

`--platform` 可选：`ont-cdna`、`ont-directrna` 或 `pacbio-isoseq`。

### CTAT 已产生结果但 IGV 打包失败时续跑

CTAT-LR-Fusion 会在准备 IGV 数据包之前写出最终预测表。如果最终表已生成，
但后续 IGV 准备失败，可以在不重做 FASTQ 比对和 CTAT 检测的情况下续跑：

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

`--resume` 要求最终表、preliminary 坐标表和 fusion GFF3 均存在且非空。恢复决策会记录在
`SAMPLE.ctat_run_status.tsv` 中。

多样本时，对每个 FASTQ 分别运行一次，并使用唯一的 `--sample`。各样本目录相互独立，
可直接用于作业数组。

## 单命令流程的处理内容

1. 自动解析 `LncFusion_v1_CTAT/ctat_genome_lib_build_dir`。
2. 使用与平台相匹配的 minimap2 splice-aware 参数，将 FASTQ 比对到该库的 `ref_genome.fa`。
3. 对 BAM 排序并建立索引。
4. 使用 BAM kickstart 模式运行 CTAT-LR-Fusion。
5. 将 CTAT 官方最终表复制到 `final/`。
6. 为最终 CTAT 结果支持的每一种精确 block chain 建立一个 entity。该无损表不会被过滤或覆盖。
7. 将 CTAT fusion-contig block 投影回 hg38，按转录顺序提取参考区间，并仅使用匹配的 LncFusion GTF 补全兼容的外侧末端 exon 边界。
8. 在独立的衍生层中，以 5 bp 断点/junction 容差、terminal-subset 折叠和忽略比对方向的规则聚类 entity，再输出逐步严格的比较层。

## 输出结构

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
│   │   └── ...与各层匹配的代表 FASTA...
│   └── SAMPLE.ctat_run_status.tsv
└── work/
    ├── SAMPLE.lncfusion_hg38.sorted.bam
    ├── SAMPLE.lncfusion_hg38.sorted.bam.bai
    ├── ctat_lr_fusion/
    └── sequence_reconstruction/
```

只有当 CTAT 生成 abridged 表时，流程才会复制该文件。BAM 和 CTAT 目录属于可重复性中间文件；
CTAT 最终表才是 fusion 检测结果。

## 检测后的条件间整合

`scripts/integrate_compare_fusion_isoforms.py` 实现了 MCF7 CTL/E2 分析中使用的整合规则。
该步骤与 fusion 检测分离，不会修改任何 CTAT 原始结果。

先合并各样本的 `01_clustered_all.tsv`，并添加 `group` 列（例如 `CTL` 或 `E2`），再运行：

```bash
python scripts/integrate_compare_fusion_isoforms.py \
  --input all_samples.clustered.tsv \
  --group-a CTL \
  --group-b E2 \
  --mapped-reads-a 7936508 \
  --mapped-reads-b 9345794 \
  --output-prefix results/MCF7_CTL_vs_E2
```

整合规则：

1. 排除非单调或内部重叠的 junction chain。
2. 只匹配相同的有方向 fusion 和 splice type。
3. 当 chain 含有超过 2 个 junction 时，仅在 core 匹配阶段去除首尾 junction。
4. 断点和 junction 均使用 5 bp 容差。
5. core 匹配允许连续的 terminal subset。
6. 保留完整代表 chain 和所有来源 cluster 的支持数。
7. 整合后检查支持读段数是否守恒。

探索性表达变化候选要求 isoform 总支持读段数 >=5，且经文库深度标准化的
`|log2FC| >= 1`。fusion 内 isoform usage 候选要求两组中的 fusion 总支持数均 >=10，
目标 isoform 总支持数 >=3，且使用比例的绝对变化 >=20 个百分点。0.5 read 伪计数只用于 log2FC 计算。

输出：

```text
PREFIX.core_isoforms.tsv
PREFIX.expression_candidates.tsv
PREFIX.isoform_usage_changes.tsv
PREFIX.dominant_switches.tsv
PREFIX.summary.json
```

这些是效应量和候选优先级结果。如果每个条件只有一个文库，就不能将其解释为统计学差异表达，
也无法估计可靠的 P 值或 FDR。正式计数模型必须保留独立的生物学重复。

## Isoform 序列的含义

FASTA 是参考引导、由 exon chain 定义的序列。每条序列根据 CTAT 观察到的 block chain
从 hg38 重建，并归一到 CTAT 的左基因 → 右基因方向。它不是最长 NR read。

该方法保留观察到的内部剪接结构和 fusion junction，但不会虚构未观察到的转录本末端。
因此，这些序列是**参考校正的 fusion isoform 核心序列**，不等同于已被实验证明的完整 5′–3′ 全长转录本。

## 可重复性记录

```bash
ctat-LR-fusion --version 2>&1 | head
minimap2 --version
samtools --version | head -1
sha256sum \
  /data3/LncFusion/Lib/LncFusion_v1_CTAT/ctat_genome_lib_build_dir/ref_genome.fa \
  /data3/LncFusion/Lib/LncFusion_v1_CTAT/ctat_genome_lib_build_dir/ref_annot.gtf
```

## 参考项目

- [LncFusion](https://github.com/CZhouLab/LncFusion)
- [CTAT-LR-Fusion](https://github.com/TrinityCTAT/CTAT-LR-fusion)
- [minimap2](https://github.com/lh3/minimap2)
