# PART III-1.2
## (1)
- 代码如下
```bash
ls
samtools flagstat COAD.ACTB.bam
```
- 结果如下
```
185650 + 0 in total (QC-passed reads + QC-failed reads)
4923 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
185650 + 0 mapped (100.00% : N/A)
0 + 0 paired in sequencing
0 + 0 read1
0 + 0 read2
0 + 0 properly paired (N/A : N/A)
0 + 0 with itself and mate mapped
0 + 0 singletons (N/A : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
```
**关注以下两行输出**
```bash
185650 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 paired in sequencing
```
两行分别表示总reads数，和双端测序的reads数。由于`paired in sequencing`的值为0，所以所有reads都是单端测序。
## (2)
- Secondary Alignment 指一个read在参考基因组上有多个可能的比对位置，但未被选为最佳（Primary）的比对结果。特点如下：
  - 该 read 在基因组其他位置也有匹配，但比对工具（如 BWA、Bowtie2）根据算法选择了其中一个作为 primary alignment。
  - 在 SAM/BAM 文件中，secondary alignment 会被标记为 FLAG 的第 9 位。
- 统计bam文件中的 secondary alignment 记录数，代码如下：
```bash
samtools view -f 256 -c COAD.ACTB.bam
```
- 结果如下：
```
4923
```
所以共有**4923**条reads属于secondary alignment
## (3)
- 脚本如下：
```bash
#!/bin/bash
awk 'BEGIN{OFS="\t"} $3 == "gene" {print $1, $4-1, $5, "ACTB_gene", ".", $7}' hg38.ACTB.gff > gene.bed
awk 'BEGIN{OFS="\t"} $3 == "exon" {print $1, $4-1, $5, "exon", ".", $7}' hg38.ACTB.gff > exon.raw.bed
sort -k1,1 -k2,2n exon.raw.bed | bedtools merge -i - > exon.merged.bed
bedtools subtract -a gene.bed -b exon.merged.bed > intron.bed
samtools view -L intron.bed -b COAD.ACTB.bam > COAD.ACTB.intron.bam
samtools fastq COAD.ACTB.intron.bam > COAD.ACTB.intron.fastq
```
- 输入如下：
```bash
ls -la hm6-3.sh
chmod u+x hm6-3.sh
ls -la hm6-3.sh
./hm6-3.sh
ls
```
- 输出如下：
```
[M::bam2fq_mainloop] discarded 0 singletons
[M::bam2fq_mainloop] processed 15132 reads
COAD.ACTB.bam           exon.merged.bed  hg38.ACTB.gff  hm6-3.sh
COAD.ACTB.intron.bam    exon.raw.bed     intron.bed
COAD.ACTB.intron.fastq  gene.bed         README.md
```
## (4)
- 代码如下
```bash
bedtools genomecov -ibam COAD.ACTB.bam -bga -split | awk 'NR==FNR {chr=$1; start=$2; end=$3; next} $1==chr && $2>=start && $3<=end' gene.bed - > ACTB_coverage.bedgraph
ls
```
- 结果如下：
```
ACTB_coverage.bedgraph  COAD.ACTB.intron.fastq  gene.bed       README.md
COAD.ACTB.bam           exon.merged.bed         hg38.ACTB.gff  hm6-3.sh
COAD.ACTB.intron.bam    exon.raw.bed            intron.bed
```


# Extra Homework
## 1.人类基因组大小与基本组成

1.1 基因组大小

|数据库	| 版本	| 总碱基数 | 更新时间 |
|--------|------|-------|-------|
|NCBI|	GRCh38.p14|	3,272,116,950 bp|	2023-03|
|Ensembl	|Release 110|	3,272,116,950 bp|	2023-07|

注：两数据库采用相同参考基因组GRCh38，差异仅在于注释更新周期

1.2 基本组成（以NCBI RefSeq为例）
  

| 成分类型         | 数量/占比              | 功能说明                     |
|------------------|-----------------------|----------------------------|
| **蛋白质编码基因** | 19,969个             | 编码蛋白质的核心功能单元     |
| **非编码RNA**     | 35,415个（GENCODE v44）| 调控基因表达（详见第2部分）  |
| **重复序列**      | 约44.3%              | 包括：                     |
| - LINE           | 20.4%               | 长散布核元件（逆转录转座子）|
| - SINE           | 13.1%               | 短散布核元件（如Alu序列）  |
| - LTR            | 8.3%                | 长末端重复逆转录病毒元件    |
| **假基因**        | 14,939个            | 丧失功能的基因遗迹          |
| **调控区域**      | >1百万个             | 启动子、增强子等调控元件    |


## 2.非编码RNA注释最新进展
- 数据来源
[RNAcentral官网](https://rnacentral.org)  
RNAcentral官网博客，“v24”版本，发布时间： 26 Mar 2024  
- 结果解释：截至2024年3月26日，人类共有759,614条ncRNA。分类、数量与功能如下表：

|类型	|数量	|主要作用简述|	
|-----|------|--------|
| lncRNA | 423,118 |调控基因表达，影响染色质修饰、转录及翻译过程，在细胞分化和疾病（如癌症）中发挥重要作用|
| rRNA | 238,884 |是核糖体的组成部分，参与蛋白质合成过程中mRNA的解读和肽链的合成|
| sncRNA | 67,404 |包括miRNA、piRNA等，调控基因表达，参与RNA干扰、基因沉默和转座子抑制|
| RNase P RNA | 1,506 |作为RNase P酶的RNA组分，催化tRNA前体的5'端剪切，参与tRNA成熟化|
| SRP RNA | 1,143 |构成信号识别颗粒（SRP），帮助分泌蛋白和膜蛋白靶向内质网膜|
| Y RNA | 723 |与Ro核糖核蛋白复合体相关，参与RNA质量控制和DNA复制启动|
| circular ncRNA | 431 |调控基因表达，可作为miRNA海绵调控mRNA稳定性，可能影响疾病发生|
| antisense RNA | 276 |与互补mRNA配对，调控其翻译和稳定性，可用于基因沉默|
| ribozyme | 33 |具有催化活性的RNA，可催化自身剪切或其他RNA的加工|
| sc RNA | 17 |与信号识别颗粒相关，可能在蛋白质转运过程中发挥作用|
| vault RNA | 16 |与vault复合体相关，可能涉及药物耐受性和细胞信号传导|
| RNase MPR RNA | 12 |与RNase MRP复合体相关，参与rRNA加工和细胞周期调控|
| telomerase RNA | 5 |作为端粒酶的模板，延长端粒DNA，维持染色体稳定性|
| guide RNA | 1 |在RNA编辑过程中引导碱基修饰，如指导trypanosome的mRNA编辑|

