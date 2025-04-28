# Part III. 2.1 
## 1) RNA-seq 中归一化基因表达值的几种基本计算方法
**1. 样本内归一化方法**

针对单个样本内基因表达的可比性，主要校正测序深度和基因长度偏差。

(1) CPM（Counts Per Million）

• 公式：


$\text{CPM} = \frac{\text{基因原始 reads 数}}{\text{样本总 reads 数}} \times 10^6$
 

• 特点：

仅校正测序深度，未考虑基因长度，适用于样本间比较（如总 RNA 量相似的实验）。

不适用于样本内不同基因的比较（长基因的 reads 数天然更高）。

(2) RPKM/FPKM（Reads/Fragments Per Kilobase per Million）

• 公式：


$\text{RPKM/FPKM} = \frac{\text{基因 reads 数}}{\text{基因长度（KB）} \times \text{样本总 reads 数（百万）}}$  

• 特点：

同时校正测序深度和基因长度，适用于单样本内不同基因的比较。

缺点：样本间总 RPKM 值差异大，跨样本比较需谨慎。

(3) TPM（Transcripts Per Million）

• 公式：


$\text{TPM} = \frac{\text{基因 reads 数}}{\text{基因长度（KB）}} \div \left( \frac{\sum (\text{所有基因 reads 数/基因长度})}{10^6} \right)$

• 特点：

改进版 RPKM，样本间总和固定为 1 百万，更适合跨样本比较。

目前被视为基因定量的“金标准”。

**2. 样本间归一化方法**

针对多样本间的技术差异（如测序深度、文库复杂度），常用于差异表达分析。

(1) TMM（Trimmed Mean of M-values）

• 原理：

假设大部分基因无表达差异，通过修剪极端值（如去除最高/低 30% 的基因）计算样本间缩放因子。

适用于差异基因分析（如 edgeR 默认方法）。

• 优点：对高表达基因和异常值稳健。

(2) DESeq 大小因子（Size Factor）

• 原理：

基于基因表达几何均值的中位数计算缩放因子，假设大多数基因表达稳定。

适用于小样本数据，常用于 DESeq2 流程。

• 公式：


$\text{Size Factor} = \text{median} \left( \frac{\text{基因 counts}}{\text{所有样本的几何均值}} \right)$  

(3) 分位数归一化（Quantile Normalization）

• 原理：

强制所有样本的表达量分布一致（均值替换为相同分位数值）。

适用于批次效应校正，但可能掩盖真实生物学差异。

**3. 高级与跨数据集归一化方法**

(1) Z-Score 标准化

• 公式：


$Z = \frac{X - \mu}{\sigma}$ 

X：基因表达值  
$\mu$ ：所有样本均值  
$\sigma$：标准差

• 应用：消除量纲差异，用于跨平台数据整合或热图可视化。

(2) 跨数据集批次校正

• 工具：如 ComBat 或 Limma，通过经验贝叶斯方法校正批次效应。

• 适用场景：整合不同实验、平台或时间点的数据。
## 2)
Answer: E; D; A
## 3)
### sequencing protocol判断
输入：
`/usr/local/bin/infer_experiment.py -r GTF/Arabidopsis_thaliana.TAIR10.34.bed -i bam/Shape02.bam`  
输出:
```bash
Reading reference gene model GTF/Arabidopsis_thaliana.TAIR10.34.bed ... Done
Loading SAM/BAM file ...  Total 200000 usable reads were sampled


This is PairEnd Data
Fraction of reads failed to determine: 0.0315
Fraction of reads explained by "1++,1--,2+-,2-+": 0.4769
Fraction of reads explained by "1+-,1-+,2++,2--": 0.4916
```
由于"1++,1--,2+-,2-+"与"1+-,1-+,2++,2--"的比例几乎相同，有很大的把握认定这个数据是由非链特异性建库得到的。
### 计算shape02的read count matrix，给出AT1G09530基因(PIF3基因)上的counts数目
输入：
```
/home/software/subread-2.0.3-source/bin/featureCounts -s 0 -p -t exon -g gene_id -a GTF/Arabidopsis_thaliana.TAIR10.34.gtf -o result/Shape02.featurecounts.exon.txt bam/Shape02.bam
```
输出：
```bash

        ==========     _____ _    _ ____  _____  ______          _____
        =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \
          =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
            ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
              ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
        ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
	  v2.0.3

//========================== featureCounts setting ===========================\\
||                                                                            ||
||             Input files : 1 BAM file                                       ||
||                                                                            ||
||                           Shape02.bam                                      ||
||                                                                            ||
||             Output file : Shape02.featurecounts.exon.txt                   ||
||                 Summary : Shape02.featurecounts.exon.txt.summary           ||
||              Paired-end : yes                                              ||
||        Count read pairs : no                                               ||
||              Annotation : Arabidopsis_thaliana.TAIR10.34.gtf (GTF)         ||
||      Dir for temp files : result                                           ||
||                                                                            ||
||                 Threads : 1                                                ||
||                   Level : meta-feature level                               ||
||      Multimapping reads : not counted                                      ||
|| Multi-overlapping reads : not counted                                      ||
||   Min overlapping bases : 1                                                ||
||                                                                            ||
\\============================================================================//

//================================= Running ==================================\\
||                                                                            ||
|| Load annotation file Arabidopsis_thaliana.TAIR10.34.gtf ...                ||
||    Features : 313952                                                       ||
||    Meta-features : 32833                                                   ||
||    Chromosomes/contigs : 7                                                 ||
||                                                                            ||
|| Process BAM file Shape02.bam...                                            ||
||    Paired-end reads are included.                                          ||
||    The reads are assigned on the single-end mode.                          ||
||    Total alignments : 2730443                                              ||
||    Successfully assigned alignments : 2559170 (93.7%)                      ||
||    Running time : 0.03 minutes                                             ||
||                                                                            ||
|| Write the final count table.                                               ||
|| Write the read assignment summary.                                         ||
||                                                                            ||
|| Summary of counting results can be found in file "result/Shape02.featurec  ||
|| ounts.exon.txt.summary"                                                    ||
||                                                                            ||
\\============================================================================//
```
输入：
```bash
cat Shape02.featurecounts.exon.txt|awk '$1 =="AT1G09530" { print $7 } '
```
输出：
```
86
```
**故AT1G09530的counts数为86**  
输入：
```bash
cat Shape02.featurecounts.exon.txt|awk '!/^#/ {print $1, $7} '>/home/test/share/Shape02.featurecounts.exon.matrix.txt
```
***输出文件表达矩阵”Shape02.featurecounts.exon.matrix.txt“另附于作业中***
## 4)
R脚本：
```r
#COAD
setwd("~/Desktop/ZhaoAnlun_Bioinfo2/tumor-transcriptome-demo/COAD")
files <- list.files(full.names = T)
COAD_df <- data.frame(c(1:2000))
for (i in files){
  COAD_vec <- read.table(i,sep = "\t",skip = 1,header = T,row.names = 1)[,6]
  COAD_df <- cbind(COAD_df,COAD_vec)
}
COAD_df <- COAD_df[,-1]  
files1 <- sub("\\.txt$","",files)
colnames(COAD_df) <- files1
rownames(COAD_df) <- read.table(files[1],sep = "\t",skip = 1,header = T)[,1]  
COAD_matrix <- as.matrix(COAD_df)
library(edgeR)
y <- DGEList(counts = COAD_matrix) 
COAD.CPM.matrix <- edgeR::cpm(y,log=F) 
log10.COAD.CPM.matrix <- log10(COAD.CPM.matrix+1) 
COAD.z.scores <- (log10.COAD.CPM.matrix - rowMeans(log10.COAD.CPM.matrix))/apply(log10.COAD.CPM.matrix,1,sd)
library(pheatmap)
pheatmap(COAD.z.scores,treeheight_row = 140,treeheight_col = 11,cluster_cols = F,cluster_rows = F,
         color = colorRampPalette(c("blue","white","red"))(20000),border_color = NA,
         fontsize = 10,fontsize_row = 10,fontsize_col = 10,display_numbers = F,
         labels_row = rep("", nrow(COAD_df)), labels_col = rep("", ncol(COAD_df))) 
  
#ESCA
setwd("~/Desktop/ZhaoAnlun_Bioinfo2/tumor-transcriptome-demo/ESCA")
files <- list.files(full.names = T)
ESCA_df <- data.frame(c(1:2000))
for (i in files){
  ESCA_vec <- read.table(i,sep = "\t",skip = 1,header = T,row.names = 1)[,6]
  ESCA_df <- cbind(ESCA_df,ESCA_vec)
}
ESCA_df <- ESCA_df[,-1]  
files1 <- sub("\\.txt$","",files)
colnames(ESCA_df) <- files1
rownames(ESCA_df) <- read.table(files[1],sep = "\t",skip = 1,header = T)[,1]  
ESCA_matrix <- as.matrix(ESCA_df)
library(edgeR)
y <- DGEList(counts = ESCA_matrix) 
ESCA.CPM.matrix <- edgeR::cpm(y,log=F) 
log10.ESCA.CPM.matrix <- log10(ESCA.CPM.matrix+1) 
ESCA.z.scores <- (log10.ESCA.CPM.matrix - rowMeans(log10.ESCA.CPM.matrix))/apply(log10.ESCA.CPM.matrix,1,sd)
library(pheatmap)
pheatmap(ESCA.z.scores,treeheight_row = 140,treeheight_col = 11,cluster_cols = F,cluster_rows = F,
         color = colorRampPalette(c("blue","white","red"))(20000),border_color = NA,
         fontsize = 10,fontsize_row = 10,fontsize_col = 10,display_numbers = F,
         labels_row = rep("", nrow(ESCA_df)), labels_col = rep("", ncol(ESCA_df))) 
  
#READ
setwd("~/Desktop/ZhaoAnlun_Bioinfo2/tumor-transcriptome-demo/READ")
files <- list.files(full.names = T)
READ_df <- data.frame(c(1:2000))
for (i in files){
  READ_vec <- read.table(i,sep = "\t",skip = 1,header = T,row.names = 1)[,6]
  READ_df <- cbind(READ_df,READ_vec)
}
READ_df <- READ_df[,-1]  
files1 <- sub("\\.txt$","",files)
colnames(READ_df) <- files1
rownames(READ_df) <- read.table(files[1],sep = "\t",skip = 1,header = T)[,1]  
READ_matrix <- as.matrix(READ_df)
library(edgeR)
y <- DGEList(counts = READ_matrix) 
READ.CPM.matrix <- edgeR::cpm(y,log=F) 
log10.READ.CPM.matrix <- log10(READ.CPM.matrix+1) 
READ.z.scores <- (log10.READ.CPM.matrix - rowMeans(log10.READ.CPM.matrix))/apply(log10.READ.CPM.matrix,1,sd)
library(pheatmap)
pheatmap(READ.z.scores,treeheight_row = 140,treeheight_col = 11,cluster_cols = F,cluster_rows = F,
         color = colorRampPalette(c("blue","white","red"))(20000),border_color = NA,
         fontsize = 10,fontsize_row = 10,fontsize_col = 10,display_numbers = F,
         labels_row = rep("", nrow(READ_df)), labels_col = rep("", ncol(READ_df)))   
```
***heatmap另附于作业中***
根据heatmap结果，**READ和COAD**的转录组最相似。
# Part III. 2.3
## 1.
- 多重检验校正（multiple testing correction）  
多重检验校正是一种统计方法，用于在同时进行多次假设检验时控制错误发现的概率（如假阳性结果）。  
当研究者对大量变量（如基因、蛋白质、代谢物等）进行独立检验时，假阳性风险会显著增加。  
多重检验校正通过调整显著性阈值或p值，确保整体错误率（如假阳性率或错误发现率）处于合理范围内。
- p值与q值（FDR）的区别  

| 特征 |	p值	| q值（FDR） |
|---------|---------|---------|
| 控制目标 |	单次检验的假阳性率 | 显著结果中的假阳性比例 |
| 适用场景 |	少量检验	| 大规模检验（如基因分析）|
| 严格性 |	未校正时假阳性率高	| 更平衡假阳性和假阴性|
| 计算方法 |	直接计算	| 基于排序p值的调整公式|

## 2.归一化（normalization）方法
### DESeq2：Median of Ratios  
通过计算样本间基因表达比例的中位数，消除测序深度差异。**假设大多数基因不差异表达**。  
步骤与公式：
1. 计算每个基因的几何均值（所有样本的计数取几何平均）：

$GM(g) = \left( \prod_{i=1}^{n} K_{i,g} \right)^{1/n}$

  ◦ $K_{i,g}$：基因g在样本i中的原始计数。

  ◦ n：样本总数。

2. 计算每个基因的比值（样本计数与几何均值的比值）：

$R_{i,g} = \frac{K_{i,g}}{GM(g)}$


  ◦ 对每个样本，仅使用非零计数的基因（避免零除问题）。

3. 取比值的中位数作为样本的归一化因子：

$s_i = \text{median}\left( R_{i,g} \right)$

  ◦ $s_i$是样本i的归一化因子。

4. 归一化后的计数（用于后续分析）：

$\hat{K}\_{i,g} = \frac{K_{i,g}}{s_i}$


### edgeR：Trimmed Mean of M-values

步骤与公式：

1. 选择参考样本：通常取所有样本的几何均值或指定一个样本。

2. 计算基因的log-fold change (M) 和绝对表达量 (A)：

$M_g = \log_2\left(\frac{K_{i,g}/N_i}{K_{r,g}/N_r}\right), \quad A_g = \frac{1}{2} \log_2\left(\frac{K_{i,g}}{N_i} \cdot \frac{K_{r,g}}{N_r}\right)$

  ◦ $K_{i,g}$：样本i中基因g的计数。

  ◦ $N_i$, $N_r$：样本i和参考样本r的原始库大小（总计数）。

  ◦ $M_g$：基因g在样本i与参考样本间的log2比值。

  ◦ $A_g$：基因g的平均表达量。

3. 修剪极端值：

  ◦ 去除M_g和A_g分布中前30%和后5%的基因（默认参数），保留中间60%的基因。

4. 计算归一化因子：

$TMM_i = 2^{\text{weighted.mean}(M_g)}$

  ◦ 加权均值（权重为1/($A_g + \epsilon$)）以降低高表达基因的影响。

5. 调整库大小：

$N_i' = N_i \cdot TMM_i$

  ◦ $N_i'$是样本i的有效库大小，用于后续CPM（Counts Per Million）计算或模型拟合。
## 3.
R脚本如下：
```r
##DESeq2
#读取数据
raw.counts.mut <- read.table("~/Desktop/ZhaoAnlun_Bioinfo2/count_exon.txt", sep='\t', header = T,row.names = 1)
#我们这里只使用突变型数据进行分析
uvr8.raw.counts <- raw.counts.mut[,c("UD1_1", "UD1_2", "UD1_3", "UD0_1", "UD0_2", "UD0_3")]
#过滤掉表达量过低的基因
uvr8.filtered.counts <- uvr8.raw.counts[rowMeans(uvr8.raw.counts) > 5, ]
# "UD1_1", "UD1_2", "UD1_3" 三个样本为control
# "UD0_1", "UD0_2", "UD0_3 三个样本对应treatment
conditions <- factor(c(rep("Control", 3), rep("Treatment", 3)),levels = c("Control","Treatment"))
colData <- data.frame(row.names = colnames(uvr8.filtered.counts),conditions=conditions)
library(DESeq2)
dds <- DESeqDataSetFromMatrix(uvr8.filtered.counts, colData, design = ~conditions)
#进行差异分析
dds <- DESeq(dds)
#获取结果
res <- results(dds)
# 筛选出有差异的基因
# 过滤标准: padj < 0.05, log2 fold change > 1
diff.table <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
# 将结果保存至"uvr8.light.vs.dark.txt
write.table(diff.table,"~/Desktop/ZhaoAnlun_Bioinfo2/deseq2.diff.uvr8.light.vs.dark.txt", sep='\t', row.names = T, quote = F)

##edgeR
#读取数据
raw.counts.mut <- read.table("~/Desktop/ZhaoAnlun_Bioinfo2/count_exon.txt", sep='\t', header = T,row.names = 1)
#我们这里只使用突变型数据进行分析
uvr8.raw.counts <- raw.counts.mut[,c("UD1_1", "UD1_2", "UD1_3", "UD0_1", "UD0_2", "UD0_3")]
#过滤掉表达量过低的基因
uvr8.filtered.counts <- uvr8.raw.counts[rowMeans(uvr8.raw.counts) > 5, ]
conditions <- factor(c(rep("Control", 3), rep("Treatment", 3)),levels = c("Control","Treatment"))
#获取design矩阵
design <- model.matrix(~conditions)
library(edgeR) 
y <- DGEList(counts = uvr8.filtered.counts) # 定义edgeR用于存储基因表达信息的DGEList对象
# TMM标准化 (TMM 实际上是edgeR的默认参数)
y <- calcNormFactors(y, method="TMM")
# 估计dispersion
y <- estimateDisp(y,design = design)
# 拟合广义线性模型
fit <- glmFit(y, design = design)
# 似然比检验
lrt <- glmLRT(fit,coef=2) 
# 默认返回10个基因，按p值排序
# 用n = nrow(y)要求返回所有基因的结果
diff.table <- topTags(lrt, n = nrow(y))$table
# 只挑选显著变化的基因
diff.table.filtered <- diff.table[abs(diff.table$logFC) > 1 & diff.table$FDR < 0.05,]
write.table(diff.table.filtered, file = '~/Desktop/ZhaoAnlun_Bioinfo2/edger.diff.uvr8.light.vs.dark.txt', sep = "\t", quote = F, row.names = T, col.names = T)
```
***文本文件“deseq2.diff.uvr8.light.vs.dark.txt”和“edger.diff.uvr8.light.vs.dark.txt”另附于作业中***
## 4.
R脚本如下：
```r
deseq2.result <- read.table("~/Desktop/大二下/生物信息学/homework7/deseq2.diff.uvr8.light.vs.dark.txt", sep='\t', header = T, row.names = NULL)
edger.result <- read.table("~/Desktop/大二下/生物信息学/homework7/edger.diff.uvr8.light.vs.dark.txt", sep='\t',header = T,row.names = NULL)
result.list <- list(DESeq2 = deseq2.result$row.names, edgeR = edger.result$row.names)
library(VennDiagram)
venn.diagram(result.list,filename = "~/Desktop/大二下/生物信息学/homework7/venn_analysis.png",
             color = "black",fill = c("red","blue"),alpha = 0.5,cat.col = c("darkred","darkblue"),
             cex=3, category.names = c("DESeq2","edgeR"),cat.cex=1.3)
```
***图片文件“venn_analysis.png”另附于作业中***
## 5.
R脚本如下：
```r
edger.result <- read.table("~/Desktop/大二下/生物信息学/homework7/edger.diff.uvr8.light.vs.dark.txt", sep='\t',header = T,row.names = 1)
edger.result1 <- edger.result[order(edger.result$logFC),]
edger.result2 <- edger.result[order(-edger.result$logFC),]
selected_genes <- row.names(rbind(edger.result1[c(1:10),],edger.result2[c(1:10),]))
raw.counts <- read.table("~/Desktop/ZhaoAnlun_Bioinfo2/count_exon.txt", sep='\t', header = T,row.names = 1)
seleted.counts <- raw.counts[selected_genes,c(7:12)]
library.size <- colSums(seleted.counts)
cpm <- t(t(seleted.counts) / library.size) * 1e6
log10.cpm.selected <- log10(cpm+1)
heatmap.matrix <- as.matrix(scale(log10.cpm.selected))
library(pheatmap)
pheatmap(heatmap.matrix,cluster_cols = T,cluster_rows = T,
         color = colorRampPalette(c("blue","white","red"))(30000),border_color = NA,
         fontsize = 10,fontsize_row = 10,fontsize_col = 10,display_numbers = F)
```
***图片文件“heatmap.png”另附于作业中***



