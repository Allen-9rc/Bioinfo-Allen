# Part III. 2.1 
## 1)

## 2)

## 3)

## 4)

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
selected_genes <- rbind(edger.result1[c(1:10),],edger.result2[c(1:10),])
selected_genes$log10CPM <- log10(2^selected_genes$logCPM)
heatmap.matrix <- as.matrix(scale(selected_genes))
library(pheatmap)
pheatmap(heatmap.matrix,treeheight_row = 140,treeheight_col = 11,cluster_cols = F,cluster_rows = F,
         color = colorRampPalette(c("yellow","white","red"))(10000),border_color = NA,
         fontsize = 10,fontsize_row = 10,fontsize_col = 10,display_numbers = F)
```
***图片文件“heatmap.png”另附于作业中***



