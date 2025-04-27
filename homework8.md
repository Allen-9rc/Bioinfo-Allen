# PARTII-3.1
## 1.
***筛选基因的R脚本”Selecting.R“、筛选出的显著上调基因”genes.txt“、GO的结果文本文件”GO analysis.txt“另附于作业中***
## 2.
**1. Fold Enrichment（富集倍数）**

公式:  

$\text{Fold Enrichment} = \frac{\text{Observed Ratio}}{\text{Expected Ratio}} = \frac{\frac{\text{Number of genes in the GO term}}{\text{Total genes in the input list}}}{\frac{\text{Number of genes in the GO term (background)}}{\text{Total genes in the background}}}$

解释:  

• Observed Ratio：输入基因列表中属于该 GO term 的基因比例。

• Expected Ratio：在整个基因组（背景基因集）中该 GO term 的基因比例。

• Fold Enrichment > 1 表示该 GO term 在输入基因列表中比随机预期更富集。

• Fold Enrichment < 1 表示该 GO term 在输入基因列表中比随机预期更少出现（可能被抑制）。

**2. P value（P 值）**

公式:  

通常使用 超几何检验（Hypergeometric Test） 计算 P value：  

$P = \sum_{k=x}^{n} \frac{\binom{K}{k} \binom{N-K}{n-k}}{\binom{N}{n}}$

其中：  

•  N  = 背景基因总数

•  K  = 背景中属于该 GO term 的基因数

•  n  = 输入基因列表的基因数

•  x  = 输入基因列表中属于该 GO term 的基因数

解释:  

• P value 衡量的是 随机情况下 观察到至少  x  个基因属于该 GO term 的概率。

• P value 越小，说明该 GO term 的富集越不可能是随机发生的（即更显著）。

**3. 为什么用 FDR 而不是 P value？**

问题：多重假设检验（Multiple Testing Problem）

• GO 富集分析通常同时测试 成百上千个 GO terms，每个 GO term 计算一个 P value。

• 如果直接使用 P value（如  P < 0.05 ），会有大量 假阳性（False Positives） 出现。例如：

测试 1000 个 GO terms，即使所有 GO terms 都不显著，仍然会有约  1000 \times 0.05 = 50  个被错误地认为是显著的。

• 解决方案：FDR（False Discovery Rate）

FDR 控制的是 被错误判定为显著的 GO terms 的比例，而不是单个 P value 的显著性。常用的方法包括 Benjamini-Hochberg (BH) 校正：

$\text{FDR} = P \times \frac{\text{Total number of tests}}{\text{Rank of the P value}}$

• FDR 的优势  :

更严格的显著性标准：FDR 考虑了多重检验的影响，减少假阳性。  

生物学解释更可靠：FDR < 0.05 意味着在所有被判定为显著的 GO terms 中，假阳性的比例不超过 5%。  

# PARTII-3.2
***KEGG分析的结果图片和文本文件”KEGG analysis.txt""KEGG analysis.png"另附于作业***  
**相同点：**  
均用于解释差异表达基因的生物学意义（特定生物学背景下的潜在作用）；  
均依赖统计显著性（p value、FDR）；  
均基于数据库中已有的生物学知识进行注释。  

**不同点：**

|分析维度	|GO	|KEGG|
|-------|------|------|
|核心意义	|功能背景（分类与定位）	|具体通路（机制解析）|
|适用阶段	|初步功能探索	|深入机制研究|
|结果特点	|层级结构、泛化、冗余性高|	通路注释特异性强、直接关联分子事件|





