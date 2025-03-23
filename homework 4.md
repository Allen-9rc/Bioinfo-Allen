# PARTIII-1 #
## 1) ##
Bowtie通过结合Burrows-Wheeler变换（BWT）的**压缩与高效搜索特性**提升运算速度，同时采用**索引结构优化和参数调整策略**降低内存需求。
### 提高运算速度 ###
1. 重复序列压缩  
BWT将基因组中重复出现的碱基序列（如ATATAT）转换为连续的相同字符（如TTTAAA），大幅减少搜索空间。这种压缩特性使得比对时无需逐字符匹配，而是通过聚类后的区块快速定位潜在匹配区域。

2. 结合FM索引实现高效查询  
BWT与FM索引（Ferragina-Manzini索引）结合，利用前缀搜索快速定位子串。FM索引通过记录字符的频率和位置偏移量，允许比对算法跳过不匹配的区域，直接锁定可能的候选位置。例如，对长度为L的读长，搜索复杂度从O(L^2)降低至接近O(L)。

3. 支持反向互补链的快速比对  
BWT生成的索引文件包含正向和反向互补链的独立索引（如.rev.1.bt2文件），避免重复计算双链匹配，直接通过索引切换实现双向搜索。

4. 局部匹配与容错性  
BWT允许部分匹配（如允许错配或插入缺失），通过动态规划算法快速筛选最优比对路径，而非全序列遍历。
### 降低内存需求 ###
1. 调整offrate参数控制索引密度  
offrate（默认值为5）决定了索引的稀疏程度。降低offrate值（如设为4）会生成更密集的索引，提高比对速度但增加内存消耗；提高offrate值（如设为6）则减少索引密度，节省内存但略微降低速度。

2. 差分覆盖采样（Difference Cover Sampling）  
构建后缀数组时，通过差分覆盖采样技术减少后缀数组的存储规模。该技术仅保留关键位置的采样点，避免存储完整的后缀数组，内存占用可降低30%-50%。

3. 分块处理与多线程优化  

  • 将大型基因组分块处理，每次仅加载部分索引到内存中，避免一次性占用全部资源。

  • 使用-p参数启用多线程，每个线程独立处理数据块，并行化降低单线程内存压力。

4. 64位版本与内存预分配  

  • 优先使用64位Bowtie版本，支持更大内存寻址空间，避免32位系统的内存溢出问题。

  • 在编译时通过make BITS=64生成优化后的可执行文件。

5. 反向链索引的独立构建  
正向链和反向互补链的索引文件（如.1.bt2和.rev.1.bt2）独立存储，比对时按需加载，减少同时占用的内存量。
## 2) ##
`bowtie -v 2 -m 10 --best --strata BowtieIndex/YeastGenome -f THA2.fa -S THA2.sam`
```
# reads processed: 1250
# reads with at least one reported alignment: 1158 (92.64%)
# reads that failed to align: 77 (6.16%)
# reads with alignments suppressed due to -m: 15 (1.20%)
Reported 1158 alignments to 1 output stream(s)
```
`awk '$3 != "*" && $0!~/^@/ {print $3}' THA2.sam | sort | uniq -c`
```
     18 chrI
     51 chrII
     15 chrIII
    194 chrIV
     25 chrIX
     12 chrmt
     33 chrV
     17 chrVI
    125 chrVII
     68 chrVIII
     71 chrX
     56 chrXI
    169 chrXII
     67 chrXIII
     58 chrXIV
    101 chrXV
     78 chrXVI
```
## 3) ##
### 3.1) ###
CIGAR string（Compact Idiosyncratic Gapped Alignment Report）是SAM/BAM文件中记录read比对到参考基因组具体细节的关键字段。它通过数字与字母组合的形式，描述read与参考序列的匹配、插入、缺失等比对关系:

1. 操作符（Operators）：

  • M：匹配或错配（Match/Mismatch）。CIGAR中的M仅表示位置对齐，不区分完全匹配或单碱基错配。

  • I：插入（Insertion），即read中存在而参考序列中缺失的碱基。

  • D：缺失（Deletion），即参考序列中存在而read中缺失的碱基。

  • S：软剪接（Soft Clip），read两端未比对的部分仍保留在SEQ字段中。

  • H：硬剪接（Hard Clip），未比对部分被完全移除且不保留在SEQ中。

  • N：跳过（Skipped Region），常用于RNA比对中表示内含子区域。

  • =和X：分别表示完全匹配和明确错配（需特定比对工具支持）。

2. 长度规则：

  • M/I/S/=/X操作符的总长度必须等于SEQ字段的碱基数。例如，CIGAR为33S117M的read，其SEQ长度必须为33+117=150 bp。
### 3.2) ###
Soft Clip（软剪接）指read的两端部分碱基未比对到参考基因组，但这些序列仍保留在SEQ字段中。其作用是处理测序错误、结构变异或重复序列导致的比对不确定性。

• CIGAR表示：使用S标记。例如，CIGAR字符串150M33S表示read前150 bp完全比对，后33 bp未比对但保留在SEQ中。

• 规则：

  S只能出现在CIGAR字符串的两端（如S100M或100MS）。

  若S出现在中间，则两侧必须为H（硬剪接）。
### 3.3) ###
Mapping Quality（MAPQ）是SAM/BAM文件中第五列的值，表示比对结果的可靠性。其定义与计算方式如下：

• 定义：MAPQ = -10 × log10(p)，其中p为比对位置错误的概率。例如，MAPQ=30表示错误概率为0.1%（即30 = -10 × log10(0.001)）。

• 信息反映：

  高MAPQ（如≥30）：比对位置唯一且置信度高，常用于变异检测中的可靠位点筛选。

  低MAPQ（如≤10）：比对到多个位置或存在高重复性区域，结果不可靠。

  特殊值255：表示MAPQ不可用或未计算。
### 3.4) ###
**不能完全推断**。原因如下：

1. CIGAR与SEQ的局限性：

  • CIGAR仅描述比对操作（如M、D），但无法明确区分匹配与错配的碱基。例如，M操作符可能对应匹配或错配，但无法通过CIGAR确定具体错配位置。

  • SEQ字段存储read的原始序列（若比对到互补链则为反向互补序列），但无法直接还原参考序列的碱基。

2. 依赖MD标签：

  • MD tag是SAM/BAM的可选字段，记录参考基因组在错配和缺失处的实际碱基信息。例如，MD:Z:50A10^TG5表示参考序列在第51位为A，第61-62位缺失TG。

  • 若未包含MD tag，则无法准确重建参考序列。

结论：仅凭CIGAR和SEQ字段无法完整还原参考基因组序列，需依赖MD tag或其他注释信息。
## 4) ##
结果见THA2-bwa.sam文件
# PARTIII-1.1 #
## 截图另附于文件夹中 ##
