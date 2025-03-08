# Part I 1.3 #
脚本以及两个输出结果文件内容见下，文件同时另附于提交的文件夹中
## bash脚本 ##
```
#!/bin/bash
echo -n "Please enter your directory: "
read dir
if [ ! -d "$dir" ];then
	echo "Such directory does not exist"
	break
else
	filename=""
	dirname=""
	cd "$dir"
	items=`ls`
	for item in $items;do
		if [ -f "$item" ];then
			filename="$filename\n$item"
		elif [ -d "$item" ];then
			dirname="$dirname\n$item"
		else
			continue
		fi
	done
	cd /home/test/linux
	echo -e "filenames:$filename\n"
	echo -e "dirnames:$dirname"
	echo  "$filename" > filename.txt
	echo  "$dirname" > dirname.txt
	echo -e `cat filename.txt | awk '{print substr($0, 3)}'` > filename.txt
	echo -e `cat dirname.txt | awk '{print substr($0, 3)}'` > dirname.txt
	echo "Already exported to /home/test/linux"
fi
exit 0
```
## filename.txt ##
```
a1.txt
a.txt
b1.txt
bam_wig.sh
b.filter_random.pl
c1.txt
chrom.size
c.txt
d1.txt
dir.txt
e1.txt
f1.txt
human_geneExp.txt
if.sh
image
insitiue.txt
mouse_geneExp.txt
name.txt
number.sh
out.bw
random.sh
read.sh
test3.sh
test4.sh
test.sh
test.txt
wigToBigWig
```
## dirname.txt ##
```
a-docker
app
backup
bin
biosoft
c1-RBPanno
datatable
db
download
e-annotation
exRNA
genome
git
highcharts
home
hub29
ibme
l-lwl
map2
mljs
module
mogproject
node_modules
perl5
postar2
postar_app
postar.docker
RBP_map
rout
script
script_backup
software
tcga
test
tmp
tmp_script
var
x-rbp
```
# Part II 1 #
## 1) ##
操作过程和结果图片见提交的文件夹，以下是对E值和P值实际意义的解释：  

1. E值  
E值表示在随机情况下，数据库中存在与当前比对得分相同或更高的匹配结果的期望次数。  
实际意义：  
E值越小，结果越显著：E值趋近于0时，表明随机匹配的可能性极低，比对结果更可能是真实的同源序列。  
判断标准：通常E值小于0.01或0.05时，认为比对结果具有统计显著性。  
与数据库大小相关：E值受数据库中序列数量的影响，数据库越大，E值可能越高，需结合实际情况分析。  
E值与序列长度相关：长序列的E值可能较高，需结合比对长度（Length）和匹配碱基数（Identities）综合判断。  

2. P值  
P值表示在随机情况下，观察到当前比对得分或更高得分的概率。  
实际意义：  
P值越小，结果越可靠：P值接近0时，表明当前得分不太可能是随机产生的，支持序列间存在真实相似性。  
与E值的关系：E值是P值的期望值。P值常用于补充E值，尤其在E值接近阈值（如0.01）时，P值可进一步验证结果的显著性。  
P值的局限性：需依赖特定的打分矩阵和比对算法，不同参数设置可能导致P值差异。  

## 2) ##
脚本、结果示例内容及对结果的解释见下，结果文件同时另附于提交的文件夹中
### bash脚本 ###
```
#!/bin/bash
cd /home/test/blast
original_seq="MSTRSVSSSSYRRMFGGPGTASRPSSSRSYVTTSTRTYSLGSALRPSTSRSLYASSPGGVYATRSSAVRL"
#创建目录存放随机序列
mkdir random_seqs
cd random_seqs
#生成随机序列
for i in `seq 1 10`;do
        #先把每个氨基酸切割在单独的一行内，再对行重排并删除换行符
        random_seq=`echo "$original_seq" | fold -w1 | shuf | tr -d '
        '`
        echo ">random_sequence_$i" > random_seq_$i.fasta
        echo "$random_seq" >> random_seq_$i.fasta
done
echo "Sequences have been generated in /home/test/blast/random_seqs"
cd ..
#创建目录存放比对结果
mkdir random_outputs
cd random_seqs
#生成两两比对的结果
for i in `seq 1 10`;do
	for j in `seq $((i+1)) 10`;do
                blastp -query random_seq_$i.fasta -subject random_seq_$j.fasta -out /home/test/blast/random_outputs/output_$i$j.txt
        done
done
cd ..
echo "Blastp aligenments are finished and kept in /home/test/blast/random_outputs"
```
### 结果示例 ###
```
BLASTP 2.6.0+


Reference: Stephen F. Altschul, Thomas L. Madden, Alejandro A.
Schaffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J.
Lipman (1997), "Gapped BLAST and PSI-BLAST: a new generation of
protein database search programs", Nucleic Acids Res. 25:3389-3402.


Reference for composition-based statistics: Alejandro A. Schaffer,
L. Aravind, Thomas L. Madden, Sergei Shavirin, John L. Spouge, Yuri
I. Wolf, Eugene V. Koonin, and Stephen F. Altschul (2001),
"Improving the accuracy of PSI-BLAST protein database searches with
composition-based statistics and other refinements", Nucleic Acids
Res. 29:2994-3005.



Database: User specified sequence set (Input: random_seq_2.fasta).
           1 sequences; 70 total letters



Query= random_sequence_1

Length=70
                                                                      Score     E
Sequences producing significant alignments:                          (Bits)  Value

  random_sequence_2                                                   11.5    1.6  


> random_sequence_2
Length=70

 Score = 11.5 bits (18),  Expect = 1.6, Method: Compositional matrix adjust.
 Identities = 4/6 (67%), Positives = 5/6 (83%), Gaps = 0/6 (0%)

Query  49  FSYSRS  54
            SY+RS
Sbjct  37  MSYTRS  42



Lambda      K        H        a         alpha
   0.307    0.116    0.309    0.792     4.96 

Gapped
Lambda      K        H        a         alpha    sigma
   0.267   0.0410    0.140     1.90     42.6     43.6 

Effective search space used: 4096


  Database: User specified sequence set (Input: random_seq_2.fasta).
    Posted date:  Unknown
  Number of letters in database: 70
  Number of sequences in database:  1



Matrix: BLOSUM62
Gap Penalties: Existence: 11, Extension: 1
Neighboring words threshold: 11
Window for multiple hits: 40
```
### 对结果的解释： ###
使用的BLASTP程序版本为2.6.0或更高版本  
比对过程中使用了BLOSUM62矩阵来计算氨基酸残基之间的相似性得分   
输入的文件（数据库）为random_seq_2.fasta，含1个序列，总共有70个氨基酸残基  
查询的序列为random_seq_1.fasta，含1个序列，70个氨基酸残基
比对得分为11.5bits，E值为1.6，P值约为0.7981  
查询序列和目标序列在比对区域中有4个相同的氨基酸残基，总比对长度为6个氨基酸残基，相同比例为67%  
查询序列和目标序列在比对区域中有5个位置上的氨基酸残基具有相似的化学性质，比例为83%  
在比对区域中有0个间隙，比例为0%  
该hit相似度较高，但由于hit长度较短，故而E值大，说明这个比对结果可能是随机产生的，不具有显著性  
整体上这两个蛋白质序列不具有显著的相似性，但这也是随机打乱生成序列的合理结果


## 3) ##
1. Seeding-and-Extending 策略

方法：

• 种子生成：将查询序列和数据库序列切分成短片段（称为“种子”或“seed words”），例如通过哈希函数生成固定长度的短序列片段。

• 索引优化：预先构建种子序列的索引（如哈希表），快速定位数据库中与种子高度匹配的候选区域，避免全局扫描。

• 局部扩展：仅对种子匹配区域进行动态规划比对，而非整个序列，显著减少计算量。

加速原理：
种子策略将全局比对问题转化为局部比对，通过索引快速缩小搜索空间，动态规划仅作用于高概率匹配区域，避免无效计算。

2. 并行化处理

方法：

• 多节点并行：将输入序列分割成多个部分，分配到不同计算节点（如高性能集群）并行处理。

• 多线程加速：在单节点内利用多核CPU并行执行多个比对任务。

加速原理：
通过分布式计算或硬件并行，将任务分解为子任务同时执行，充分利用计算资源，降低总耗时。

3. 参数调整与启发式算法

方法：

• 调整阈值参数：例如增大“word-size”（种子长度阈值）或放宽E值阈值，减少需要比对的种子数量。

• Fast模式：在BLASTX等工具中，通过简化算法（如减少比对长度或调整打分矩阵）加速计算，牺牲部分灵敏度换取速度。

加速原理：
通过参数调整过滤低概率匹配，减少动态规划的计算量；启发式算法（如贪婪算法）优先选择高得分路径，避免全局最优解的复杂计算。

4. 低复杂度屏蔽

方法：

• 屏蔽低复杂度区域：例如通过K值（信息含量）过滤重复序列（如微卫星重复），避免产生大量假阳性种子匹配。

加速原理：
低复杂度区域（如纯碱基重复）通常缺乏生物学意义，屏蔽后可减少无效种子生成和比对，提升整体效率。

5. 数据库预处理

方法：

• 本地化数据库：将数据库文件复制到计算节点本地磁盘，减少网络I/O延迟。

• 条带化存储：在分布式文件系统中对数据库文件进行条带化存储，提升并行读写速度。

加速原理：
优化数据访问路径，减少I/O瓶颈，使计算节点能更快读取所需数据。
## 4) ##
1. 对称PAM250矩阵

• 对称性原因：基于双向替换概率，即假设从氨基酸A到B的替换概率与从B到A的替换概率相等。

应用

• 常规序列比对：如BLAST、Clustal Omega等工具默认使用对称PAM矩阵，假设进化过程中替换是双向等概率的。

• 进化树构建：在对称模型下，进化路径的方向不影响替换概率，适合无根树的拓扑分析。

• 密码学与随机序列生成：对称性简化了替换规则的建模。

2. 非对称PAM250矩阵

• 非对称性原因：基于单向替换概率，即允许  M_{AB} \neq M_{BA} 。

应用

• 方向敏感型比对：如RNA结构预测中，考虑密码子翻译方向（5'→3'）的替换倾向。

• 人工序列设计：在密码学中，非对称矩阵可增强加密安全性（如替换规则不可逆）。

• 定向进化分析：在定向进化实验中，非对称矩阵可模拟特定方向的氨基酸替换压力。
