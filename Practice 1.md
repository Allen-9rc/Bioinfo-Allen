step 0  
```
cd linux
gunzip 1.gtf.gz
ls
```
`1.gtf  bash_homework  file  test_command.gtf`

step 1  
`cat 1.gtf | head`
```
#!genome-build R64-1-1
#!genome-version R64-1-1
#!genome-date 2011-09
#!genome-build-accession :GCA_000146045.2
#!genebuild-last-updated 2011-12
IV	ensembl	gene	1802	2953	.	+	.	gene_id "YDL248W"; gene_version "1"; gene_name "COS7"; gene_source "ensembl"; gene_biotype "protein_coding";
IV	ensembl	transcript	1802	2953	.	+	.	gene_id "YDL248W"; gene_version "1"; transcript_id "YDL248W"; transcript_version "1"; gene_name "COS7"; gene_source "ensembl"; gene_biotype "protein_coding"; transcript_name "COS7"; transcript_source "ensembl"; transcript_biotype "protein_coding";
IV	ensembl	exon	1802	2953	.	+	.	gene_id "YDL248W"; gene_version "1"; transcript_id "YDL248W"; transcript_version "1"; exon_number "1"; gene_name "COS7"; gene_source "ensembl"; gene_biotype "protein_coding"; transcript_name "COS7"; transcript_source "ensembl"; transcript_biotype "protein_coding"; exon_id "YDL248W.1"; exon_version "1";
IV	ensembl	CDS	1802	2950	.	+	0	gene_id "YDL248W"; gene_version "1"; transcript_id "YDL248W"; transcript_version "1"; exon_number "1"; gene_name "COS7"; gene_source "ensembl"; gene_biotype "protein_coding"; transcript_name "COS7"; transcript_source "ensembl"; transcript_biotype "protein_coding"; protein_id "YDL248W"; protein_version "1";
IV	ensembl	start_codon	1802	1804	.	+	0	gene_id "YDL248W"; gene_version "1"; transcript_id "YDL248W"; transcript_version "1"; exon_number "1"; gene_name "COS7"; gene_source "ensembl"; gene_biotype "protein_coding"; transcript_name "COS7"; transcript_source "ensembl"; transcript_biotype "protein_coding";
```
`cat 1.gtf | tail`
```
Mito	ensembl	exon	85035	85112	.	+	.	gene_id "tM(CAU)Q2"; gene_version "1"; transcript_id "tM(CAU)Q2"; transcript_version "1"; exon_number "1"; gene_source "ensembl"; gene_biotype "tRNA"; transcript_name "tM(CAU)Q2"; transcript_source "ensembl"; transcript_biotype "tRNA"; exon_id "tM(CAU)Q2.1"; exon_version "1";
Mito	ensembl	gene	85295	85777	.	+	.	gene_id "RPM1"; gene_version "1"; gene_source "ensembl"; gene_biotype "ncRNA";
Mito	ensembl	transcript	85295	85777	.	+	.	gene_id "RPM1"; gene_version "1"; transcript_id "RPM1"; transcript_version "1"; gene_source "ensembl"; gene_biotype "ncRNA"; transcript_name "RPM1"; transcript_source "ensembl"; transcript_biotype "ncRNA";
Mito	ensembl	exon	85295	85777	.	+	.	gene_id "RPM1"; gene_version "1"; transcript_id "RPM1"; transcript_version "1"; exon_number "1"; gene_source "ensembl"; gene_biotype "ncRNA"; transcript_name "RPM1"; transcript_source "ensembl"; transcript_biotype "ncRNA"; exon_id "RPM1.1"; exon_version "1";
Mito	ensembl	gene	85554	85709	.	+	.	gene_id "Q0297"; gene_version "1"; gene_source "ensembl"; gene_biotype "protein_coding";
Mito	ensembl	transcript	85554	85709	.	+	.	gene_id "Q0297"; gene_version "1"; transcript_id "Q0297"; transcript_version "1"; gene_source "ensembl"; gene_biotype "protein_coding"; transcript_name "Q0297"; transcript_source "ensembl"; transcript_biotype "protein_coding";
Mito	ensembl	exon	85554	85709	.	+	.	gene_id "Q0297"; gene_version "1"; transcript_id "Q0297"; transcript_version "1"; exon_number "1"; gene_source "ensembl"; gene_biotype "protein_coding"; transcript_name "Q0297"; transcript_source "ensembl"; transcript_biotype "protein_coding"; exon_id "Q0297.1"; exon_version "1";
Mito	ensembl	CDS	85554	85706	.	+	0	gene_id "Q0297"; gene_version "1"; transcript_id "Q0297"; transcript_version "1"; exon_number "1"; gene_source "ensembl"; gene_biotype "protein_coding"; transcript_name "Q0297"; transcript_source "ensembl"; transcript_biotype "protein_coding"; protein_id "Q0297"; protein_version "1";
Mito	ensembl	start_codon	85554	85556	.	+	0	gene_id "Q0297"; gene_version "1"; transcript_id "Q0297"; transcript_version "1"; exon_number "1"; gene_source "ensembl"; gene_biotype "protein_coding"; transcript_name "Q0297"; transcript_source "ensembl"; transcript_biotype "protein_coding";
Mito	ensembl	stop_codon	85707	85709	.	+	0	gene_id "Q0297"; gene_version "1"; transcript_id "Q0297"; transcript_version "1"; exon_number "1"; gene_source "ensembl"; gene_biotype "protein_coding"; transcript_name "Q0297"; transcript_source "ensembl"; transcript_biotype "protein_coding";
```
`cat 1.gtf | head -2`
```
#!genome-build R64-1-1
#!genome-version R64-1-1
```
`ls -lh 1.gtf`  
`-rw-rw-r-- 1 test test 12M Sep 11  2018 1.gtf`  
`wc -l 1.gtf`  
`42252 1.gtf`  
`grep -v "^#" 1.gtf | grep -v "^$" | wc -l`  
`42247`  
`cat 1.gtf | awk '$0!~/^\s*$/{print}' | head -10`或`grep -v '^\s*$' 1.gtf | head -10`  
```
#!genome-build R64-1-1
#!genome-version R64-1-1
#!genome-date 2011-09
#!genome-build-accession :GCA_000146045.2
#!genebuild-last-updated 2011-12
IV	ensembl	gene	1802	2953	.	+	.	gene_id "YDL248W"; gene_version "1"; gene_name "COS7"; gene_source "ensembl"; gene_biotype "protein_coding";
IV	ensembl	transcript	1802	2953	.	+	.	gene_id "YDL248W"; gene_version "1"; transcript_id "YDL248W"; transcript_version "1"; gene_name "COS7"; gene_source "ensembl"; gene_biotype "protein_coding"; transcript_name "COS7"; transcript_source "ensembl"; transcript_biotype "protein_coding";
IV	ensembl	exon	1802	2953	.	+	.	gene_id "YDL248W"; gene_version "1"; transcript_id "YDL248W"; transcript_version "1"; exon_number "1"; gene_name "COS7"; gene_source "ensembl"; gene_biotype "protein_coding"; transcript_name "COS7"; transcript_source "ensembl"; transcript_biotype "protein_coding"; exon_id "YDL248W.1"; exon_version "1";
IV	ensembl	CDS	1802	2950	.	+	0	gene_id "YDL248W"; gene_version "1"; transcript_id "YDL248W"; transcript_version "1"; exon_number "1"; gene_name "COS7"; gene_source "ensembl"; gene_biotype "protein_coding"; transcript_name "COS7"; transcript_source "ensembl"; transcript_biotype "protein_coding"; protein_id "YDL248W"; protein_version "1";
IV	ensembl	start_codon	1802	1804	.	+	0	gene_id "YDL248W"; gene_version "1"; transcript_id "YDL248W"; transcript_version "1"; exon_number "1"; gene_name "COS7"; gene_source "ensembl"; gene_biotype "protein_coding"; transcript_name "COS7"; transcript_source "ensembl"; transcript_biotype "protein_coding";
```
step 2  
`cat 1.gtf | awk ' { print $1, $2, $3 } ' | head`或`cat 1.gtf | cut -f 1,2,3 | head`  
```
#!genome-build R64-1-1
#!genome-version R64-1-1
#!genome-date 2011-09
#!genome-build-accession :GCA_000146045.2
#!genebuild-last-updated 2011-12
IV ensembl gene
IV ensembl transcript
IV ensembl exon
IV ensembl CDS
IV ensembl start_codon
```
`cat 1.gtf | awk '$3 =="gene" { print $1, $3, $9 } ' | head`
```
IV gene gene_id
IV gene gene_id
IV gene gene_id
IV gene gene_id
IV gene gene_id
IV gene gene_id
IV gene gene_id
IV gene gene_id
IV gene gene_id
IV gene gene_id
```
step 3  
`grep -v '^#' 1.gtf |awk '{print $3}'| sort | uniq -c`
```
   7050 CDS
   7553 exon
   7126 gene
   6700 start_codon
   6692 stop_codon
   7126 transcript
```
`cat 1.gtf | awk 'BEGIN{L=0;}$3 =="CDS"{L+=$5-$4 + 1;}END{print L;}'`  
`9030648`  
`awk 'BEGIN  {s = 0;line = 0;}$3 =="CDS" && $1 =="I"{ s += $5-$4+1;line += 1}END {print "mean="s/line}' 1.gtf`  
`mean=1239.52`  
`cat 1.gtf | awk -F"\t" '$3 == "gene"{split($9,x,";");name = x[1];gsub("\"", "", name);gsub("gene_id", "", name);print "gene_id:",name,"length=",$5-$4+1}' | head` 
```
gene_id:  YDL248W length= 1152
gene_id:  YDL247W-A length= 75
gene_id:  YDL247W length= 1830
gene_id:  YDL246C length= 1074
gene_id:  YDL245C length= 1704
gene_id:  YDL244W length= 1023
gene_id:  YDL243C length= 990
gene_id:  YDL242W length= 354
gene_id:  YDL241W length= 372
gene_id:  YDL240C-A length= 138
```
step 4  
```
grep exon 1.gtf | awk '{print $5-$4+1}' | sort -n | tail -3 > 1.txt
cat 1.txt
mv 1.txt /home/test/share
```
或  
```
vi run.sh
#!/bin/bash   
grep exon *.gtf | awk '{print $5-$4+1}' | sort -n | tail -3
chmod u+x run.sh
./run.sh
```
```
12279
14730
14733
```
