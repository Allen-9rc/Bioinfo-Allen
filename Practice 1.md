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
`cat 1.gtf | awk '$0!~/^\s*$/{print}' | head -10`  
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
