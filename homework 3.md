# Part I 1.3 #
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
