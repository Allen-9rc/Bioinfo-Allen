# 1. Linux #
## 1.1 Basic Command ##
### 1. ###
`wc -l test_command.gtf`   
`wc -c test_command.gtf`  
### 2. ###
`grep 'chr_' test_command.gtf | grep 'YDL248W'`  
### 3. ###
`sed 's/chr_/chromosome_/g' test_command.gtf | awk '{print $1 $3 $4 $5}'`
### 4. ###
