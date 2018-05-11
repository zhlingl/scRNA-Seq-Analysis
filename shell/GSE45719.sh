#!/bin/bash

# download data
wget -O data.tar \
'http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE45719&format=file'
mkdir data
tar -C data -xvf data.tar
gunzip data/*
 ls data -lh| wc -l
318
#解压后共317个文件(ls会多出一个total的)

#zy
paste data/GSM111*_zy* | \
awk '{for (i = 4; i <= NF; i += 6) printf ("%s%c", $i, i + 6 <= NF ? "\t" : "\n");}' > \
reads_zy.txt
paste data/GSM111*_early2cell_* | \
awk '{for (i = 4; i <= NF; i += 6) printf ("%s%c", $i, i + 6 <= NF ? "\t" : "\n");}' > \
reads_early2cell.txt
paste data/GSM111*_mid2cell_* | \
awk '{for (i = 4; i <= NF; i += 6) printf ("%s%c", $i, i + 6 <= NF ? "\t" : "\n");}' > \
reads_mid2cell.txt
paste data/GSM111*_late2cell_* | \
awk '{for (i = 4; i <= NF; i += 6) printf ("%s%c", $i, i + 6 <= NF ? "\t" : "\n");}' > \
reads_late2cell.txt
paste data/GSM111*_4cell_* | \
awk '{for (i = 4; i <= NF; i += 6) printf ("%s%c", $i, i + 6 <= NF ? "\t" : "\n");}' > \
reads_4cell.txt
# 8cell files have a bit different notation
paste data/*_8cell_*-* | \
awk '{for (i = 4; i <= NF; i += 6) printf ("%s%c", $i, i + 6 <= NF ? "\t" : "\n");}' > \
reads_8cell.txt
paste data/GSM111*_16cell_* | \
awk '{for (i = 4; i <= NF; i += 6) printf ("%s%c", $i, i + 6 <= NF ? "\t" : "\n");}' > \
reads_16cell.txt
paste data/GSM111*_earlyblast_* | \
awk '{for (i = 4; i <= NF; i += 6) printf ("%s%c", $i, i + 6 <= NF ? "\t" : "\n");}' > \
reads_earlyblast.txt
paste data/GSM111*_midblast_* | \
awk '{for (i = 4; i <= NF; i += 6) printf ("%s%c", $i, i + 6 <= NF ? "\t" : "\n");}' > \
reads_midblast.txt
paste data/GSM111*_lateblast_* | \
awk '{for (i = 4; i <= NF; i += 6) printf ("%s%c", $i, i + 6 <= NF ? "\t" : "\n");}' > \
reads_lateblast.txt

paste reads_zy.txt reads_early2cell.txt reads_mid2cell.txt reads_late2cell.txt reads_4cell.txt reads_8cell.txt reads_16cell.txt reads_earlyblast.txt reads_midblast.txt reads_lateblast.txt > deng.txt
sed -i '1d' deng.txt #删除第一行
head -n 1 deng.txt |awk '{print NF;}' 
268  
#共268个细胞系

cd data
for i in `ls GSM111*_zy*`; do  echo -e "$i\t\c" >> ../filename.txt ;  done;
for i in `ls GSM111*_early2cell_*`; do  echo -e "$i\t\c" >> ../filename.txt ;  done;
for i in `ls GSM111*_mid2cell_*`; do  echo -e "$i\t\c" >> ../filename.txt ;  done;
for i in `ls GSM111*_late2cell_*`; do  echo -e "$i\t\c" >> ../filename.txt ;  done;
for i in `ls GSM111*_4cell_*`; do  echo -e "$i\t\c" >> ../filename.txt ;  done;
for i in `ls *_8cell_*-*`; do  echo -e "$i\t\c" >> ../filename.txt ;  done;
for i in `ls GSM111*_16cell_*`; do  echo -e "$i\t\c" >> ../filename.txt ;  done;
for i in `ls GSM111*_earlyblast_*`; do  echo -e "$i\t\c" >> ../filename.txt ;  done;
for i in `ls GSM111*_midblast_*`; do  echo -e "$i\t\c" >> ../filename.txt ;  done;
for i in `ls GSM111*_lateblast_*`; do  echo -e "$i\t\c" >> ../filename.txt ;  done;
echo  >> ../filename.txt
cd ../data
head -n 1 filename.txt |awk '{print NF;}' 
268
#共268个文件名
awk 'BEGIN{FS="\t"; OFS="\t";}{for(i=1; i<NF; i+=1){split($i,a,"_"); printf("%s\t",a[1]);}}END{print "";}' filename.txt > GSM_num.txt

cat GSM_num.txt deng.txt >deng_name.txt
wc -l deng.txt 
22958 deng.txt
#共22958个基因
wc -l deng_name.txt 
22959 deng_name.txt

awk -F"\t" '{if ($1) print $1}' data/GSM1112767_zy2_expression.txt > gene-names.txt
wc -l gene-names.txt 
22959 gene-names.txt

paste gene-names.txt deng_name.txt > GSE45719_name.txt
sed -i '1s/^#Gene_symbol/id/' GSE45719_name.txt

awk 'NR==FNR{hash[$5]=$4;} NR!=FNR{ if($1 == "id") {print $0;} else {if(hash[$1]){gsub($1,hash[$1]); print $0;}}}' /data2/zll/project/deepBaseV3/genome/mm10/Mus_musculus.GRCm38.92.chr.bed GSE45719_name.txt > GSE45719.txt
