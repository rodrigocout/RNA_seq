#!/bin/bash
#$ -S /bin/bash
#$ -q all.q
#$ -N processtable

cd /exports/molepi/RAAK_GARP_RNA_Seq/RNA_seq_BonesCart

for FILE1 in *1.fastq.gz
do
NEWFILE1=sample_$FILE1 

for FILE2 in *2.fastq.gz
NEWFILE2=sample_$FILE2

paste $FILE1 $FILE2 | column -s $'\t' -t > rnaseqsamples1.tsv
done

for FILE in *1.fastq.gz.txt
do
NEWFILE3=md5$FILE

for FILE in *2.fastq.gz.txt
do
NEWFILE4=md5$FILE

cat NEWFILE3 NEWFILE4 | column -s $'\t' -t > md5files_info.tsv
sort -k2 md5files_info.tsv > md5files_info.txt
done     