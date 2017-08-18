#!/bin/bash
#author: Jennifer Havens
#makes pheno_sheet.txt for monocle analysis
#run after split.sh
#once done can run monocle analysis


scFILES=/home/exacloud/lustre1/SpellmanLab/havens/ATAC/scAlign/* #gives directory where aligned files split by cells are put, end with *
OUTfile=/home/exacloud/lustre1/SpellmanLab/havens/ATAC/pheno_sheet.txt


# rows cells, columns are attibutes (sample and sites (counted unique reads)) 

echo cells    sample    sites > $OUTfile

#pattern of sc alignment file names $OUTdir/$samp.$cellID.bam

for cellFILE in $scFILES
do
	mapped=$(samtools view -F 0x904 -c $cellFILE)
	fileName=$(echo $cellFile | awk -F'/' '{print $(NF)}')
	cell=$(echo $fileName | awk -F'.' '{print $2}')
	sample=$(echo $fileName | awk -F'.' '{print $1}')
	echo $cell    $sample    $mapped >> pheno_sheet.txt
done



