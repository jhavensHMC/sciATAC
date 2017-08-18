#!/bin/bash
#author: Jennifer Havens

#splits alignment.bam files which are seperated by samples into alignment.bam files seperated by cells
#run after Snakefile.preprocessing and SCIseq_MakeCellIDList.pl
#once done can do SCRAT analysis and count cell reads in peaks

OUTdir=/home/exacloud/lustre1/SpellmanLab/havens/ATAC/scAlign #gives directory where files split by cell will be output, end without /
cellIDList=/home/exacloud/lustre1/SpellmanLab/havens/ATAC/cellIDList.txt #generate with SCIseq_MakeCellIDList.pl
sampFILES=/home/exacloud/lustre1/SpellmanLab/havens/ATAC/sampAlign/* #gives directory where aligned files split by samples are, end with /*
scFILES=$OUTdir/* #gives directory where aligned files split by cells are put, end with *



for samp in $sampFILES
do
	for cellID in $cellIDList
	do
		perl SCIseq_FilterBamToCellIDList.pl $samp $cellID $OUTdir/$samp.$cellID.bam
	done
done

for cell in $scFILES
do
	mapped=$(samtools view -F 0x904 -c $cell)
	if [0 -eq mapped]
	then
		rm $cell
	fi
done


