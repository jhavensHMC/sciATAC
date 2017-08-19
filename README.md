
ATACnotes1.2
#author: Jennifer Havens

AIM
	to identify characteristics of chromatin state which are associated with metastatic ability,
	and to do this by comparing chromatin accessibility in primary, low burden, and high burden tumors

Note: This pipeline was built with a focus on ease of use and getting results, it has not been optimized for efficiency 

Overview of Pipeline:
	Prepare read files - bcl2fast2q and Vitak
	Align - BWA and samtools
	Clean up sam/bam files - samtools Vitak and Picardtools
	Peak calling on aggregated samples - MACS2
	Look at aggregated samples in peaks - chromVAR analysis (GC bias correction)
	Split aggregated alignments to single cell bam files - Vitak
	Look at single cells in peaks
		SCRAT
		Monocle (TF-IDF normalization)
		TSCAN (TF-IDF normalization)
	look at expression along genome - custom (TF-IDF normalization)

Running Files:
	Snakefile.preprocessing
	Snakefile.align
	split.sh
	Snakefile.count
	tfIDF.r
	phenoSheet.sh
	monocle.r
	chromVarOutline.r
	Snakefile.linearCount
	putTogether.sh
	agCells.r


More detailed Pipeline:

	Prepare read files 
		#unfortunatly beacause I am unsure what format the files will come in I cant be more specific
		#this process should be identical to whatever is being done to process the scRNAseq files by the rest of the lab
		bcl2fast2q 
			with --with-failed-reads and --create-fastq-for-index-reads
	 	run Vitak's SCIseq_NextSeqFastq_to_SCIseqFastq.pl to remove indexes from reads and lable read 

	Align - run snakemake -s Snakefile.preprocessing makeIndex and snakemake -s Snakefile.preprocessing BWAalign 
		(put Snakefile.align in align/Snakefile)
		#can use something other than BWA, just chosen because it was used by Cusanovich et al 2017
		export PATH=$PATH:/opt/installed/bwa-0.7.12 
		bwa index -p hg19genome -a bwtsw hg19genome.fa #make index for refernace genome
		bwa mem -t {threads} hg19genome.fa read1.fq read2.fq > sampleA_full.sam
			#use -M for picard and GATK
		
	Convert to and clean up bam files - run snakemake -s Snakefile.preprocessing (does not include EstimateLibraryComplexity)
		convert sam to bam files
			samtools view -b sampleA_full.sam -o sampleA_full.bam
		sort bam files with samtools
			samtools sort sampleA_full.bam sampleA_full.sorted.bam
		RemoveDuplicates - Vitak
			SCIseq_RemoveDuplicates.pl sampleA_full.sorted.bam sampleA_noDup.bam
				#optional add path to SCIseq_RemoveDuplicatesPlot.r to get complexity figures
		FilterToReadThreshold - Vitak
			#optional: define Gaussian mixture model to the number of reads for each cell to choose the minimum depth (Cusanovich et al 2017)
			SCIseq_FilterBamToReadThreshold.pl sampleA_noDup.bam {threshold} sampleA.bam
		AddRGtoBam - Vitak
			SCIseq_AddRGtoBam.pl sampleA.bam sampleA_rg.bam
		EstimateLibraryComplexity - picard tools #optional

	Peak calling - MACS2 - in Snakefile.preprocessing
		macs2 callpeak -t sampleA_rg.bam -f BAM -n sampleA -B
	
	look at variation in peaks with RG bam files  
		chromVAR analysis parts - see chromVarOutline.R for detials and execution outline
			set up chromVar
			motif/annotate
			deviations
			diffential access
			variabilty
			cluster

	Split to single cell bam files - Vitak scripts -  
		split bam files in sampAlign to single cell bam files in scAlign folder:
		run perl SCIseq_MakeCellIDList.pl with correct indexes
		run chmod +x split.sh and ./split.sh
			
	Look at variation in features between single cells
		SCRAT
			#note: refreshing will erase completed work		
			load sc bam files into SCRAT
			add features - ENCODE, gene, custom (BED from MACS2)
				can try use motif and gene-set but not currently working
			similarity to existing cell types (makes heatmap)
			cluster w/ t-SNE
			can run on web app or download and run through shiny:
				devtools::install_github("zji90/SCRAT")
				library("SCRAT")
				install_github("SCRATexample", "zji90")
			to run: SCRATui()

		count reads for sc in features - using cufflinks by running: make sample_sheet.txt before run snakemake -s Snakefile.count
 
 			#before running please set up sample_sheet.txt with header:
			#sample_name                          group


			cuffquant -o CellA_cuffquant ucsc.hg.annotations.gtf cellA.bam
			cuffnorm --use-sample-sheet -o fpkm_matrix ucsc.hg.annotations.gtf sample_sheet.txt

			#sample_sheet.txt
			sample_name                          group
			cellA_cuffquant/abundances.cxb     cellA
			cellB_cuffquant/abundances.cxb     cellB

			#make into "fpkm_matrix.txt" #rows are gene names, columns are cells 

			apply TF-IDF normalization

		monocle.r 
			make pheno_sheet
				must split into single cell bam files first
				run chmod +x phenoSheet.sh
				./phenoSheet.sh
			install package
			newCellDataSet
			reduceDimension
			orderCells
			minSpaningTree -> plot_spanning_tree
			fitModel - with total # sites accesable as covariate
			differntialGeneTest 

		TSCAN
			https://zhiji.shinyapps.io/TSCAN/
			load fpkm_matrix.normalized.txt 
			recomend working on local R session with multithreading, rather than web, to speed up analysis
			can do pseudotime analysis and expression which is differntially expressed in pseudotime
			limited: cannot mark cells based on their samples with in this ordering 
				though the Miscellaneous option does allow for some sample based analysis

	make raw count matrix, aggregate, and analyze - custom  
		Make 5k genomic windows - bedtools
		count reads from sc w/ samtools and concatonate into 1 matrix
		
		Run:
		snakemake -s Snakefile.linearCount 
		printf "'location' \n" > exprssRawCount
		awk '{print $1"_"$2}' cellA_count.txt >> exprssRawCount
		chmod +x putTogether.sh
		snakemake -s Snakefile.linearCount putTogether

		Rscript agCells.r {Input_rawCountMatrix} {Out_nomralizedMatrix} {Out_aggregateMatrix} {Out_aggregateingInfoMatrix}

			apply TF-IDF normalization - this reduces the importance of the most common peaks - used by Cusanovich et al 2017
			use clustering based on distance to make representative cell profiles 
			#25-50 cells in aggregation is reasonable based on Zhou et al 2016
			#clustering may amplify the clustering of low read depth cells, check variablity of avg reads per cell
			cluster method:
				while meadian cells in cluster < 30:
					#use median rather than average so as to not overweight superclusters which have only power of 1 cell in downstream analysis
					make tracking table: contains list of sampleName/cellName and number of cells  
					measure distance between all pairs
						sum((differance(cellAsite*, cellBsite*))^2)
					combine 2 closest into cluster (avg weighted by number of cells)
					update tracking table
		
		measure variablity withen samples and compare between for sample heterogenity #to check for batch effects compare withen seperate replicates to between sample types
		cluster w/ t-SNE and plot: install.packages("tsne") and in R environment plot(tsne(matrix, k=2))

		replicate and compare monocle.r and TSCAN analysis on clusteredMatrix
