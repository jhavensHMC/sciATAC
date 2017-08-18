#author: Jennifer Havens
#Motif analysis of peaks from aggregated sci-ATAC-seq using chromVar	
#run after Snakefile.preprocessing 		

devtools::install_github("GreenleafLab/chromVAR")
source("https://bioconductor.org/biocLite.R")
biocLite("JASPAR2016")
devtools::install_github("GreenleafLab/motifmatchr")
library(BiocParallel)
register(MulticoreParam(16)) 
library(BSgenome.Hsapiens.UCSC.hg19)
library(pheatmap)


#note need to change smaple names, alignment file name, and sampleList  

#combine all peak files into one allPeaks.bed
peakfile <- "allPeaks.bed"
peaks <- getPeaks(peakfile)

bamfiles <- c("sampleA.bam","sampleB.bam") #list of all sample alignment files
sampleList <- c("A","B")
fragment_counts <- getCounts(bamfiles, peaks, 
                              paired =  TRUE, 
                              by_rg = TRUE, 
                              format = "bam",
                              colData = DataFrame(Sample = sampleList))



fragment_counts_GC <- addGCBias(fragment_counts, 
                              genome = BSgenome.Hsapiens.UCSC.hg19)
counts_filtered <- filterSamples(fragment_counts_GC, min_depth = 500,
                                  min_in_peaks = 0.15)
counts_filtered <- filterPeaks(counts_filtered)

motifs <- getJasparMotifs()
motif_ix <- matchMotifs(motifs, counts_filtered,
                         genome = BSgenome.Hsapiens.UCSC.hg19)


#Annotation process

jaspar_motifs <- getJasparMotifs()

motif_ix <- matchMotifs(jaspar_motifs, counts_filtered, genome = BSgenome.Hsapiens.UCSC.hg19)

kmer_ix <- matchKmers(6, counts_filtered, genome = BSgenome.Hsapiens.UCSC.hg19)

cis_ix <- getCisGroups(counts_filtered, grpsize = 25, stepsize = 10) 


#Compute Deviations

dev <- computeDeviations(object = counts_filtered, annotations = motif_ix)

bg <- getBackgroundPeaks(object = counts_filtered)

dev <- computeDeviations(object = counts_filtered, annotations = motif_ix,
                         background_peaks = bg)

#Variablity
variability <- computeVariability(dev)

plotVariability(variability, use_plotly = FALSE) 

#Visulaizing Deviations t-SNE
tsne_results <- deviationsTsne(dev, threshold = 1.5, perplexity = 10)

tsne_plots <- plotDeviationsTsne(dev, tsne_results, annotation = "TEAD3", 
                                   sample_column = "Sample", shiny = FALSE)


#with heatmap
sample_cor <- getSampleCorrelation(dev)

pheatmap(as.dist(sample_cor), 
         annotation_row = colData(dev), 
         clustering_distance_rows = as.dist(1-sample_cor), 
         clustering_distance_cols = as.dist(1-sample_cor))


diff_acc <- differentialDeviations(dev, "Sample")
diff_var <- differentialVariability(dev, "Sample")
