#author: Jennifer Havens
#Pseudotime analysis of sci-ATAC-seq data using monocle
#Run after Snakefile.count and phenoSheet.sh


source("http://bioconductor.org/biocLite.R")
biocLite("monocle")
library(monocle)
biocLite(c("DDRTree", "pheatmap"))

#load data 
#they recomend FPKM or TPM if not using UMI, but we dont need that we can use raw read counts or FPKM/TPM


#exprs, a numeric matrix of expression values, where rows are genes, and columns are cells
#phenoData, an AnnotatedDataFrame object, where rows are cells, and columns are cell attributes (such as cell type, culture condition, day captured, etc.)
#featureData, an AnnotatedDataFrame object, where rows are features (e.g. genes), and columns are gene attributes, such as biotype, gc content, etc.

expr_matrix <- read.table("fpkm_matrix.normalized.txt") #rows are gene names, columns are cells: make with tfIDF.r
pheno_sheet <- read.delim("pheno_sheet.txt") # rows cells, columns are attibutes (sample and sites (counted unique reads)) #make with phenoSheet.sh
#gene_annotation <- read.delim("gene_annotations.txt") #suggested to add if there are gene atributes you want to look at

pd <- new("AnnotatedDataFrame", data = pheno_sheet) #
#fd <- new("AnnotatedDataFrame", data = gene_annotation)
cds <- newCellDataSet(as(expr_matrix, "sparseMatrix"), phenoData = pd, expressionFamily=negbinomial.size())


#classify cells - skipping, look into projects where people are identifing 


#order cells in pseudotime - currently unsupervised

#finding relevant genes
expressed_genes <- row.names(subset(fData(cds), num_cells_expressed >= 10))

#not use sites in model to  bias created by different assay efficiency in different cells - Pliner et al 2017
full_model_fits <- fitModel(cds[expressed_genes,],  modelFormulaStr = "~sample + sites")
reduced_model_fits <- fitModel(cds[expressed_genes,], modelFormulaStr = "~sites")
diff_test_res <- compareModels(full_model_fits, reduced_model_fits)


ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))

cdsOrdering <- setOrderingFilter(cds, ordering_genes)
plot_ordering_genes(cdsOrdering)


#reduce dimention
cdsOrdering <- reduceDimension(cdsOrdering, max_components = 2, method = 'DDRTree')

#order and plot

cdsOrdering <- orderCells(cdsOrdering)

plot_cell_trajectory(cdsOrdering, color_by = "sample")

plot_cell_trajectory(cdsOrdering, color_by = "Pseudotime")


#cluster cells
cds <- clusterCells(cds)

#test diff expression, between samples
full_model_fits <- fitModel(cds_subset,  modelFormulaStr = "~sample + sites")
reduced_model_fits <- fitModel(cds_subset, modelFormulaStr = "~sites")
diff_test_res_sample <- compareModels(full_model_fits, reduced_model_fits) #this will take a while to test everything
sig_genes_sample <- subset(diff_test_res, qval < 0.1)


#test diff expression, with Pseudotime
full_model_fits <- fitModel(cds_subset,  modelFormulaStr = "~sm.ns(Pseudotime) + sites") 
reduced_model_fits <- fitModel(cds_subset, modelFormulaStr = "~sites")
diff_test_res_time <- compareModels(full_model_fits, reduced_model_fits)
sig_genes_time <- subset(diff_test_res_time, qval < 0.1)

sig_names_time <- row.names(sig_genes_time)
cds_subset <- cds[sig_names_time] #not sure about this syntax
plot_genes_in_pseudotime(cds_subset, color_by = "sample")

plot_pseudotime_heatmap(cds[sig_names_time,],
                        num_clusters = 3,
                        cores = 1,
                        show_rownames = T)