#author: Jennifer Havens
#dependancy of Snakemake.count
#normalizes expression matrix with TF-IDF


#first argument is file name of input matrix 
#second argument is file name of normalized by output matrix
#normalizaton is done by term frequency–inverse document frequency alg, used by Cusanovich et al 2017
	#demphsizes signals which are overlly common
#rows are genes or peak sites, columns are cells, should have names in input matrix



args<-(commandArgs(TRUE))


tbl <- read.table(args[1], head=TRUE, row.names=1)
sink(args[2])

normTbl <- tbl

for (site in row.names(tbl)){
    for(cell in colnames(tbl)){
        peakFreq <- length(tbl[site,])/sum(tbl[site,])
        if (is.infinite(peakFreq)) {peakFreq <- 0}
        normTbl[cell][site,] <- (tbl[cell][site,])*1/(sum(tbl[cell]))*log(1+peakFreq)

    }
}

normTbl

