#author: Jennifer Havens
#dependancy of Snakemake.count
#normalizes expression matrix with TF-IDF and filters based on threshold


#first argument is file name of input matrix 
#second argument is file name of normalized by output matrix
#normalizaton is done by term frequency inverse document frequency alg, used by Cusanovich et al 2017
	#demphsizes signals which are very common
#rows are genes or peak sites, columns are cells, should have names in input matrix

#then filteres cells based on low accesible sites and filters sites based on their accesablity across cells
#can not do by returning normTbl instead of filterNormTbl


args<-(commandArgs(TRUE))


tbl <- read.table(args[1], head=TRUE, row.names=1)
sink(args[2])

#can change these thresholds as needed
ACCESSTHRESHOLD <- 500
FREQTHRESHOLD <- 0.1

normTbl <- tbl

for (site in row.names(tbl)){
    for(cell in colnames(tbl)){
        peakFreq <- length(tbl[site,])/sum(tbl[site,])
        if (is.infinite(peakFreq)) {peakFreq <- 0}
        normTbl[cell][site,] <- (tbl[cell][site,])*1/(sum(tbl[cell]))*log(1+peakFreq)

    }
}




#filters out cells with too few reads (ACCESSTHRESHOLD)
filterLowRead <- function(expressMat){
  cellList <- colnames(expressMat)
  saveList <- c()
  for (cell in cellList){
    if (sum(expressMat[,cell]) > ACCESSTHRESHOLD){
      saveList <- c(saveList, cell)
    } 
  }
  return(subset(expressMat, select = saveList))
}


#filters the sites with amount of observation below threshold across cells (FREQTHRESHOLD)
filterLowFreq <- function(expressMat){
  useMat <- t(expressMat)
  siteList <- colnames(useMat)
  cellNum <- length(row.names(useMat))
  saveList <- c()
  for (site in siteList){
    if (sum(useMat[,site])/cellNum > FREQTHRESHOLD){
      saveList <- c(saveList, site)
    }
  }
  outMat <- subset(useMat, select = saveList)
  return(t(outMat))
}



filterNormTbl <- filterLowRead(filterLowFreq(normTbl))

filterNormTbl

