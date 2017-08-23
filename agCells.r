#author: Jennifer Havens
#takes matrix of raw counts and returns matrixes of filterend and TF-IDF normalized and a matrix with filtered and normalized cells aggregated into groups
#run after making raw matrix with Snakefile.linearCount
#once done can re-do monocle.r and TSCAN analysis and compare to unag results
#filtering cells based on low accesible sites would have been done already based on outlined pipeline but if alternative pipeline used option to filter is here


args<-(commandArgs(TRUE))
#raw count matrix 1st input, cells in cols, names should be in sample.cell format, will chage to seperate by "/" for interal work, outputs will change back
exprssM <- read.table(args[1], head=TRUE, row.names=1)
colnames(exprssM) <- gsub(".", "/", colnames(exprssM))
#second argument is output for normalized matrix
outNomrMName <- args[2]

#third argument is output for aggregated matrix
outClustMName <- args[3]

#forth argument is output for aggregating information matrix
outClustInfoMName <- args[4]


medClustSize <- 30
ACCESSTHRESHOLD <- 500
FREQTHRESHOLD <- 0.1


#normalizes raw count expression matrix with TF-IDF alg
#returns normalized matrix and saves into outNormM file
tfIDF <- function(expressM){

  normM <- expressM
  
  for (site in row.names(expressM)){
    for(cell in colnames(expressM)){
      peakFreq <- length(expressM[site,])/sum(expressM[site,])
      if (is.infinite(peakFreq)) {peakFreq <- 0}
      normM[site,][cell] <- (expressM[site,][cell])*1/(sum(expressM[,cell]))*log(1+peakFreq)
      
    }
  }
  return(normM)
  
  #outputs normalized matrix with sample.cell naming format
  outNormM <- normM
  colnames(outNormM) <- gsub("/", ".", colnames(outNormM))
  write.table(outNormM, outNomrMName, sep="\t")
  
}

#returns a matrix with clusterNames in col and rows with total number of cells in cluster (all 1) and the fraction of cells from each sample (all 0 or 1)
makeInfoMat <- function(mat){
  #set up matrix structure and names
  names <- strsplit(colnames(mat), "/")
  sampleL <- c()
  for (i in 1: length(names)){
    sampleL<- c(sampleL, names[[i]][1])
  }
  samplesUniq <- unique(sampleL)
  infoMat <- mat.or.vec(1+length(samplesUniq),length(colnames(mat)))
  colnames(infoMat) <- colnames(mat)
  row.names(infoMat) <- c("numberOfCells", samplesUniq)
  
  #populate matrix with correct values
  infoMat["numberOfCells",] <- seq(from=1, to =1, length.out=length(colnames(mat)))
  for (samp in samplesUniq){
    for (i in 1:length(colnames(infoMat))){
      if(sampleL[i]==samp){
        infoMat[samp,][i] <- 1
      }
    }
  }
  return(infoMat)
}

#finds the square of the differance between 2 cells across all sties, returns distance
findSqDis <- function(mat, clustA, clustB){
  disScore <- 0
  for (site in row.names(mat)){
    localDis <- (mat[,clustA][site] - mat[,clustB][site])^2
    disScore <- disScore + localDis
  }
  return(disScore)
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



  
#returns clustered matrix and outputs clustered matrix and information matrix
cluster <- function(clustMat){
  clustMat <- filterLowFreq(clustMat) #optional filtering step, can adjust threshold at top of document
  clustMat <- filterLowRead(clustMat) #optional filtering step, can adjust threshold at top of document
  clustInfoMat <- makeInfoMat(clustMat)
  #for(k in 1:2){
  while(median(clustInfoMat["numberOfCells",])<medClustSize){
    clustNames <- colnames(clustMat)
    totalClusters <- length(clustNames)
    minScore <- Inf
    #comparing distances between all cells
    for (i in 1:totalClusters){
      for (j in i:totalClusters){
        if (j>i){
          score <- findSqDis(clustMat, clustNames[i], clustNames[j])
          if (score < minScore){
            clustA <- clustNames[i]
            clustB <- clustNames[j]
          }
        }
      }
    }
   
    #combinding 2 nearest clusters in info mat
    newCellNumber <- clustInfoMat[,clustA]["numberOfCells"] + clustInfoMat[,clustB]["numberOfCells"]
    clustACellNumber <- clustInfoMat[,clustA]["numberOfCells"]
    clustBCellNumber <- clustInfoMat[,clustB]["numberOfCells"]
    newInfoClust <- clustInfoMat[,clustA]*(clustACellNumber)/newCellNumber +  clustInfoMat[,clustB]*(clustBCellNumber)/newCellNumber
    clustInfoMat <- cbind(clustInfoMat, newInfoClust)
    newName <- paste(clustA, clustB, sep = "_")
    fullNames <- c(clustNames, newName)
    colnames(clustInfoMat) <- fullNames
    clustInfoMat[,newName]["numberOfCells"] <- newCellNumber

    outL <- c()
    for (n in fullNames){if (n != clustB & n != clustA){outL <- c(outL, n)} }

    clustInfoMat <- subset(clustInfoMat, select = outL)
    
    #combinding 2 nearest clusters in expression mat
    newClust <- clustMat[,clustA]*clustACellNumber/newCellNumber +  clustMat[,clustB]*clustBCellNumber/newCellNumber
    clustMat <- cbind(clustMat, newClust)
    colnames(clustMat) <-  fullNames

    clustMat <- subset(clustMat, select = outL)
    
  }
  
  return(clustMat)
  
  #outputs clustered matrix and its info matrix with sample.cell naming format
  outClustM <- clustMat
  colnames(outClustM) <- gsub("/", ".", colnames(outClustM))
  write.table(outClustM, outClustMName, sep="\t")
  outclustInfoMat<- clustInfoMat
  colnames(outclustInfoMat) <- gsub("/", ".", colnames(outclustInfoMat))
  write.table(outclustInfoMat, outClustInfoMName, sep="\t")
}


cluster(exprssM) 

