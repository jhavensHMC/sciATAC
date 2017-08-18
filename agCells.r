#author: Jennifer Havens
#takes matrix of raw counts and returns matrixes TF-IDF normalized and a matrix with cells aggregated into groups
#run after making raw matrix with Snakefile.linearCount
#once done can re-do monocle.r and TSCAN analysis and compare to unag results

args<-(commandArgs(TRUE))
#raw count matrix 1st input, cells in cols, names should be in sample.cell format, will chage to seperate by "/" for interal work, outputs will change back
exprssM <- read.table(args[1], head=TRUE, row.names=1)
colnames(exprssM) <- gsub("[.]", "/", colnames(exprssM))
#second argument is output for normalized matrix
outNomrMName <- args[2]

#third argument is output for aggregated matrix
outClustMName <- args[3]

#forth argument is output for aggregating information matrix
outClustInfoMName <- args[4]


medClustSize <- 30

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
  
  outputs normalized matrix with sample.cell naming format
  outNormM <- normM
  colnames(outNormM) <- gsub("/", "[.]", colnames(outNormM))
  write.table(outNormM, outNomrMName, sep="\t")
  
}

#returns a matrix with clusterNames in col and rows with total number of cells in cluster (all 1) and the fraction of cells from each sample (all 0 or 1)
makeInfoMat <- function(mat){
  #set up matrix structure and names
  names <- strsplit(colnames(mat), "_")
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

#finds the differance between 2 cells across all sties, returns distance
findDis <- function(mat, clustA, clustB){
  disScore <- 0
  for (site in row.names(mat)){
    localDis <- (mat[,clustA][site] - mat[,clustB][site])^2
    disScore <- disScore + localDis
  }
  return(disScore)
}


  
#returns clustered matrix and outputs clustered matrix and information matrix
cluster <- function(clustMat){
  clustInfoMat <- makeInfoMat(clustMat)
  while(median(clustInfoMat["numberOfCells",])<medClustSize){
    clustNames <- colnames(clustMat)
    totalClusters <- length(clustNames)
    minScore <- Inf
    #comparing distances between all cells
    for (i in 1:totalClusters){
      for (j in i:totalClusters){
        if (j>i){
          score <- findDis(clustMat, clustNames[i], clustNames[j])
          if (score < minScore){
            clustA <- clustNames[i]
            clustB <- clustNames[j]
          }
        }
      }
    }
    
    
    #combinding 2 nearest clusters
    newClust <- (clustMat[,clustA]*clustInfoMat[,clustA]["numberOfCells"] +  clustMat[,clustB]*clustInfoMat[,clustB]["numberOfCells"])*0.5
    newNames <- c(clustNames, paste(clustA, clustB, sep = "_"))
    clustMat <- cbind(clustMat, newClust)
    colnames(clustMat) <- newNames

    outL <- c()
    for (n in newNames){if (n != clustB & n != clustA){outL <- c(outL, n)} }
    clustMat <- subset(clustMat, select = outL)


    newClust <- clustInfoMat[,clustA]*0.5 +  clustInfoMat[,clustB]*0.5
    clustInfoMat <- cbind(clustInfoMat, newClust)
    newName <- paste(clustA, clustB, sep = "_")
    colnames(clustInfoMat) <- c(clustNames, newName)
    clustInfoMat[,newName]["numberOfCells"] <- 2*clustInfoMat[,newName]["numberOfCells"]
    clustInfoMat <- subset(clustInfoMat, select = outL)

  }
  return(clustMat)
  
  outputs clustered matrix and its info matrix with sample.cell naming format
  outClustM <- clustMat
  colnames(outClustM) <- gsub("/", "[.]", colnames(outClustM))
  write.table(outClustM, outClustMName, sep="\t")
  outclustInfoMat<- clustInfoMat
  colnames(outclustInfoMat) <- gsub("/", "[.]", colnames(outclustInfoMat))
  write.table(outclustInfoMat, outClustInfoMName, sep="\t")
}


cluster(exprssM)

