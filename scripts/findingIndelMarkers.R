#Setup analysis
setwd("~/Experiments/Tspicata/")
indelFile <- "./data/data_ignored/Tspicata_both.d8.merged_indelSubset.csv"

#Load data
df1 <- read.csv(indelFile)
colnames(df1)[colnames(df1)=="X.CHROM"]<-"CHROM"
dfOrig<-df1

#Recalculate alt allele length
maxAlleleDiff <- function(x){
  y <-nchar(unlist(strsplit(x,split = ",")))
  maxInd <- which.max(y)
  zMax <- y[maxInd]
  zMin <- max(y[-maxInd])
  if(length(y)==1){
    stop("Error!")
  }else{
    return(zMax-zMin)
  }
}
alleleVec <- paste0(df1$REF,",",df1$ALT)
df1$alleleDiff <- sapply(alleleVec,maxAlleleDiff)

#Set filters
opt<-list()
opt$notSamples <- c('X.CHROM','CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','altLen','kept')
opt$startCol <- 10
opt$maxMiss <- 0.5
opt$minHomo <- 0.8
opt$minAlleleDiff <- 10
opt$outDir  <- "./data/"
opt$outName <- file.path(opt$outDir,paste0(
  gsub("\\.csv$","",basename(indelFile)),
  "_lt",opt$maxMiss,"Miss",
  "_gt",opt$minAlleleDiff-1,"bpDiff",
  ".csv"
))

##Filtering calculations
dfCalls <- df1[,!colnames(dfOrig)%in%opt$notSamples]
dfCalls[dfCalls=="./."]<-NA
df1$missFreq  <- rowMeans(is.na(dfCalls))
hetGenos     <- apply(combn(0:4,2),2,paste,collapse="/")
df1$hetFreq <- 1-rowMeans(!apply(dfCalls,2,`%in%`,hetGenos))

updateKept <- function(prevLogic,logicFun){
  afterLogic <- prevLogic & logicFun
  print(sum(afterLogic))
  return(afterLogic)
}

## Reset
df1$kept1 <- T
## Checking missing freq
df1$missLogic <- df1$missFreq<opt$maxMiss
df1$kept1 <- updateKept(df1$kept1,logicFun = df1$missLogic)
## Check min allele length
df1$minLenLogic <- df1$altLen>=opt$minAlleleDiff
df1$kept1 <- updateKept(df1$kept1,logicFun = df1$minLenLogic)

## Check for biallelic
df1$biallelic <- !grepl(",",df1$ALT)
df1$kept1 <- updateKept(df1$kept1,logicFun = df1$biallelic)

#Out name
write.csv(x = df1[df1$kept1,],file = opt$outName)

