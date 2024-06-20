


#### Setup section ####
#Setup variables
pathToOrig <- "~/Experiments/Tripogon/TS_all_filtered.recode_wide1.tsv"

#### Prepare for recoding ####
#Load data
dfOrig <- read.delim(file = pathToOrig)
head(dfOrig)

#Remove empty columns
dfMiss <- dfOrig=="./."|dfOrig=="."
hist(rowSums(dfMiss),1000)
dfOrig <- dfOrig[,colMeans(dfMiss)<1]

#Create subset of data for recoding
colsToAjhust <- grep("SAMPLE",colnames(dfOrig))+1
dfGenotypes  <- as.matrix(dfOrig[,colsToAjhust])
dfRecoded    <- dfGenotypes
dfRecoded[,] <-"NN"

#Generate logic for controlling recoding
dfHomo1 <- dfGenotypes=="0/0" ; hist(rowSums(dfHomo1),1000)
dfHet   <- dfGenotypes=="0/1" ; hist(rowSums(dfHet  ),1000)
dfHomo2 <- dfGenotypes=="1/1" ; hist(rowSums(dfHomo2),1000)
homo1 <- paste0(dfOrig$REF,dfOrig$REF)
het <- paste0(dfOrig$REF,dfOrig$ALT)
homo2<-paste0(dfOrig$ALT,dfOrig$ALT)

#Alphabetize het
unihet <- sort(unique(het))
sortedHet <- unlist(lapply(strsplit(unihet,split = ""),function(x){paste0(sort(x),collapse = "")}))
for(i in which(unihet!=sortedHet)){
  het[het==unihet[i]]<-sortedHet[i]
}

#### Recode genotypes #####
unihomo1 <- unique(homo1)
unihet   <- unique(het)
unihomo2 <- unique(homo2)
for(i in 1:length(unihomo1)){
  currRows <- matrix(homo1==unihomo1[i],nrow = nrow(dfRecoded),ncol=ncol(dfRecoded))
  dfRecoded[currRows&dfHomo1]<-unihomo1[i]
  print(i)
}
for(i in 1:length(unihet)){
  currRows <- matrix(het==unihet[i],nrow = nrow(dfRecoded),ncol=ncol(dfRecoded))
  dfRecoded[currRows&dfHet]<-unihet[i]
  print(i)
}
for(i in 1:length(unihomo2)){
  currRows <- matrix(homo2==unihomo2[i],nrow = nrow(dfRecoded),ncol=ncol(dfRecoded))
  dfRecoded[currRows&dfHomo2]<-unihomo2[i]
  print(i)
}

#### Generate a subset of loci for concatentation ####
## Subset by missing call frequency. 
missCheck <- rowMeans(dfRecoded=="NN")==0
## Subset by proximity
minDiff <- 30000
dfSub <- dfOrig[missCheck,]
dfSub$posBin1 <- paste0(dfSub$CHROM,floor(dfSub$POS/minDiff))
dfSub <- dfSub[!duplicated(dfSub$posBin1),]
dfSub$posBin2 <- paste0(dfSub$CHROM,floor((dfSub$POS+minDiff/2)/minDiff))
dfSub <- dfSub[!duplicated(dfSub$posBin2),]
posCheck2 <- c(TRUE,(diff(dfSub$POS)>minDiff)|diff(as.numeric(as.factor(dfSub$CHROM)))!=0)
dfSub <- dfSub[posCheck2,]
hist(diff(dfSub$POS)[diff(as.numeric(as.factor(dfSub$CHROM)))==0],100000,xlim=c(0,minDiff*10))
## Subset by number of markers on chr
chrCheck <- dfSub$CHROM%in%names(table(dfSub$CHROM))[table(dfSub$CHROM)>100]
dfSub <- dfSub[chrCheck,]
##Isolate the recoded values
dfRecodedSub <- dfRecoded[(1:nrow(dfRecoded))%in%rownames(dfSub),]

#Use the recoded genotpyes to create alignable sequence
sequences <-apply(dfRecodedSub,2,paste0,collapse=""); nchar(sequences)
sequences <- unlist(as.vector(rbind(paste0(">",names(sequences)),sequences))); nchar(sequences)
write(sequences,file = paste0(gsub("\\.tsv$","",pathToOrig),".fasta"))

#### Phylogentic analysis  ####
# I downloaded Mega11
# I imported the file as a nucleotide sequence with "N" as the missing indicator
# It is non-protein coding sequence
# It is already aligned, so I jumped to calculating a phylogeny
# I used maximum likelihood with bootstrapping x1000
# Nei-tamura method
# Uniform rates
# Use all sites
# Tree building:  NNI
# Make initital tree automatically
#### Setup section ####
#Setup variables
pathToOrig <- "~/Experiments/Tripogon/TS_all_filtered.recode_wide1.tsv"

#### Prepare for recoding ####
#Load data
dfOrig <- read.delim(file = pathToOrig)

#Remove empty columns
dfMiss <- dfOrig=="./."|dfOrig=="."
hist(rowSums(dfMiss),1000)
dfOrig <- dfOrig[,colMeans(dfMiss)<1]

#Create subset of data for recoding
colsToAjhust <- grep("SAMPLE",colnames(dfOrig))+1
dfGenotypes  <- as.matrix(dfOrig[,colsToAjhust])
dfRecoded    <- dfGenotypes
dfRecoded[,] <-"NN"

#Generate logic for controlling recoding
dfHomo1 <- dfGenotypes=="0/0" ; hist(rowSums(dfHomo1),1000)
dfHet   <- dfGenotypes=="0/1" ; hist(rowSums(dfHet  ),1000)
dfHomo2 <- dfGenotypes=="1/1" ; hist(rowSums(dfHomo2),1000)
homo1 <- paste0(dfOrig$REF,dfOrig$REF)
het <- paste0(dfOrig$REF,dfOrig$ALT)
homo2<-paste0(dfOrig$ALT,dfOrig$ALT)

#Alphabetize het
unihet <- sort(unique(het))
sortedHet <- unlist(lapply(strsplit(unihet,split = ""),function(x){paste0(sort(x),collapse = "")}))
for(i in which(unihet!=sortedHet)){
  het[het==unihet[i]]<-sortedHet[i]
}

#### Recode genotypes #####
unihomo1 <- unique(homo1)
unihet   <- unique(het)
unihomo2 <- unique(homo2)
for(i in 1:length(unihomo1)){
  currRows <- matrix(homo1==unihomo1[i],nrow = nrow(dfRecoded),ncol=ncol(dfRecoded))
  dfRecoded[currRows&dfHomo1]<-unihomo1[i]
  print(i)
}
for(i in 1:length(unihet)){
  currRows <- matrix(het==unihet[i],nrow = nrow(dfRecoded),ncol=ncol(dfRecoded))
  dfRecoded[currRows&dfHet]<-unihet[i]
  print(i)
}
for(i in 1:length(unihomo2)){
  currRows <- matrix(homo2==unihomo2[i],nrow = nrow(dfRecoded),ncol=ncol(dfRecoded))
  dfRecoded[currRows&dfHomo2]<-unihomo2[i]
  print(i)
}

#### Generate a subset of loci for concatentation ####
## Subset by missing call frequency. 
missCheck <- rowMeans(dfRecoded=="NN")==0
## Subset by proximity
minDiff <- 30000
dfSub <- dfOrig[missCheck,]
dfSub$posBin1 <- paste0(dfSub$CHROM,floor(dfSub$POS/minDiff))
dfSub <- dfSub[!duplicated(dfSub$posBin1),]
dfSub$posBin2 <- paste0(dfSub$CHROM,floor((dfSub$POS+minDiff/2)/minDiff))
dfSub <- dfSub[!duplicated(dfSub$posBin2),]
posCheck2 <- c(TRUE,(diff(dfSub$POS)>minDiff)|diff(as.numeric(as.factor(dfSub$CHROM)))!=0)
dfSub <- dfSub[posCheck2,]
hist(diff(dfSub$POS)[diff(as.numeric(as.factor(dfSub$CHROM)))==0],100000,xlim=c(0,minDiff*10))
## Subset by number of markers on chr
chrCheck <- dfSub$CHROM%in%names(table(dfSub$CHROM))[table(dfSub$CHROM)>100]
dfSub <- dfSub[chrCheck,]
##Isolate the recoded values
dfRecodedSub <- dfRecoded[(1:nrow(dfRecoded))%in%rownames(dfSub),]

#Use the recoded genotpyes to create alignable sequence
sequences <-apply(dfRecodedSub,2,paste0,collapse=""); nchar(sequences)
sequences <- unlist(as.vector(rbind(paste0(">",names(sequences)),sequences))); nchar(sequences)
write(sequences,file = paste0(gsub("\\.tsv$","",pathToOrig),".fasta"))

#### Phylogentic analysis  ####
# I downloaded Mega11
# I imported the file as a nucleotide sequence with "N" as the missing indicator
# It is non-protein coding sequence
# It is already aligned, so I jumped to calculating a phylogeny
# I used maximum likelihood with bootstrapping x1000
# Nei-tamura method
# Uniform rates
# Use all sites
# Tree building:  NNI
# Make initital tree automatically
