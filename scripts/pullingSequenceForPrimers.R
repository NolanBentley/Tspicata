BiocManager::install("Biostrings") 
library("Biostrings") #

ref <- Biostrings::readDNAStringSet(filepath = "./data/data_ignored/T_spicatus_V2.2.fasta")

chrVec <- c('Contig_2','Contig_4','Contig_6','Contig_3','Contig_3','Contig_7','Contig_8')
posVec <- c(20598628,6451385,13247919,1507108,3618088,14456053,36153)
buffer <- 150
for(i in 1:length(chrVec)){
  currChr <-chrVec[i]
  currPos <- posVec[i]
  currSeq <- as.character(subseq(x = ref[[currChr]],currPos-buffer,currPos+buffer))
  outSeq <- paste0(substr(currSeq,1,buffer),
                   "[",substr(currSeq,buffer+1,buffer +1),"]",
                   substr(currSeq,buffer+2,buffer*2+1))
  if(i==1){outVec <- NULL}
  outVec <- c(outVec,paste0(">",currChr,"_",currPos),outSeq)
}

write(outVec,file="./data/indelSeqs_targettingRef_20240620.fasta")