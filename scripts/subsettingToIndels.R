#Setup analysis
setwd("~/Experiments/Tspicata/")
vcfFile      <- "./data/data_ignored/Tspicata_both.d8.merged.vcf"
indelFileOut <- paste0(gsub("\\.vcf","",vcfFile),"_indelSubset.csv")

#Load data
library(data.table)
headerRow <- grep("#CHROM",readLines(vcfFile,n = 1000))
vcf <- fread(vcfFile,skip = headerRow-1)

#Remove SNPs
vcf$refLen <- nchar(vcf$REF)
vcf$altLen <- nchar(vcf$ALT)
vcf<-vcf[vcf$altLen!=vcf$refLen,]

#Save file
write.table(x = vcf,indelFileOut,row.names = F,sep = ",")

