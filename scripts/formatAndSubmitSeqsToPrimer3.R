#Setup
setwd("~/Experiments/Tspicata/")
faFile <- "data/indelSeqs_targettingRef_20240620.fasta"
p3CorePath <- "~/primer3/src/primer3_core"
p3InputDir <- "data/primer3/inputs"
p3OutputDir <- "data/primer3/outputs"

#Read into memory
fa <- readLines(faFile)

#Get headers
fa_headers_rows <- grep("^>",fa)
fa_headers <- fa[fa_headers_rows]

#Get and collapse seqs
seqPosList <- mapply(seq,from=fa_headers_rows+1,to=c(fa_headers_rows[-1]-1,length(fa)))
fa_Seqs <- lapply(seqPosList,function(x,f){paste0(f[x],collapse = "")},f=fa)
fa_targets <- data.frame(
  start=unlist(lapply(gregexpr("\\[",fa_Seqs),min)),
  stop=unlist(lapply(gregexpr("\\]",fa_Seqs),max))-2
)
fa_targets$len <- fa_targets$stop-fa_targets$start+1
fa_Seqs <- gsub("\\[|\\]","",fa_Seqs)


#Make primer 3 input
primer3File <- paste0("SEQUENCE_ID=",fa_headers,"
SEQUENCE_TEMPLATE=",fa_Seqs,"
SEQUENCE_TARGET=",fa_targets$start,",",fa_targets$len,"
PRIMER_TASK=generic
PRIMER_PICK_LEFT_PRIMER=1
PRIMER_PICK_RIGHT_PRIMER=1
PRIMER_OPT_SIZE=20
PRIMER_MIN_SIZE=18
PRIMER_MAX_SIZE=23
PRIMER_PRODUCT_SIZE_RANGE=150-250 100-300 75-350
PRIMER_EXPLAIN_FLAG=0
="
)
dir.create(p3InputDir,recursive = T)
inp_filenames <- file.path(p3InputDir,paste0(gsub(">","",fa_headers),"_",1:length(primer3File),".inp"))
for(i in 1:length(primer3File)){
  write(primer3File[i], file = inp_filenames[i])
}

dir.create(p3OutputDir,recursive = T)
out_filenames <- file.path(p3OutputDir,paste0(gsub(">","",fa_headers),"_",1:length(primer3File),".out"))
outputCode <- paste0(p3CorePath," --output=",out_filenames," ",inp_filenames)
sapply(outputCode,system)
