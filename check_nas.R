library(tidyverse)
library(data.table)

lf <- list.files('/BCGLAB/2020_signatures/stats',pattern = '_Pstep_100000',full.names = TRUE)

alldata <- list()

for(i in seq_len(length(lf))){
  
  message(lf[i])
  
  mat <- fread(file = file.path(lf[i],'data.tsv'),data.table = FALSE,nThread = 2)
  
  pos <- read.delim(file = file.path(lf[i],'positions.tsv'),stringsAsFactors = F,header = TRUE)
  
  ids <- readLines(file.path(lf[i],'samples.txt'))
  
  tumor  <- grep(ids,pattern = 'tumor$')
  normal <- grep(ids,pattern = 'normal$')
  
  pos <- pos %>% 
    mutate(na_tumor_count  = colSums(is.na(mat[tumor,])),
           na_normal_count = colSums(is.na(mat[normal,]))) %>% 
    mutate(na_tumor_fraction  = round(na_tumor_count/length(tumor),2),
           na_normal_fraction = round(na_normal_count/length(normal),2))
  
  alldata[basename(lf[i])] <- pos
  
}

save(alldata,file = 'alldata_nas.RData',compress = TRUE)