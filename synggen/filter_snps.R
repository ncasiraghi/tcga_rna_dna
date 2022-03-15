library(data.table)
library(tidyverse)

setwd('/BCGLAB/2020_signatures/synggen_inputs/pacbam')

pb <- list.files('/BCGLAB/2020_signatures/synggen_inputs/pacbam/snps_original/',pattern = '\\.snps$',full.names = TRUE)

for(snp in pb){
  
  message(snp)
  
  m <- fread(snp,data.table = FALSE,nThread = 20)
  
  mf <- m %>% filter(cov > 0 & af > 0)
  
  write.table(x = mf,file = basename(snp),quote = FALSE,sep = '\t',col.names = TRUE,row.names = FALSE)
  
}
