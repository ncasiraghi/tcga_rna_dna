library(tidyverse)

setwd('/BCGLAB/ncasiraghi/tcga_rna_dna/data')

datasource <- list.files('/BCGLAB/rscandino/2020_signatures/synggen_BRCA/Workflow/results',full.names = TRUE)

input_data_folder <- '/BCGLAB/2020_signatures/input_data'

df <- c()

for (dt in datasource) {
  
  message(dt)
  
  snps <- list.files(dt,full.names = TRUE,pattern = '\\.snps$',recursive = TRUE)
  
  m <- data.frame(data.source = paste0('synggen_',basename(dt)),
                  file.snps = snps,
                  sample.id = paste0(basename(dirname(snps)),'_',str_remove(basename(snps),pattern = '\\.sorted\\.snps$')),
                  ref.genome = 'GRCh38',
                  data.type = 'DNA',
                  out.folder = file.path(input_data_folder,paste0('synggen_',basename(dt))),
                  stringsAsFactors = FALSE)
  
  df <- rbind(df,m)
  
}

write.table(df,file = 'synggen_dataset.tsv',quote = FALSE,col.names = TRUE,row.names = FALSE,sep = '\t')