library(tidyverse)

setwd('/BCGLAB/ncasiraghi/tcga_rna_dna/data')

datasource <- list.files('/BCGLAB/2020_signatures/pileup',full.names = TRUE)

input_data_folder <- '/BCGLAB/2020_signatures/input_data'

df <- c()

for (dt in datasource) {
  
  message(dt)
  
  snps <- list.files(dt,full.names = TRUE,pattern = '\\.snps$',recursive = TRUE)
  
  if(str_detect(dt,pattern = 'TCGA')){
    
    m <- data.frame(data.source = basename(dt),
                    file.snps = snps,
                    sample.id = basename(dirname(snps)),
                    stringsAsFactors = FALSE) %>% 
      mutate(ref.genome = 'GRCh38') %>% 
      mutate(data.type = case_when(str_detect(sample.id,pattern = '_DNA_') ~ 'DNA',
                                   str_detect(sample.id,pattern = '_RNA_') ~ 'RNA')) %>% 
      mutate(out.folder = file.path(input_data_folder,basename(dt)))
    
  }
  
  if(str_detect(dt,pattern = 'SU2C')){
    
    m <- data.frame(data.source = basename(dt),
                    file.snps = snps,
                    sample.id = gsub(basename(snps),pattern = '\\.snps$',replacement = ''),
                    stringsAsFactors = FALSE) %>% 
      mutate(ref.genome = 'GRCh37') %>% 
      mutate(data.type = 'DNA') %>% 
      mutate(out.folder = file.path(input_data_folder,basename(dt)))
    
  }
  
  if(str_detect(dt,pattern = 'WES-Plasma')){
    
    m <- data.frame(data.source = paste0(basename(dt),'_',basename(dirname(snps))),
                    file.snps = snps,
                    sample.id = gsub(basename(snps),pattern = '\\.snps$',replacement = ''),
                    stringsAsFactors = FALSE) %>% 
      mutate(ref.genome = 'GRCh37') %>% 
      mutate(data.type = 'DNA') %>% 
      mutate(out.folder = file.path(input_data_folder,paste0(basename(dt),'_',basename(dirname(file.snps)))))
    
  }

  df <- rbind(df,m)
  
}

write.table(df,file = 'datasetSNPs.tsv',quote = FALSE,sep = '\t',col.names = TRUE,row.names = FALSE)

