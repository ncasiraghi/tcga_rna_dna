library(dplyr)

setwd('/BCGLAB/ncasiraghi/tcga_rna_dna/Pstep_modified_weighted')

pileup <- list.files('/BCGLAB/2020_signatures/pileup',full.names = TRUE)

rscript <- 'Rscript /BCGLAB/ncasiraghi/tcga_rna_dna/Pstep_modified_weighted/computeData.R'
pileup <- grep(pileup,pattern = 'TCGA',value = TRUE)
type <- c('DNA','RNA')
aggregate_as <- c('mean','median')
min_vaf <- 0.1
max_vaf <- 0.9
min_cov <- c(0,10)
sw <- c(25,50,100)
pstep <- c(0.25,0.75,1)
n <- 50

df <- expand.grid(rscript,pileup,type,aggregate_as,min_vaf,max_vaf,min_cov,sw,pstep,n,stringsAsFactors = FALSE)

colnames(df) <- c('rscript','pileup','type','aggregate_as','min_vaf','max_vaf','min_cov','sw','pstep','n')

m <- df %>% 
  filter(!grepl('-BRCA',pileup)) %>% 
  filter(!grepl('-OV',pileup)) %>% 
  filter(type == 'DNA') %>% 
  filter(pstep == 0.75) %>% 
  filter(aggregate_as == 'mean') %>%
  filter(sw == 50) %>% 
  arrange(pileup)

m$out <- file.path('/BCGLAB/2020_signatures/input_data',basename(m$pileup))

write.table(m, file = 'auto_computeData.sh',quote = FALSE,sep = ' ',col.names = FALSE,row.names = FALSE)