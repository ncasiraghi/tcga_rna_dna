library(tidyverse)

# build up the sh script to run computeData.R

setwd('/BCGLAB/ncasiraghi/tcga_rna_dna/computeData')

rscript <- 'Rscript /BCGLAB/ncasiraghi/tcga_rna_dna/computeData/computeData.R'

# select samples of interest

sif.file <- '/BCGLAB/ncasiraghi/tcga_rna_dna/data/synggen_dataset.tsv'

sif <- read.delim(file = sif.file,header = TRUE,stringsAsFactors = FALSE)

cat(unique(sif$data.source),sep = '\n')

dataset <- unique(sif$data.source)
omic <- 'DNA'
aggregate_as <- 'mean'
min_vaf <- 0.1
max_vaf <- 0.9
min_cov <- 10
sw <- 50
pstep <- 0.75
cores <- 50
custom_borders <- c(FALSE)

df <- expand.grid(rscript,sif.file,omic,dataset,aggregate_as,min_vaf,max_vaf,min_cov,sw,pstep,custom_borders,cores,stringsAsFactors = FALSE)

colnames(df) <- c('rscript','sif.file','omic','dataset','aggregate_as','min_vaf','max_vaf','min_cov','sw','pstep','custom_borders','cores')

write.table(df, file = 'auto_computeData.sh',quote = FALSE,sep = ' ',col.names = FALSE,row.names = FALSE)
