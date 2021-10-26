library(dplyr)

setwd('/BCGLAB/ncasiraghi/tcga_rna_dna/Pstep_modified_weighted')

rscript <- 'Rscript /BCGLAB/ncasiraghi/tcga_rna_dna/Pstep_modified_weighted/computeData.R'

sif.file <- '/BCGLAB/ncasiraghi/tcga_rna_dna/data/full_dataset.tsv'

sif <- read.delim(file = sif.file,header = TRUE,stringsAsFactors = FALSE)

cat(unique(sif$data.source),sep = '\n')

dataset <- grep(unique(sif$data.source),pattern = 'WES-Plasma',value = TRUE)
omic <- 'DNA'
aggregate_as <- 'mean'
min_vaf <- 0.3
max_vaf <- 0.7
min_cov <- 10
sw <- 50
pstep <- 0.75
cores <- 50
custom_borders <- c(TRUE,FALSE)

df <- expand.grid(rscript,sif.file,omic,dataset,aggregate_as,min_vaf,max_vaf,min_cov,sw,pstep,custom_borders,cores,stringsAsFactors = FALSE)

colnames(df) <- c('rscript','sif.file','omic','dataset','aggregate_as','min_vaf','max_vaf','min_cov','sw','pstep','custom_borders','cores')

write.table(df, file = 'auto_computeData.sh',quote = FALSE,sep = ' ',col.names = FALSE,row.names = FALSE)
