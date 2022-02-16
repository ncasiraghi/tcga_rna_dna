library(tidyverse)
library(data.table)

setwd('/BCGLAB/ncasiraghi/tcga_rna_dna/data')

sif <- read.delim('/BCGLAB/2020_signatures/synthetic_dilutions/input_sif_synthetic_dilutions.tsv',stringsAsFactors = FALSE)

# pacbam

snps_data <- list.files('/BCGLAB/2020_signatures/synthetic_dilutions/pacbam/insilico_data',pattern = '\\.snps$',full.names = TRUE)

m <- c()
for(fs in snps_data){
  
  message(basename(fs))
  
  sample.id <- sif$plasma[grep(sif$plasma.bam,pattern = unlist(str_split(string = basename(fs),pattern = '_adm'))[1])]
  
  tc <- unlist(str_split(string = basename(fs),pattern = '_adm'))[2] %>% 
    str_remove(pattern = '\\.sorted\\.snps$') %>%
    str_replace(pattern = '_',replacement = '.') %>% 
    as.numeric()
  
  out <- data.frame(data.source='TCGA-BRCA-insilico',
                    file.snps=fs,
                    sample.id=paste(sample.id,1-(tc/100),sep = '_'),
                    ref.genome='GRCh38',
                    data.type='DNA',
                    out.folder='/BCGLAB/2020_signatures/input_data/TCGA-BRCA-insilico',
                    tumor.content=1-(tc/100),
                    stringsAsFactors = FALSE)
  m <- rbind(m,out)
  
}

write.table(m,file = 'insilico_dataset.tsv',col.names = TRUE,row.names = FALSE,quote = FALSE,sep = '\t')