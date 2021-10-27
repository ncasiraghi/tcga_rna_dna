library(tidyverse)
library(data.table)

setwd('/BCGLAB/2020_signatures/synthetic_dilutions/pacbam/outs/tcga_cbioportal_snvs')

sif <- read.delim('/BCGLAB/2020_signatures/synthetic_dilutions/input_sif_synthetic_dilutions.tsv',stringsAsFactors = FALSE) %>% mutate(id = str_remove(plasma,pattern = '_DNA_tumor'))

# SNVs from TCGA

snvs <- read.delim('/BCGLAB/2020_signatures/synthetic_dilutions/pacbam/outs/tcga_cbioportal_snvs/data_mutations_extended.txt',stringsAsFactors = FALSE) %>% 
  filter(Tumor_Sample_Barcode %in% unique(sif$id)) %>% 
  filter(Variant_Type == 'SNP') %>% 
  filter(grepl(Variant_Classification,pattern = '_Mutation')) %>% 
  select(Chromosome,Start_Position,End_Position,Reference_Allele,Tumor_Seq_Allele2,Tumor_Sample_Barcode,NCBI_Build) %>% 
  mutate(Chromosome = paste0('chr',Chromosome)) %>% 
  unite("id",Chromosome,Start_Position,End_Position,remove = FALSE)

bed <- snvs %>% 
  select(Chromosome:End_Position,Chromosome,id) %>% 
  mutate(Start_Position = Start_Position - 1) 

write.table(bed,file = 'tcga_brca_snvs_GRCh37.bed',quote = FALSE,col.names = FALSE,row.names = FALSE,sep = '\t')

tcga_snvs <- read.delim("tcga_brca_snvs_GRCh38.bed",stringsAsFactors = FALSE,col.names = c('Chromosome','Start_Position_GRCh38','End_Position_GRCh38','id'),header = FALSE)
  
tcga_snvs <- left_join(tcga_snvs,snvs,by = c('Chromosome','id'))

VCF <- data.frame(chr = tcga_snvs$Chromosome,
                  POS = tcga_snvs$End_Position_GRCh38,
                  ID = '.',
                  REF = tcga_snvs$Reference_Allele,
                  ALT = tcga_snvs$Tumor_Seq_Allele2,
                  QUAL = '.',
                  FILTER = '.',
                  INFO = '.',
                  stringsAsFactors = FALSE) 

colnames(VCF)[1] <- '#CHROM'

write.table(VCF,file = 'tcga_brca_snvs_GRCh38.vcf',quote = FALSE,col.names = TRUE,row.names = FALSE,sep = '\t')

# perform pacbam on in silico dilutions