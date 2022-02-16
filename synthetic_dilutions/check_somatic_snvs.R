library(tidyverse)
library(data.table)

setwd('/BCGLAB/2020_signatures/synthetic_dilutions/pacbam/outs/tcga_cbioportal_snvs')

sif <- read.delim('/BCGLAB/2020_signatures/synthetic_dilutions/input_sif_synthetic_dilutions.tsv',stringsAsFactors = FALSE) %>% mutate(id = str_remove(plasma,pattern = '_DNA_tumor'))

# cat(sif$plasma.bam,file = '/BCGLAB/2020_signatures/synthetic_dilutions/pacbam/original_data/tcga_cbioportal_snvs/bamlist.txt',sep = '\n')

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
  mutate(Start_Position = Start_Position - 1) %>% 
  arrange(Chromosome,Start_Position,End_Position)

# write.table(bed,file = 'tcga_brca_snvs_GRCh37.bed',quote = FALSE,col.names = FALSE,row.names = FALSE,sep = '\t') 

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
                  stringsAsFactors = FALSE) %>% 
  arrange(chr,POS)

colnames(VCF)[1] <- '#CHROM'

# write.table(VCF,file = 'tcga_brca_snvs_GRCh38.vcf',quote = FALSE,col.names = TRUE,row.names = FALSE,sep = '\t')

# perform pacbam on in silico dilutions

file.snps <- c(list.files('/BCGLAB/2020_signatures/synthetic_dilutions/pacbam/insilico_data/tcga_cbioportal_snvs',pattern = '\\.snps$',full.names = TRUE),
               list.files('/BCGLAB/2020_signatures/synthetic_dilutions/pacbam/original_data/tcga_cbioportal_snvs',pattern = '\\.snps$',full.names = TRUE))
               
m <- c()
for(pp in file.snps){
  
  h <- unlist(strsplit(basename(pp),split = '_adm'))
  
  if(length(h) == 1){
    pattern <- str_remove(h,pattern = '.snps')
    tc <- NA
  } else {
    pattern <- str_remove(h[1],pattern = '.snps')
    tc <- 1 - as.numeric(str_replace(str_remove(h[2],pattern='.sorted.snps'),pattern = '_',replacement = '.'))/100 
  }

  message(pattern)
  
  tt <- tcga_snvs %>% 
    filter(Tumor_Sample_Barcode == sif$id[grep(sif$plasma.bam,pattern = pattern)]) %>% 
    unite("vcf",Chromosome,End_Position_GRCh38)
  
  if(nrow(tt) > 0){
    
    if(is.na(tc)){
      tc <- sif %>% filter(id == unique(tt$Tumor_Sample_Barcode)) %>% pull(TC)
    }
    
    x <- read.delim(pp,header = TRUE,stringsAsFactors = FALSE) %>% 
      unite("vcf",chr,pos) %>% 
      filter(vcf %in% tt$vcf) %>% 
      mutate(tumor.sample = unique(tt$Tumor_Sample_Barcode)) %>% 
      mutate(tumor.content = tc)
    
    m <- rbind(m,x)

  }
  
}

# plot

focus <- m %>% 
  group_by(tumor.sample) %>% 
  summarise(ntc = n_distinct(tumor.content)) %>% 
  filter(ntc == 4) %>% 
  pull(tumor.sample)

p <- ggplot(m %>% filter(tumor.sample %in% focus), aes(x=af,y = as.factor(tumor.content))) +
  # geom_violin() +
  # geom_boxplot(width=0.1) +
  geom_boxplot() + 
  facet_wrap(~tumor.sample,scales = 'free_y',nrow = 1) + ylab('Tumor content') + theme_bw() + ggtitle('somatic SNVs')

ggsave(filename = 'somaticsnvs.pdf',plot = p, width = 300,height = 150,dpi = 300,units = 'mm',device = 'pdf')




