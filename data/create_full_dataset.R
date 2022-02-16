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
      mutate(out.folder = file.path(input_data_folder,basename(dt))) %>% 
      mutate(is.tumor = case_when(str_detect(sample.id,pattern = 'tumor') ~ 1,
                                  str_detect(sample.id,pattern = 'normal') ~ 0))
  }
  
  if(str_detect(dt,pattern = 'SU2C')){
    
    tumr.bams <- readLines("/BCGLAB/2020_signatures/pileup/SU2C/BAMs_Tumor.txt") %>% 
      unique() %>% 
      basename() %>% 
      str_remove(pattern = '\\.bam$')
      
    # ctrl.bams <- readLines("/BCGLAB/2020_signatures/pileup/SU2C/BAMs_Normal.txt") %>%
    #   unique() %>% 
    #   basename() %>% 
    #   str_remove(pattern = '\\.bam$')
    
    m <- data.frame(data.source = basename(dt),
                    file.snps = snps,
                    sample.id = gsub(basename(snps),pattern = '\\.snps$',replacement = ''),
                    stringsAsFactors = FALSE) %>% 
      mutate(ref.genome = 'GRCh37') %>% 
      mutate(data.type = 'DNA') %>% 
      mutate(out.folder = file.path(input_data_folder,basename(dt))) %>% 
      mutate(is.tumor = if_else( sample.id %in% tumr.bams, 1, 0))
    
  }
  
  if(str_detect(dt,pattern = 'WES-Plasma')){
    
    m <- data.frame(data.source = paste0(basename(dt),'_',basename(dirname(snps))),
                    file.snps = snps,
                    sample.id = gsub(basename(snps),pattern = '\\.snps$',replacement = ''),
                    stringsAsFactors = FALSE) %>% 
      mutate(ref.genome = 'GRCh37') %>% 
      mutate(data.type = 'DNA') %>% 
      mutate(out.folder = file.path(input_data_folder,paste0(basename(dt),'_',basename(dirname(file.snps))))) %>% 
      mutate(is.tumor = case_when(str_detect(sample.id,pattern = 'germline') ~ 0,
                                  str_detect(sample.id,pattern = 'Ctrl') ~ 0)) %>% 
      mutate(is.tumor = replace_na(is.tumor, replace = 1))
    
  }

  df <- rbind(df,m)
  
}

write.table(df,file = 'datasetSNPs.tsv',quote = FALSE,sep = '\t',col.names = TRUE,row.names = FALSE)

#  stats 

outdir <- '/BCGLAB/2020_signatures/stats/samples_count'

bp <- df %>%
  filter(data.type == "DNA") %>% 
  group_by(data.source, ref.genome, is.tumor) %>% 
  summarise(n = n()) 

write.table(bp,file = file.path(outdir,'table_nsamples.tsv'),quote = FALSE,col.names = TRUE,row.names = FALSE,sep = '\t')

p <- ggplot(data=bp, aes(x=data.source, y=n, fill=as.factor(is.tumor))) +
  geom_bar(stat="identity", position=position_dodge()) +
  coord_flip() + theme_bw() 
  
ggsave(filename = file.path(outdir,'barplot_nsamples.pdf'),plot = p,width = 300,height = 200,dpi = 300,units = 'mm',device = 'pdf')
