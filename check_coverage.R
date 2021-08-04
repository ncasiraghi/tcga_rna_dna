library(tidyverse)
library(data.table)

setwd('/BCGLAB/2020_signatures/synthetic_dilutions/coverage')

sif <- read.delim('/BCGLAB/2020_signatures/annotations/TCGA-BRCA/BRCA_samples_info_file.txt') 

dil <- sif %>%
  filter(note == 'dilution' & batch_approved == 1)

pileup_data <- '/BCGLAB/2020_signatures/synthetic_dilutions/pacbam'

df <- c()

for(i in seq_len(nrow(dil))){
  
  message(dil$filename[i])
  
  snps.file <- file.path(pileup_data,gsub(dil$filename[i],pattern = '\\.bam$',replacement = '.snps'))
  
  snps <- fread(snps.file,data.table = FALSE,showProgress = FALSE,nThread = 5,nrows = 100)
    
  snps <- snps %>% 
    select(af,cov) %>% 
    mutate(file = 'snps', sample = dil$sample_bcglab[i],type = dil$type[i])
  
  df <- rbind(df,snps)
  
  pabs.file <- file.path(pileup_data,gsub(dil$filename[i],pattern = '\\.bam$',replacement = '.pabs'))
  
  pabs <- fread(pabs.file,data.table = FALSE,showProgress = FALSE,nThread = 5,nrows = 100)
  
  pabs <- pabs %>% 
    select(af,cov) %>% 
    mutate(file = 'pabs', sample = dil$sample_bcglab[i],type = dil$type[i])
  
  df <- rbind(df,pabs)
  
}

df <- df %>% mutate(type = ifelse(type == 1, 'tumor', 'normal'))

# stats

df %>% 
  filter(cov > 10) %>% 
  filter(af >= 0.1 & af <= 0.9) %>% 
  group_by(type,file) %>% 
  summarise(mean_cov = mean(cov), median_cov = median(cov)) 

tab <- df %>% 
  filter(file == 'snps') %>% 
  filter(cov > 10) %>% 
  filter(af >= 0.1 & af <= 0.9) %>% 
  group_by(sample,file,type) %>% 
  summarise(mean_cov = mean(cov), median_cov = median(cov))

write.table(tab,file = 'coverage_per_sample.tsv',quote = FALSE,row.names = FALSE,col.names = TRUE,sep = '\t')

p <- ggplot(tab,aes(x = sample,y = mean_cov, fill = type)) +
  geom_bar(stat="identity") +
  facet_wrap(~type,scales = 'free_x') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(filename = 'coverage_per_sample.pdf', plot = p, width = 210,height = 150,dpi = 300,units = 'mm',device = 'pdf')
