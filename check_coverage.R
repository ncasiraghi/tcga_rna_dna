library(tidyverse)
library(data.table)

sif <- read.delim('/BCGLAB/2020_signatures/stats/annotations/TCGA-BRCA/BRCA_samples_info_file.txt') 

dil <- sif %>%
  filter(note == 'dilution' & batch_approved == 1)

pileup_data <- '/BCGLAB/2020_signatures/pileup/TCGA-BRCA/'

df <- c()

for(i in seq_len(nrow(dil))){
  
  message(dil$sample_bcglab[i])
  
  snps.file <- list.files(path = file.path(pileup_data,dil$sample_bcglab[i]),pattern = '\\.snps$',full.names = TRUE)
  
  if(length(snps.file) > 0){
    
    snps <- fread(snps.file,data.table = FALSE,showProgress = FALSE,nThread = 5)
    
    snps <- snps %>% 
      select(af,cov) %>% 
      mutate(file = 'snps', sample = dil$sample_bcglab[i],type = dil$type[i])
    
    df <- rbind(df,snps)
    
  }
  
  
  pabs.file <- list.files(path = file.path(pileup_data,dil$sample_bcglab[i]),pattern = '\\.pabs$',full.names = TRUE)
  
  if(length(pabs.file) > 0){
    
    pabs <- fread(pabs.file,data.table = FALSE,showProgress = FALSE,nThread = 5)
    
    pabs <- pabs %>% 
      select(af,cov) %>% 
      mutate(file = 'pabs', sample = dil$sample_bcglab[i],type = dil$type[i])
    
    df <- rbind(df,pabs)
    
  }
  

}

df <- df %>% mutate(type = ifelse(type == 1, 'tumor', 'normal'))

# stats

df %>% 
  filter(cov > 10) %>% 
  filter(af >= 0.1 & af <= 0.9) %>% 
  group_by(type,file) %>% 
  summarise(mean_cov = mean(cov), median_cov = median(cov)) 

df %>% 
  filter(file == 'snps') %>% 
  filter(cov > 10) %>% 
  filter(af >= 0.1 & af <= 0.9) %>% 
  group_by(sample,file,type) %>% 
  summarise(mean_cov = mean(cov), median_cov = median(cov)) %>% 
  ggplot(.,aes(x = sample,y = mean_cov, fill = type)) +
  geom_bar(stat="identity") +
  facet_wrap(~type,scales = 'free_x') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

tp <- df %>% 
  filter(cov > 10) %>% 
  filter(af >= 0.1 & af <= 0.9)

ggplot(tp, aes(x=file, y=cov, fill = type)) + 
  geom_boxplot() + coord_cartesian(ylim = c(0,200))

ggplot(tp, aes(x=file, y=af, fill = type)) + 
  geom_boxplot() 

