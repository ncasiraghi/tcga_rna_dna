library(tidyverse)
library(data.table)
library(parallel)

mc.cores <- 10

setwd('/BCGLAB/2020_signatures/synthetic_dilutions/pacbam/outs/')

sif <- read.delim('/BCGLAB/2020_signatures/synthetic_dilutions/input_sif_synthetic_dilutions.tsv',stringsAsFactors = FALSE)

# pacbam

snps_data <- c(list.files('/BCGLAB/2020_signatures/synthetic_dilutions/pacbam/original_data',pattern = '\\.snps$',full.names = TRUE),
               list.files('/BCGLAB/2020_signatures/synthetic_dilutions/pacbam/insilico_data',pattern = '\\.snps$',full.names = TRUE))


df <- c()
for(i in seq_len(nrow(sif))){
  
  case.bam.id <- str_remove(string = basename(sif$plasma.bam[i]),pattern = '\\.bam$')   
  ctrl.bam.id <- str_remove(string = basename(sif$germline.bam[i]),pattern = '\\.bam$')
  
  case.snps <- grep(snps_data, pattern = case.bam.id, value = TRUE)
  ctrl.snps <- grep(snps_data, pattern = ctrl.bam.id, value = TRUE)
  
  tc <- str_remove(basename(case.snps),pattern = case.bam.id) 
  
  tc[1] <- sif$TC[i]
  
  if(length(tc) > 1){
    tc <- tc %>% 
      str_remove(pattern = '^_adm') %>% 
      str_remove(pattern = '\\.sorted\\.snps$') %>%
      str_replace(pattern = '_',replacement = '.') %>% 
      as.numeric() 
  }
  
  h <- data.frame(patient = sif$patient[i],
                   case.bam.id = case.bam.id,
                   ctrl.bam.id = case.bam.id,
                   snps = c(ctrl.snps, case.snps),
                   tc = as.numeric(c(NA,tc)),
                   stringsAsFactors = FALSE) %>%
    mutate(tc = if_else(tc > 1, true = 1-(tc/100),false = tc))
    
  df <- rbind(df,h)
  
}

focus <- df %>% filter(tc < 0.05) %>% pull(patient) %>% unique()

m <- df %>% 
  filter(patient %in% focus) %>% 
  group_by(patient) %>% 
  group_split()

# extract informative snps

getInfoSNPs <- function(i,case.snps,ctrl.snps,patient.id){
  
  snps <- fread(case.snps$snps[i],data.table = FALSE,showProgress = FALSE,nThread = mc.cores) %>% 
    filter(cov >= 10) %>% 
    filter(rsid %in% unique(ctrl.snps$rsid)) %>% 
    select(rsid,af,cov) %>% 
    mutate(af.mir = if_else(condition = af < 0.5,true = 1 - af,false = af)) %>% 
    mutate(patient = patient.id, tc = case.snps$tc[i])
  
  return(snps)
  
}

mt <- c()

for(i in seq_len(length(m))){
  
  ptnt <- m[[i]]
  cat(paste("[",Sys.time(),"]",unique(ptnt$patient),"\n"))
  
  cat(paste("[",Sys.time(),"]","Reading SNPs of the NORMAL sample","\n"))
  ctrl.snps <- ptnt %>% 
    filter(is.na(tc)) %>%
    pull(snps) %>% 
    fread(.,data.table = FALSE,showProgress = FALSE,nThread = mc.cores) %>% 
    filter(cov >= 10) %>% 
    select(rsid,af,cov) %>% 
    filter(af >= 0.1 & af <= 0.9) %>%
    mutate(af.mir = if_else(condition = af < 0.5,true = 1 - af,false = af)) %>% 
    mutate(patient = unique(ptnt$patient), tc = NA)
  
  case.snps <- ptnt %>% filter(!is.na(tc))
  
  cat(paste("[",Sys.time(),"]","Reading SNPs of the TUMOR samples","\n"))
  
  out <- lapply(seq_len(nrow(case.snps)),getInfoSNPs, case.snps=case.snps, ctrl.snps=ctrl.snps, patient.id = unique(ptnt$patient))
  
  cat(paste("[",Sys.time(),"]","Merging together NORMAL and TUMOR data","\n"))
  
  out <- do.call(rbind,out)
  
  mt <- rbind(mt,ctrl.snps,out)
  
}

save(mt,file = 'infosnps.RData',compress = TRUE)

# load('infosnps.RData')

# plots

med <- mt %>% 
  group_by(patient,tc) %>% 
  summarise(median = median(af,na.rm = TRUE))

p <- ggplot(mt, aes(x=af,y = as.factor(tc))) +
  geom_violin() +
  geom_boxplot(width=0.1) +
  facet_wrap(~patient,scales = 'free_y',nrow = 1) + ylab('Tumor content') + theme_bw()

ggsave(filename = 'infosnps.pdf',plot = p, width = 300,height = 150,dpi = 300,units = 'mm',device = 'pdf')

p <- ggplot(mt, aes(x=af.mir,y = as.factor(tc))) +
  geom_violin() +
  geom_boxplot(width=0.1) +
  facet_wrap(~patient,scales = 'free_y',nrow = 1) + ylab('Tumor content') + theme_bw()

ggsave(filename = 'infosnps_AFmir.pdf',plot = p, width = 300,height = 150,dpi = 300,units = 'mm',device = 'pdf')
