library(data.table)
library(tidyverse)

setwd('/BCGLAB/2020_signatures/synggen_inputs')

ids <- read.delim('/BCGLAB/2020_signatures/synggen_inputs/ids_spice_tcga_brca.txt',stringsAsFactors = F,header = F) %>% select(c(1,3,6,7))

colnames(ids) <- c('spice_id','tcga_id','clonet_purity','clonet_ploidy')

head(ids)

sif <- read.delim('/BCGLAB/2020_signatures/synthetic_dilutions/input_sif_synthetic_dilutions.tsv') %>% pull(patient) %>% unique()

ids <- ids %>% filter(tcga_id %in% sif)

seg <- fread(input = '/BCGLAB/2020_signatures/synggen_inputs/brca.seg',data.table = FALSE,stringsAsFactors = FALSE,colClasses = list(numeric=2),nThread = 4) %>% distinct()

head(seg)

m <- seg %>% 
  filter(!is.na(log2_corr)) %>% 
  filter(sample_id %in% unique(ids$spice_id)) %>%
  select(-sample_id) %>%
  separate(segment_id,into = c('spice_id','chrom','start','end')) %>% 
  select(c('spice_id','chrom','start','end','as_cn_disc','log2_corr')) %>% 
  mutate(start = as.numeric(start),end = as.numeric(end)) %>% 
  arrange(spice_id,chrom,start,end)

df <- merge(x = m,y = ids,by = 'spice_id',all.x = TRUE)

p <- ggplot(df, aes(x=log2_corr)) + 
  geom_histogram(bins = 100) +
  facet_wrap(~as_cn_disc,scales = 'free')

ggsave(filename = 'log2corr_by_cna.pdf',plot = p,width = 300,height = 200,dpi = 200,units = 'mm',device = 'pdf')
