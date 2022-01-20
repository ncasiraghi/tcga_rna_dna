library(tidyverse)
library(data.table)

setwd('/Users/ncasiraghi/Desktop')

file.seg <- '/BCGLAB/Nicolo/Tesi/point_mutation/brca_tcga_pan_can_atlas_2018/data_cna_hg19.seg'

ids <- c('TCGA-A2-A04R-01',
         'TCGA-BH-A1FL-01',
         'TCGA-D8-A1XA-01',
         'TCGA-E2-A2P6-01')

seg <- read.delim(file.seg,stringsAsFactors = FALSE) %>% filter(ID %in% ids)

th.ampl <- 0.6
th.left.del <- -1.5
th.right.del <- -0.3
th.left.diploid <- -0.1
th.right.diploid <- 0.2

min.cov <- 10

p <- ggplot(seg, aes(x=seg.mean)) + 
  geom_histogram(bins = 100) +
  facet_wrap(~ID,nrow = 4) +
  geom_vline(xintercept = th.ampl, color = "red") +
  geom_vline(xintercept = c(th.left.del,th.right.del), color = "blue") +
  geom_vline(xintercept = c(th.left.diploid,th.right.diploid), color = "forestgreen")
  
ggsave(filename = 'hist_seg.mean.pdf',plot = p,width = 300,height = 200,dpi = 200,units = 'mm',device = 'pdf')

# convert SEG hg19 into hg38

bed <- seg %>% 
  unite(info,ID:seg.mean,remove = FALSE,sep = '_') %>% 
  select(chrom,loc.start,loc.end,info) %>% 
  mutate(chrom = paste0('chr',chrom))

write.table(bed,file = 'hg19_seg.bed',quote = FALSE,sep = '\t',col.names = FALSE,row.names = FALSE)

# BED hg38

bed <- read.delim('hg38_seg.bed',stringsAsFactors = FALSE,header = FALSE,col.names = c('chrom','loc.start','loc.end','info'))

seg <- bed %>% 
  separate(col = info,into = c('ID','chr','start','end','num.mark','seg.mean'),sep = '_') %>% 
  select(ID,chrom,loc.start,loc.end,num.mark,seg.mean) %>% 
  mutate(seg.mean = as.numeric(seg.mean))

write.table(seg,file = 'hg38.seg',quote = FALSE,sep = '\t',col.names = FALSE,row.names = FALSE)

# pileups

ps <- list.files('/BCGLAB/2020_signatures/pileup/TCGA-BRCA',full.names = TRUE)

for(smpl in ids){
  message(smpl)
  
  h <- grep(ps,pattern = str_remove(smpl,pattern = '-01'),value = TRUE) %>% grep(pattern = 'DNA',value = TRUE)
  
  ctrl <- fread(list.files(grep(h,pattern = 'normal$',value = TRUE),full.names = TRUE,pattern = '\\.snps$'),data.table = FALSE,stringsAsFactors = FALSE,header = TRUE)
  case <- fread(list.files(grep(h,pattern = 'tumor$',value = TRUE),full.names = TRUE,pattern = '\\.snps$'),data.table = FALSE,stringsAsFactors = FALSE,header = TRUE)

  ctrl <- ctrl %>% 
    filter(af > 0.2 & af < 0.8)
  
  case <- case %>% 
    filter(rsid %in% unique(ctrl$rsid))
  
  # select based on segments
  
  diploid <- seg %>% 
    filter(ID == smpl) %>% 
    filter(seg.mean < th.right.diploid & seg.mean > th.left.diploid)
  
  del <- seg %>% 
    filter(ID == smpl) %>% 
    filter(seg.mean > th.left.del & seg.mean < th.right.del)
  
  ampl <- seg %>% 
    filter(ID == smpl) %>% 
    filter(seg.mean > th.ampl)

  int_seeg_snp <- function(i,seg,snps){
    
    m <- snps %>% 
      filter(chr == seg[i,'chrom']) %>% 
      filter(pos >= seg[i,'loc.start'] & pos < seg[i,'loc.end'])
    
    return(m)
    
  }
  
  snps.ctrl.diploid <- do.call(rbind,lapply(seq_len(nrow(diploid)),int_seeg_snp,seg = diploid, snps = ctrl)) %>% mutate(class = 'snps.ctrl.diploid')
  snps.case.diploid <- do.call(rbind,lapply(seq_len(nrow(diploid)),int_seeg_snp,seg = diploid, snps = case)) %>% mutate(class = 'snps.case.diploid')
  snps.case.del <- do.call(rbind,lapply(seq_len(nrow(del)),int_seeg_snp,seg = del, snps = case)) %>% mutate(class = 'snps.case.del')
  snps.case.ampl <- do.call(rbind,lapply(seq_len(nrow(ampl)),int_seeg_snp,seg = ampl, snps = case)) %>% mutate(class = 'snps.case.ampl')
  
  df <- rbind(snps.ctrl.diploid,snps.case.diploid,snps.case.del,snps.case.ampl) %>% 
    filter(cov > min.cov) %>% 
    mutate(af = if_else(af < 0.5,true = 1 - af,false = af))
  
  p <- ggplot(df, aes(x=class, y=af)) + 
    geom_boxplot(varwidth = TRUE) +
    ggtitle(smpl)
  
  ggsave(filename = paste0(smpl,'.pdf'),plot = p,width = 300,height = 200,dpi = 200,units = 'mm',device = 'pdf')
  
}
