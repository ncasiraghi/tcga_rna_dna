library(data.table)
library(tidyverse)
library(parallel)

wd <- '/BCGLAB/2020_signatures/stats/modified_Pstep'

setwd(wd)

samples.in.parallel <- 10

lf <- list.files('/BCGLAB/2020_signatures/pileup/TCGA-BRCA',pattern = '_DNA',full.names = TRUE)

bands <- read.delim('/BCGLAB/ncasiraghi/gene_panel/input/cytobands/ucsc_cytobands_GRCh38p13.tsv',stringsAsFactors = F)

bands$arm <- 'p'
bands$arm[grep(bands$name,pattern = 'q')] <- 'q'

bands <- bands %>% 
  filter(!chrom %in% c('chrY','chrM')) %>% 
  group_by(chrom,arm) %>% 
  summarise(start = min(chromStart), end = max(chromEnd))

SnpSelect <- function(i, sl, bands){
  
  df <- sl[[i]]
  df$arm <- NA
  
  chr <- unique(df$chr)
  
  p.coord <- bands %>% filter(chrom == chr, arm == 'p')
  df$arm[which( df$pos >= p.coord$start & df$pos <= p.coord$end )] <- 'p'
  
  q.coord <- bands %>% filter(chrom == chr, arm == 'q')
  df$arm[which( df$pos >= q.coord$start & df$pos <= q.coord$end )] <- 'q'
  
  out <- df %>%
    group_by(arm) %>%
    slice(which.min(pos),which.max(pos)) %>% 
    mutate(check = min(pos))
  
  out$flag <- 'last'
  out$flag[which(out$pos == out$check)] <- 'first'
  
  return(out %>% select(-check))
  
}

getArmsMargin <- function(i,lf,bands,min_vaf,max_vaf,min_cov){
  
  file <- lf[i]
  
  snps <- fread(file,data.table = FALSE,nThread = 5,verbose = FALSE) %>% 
    filter(af >= min_vaf, af <= max_vaf, cov >= min_cov)
  
  if(str_detect(snps$chr[1],'chr',negate = TRUE)){
    snps <- snps %>% 
      mutate(chr = paste0('chr',chr)) %>% 
      filter(chr %in% paste0('chr',1:22))
  }else{
    snps <- snps %>% 
      filter(chr %in% paste0('chr',1:22))
  }
  
  tomi <- which(snps$af < 0.5)
  snps$af[tomi] <- (1 - snps$af[tomi])
  
  sl <- snps %>% 
    group_by(chr) %>% 
    group_split()
  
  deck <- do.call(rbind,lapply(seq_len(length(sl)), SnpSelect, sl, bands))
  
  return(deck)
  
}

listmat <- mclapply(seq_len(length(lf)),getArmsMargin,lf=lf,bands=bands,min_vaf=min_vaf,max_vaf=max_vaf,min_cov=min_cov,mc.preschedule = TRUE,mc.cores = samples.in.parallel)

snps_borders <- do.call(rbind,listmat)

save(snps_borders,file = 'snps_borders.RData')

ab_first <- snps_borders %>%
  filter(flag == 'first') %>% 
  group_by(chr, arm) %>% 
  slice(which.max(pos))

ab_last <- snps_borders %>%
  filter(flag == 'last') %>% 
  group_by(chr, arm) %>% 
  slice(which.min(pos))

arm_borders <- rbind(ab_first,ab_last) %>% arrange(chr,arm,pos)

write.table(arm_borders,file = 'arm_borders.tsv',sep = '\t',row.names = F,quote = F,col.names = TRUE)
