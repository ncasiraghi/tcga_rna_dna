library(data.table)
library(tidyverse)
library(parallel)

setwd('/BCGLAB/2020_signatures/stats/')

# pileup on snps

lf <- list.files('/BCGLAB/2020_signatures/pileup/TCGA-BRCA/',pattern = '_DNA_tumor$',full.names = TRUE)

# cytobands

bands <- read.delim('/BCGLAB/ncasiraghi/gene_panel/input/cytobands/ucsc_cytobands_GRCh38p13.tsv',stringsAsFactors = F)

bands$arm <- 'p'
bands$arm[grep(bands$name,pattern = 'q')] <- 'q'

bands <- bands %>% 
  group_by(chrom,arm) %>% 
  summarise(start = min(chromStart), end = max(chromEnd))

# run sliding window

main_out <- list()

for (id in lf[1:10]) {
  message(id)
  file <- list.files(id,full.names = TRUE,pattern = '\\.snps$')
  
  snps <- fread(file,data.table = FALSE,nThread = 4) 
  
  sl <- snps %>% 
    filter(af >= 0.2, af <= 0.8) %>% 
    group_by(chr) %>% 
    group_split()
  
  runsw <- function(i, sl, bands, length.sw){
    
    df <- sl[[i]]
    df$arm <- NA
    
    chr <- unique(df$chr)
    
    p.coord <- bands %>% filter(chrom == chr, arm == 'p')
    df$arm[which( df$pos >= p.coord$start & df$pos <= p.coord$end )] <- 'p'
    
    q.coord <- bands %>% filter(chrom == chr, arm == 'q')
    df$arm[which( df$pos >= q.coord$start & df$pos <= q.coord$end )] <- 'q'
    
    stats_by_arm <- function(df,which.arm,n){
      
      adf <- df %>% filter(arm == which.arm)
      
      if(nrow(adf)==0){
        return(NA)
      }
      
      stop <- nrow(adf) - (n-1)
      
      if(stop < 0){
        stop <- 1
        n <- nrow(adf)
      }
      
      stats <- c()
      for(indx in seq_len(stop)){
        
        m <- adf[indx:( indx + (n-1) ),]
        
        out <- data.frame(chr=chr,
                          arm=unique(adf$arm),
                          bp = max(m$pos) - min(m$pos),
                          median_af = median(m$af,na.rm = TRUE),
                          median_cov = median(m$cov,na.rm = TRUE),
                          stringsAsFactors = FALSE)
        
        stats <- rbind(stats,out)
        
      }
      
      return(stats)
      
    }
    
    return(rbind(stats_by_arm(df,which.arm = 'p',n = length.sw),
                 stats_by_arm(df,which.arm = 'q',n = length.sw)))

  }
  
  deck <- mclapply(seq_len(length(sl)), runsw, sl, bands, length.sw = 10,mc.cores = 4)
  
  main_out[[basename(id)]] <- deck

}

save(main_out,file = 'main_out.RData',compress = TRUE)