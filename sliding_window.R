library(data.table)
library(tidyverse)
library(parallel)

options(dplyr.summarise.inform = FALSE)

setwd('/BCGLAB/2020_signatures/stats/')

# params

mc.cores <- 24

length.sw <- 10

sample_type <- 'DNA_tumor'
  
# pileup on snps

lf <- list.files('/BCGLAB/2020_signatures/pileup/TCGA-BRCA',pattern = paste0('_',sample_type,'$'),full.names = TRUE)

# cytobands

bands <- read.delim('/BCGLAB/ncasiraghi/gene_panel/input/cytobands/ucsc_cytobands_GRCh38p13.tsv',stringsAsFactors = F)

bands$arm <- 'p'
bands$arm[grep(bands$name,pattern = 'q')] <- 'q'

bands <- bands %>% 
  group_by(chrom,arm) %>% 
  summarise(start = min(chromStart), end = max(chromEnd))

# run sliding window

main_out <- list()

for (id in lf) {
  message(paste("[",Sys.time(),"]\t",id))
  
  file <- list.files(id,full.names = TRUE,pattern = '\\.snps$')
  
  snps <- fread(file,data.table = FALSE,nThread = 5) 
  
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
    
    stats_by_arm <- function(df,chr,which.arm,n){
      
      adf <- df %>% filter(arm == which.arm)
      
      if(nrow(adf)==0){
        out <- data.frame(chr=chr,
                          arm=which.arm,
                          bp = NA,
                          stringsAsFactors = FALSE)
        return(out)
      }
      
      stop <- nrow(adf) - (n-1)
      
      if(stop <= 0){
        stop <- 1
        n <- nrow(adf)
      }
      
      stats <- c()
      for(indx in seq_len(stop)){
        
        m <- adf[indx:( indx + (n-1) ),]
        
        out <- data.frame(chr=chr,
                          arm=unique(adf$arm),
                          bp = max(m$pos) - min(m$pos),
                          stringsAsFactors = FALSE)
        
        stats <- rbind(stats,out)
        
      }
      
      return(stats)
      
    }
    
    tmp <- rbind(stats_by_arm(df,chr = chr, which.arm = 'p',n = length.sw), 
                 stats_by_arm(df,chr = chr, which.arm = 'q',n = length.sw))
    
    out <- tmp %>% 
      filter(!is.na(chr)) %>% 
      group_by(chr,arm) %>% 
      summarise(n_windows=n(),
                median_bp = median(bp,na.rm = T),
                mean_bp = mean(bp,na.rm = T),
                sd_bp = sd(bp,na.rm = T))
    
    return(out)

  }
  
  deck <- mclapply(seq_len(length(sl)), runsw, sl, bands, length.sw = length.sw, mc.cores = mc.cores)
  
  main_out[[basename(id)]] <- do.call(rbind,deck)

}

save(main_out,file = paste0('main_out_sw',length.sw,'_',sample_type,'.RData'),compress = TRUE)