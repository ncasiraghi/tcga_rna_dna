library(data.table)
library(tidyverse)
library(parallel)

options(dplyr.summarise.inform = FALSE)
options(scipen = 9999)

wd <- '/BCGLAB/2020_signatures/stats/'

setwd(wd)

# params

mc.cores <- 24

length.sw <- 10

sample_type <- 'DNA_tumor'

pos_step <- 500000
  
# outs

outdir <- paste0('main_out_afs_sw',length.sw,'_',sample_type,'_Pstep_',format(pos_step,scientific = FALSE))

if(!file.exists(file.path(wd,outdir))){
  dir.create(file.path(wd, outdir), showWarnings = FALSE)
} else{
  unlink(file.path(wd,outdir),recursive = TRUE)
  dir.create(file.path(wd, outdir), showWarnings = FALSE)
}

# pileup on snps

lf <- list.files('/BCGLAB/2020_signatures/pileup/TCGA-BRCA',pattern = paste0('_',sample_type,'$'),full.names = TRUE)

# cytobands

bands <- read.delim('/BCGLAB/ncasiraghi/gene_panel/input/cytobands/ucsc_cytobands_GRCh38p13.tsv',stringsAsFactors = F)

bands$arm <- 'p'
bands$arm[grep(bands$name,pattern = 'q')] <- 'q'

bands <- bands %>% 
  group_by(chrom,arm) %>% 
  summarise(start = min(chromStart), end = max(chromEnd))

# Positions to aggregate data

positions <- c()

for(i in seq_len(nrow(bands))){
  
  x <- data.frame(chrom=bands$chrom[i],
                  arm=bands$arm[i],
                  pos=seq(from=bands$start[i],bands$end[i],by = pos_step),
                  stringsAsFactors = FALSE)
  
  positions <- rbind(positions,x)
  
}

write.table(positions,file = file.path(wd,outdir,'positions.tsv'),col.names = TRUE,quote = FALSE,row.names = FALSE,sep = '\t')

# run sliding window

for (id in lf) {
  message(paste("[",Sys.time(),"]\t",id))
  
  file <- list.files(id,full.names = TRUE,pattern = '\\.snps$')
  
  snps <- fread(file,data.table = FALSE,nThread = mc.cores) 
  
  sl <- snps %>% 
    filter(af >= 0.2, af <= 0.8) %>% 
    group_by(chr) %>% 
    group_split()
  
  runsw <- function(i, sl, bands, length.sw, positions){
    
    df <- sl[[i]]
    df$arm <- NA
    
    chr <- unique(df$chr)
    
    p.coord <- bands %>% filter(chrom == chr, arm == 'p')
    df$arm[which( df$pos >= p.coord$start & df$pos <= p.coord$end )] <- 'p'
    
    q.coord <- bands %>% filter(chrom == chr, arm == 'q')
    df$arm[which( df$pos >= q.coord$start & df$pos <= q.coord$end )] <- 'q'
    
    stats_by_arm <- function(df,chr,which.arm,n,positions){
      
      adf <- df %>% filter(arm == which.arm)
      
      if( nrow(adf) == 0 ){
        
        stats <- data.frame(chr=chr,
                            arm=which.arm,
                            wnd_start=NA,
                            wnd_end=NA,
                            median_af=NA,
                            mean_af=NA,
                            stringsAsFactors = FALSE)
        
      } else {
        
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
                            wnd_start=min(m$pos),
                            wnd_end=max(m$pos),
                            median_af=median(m$af,na.rm = TRUE),
                            mean_af=mean(m$af,na.rm = TRUE),
                            stringsAsFactors = FALSE)
          
          stats <- rbind(stats,out)
          
        }
        
      }
      
      afs <- c()
      
      pos.arm <- as.numeric(positions %>% filter(arm == which.arm, chrom == chr ) %>% pull(pos))
      
      for(px in pos.arm){
        
        x <- stats %>% 
          filter(wnd_start <= px & wnd_end >=px) %>% 
          pull(median_af)
        
        afs <- c(afs,median(x,na.rm = TRUE))
        
      }
      
      return(afs)
      
    }
    
    out <- c(stats_by_arm(df,chr = chr, which.arm = 'p',n = length.sw, positions = positions),
             stats_by_arm(df,chr = chr, which.arm = 'q',n = length.sw, positions = positions))
    
    return(out)

  }
  
  deck <- mclapply(seq_len(length(sl)), runsw, sl, bands, length.sw = length.sw, positions = positions, mc.cores = mc.cores)
  
  write(as.numeric(unlist(deck)),file = file.path(wd,outdir,'data.tsv'),sep = '\t',append = TRUE,ncolumns = length(unlist(deck))) 
  
  write(basename(id),file = file.path(wd,outdir,'samples.txt'),sep = '\t',append = TRUE)

}
