library(data.table)
library(tidyverse)
library(parallel)

options(dplyr.summarise.inform = FALSE)
options(scipen = 9999)

source('/BCGLAB/ncasiraghi/tcga_rna_dna/functions_sliding_window_stats.R')

args <- commandArgs(trailingOnly = TRUE)

wd <- args[1]

length.sw <- as.numeric(args[2])

samples.in.parallel <- as.numeric(args[3])

setwd(wd)

# outs

outdir <- paste0('sw_stats_wnd_',length.sw)

if(!file.exists(file.path(wd,outdir))){
  dir.create(file.path(wd, outdir), showWarnings = FALSE)
} else{
  unlink(file.path(wd,outdir),recursive = TRUE)
  dir.create(file.path(wd, outdir), showWarnings = FALSE)
}

# pileup on snps

lf <- list.files('/BCGLAB/2020_signatures/pileup/TCGA-BRCA',pattern = '_DNA',full.names = TRUE)

# cytobands

bands <- read.delim('/BCGLAB/ncasiraghi/gene_panel/input/cytobands/ucsc_cytobands_GRCh38p13.tsv',stringsAsFactors = F)

bands$arm <- 'p'
bands$arm[grep(bands$name,pattern = 'q')] <- 'q'

bands <- bands %>% 
  filter(!chrom %in% c('chrY','chrM')) %>% 
  group_by(chrom,arm) %>% 
  summarise(start = min(chromStart), end = max(chromEnd))

# run sliding window

getMatrix <- function(i,lf,bands,length.sw){
  
  id <- lf[i]
  
  file <- list.files(id,full.names = TRUE,pattern = '\\.snps$')
  
  snps <- fread(file,data.table = FALSE,verbose = FALSE) %>% 
    filter(af >= 0.1, af <= 0.9) %>% 
    filter(!chr %in% c('chrY','chrM'))
  
  tomi <- which(snps$af < 0.5)
  snps$af[tomi] <- (1 - snps$af[tomi])
  
  sl <- snps %>% 
    group_by(chr) %>% 
    group_split()
  
  deck <- lapply(seq_len(length(sl)), runsw, sl = sl, bands = bands, length.sw = length.sw)
  deck <- do.call(rbind,deck) %>% mutate(sample_id = basename(id))
  
  return(deck)
  
}

listmat <- mclapply(seq_len(length(lf)),getMatrix,lf=lf,bands = bands,length.sw=length.sw,mc.cores = samples.in.parallel)

mat <- do.call(rbind,listmat)

mat$type <- 'normal'
mat$type[grep(mat$sample_id,pattern = '_tumor')] <- 'tumor'

save(mat,file = file.path(wd, outdir,'alldata.RData'),compress = TRUE)

chrOrder <- paste0('chr',c(1:22,'X'))

mat$chr <- factor(mat$chr, levels = chrOrder)

p <- ggplot(mat, aes(x=arm, y=n_wnds, fill=type)) +
  geom_boxplot() +
  facet_wrap(~chr,ncol = 6,scales = 'free_y') +
  ggtitle(paste('window length =',length.sw,'SNPs'))

ggsave(filename = file.path(wd, outdir,'boxplot_n_wnds.pdf'),width = 300,height = 200,dpi = 300,units = 'mm',device = 'pdf')

p <- ggplot(mat, aes(x=arm, y=median_wnd_length, fill=type)) +
  geom_boxplot() +
  facet_wrap(~chr,ncol = 6,scales = 'free_y') +
  ggtitle(paste('window length =',length.sw,'SNPs'))

ggsave(filename = file.path(wd, outdir,'boxplot_length_wnds.pdf'),width = 300,height = 200,dpi = 300,units = 'mm',device = 'pdf')