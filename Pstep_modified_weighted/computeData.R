library(data.table)
library(tidyverse)
library(parallel)

options(dplyr.summarise.inform = FALSE)
options(scipen = 9999)

source('/BCGLAB/ncasiraghi/tcga_rna_dna/Pstep_modified_weighted/computeData_functions.R')

args <- commandArgs(trailingOnly = TRUE)

configure_file <- args[1]

source(configure_file)

# set working dir

setwd(wd)

# outs

outdir <- paste0('main_out_afs_sw_',length.sw,'_Pstep_',format(pos_step,scientific = FALSE))

if(!file.exists(file.path(wd,outdir))){
  dir.create(file.path(wd, outdir), showWarnings = FALSE)
} else{
  unlink(file.path(wd,outdir),recursive = TRUE)
  dir.create(file.path(wd, outdir), showWarnings = FALSE)
}

# pileup on snps

lf <- list.files(data_folder,pattern = paste0('_',data_type),full.names = TRUE)

write(x = basename(lf),file = file.path(wd,outdir,'samples.txt'),ncolumns = 1)

# cytobands

bands <- read.delim('/BCGLAB/ncasiraghi/gene_panel/input/cytobands/ucsc_cytobands_GRCh38p13.tsv',stringsAsFactors = F)

bands$arm <- 'p'
bands$arm[grep(bands$name,pattern = 'q')] <- 'q'

bands <- bands %>% 
  filter(!chrom %in% c('chrY','chrM')) %>% 
  group_by(chrom,arm) %>% 
  summarise(start = min(chromStart), end = max(chromEnd))

# Positions to aggregate data

arm_borders <- read.delim('/BCGLAB/2020_signatures/stats/modified_Pstep/arm_borders.tsv',stringsAsFactors = F)

arm_borders <- arm_borders %>% 
  group_by(chr,arm) %>% 
  summarise(start = min(pos), end = max(pos))

positions <- c()

for(j in seq_len(nrow(arm_borders))){
  
  pos <- c(seq(from=arm_borders$start[j],arm_borders$end[j],by = pos_step * 1e6), arm_borders$end[j])
  
  x <- data.frame(chrom=arm_borders$chr[j],
                  arm=arm_borders$arm[j],
                  pos=unique(pos),
                  stringsAsFactors = FALSE)
  
  positions <- rbind(positions,x)
  
}

# genes of interest

goi <- c(readLines('/BCGLAB/ncasiraghi/gene_panel/input/oncogenes.txt'),
         readLines('/BCGLAB/ncasiraghi/gene_panel/input/cancer_genes.txt'),
         readLines('/BCGLAB/ncasiraghi/gene_panel/input/tumor_suppressors.txt')) %>% 
  unique() %>% 
  sort()

ensembl <- read.delim('/BCGLAB/ncasiraghi/gene_panel/input/mart_export_GRCh38p13.tsv',check.names = F,stringsAsFactors = F,col.names = c('ensg','start','end','chr','band','Hugo_Symbol'))

ensembl <- ensembl %>% mutate(gene_center = (start + end) %/% 2)

gc <- ensembl %>% filter(Hugo_Symbol %in% goi)

gc$chr <- paste0('chr',gc$chr)
gc$arm <- 'p'
gc$arm[grep(gc$band,pattern = 'q')] <- 'q'

positions$dist_cancer_gene <- NA

for(idx in seq_len(nrow(positions))){
  
  u <- gc %>% 
    filter(chr == positions$chrom[idx] & arm == positions$arm[idx])
  
  if(nrow(u) > 0){
    positions$dist_cancer_gene[idx] <- min(abs(positions$pos[idx] - u$gene_center))
  }
  
}

write.table(positions,file = file.path(wd,outdir,'positions.tsv'),col.names = TRUE,quote = FALSE,row.names = FALSE,sep = '\t')

# run sliding window

getMatrix <- function(i,lf,bands,length.sw,positions,aggregate_as){
  
  id <- lf[i]
  
  file <- list.files(id,full.names = TRUE,pattern = '\\.snps$')
  
  snps <- fread(file,data.table = FALSE,verbose = FALSE,stringsAsFactors = FALSE) %>% 
    filter(af >= min_vaf, af <= max_vaf) %>% 
    filter(cov >= min_cov) %>% 
    filter(!chr %in% exclude_chromosomes)
  
  tomi <- which(snps$af < 0.5)
  snps$af[tomi] <- (1 - snps$af[tomi])
  
  sl <- snps %>% 
    group_by(chr) %>% 
    group_split()
  
  deck <- lapply(seq_len(length(sl)), runsw, sl, bands, length.sw = length.sw, positions = positions, aggregate_as = aggregate_as)
  
  return(as.numeric(unlist(deck)))
  
}

listmat <- mclapply(seq_len(length(lf)),getMatrix,lf=lf,bands=bands,length.sw=length.sw,positions=positions,aggregate_as = aggregate_as,mc.cores = samples.in.parallel)

mat <- do.call(rbind,listmat)

write.table(x = mat,file = file.path(wd,outdir,'data.tsv'),sep = '\t',quote = FALSE,col.names = FALSE,row.names = FALSE)
