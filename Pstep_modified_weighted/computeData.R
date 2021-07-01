library(data.table)
library(tidyverse)
library(parallel)

options(dplyr.summarise.inform = FALSE)
options(scipen = 9999)

source('/BCGLAB/ncasiraghi/tcga_rna_dna/Pstep_modified_weighted/computeData_functions.R')

args <- commandArgs(trailingOnly = TRUE)

data_folder <- args[1]

data_type <- args[2]

aggregate_as  <-  args[3]

min_vaf <- as.numeric(args[4])

max_vaf <- as.numeric(args[5])

min_cov <- as.numeric(args[6])

length.sw <- as.numeric(args[7])

pos_step <- as.numeric(args[8])

samples.in.parallel <- as.numeric(args[9])

wd <- args[10]

exclude_chromosomes <- c('chrM')

# set working dir

# outs

outdir_level_1 <- file.path(wd,paste('minvaf',min_vaf,'maxvaf',max_vaf,'mincov',min_cov,sep = '_'))

if(!file.exists(file.path(outdir_level_1))){
  dir.create(file.path(outdir_level_1), showWarnings = FALSE)
}

outdir_level_2 <- file.path(outdir_level_1,paste('aggregate',aggregate_as,sep = '_'))

if(!file.exists(file.path(outdir_level_2))){
  dir.create(file.path(outdir_level_2), showWarnings = FALSE)
}

outdir_level_3 <- file.path(outdir_level_2,paste('sw',length.sw,sep = '_'))

if(!file.exists(file.path(outdir_level_3))){
  dir.create(file.path(outdir_level_3), showWarnings = FALSE)
}

outdir <- file.path(outdir_level_3,paste('pstep',pos_step,sep = '_'))

if(!file.exists(file.path(outdir))){
  dir.create(file.path(outdir), showWarnings = FALSE)
} else{
  unlink(file.path(outdir),recursive = TRUE)
  dir.create(file.path(outdir), showWarnings = FALSE)
}

setwd(outdir)

# pileup on snps

snps_list <- list.files(data_folder,pattern = '\\.snps$',recursive = TRUE,full.names = TRUE)

datatab_samples <- data.frame(file.snps = snps_list,
                              sample.id = basename(dirname(snps_list)),
                              stringsAsFactors = FALSE)

datatab_samples <- filter(datatab_samples,grepl(pattern = paste0('_',data_type,'_'),sample.id))

write.table(datatab_samples,file = file.path(outdir,'samples.tsv'),sep = '\t',col.names = TRUE,row.names = FALSE,quote = FALSE)

lf <- unique(datatab_samples$file.snps)

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

write.table(positions,file = file.path(outdir,'positions.tsv'),col.names = TRUE,quote = FALSE,row.names = FALSE,sep = '\t')

# run sliding window

getMatrix <- function(i,lf,bands,length.sw,positions,aggregate_as){
  
  file <- lf[i]
  
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

write.table(x = mat,file = file.path(outdir,'data.tsv'),sep = '\t',quote = FALSE,col.names = FALSE,row.names = FALSE)
