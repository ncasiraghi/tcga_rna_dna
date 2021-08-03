library(tidyverse)

setwd('/BCGLAB/2020_signatures/annotations')

sif.file <- '/BCGLAB/ncasiraghi/tcga_rna_dna/data/datasetSNPs.tsv'

sif <- read.delim(file = sif.file,header = TRUE,stringsAsFactors = FALSE)

# SU2C

ctrls <- readLines('/BCGLAB/2020_signatures/pileup/SU2C/BAMs_Normal.txt')
cases <- readLines('/BCGLAB/2020_signatures/pileup/SU2C/BAMs_Tumor.txt')

su2c <- sif %>% filter(data.source == 'SU2C')

su2c$is.tumor <- 0
su2c$is.tumor[which(su2c$sample.id %in% gsub(basename(cases),pattern = '\\.bam',replacement = ''))] <- 1

su2c$is.tissue <- 1

su2c <- su2c %>% select(sample.id,is.tumor,is.tissue,data.source)

# WES-Plasma

wesp <- sif %>% filter(data.source %in% c("WES-Plasma_DEDUP","WES-Plasma_DUP"))

dedup <- wesp %>% filter(data.source == "WES-Plasma_DEDUP") %>% pull(sample.id) %>% unique()
dup <- wesp %>% filter(data.source == "WES-Plasma_DUP") %>% pull(sample.id) %>% unique()

setdiff(dedup,dup)
setdiff(dup,dedup)

wesp$is.tumor <- 1
wesp$is.tumor[grep(wesp$sample.id,pattern = '-germline$')] <- 0

wesp$is.tissue <- 0

wesp <- wesp %>% select(sample.id,is.tumor,is.tissue,data.source)

# WES-Tissue

west <- sif %>% filter(data.source == 'WES-Plasma_TISSUES')

west$is.tumor <- 1
west$is.tumor[str_which(string = west$sample.id,fixed(pattern = 'ctrl',ignore_case = TRUE))] <- 0

west$is.tissue <- 1

west <- west %>% select(sample.id,is.tumor,is.tissue,data.source)

m <- rbind(su2c,wesp,west)

write.table(m,file = 'WES_plasma_tissues_SU2C.tsv',quote = FALSE,sep = '\t',row.names = F,col.names = TRUE)



