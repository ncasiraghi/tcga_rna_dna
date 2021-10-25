library(tidyverse)

setwd('/BCGLAB/2020_signatures/annotations/WES-Plasma/')

sif.file <- '/BCGLAB/ncasiraghi/tcga_rna_dna/data/datasetSNPs.tsv'

sif <- read.delim(file = sif.file,header = TRUE,stringsAsFactors = FALSE)

# WES-Plasma

wesp <- sif %>% filter(data.source %in% c("WES-Plasma_DEDUP","WES-Plasma_DUP"))

wesp$is.tumor <- 1
wesp$is.tumor[grep(wesp$sample.id,pattern = '-germline$')] <- 0

wesp$is.tissue <- 0

wesp <- wesp %>% 
  select(sample.id,is.tumor,is.tissue) %>% 
  unique()

head(wesp)

# WES-Tissue

west <- sif %>% filter(data.source == 'WES-Plasma_TISSUES')

west$is.tumor <- 1
west$is.tumor[str_which(string = west$sample.id,fixed(pattern = 'ctrl',ignore_case = TRUE))] <- 0

west$is.tissue <- 1

west <- west %>% select(sample.id,is.tumor,is.tissue) %>% unique()

m <- rbind(wesp,west)

m %>% 
  group_by(is.tumor,is.tissue) %>% 
  summarise(n = n())

# TC annotations
# from https://pubmed.ncbi.nlm.nih.gov/32091413

tab <- read.delim('PMID_32091413_TabS4.txt',stringsAsFactors = FALSE) %>% 
  mutate(Sample = gsub(Sample,pattern = 'WCM',replacement = 'PM')) %>% 
  filter(Type == 'Tissue')

mt <- m %>% filter(is.tissue == 1, is.tumor == 1)

mt$Sample <- gsub(mt$sample.id,pattern = 'Sample_',replacement = '')
mt$Sample <- gsub(mt$Sample,pattern = '_Case_HALO.md|_CASE_HALO.md|_Case_HALO.md.filt',replacement = '')
mt$Sample <- gsub(mt$Sample,pattern = '_Case_SS.sorted.dedup.realigned.recalibrated.md.filt|_Case.bwa.sorted.dedup.realigned.recalibrated.md.filt|.bwa.sorted.dedup.realigned.recalibrated.md.filt',replacement = '')
mt$Sample <- gsub(mt$Sample,pattern = '_Case.md.filt|_Case_Halo.md.filt',replacement = '')

unique(mt$Sample)
unique(tab$Sample)

mt <- left_join(x = mt, y = tab, by = 'Sample')

# plasma

mp <- m %>% filter(is.tissue == 0, is.tumor == 1) %>% mutate(Sample = sample.id)

tab <- read.delim('PMID_32091413_TabS4.txt',stringsAsFactors = FALSE) %>% 
  mutate(Sample = gsub(Sample,pattern = 'WCM',replacement = 'PM')) %>% 
  filter(Type == 'Plasma') %>% 
  mutate(Sample = paste0(Sample,'-1st'))

tab$Sample[grep(tab$Sample,pattern = '-TP')] <- gsub(tab$Sample[grep(tab$Sample,pattern = '-TP')],pattern = '-1st',replacement = '')

mp <- left_join(x = mp, y = tab, by = 'Sample')

out <- rbind(mp,mt) %>% select(-Type,-Sample)

ctrls <- m %>% 
  filter(is.tumor == 0) %>% 
  mutate(Class = NA, TC = NA, Ploidy = NA, Site = NA)

final <- rbind(out,ctrls) %>% arrange(desc(is.tissue),desc(is.tumor))

write.table(final,file = '/BCGLAB/2020_signatures/annotations/WES-Plasma/wes-plasma_sif.tsv',quote = F,col.names = T,row.names = F,sep = '\t')