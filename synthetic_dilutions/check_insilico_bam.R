library(tidyverse)
library(data.table)
library(parallel)

folder <- '/BCGLAB/2020_signatures/synthetic_dilutions/bams_insilico'

patients <- list.files(folder,pattern = 'TCGA',full.names = TRUE)

m <- c()

for(pt in patients){
  message(basename(pt))
  bams <- list.files(pt,full.names = FALSE,pattern = '\\.bam$')
  adm <- bams %>% 
    gsub(.,pattern = '\\.bam$',replacement = "") %>% 
    strsplit(.,split = '_adm') %>% 
    lapply(., `[`,2) %>% 
    unlist() %>% 
    str_remove(.,pattern = '.sorted') %>%
    str_replace(.,pattern = '_',replacement = '.') %>% 
    as.numeric() 
  
  tc <- 1 - (adm/100)
  
  df <- data.frame(id = basename(pt),
                   bam = basename(bams),
                   tc = tc,
                   stringsAsFactors = FALSE)
  
  m <- rbind(m,df)
  
}

head(m)

# focus <- m %>% 
#   group_by(id) %>% 
#   summarise(n = n_distinct(tc)) %>% 
#   filter(n == 7) %>% 
#   pull(id)
# 
# ggplot(data = m %>% filter(id %in% focus), aes(x = as.factor(tc), y = tc)) +
#   geom_bar(stat="identity") +
#   facet_wrap(~id,)

ggplot(data = m, aes(x = as.factor(tc), y = tc)) +
  geom_bar(stat="identity") +
  facet_wrap(~id,)

