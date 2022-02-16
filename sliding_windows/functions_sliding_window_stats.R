compute_stats <- function(indx, chr, adf, tab){
  
  m <- adf[tab$start[indx]:tab$stop[indx],]
  
  out <- data.frame(chr=chr,
                    arm=unique(adf$arm),
                    wnd_length=max(m$pos)-min(m$pos),
                    stringsAsFactors = FALSE)
  return(out)
}

stats_by_arm <- function(df,chr,which.arm,wnd){
  
  adf <- df %>% filter(arm == which.arm)
  
  if( nrow(adf) == 0 ){
    
    stats <- data.frame(chr=chr,
                        arm=which.arm,
                        wnd_length=NA,
                        stringsAsFactors = FALSE)
    
    tab <- data.frame()
    
  } else {
    
    if( wnd <= nrow(adf) ){
      start <- 1:(nrow(adf)-(wnd-1))
      
      stop <- start+(wnd-1)
      
      tab <- data.frame(start = start, stop = stop, stringsAsFactors = FALSE)
    } else {
      tab <- data.frame(start = 1, stop = nrow(adf))
    }
    
    stats <- lapply(seq_len(nrow(tab)), compute_stats, chr = chr, adf = adf, tab = tab)
    stats <- do.call(rbind,stats)
  }
  
  stats <- stats %>% 
    group_by(chr,arm) %>% 
    summarise(median_wnd_length = median(wnd_length, na.rm = TRUE)) %>% 
    mutate(n_wnds = nrow(tab)) %>% 
    ungroup()
  
  return(stats)
  
}

runsw <- function(i, sl, bands, length.sw){
  
  df <- sl[[i]]
  df$arm <- NA
  
  chr <- unique(df$chr)
  
  p.coord <- bands %>% filter(chrom == chr, arm == 'p')
  df$arm[which( df$pos >= p.coord$start & df$pos <= p.coord$end )] <- 'p'
  
  q.coord <- bands %>% filter(chrom == chr, arm == 'q')
  df$arm[which( df$pos >= q.coord$start & df$pos <= q.coord$end )] <- 'q'
  
  out <- rbind(stats_by_arm(df,chr = chr, which.arm = 'p',wnd = length.sw),
               stats_by_arm(df,chr = chr, which.arm = 'q',wnd = length.sw))
  
  return(out)
  
}
