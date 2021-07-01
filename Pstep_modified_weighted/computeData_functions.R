compute_stats <- function(indx, chr, adf, tab){
  
  m <- adf[tab$start[indx]:tab$stop[indx],]
  
  out <- data.frame(chr=chr,
                    arm=unique(adf$arm),
                    wnd_start=min(m$pos),
                    wnd_end=max(m$pos),
                    snp_pos=paste(m$pos,collapse = ';'),
                    snp_af=paste(m$af,collapse = ';'),
                    stringsAsFactors = FALSE)
  return(out)
}

get_afs <- function(px, stats, aggregate_as){
  
  x <- stats %>% 
    filter(wnd_start <= px & wnd_end >= px) 
  
  compute_wafs <- function(j,x,px){
    
    dst <- abs(as.numeric(unlist(strsplit(x$snp_pos[j],split = ';'))) - px) + 1
    
    waf <- weighted.mean(x = as.numeric(unlist(strsplit(x$snp_af[j],split = ';'))), w = 1/dst)
    
    return(waf)
  }
  
  if(nrow(x) > 0){
    
    wout <- unlist(lapply(seq_len(nrow(x)),compute_wafs,x=x,px=px))
    
    if(aggregate_as == 'median'){
      return( median(wout,na.rm = TRUE) )
    }
    
    if(aggregate_as == 'mean'){
      return( mean(wout,na.rm = TRUE) )
    }
    
  } else {
    
    return( NA )
    
  }

}

stats_by_arm <- function(df,chr,which.arm,wnd,positions,aggregate_as){
  
  adf <- df %>% filter(arm == which.arm)
  
  if( nrow(adf) == 0 ){
    
    stats <- data.frame(chr=chr,
                        arm=which.arm,
                        wnd_start=NA,
                        wnd_end=NA,
                        snp_pos=NA,
                        snp_af=NA,
                        stringsAsFactors = FALSE)
    
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
  
  pos.arm <- as.numeric(positions %>% filter(arm == which.arm, chrom == chr ) %>% pull(pos))
  
  afs <- lapply(pos.arm, get_afs, stats = stats, aggregate_as = aggregate_as)
  afs <- as.numeric(unlist(afs))
  
  return(afs)
  
}

runsw <- function(i, sl, bands, length.sw, positions, aggregate_as){
  
  df <- sl[[i]]
  df$arm <- NA
  
  chr <- unique(df$chr)
  
  p.coord <- bands %>% filter(chrom == chr, arm == 'p')
  df$arm[which( df$pos >= p.coord$start & df$pos <= p.coord$end )] <- 'p'
  
  q.coord <- bands %>% filter(chrom == chr, arm == 'q')
  df$arm[which( df$pos >= q.coord$start & df$pos <= q.coord$end )] <- 'q'
  
  out <- c(stats_by_arm(df,chr = chr, which.arm = 'p',wnd = length.sw, positions = positions, aggregate_as = aggregate_as),
           stats_by_arm(df,chr = chr, which.arm = 'q',wnd = length.sw, positions = positions, aggregate_as = aggregate_as))
  
  return(out)
  
}
