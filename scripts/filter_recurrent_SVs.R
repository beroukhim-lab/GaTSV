library(dplyr)
library(GenomicRanges)
library(igraph)
library(data.table)

sv_bedpe <- classified_bedpe[predicted_class=='SOMATIC',] #load the classified bedpe result from GaTSV and filter for somatic SVs
slop <- 10 #the window within which an SV is considered 'matching' to another SV

# Main loop
max_iter <- 10
iter <- 1
converged <- F

repeat {
  message("Iteration = ", iter )
  # Construct GRanges
  gr1 <- GRanges(
    seqnames = Rle(paste0("chr", sv_bedpe$chrom1)),
    ranges   = IRanges(start = sv_bedpe$start1, end = sv_bedpe$end1)
  )
  gr2 <- GRanges(
    seqnames = Rle(paste0("chr", sv_bedpe$chrom2)),
    ranges   = IRanges(start = sv_bedpe$start2, end = sv_bedpe$end2)
  )
  
  sv_bedpe$sv_id <- seq_len(nrow(sv_bedpe))
  # Direct orientation overlaps
  o11 <- findOverlaps(gr1, gr1, maxgap = slop)
  o22 <- findOverlaps(gr2, gr2, maxgap = slop)
  df11 <- as.data.frame(o11)[, c("queryHits","subjectHits")]
  colnames(df11) <- c("i","j")
  df22 <- as.data.frame(o22)[, c("queryHits","subjectHits")]
  colnames(df22) <- c("i","j")
  direct_pairs <- inner_join(df11, df22, by = c("i","j")) %>%
    filter(i < j)
  single_direct <- bind_rows(df11, df22) %>%
    filter(i != j) %>%
    distinct()
  
  # Swapped orientation overlaps
  o12 <- findOverlaps(gr1, gr2, maxgap = slop)
  o21 <- findOverlaps(gr2, gr1, maxgap = slop)
  df12 <- as.data.frame(o12)[, c("queryHits","subjectHits")]
  colnames(df12) <- c("i","j")
  df21 <- as.data.frame(o21)[, c("queryHits","subjectHits")]
  colnames(df21) <- c("i","j")
  swapped_pairs <- inner_join(df12, df21, by = c("i","j")) %>%
    filter(i != j)
  single_swapped <- bind_rows(df12, df21) %>%
    filter(i != j) %>%
    distinct()
  
  # Combine edges
  edges <- bind_rows(
    direct_pairs,
    swapped_pairs,
    single_direct,
    single_swapped
  ) %>%
    distinct() %>%
    filter(i < j)
  
  if (nrow(edges) == 0) {
    warning("No overlapping SVs found at slop ", slop,'. Stopping')
    sv_bedpe$SV_DUP_ID <- seq_len(nrow(sv_bedpe))
    break
  }
  
  g_dup <- graph_from_data_frame(
    edges,
    directed = FALSE,
    vertices = data.frame(sv_id = sv_bedpe$sv_id)
  )
  
  comp_dup <- components(g_dup)
  
  sv_bedpe$DUP_GROUP <- comp_dup$membership
  
  dup_stats <- sv_bedpe %>%
    filter(sv_id %in% unique(c(edges$i, edges$j))) %>%
    group_by(DUP_GROUP) %>%
    summarise(
      N_OCCURRENCES  = n(),
      SV_NUM_SAMPLES = n_distinct(sample),
      .groups = "drop"
    )
  percent_recurrent <- round((sum(dup_stats$N_OCCURRENCES)/nrow(sv_bedpe))*100,2)
  if (percent_recurrent <= 2){
    message('≤ 2% duplicated SVs. Stopping deduplication.')
    converged <- T
    break
  }
  
  message(paste("Currently, ",percent_recurrent,'% of your SVs are recurrent'))
  
  target_recurrence <- 0.02*nrow(sv_bedpe) #2% recurrence threshold
  dup_stats<-dup_stats[order(-dup_stats$N_OCCURRENCES,-dup_stats$SV_NUM_SAMPLES),]
  dup_stats$cumsum_N <- cumsum(dup_stats$N_OCCURRENCES)
  
  dup_groups_to_lose <- as.vector(dup_stats[dup_stats$cumsum_N<target_recurrence,]$DUP_GROUP)
  sv_bedpe_dedup <- sv_bedpe %>%
    filter(!(DUP_GROUP %in% dup_groups_to_lose))
  
  iter <- iter+1
  sv_bedpe <- sv_bedpe_dedup #updates sv_bedpe
  
  if (iter>max_iter){
    warning("Maximum iterations reached. Deduplication may be incomplete.\nEither increase the maximum iterations threshold or manually evaluate recurrent call set")
    break
  }
}

#This stores the best available de-duplicated call set in sv_bedpe, and whether the solution converged (T) within the rounds of de-duplication
saveRDS(
  list(
    sv_bedpe = sv_bedpe,
    converged = converged
  ),
  file.path(
    paste0("classified_bedpe_dedup_slop_", slop, ".rds")
  )
)

