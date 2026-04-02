## This script identifies SOMATIC-labeled SVs where both breakends are within a specified distance (slop) of another SV's breakends
## It then arranges SVs based on how many times they occur and removes the most recurrent ones
## until the target recurrence rate is met. "Removed" SVs are re-labeled as GERMLINE in the predicted_class column.
## For details, see comment at the end of the script.

## Input requirements:
## - A BEDPE data frame with columns: chrom1, start1, end1, chrom2, start2, end2, name, and predicted_class
##   - Each SV should have a unique name in the "name" column
##   - Chromosomes should be in the format "1", "2", ..., "X", "Y"
##   - predicted_class should have values "SOMATIC" or "GERMLINE"

## Output:
## - A BEDPE data frame with the same columns as input, but with some SVs re-labeled as GERMLINE
## - The output is saved in a file called "classified_bedpe_dedup_slop_{slop}_target_recurrence_{target_recurrence}.rds" 

library(dplyr)
library(GenomicRanges)
library(igraph)
library(parallel)

# Load your SV BEDPE data frame here

# sv_bedpe <- readRDS("../data//CCLE356_classified.rds")

# Parameters:

slop <- 10                         # Window within which an SV is considered 'matching' to another SV
target_recurrence <- 1734 / 114656 # Defaults to recurrence-rate observed in TCGA (~1.5%)
n_cores <- 8                       # Number of cores to use for parallel processing

# Helper Funtctions:

# Add column to sv_bedpe indicating component ID of each SV,
# where components are defined as sets of SVs where each SV is within a given slop of at least one other SV in the component
#
# Note that for this function, each breakpoint is allowed to be offset by slop bp
# `calc_sv_frequency` uses a more stringent distance calculation, where the sum
# of the offset distances of the two breakpoints of an SV must be less than or equal to slop
#
# Parameters:
# - sv_bedpe: bedpe dataframe with chrom1, start1, end1, chrom2, start2, and end2 columns
# - slop: integer indicating number of base pairs of offset to allow when determining if two SVs are recurrent
#
# Returns: original sv_bedpe with additional column SV_DUP_ID (component ID)
mark_recurrence <- function(sv_bedpe, slop) {
  # maxgap is offset by 1 from slop, e.g. maxgap of -1 corresponds to exact matches only (slop of 0)
  maxgap <- slop - 1

  # Construct GRanges
  gr1 <- GRanges(
    seqnames = Rle(paste0("chr", sv_bedpe$chrom1)),
    ranges   = IRanges(start = sv_bedpe$start1, end = sv_bedpe$end1)
  )
  gr2 <- GRanges(
    seqnames = Rle(paste0("chr", sv_bedpe$chrom2)),
    ranges   = IRanges(start = sv_bedpe$start2, end = sv_bedpe$end2)
  )

  # Label each SV with a unique ID
  sv_bedpe$sv_id <- seq_len(nrow(sv_bedpe))

  ## Overlaps where both ends match in the same orientation (1-1 and 2-2)
  ## Remove self-matches (e.g. (1,1))
  ## Remove double-counting of pairs (e.g. (1,2) and (2,1))
  o11 <- findOverlaps(gr1, gr1, maxgap = maxgap, type = "any")
  o22 <- findOverlaps(gr2, gr2, maxgap = maxgap, type = "any")
  df11 <- as.data.frame(o11)[, c("queryHits","subjectHits")]
  colnames(df11) <- c("i","j")
  df11 <- df11 %>% filter(i < j)
  df22 <- as.data.frame(o22)[, c("queryHits","subjectHits")]
  colnames(df22) <- c("i","j")
  df22 <- df22 %>% filter(i < j)

  ## Inner join to only keep pairs that match on both ends
  direct_pairs <- inner_join(df11, df22, by = c("i","j"))

  # Overlaps where both ends match in swapped orientation (1-2 and 2-1)
  o12 <- findOverlaps(gr1, gr2, maxgap = maxgap, type = "any")
  o21 <- findOverlaps(gr2, gr1, maxgap = maxgap, type = "any")
  df12 <- as.data.frame(o12)[, c("queryHits","subjectHits")]
  colnames(df12) <- c("i","j")
  df12 <- df12 %>% filter(i < j)
  df21 <- as.data.frame(o21)[, c("queryHits","subjectHits")]
  colnames(df21) <- c("i","j")
  df21 <- df21 %>% filter(i < j)

  swapped_pairs <- inner_join(df12, df21, by = c("i","j"))

  # Combine edges
  edges <- bind_rows(direct_pairs, swapped_pairs) %>%
    dplyr::distinct()

  if (nrow(edges) == 0) {
    sv_bedpe <- sv_bedpe %>%
      mutate(SV_DUP_ID = row_number()) %>%
      select(-sv_id)
    return(sv_bedpe)
  }

  g_dup <- graph_from_data_frame(
    edges,
    directed = FALSE,
    vertices = data.frame(sv_id = sv_bedpe$sv_id)
  )

  comp <- igraph::components(g_dup)
  sv_bedpe <- sv_bedpe %>%
      mutate(SV_DUP_ID = comp$membership) %>%
      select(-sv_id)
  return(sv_bedpe)
}

# Calculate number of SVs up to `slop` distance away.
# Distance is calculated in a more stringent way than in `mark_recurrence`
# - For `mark_recurrence`, each breakpoint is allowed to be offset by `slop` bp
# - Here, the sum of the offset distances of the two breakpoints of an SV must be less than or equal to `slop`
#
# Calculation is sped up by restricting to SVs in the same SV_DUP_ID group, since SVs in different groups
# are guaranteed to be more than `slop` distance away from each other. 
#
# Parameters:
# - sv_bedpe: bedpe with SV_DUP_ID column indicating SV groups within a given slop of each other 
# - slop: Slop value that was used for `mark_recurrence`
# - n_cores: Number of cores to use for parallel processing
# 
# Returns:
# - A tibble with columns "name" and "frequency", where "name" is the name of each SV and "frequency" is the number of SVs within `slop` distance of that SV
calc_sv_frequency <- function(sv_bedpe, slop, n_cores = 1) {

  dists <- sv_bedpe %>%
    # Arrange so chrom1 is always "smaller" chromosome and start1 is always before start2
    mutate(
      chrom1 = case_when(
        chrom1 == "X" ~ "23",
        chrom1 == "Y" ~ "24",
        TRUE          ~ chrom1
      ),
      chrom2 = case_when(
        chrom2 == "X" ~ "23",
        chrom2 == "Y" ~ "24",
        TRUE          ~ chrom2
      ),
      chrom1 = as.numeric(chrom1),
      chrom2 = as.numeric(chrom2),
      chrom1_ord  = if_else(chrom1 < chrom2, chrom1, chrom2),
      chrom2_ord  = if_else(chrom2 > chrom1, chrom2, chrom1),
      start1_ord      = case_when(
        chrom1 < chrom2 ~ start1,
        chrom1 > chrom2 ~ start2,
        TRUE         ~ pmin(start1, start2)
      ),
      start2_ord      = case_when(
        chrom1 < chrom2 ~ start2,
        chrom1 > chrom2 ~ start1,
        TRUE            ~ pmax(start1, start2)
    )) %>%
    mutate(start1 = start1_ord, start2 = start2_ord, chrom1 = chrom1_ord, chrom2 = chrom2_ord, group_ID = SV_DUP_ID, ID = name, .keep = "none") 

    dists_list <- split(dists, dists$group_ID)
    res <- do.call(bind_rows, parallel::mclapply(dists_list, calc_sv_frequency_helper, slop = slop, mc.cores = n_cores))

    return(res)
}

# Helper function to calculate number of SVs within `slop` distance of each SV in a given group of SVs (i.e. SVs with the same SV_DUP_ID)
calc_sv_frequency_helper <- function(df, slop) {
  if (nrow(df) == 1) return(tibble(name = df$ID, frequency = 1))

  freqs <- sapply(1:nrow(df), function(i) {
    sv_i <- df[i,]
    other_svs <- df[-i,]

    # Calculate distance to each other SV in group
    dist_to_others <- abs(sv_i$start1 - other_svs$start1) + abs(sv_i$start2 - other_svs$start2)

    res <- sum(dist_to_others <= slop) + 1 # Add 1 to include the SV itself
    return(res)
  })

  return(tibble(name = df$ID, frequency = freqs))
}

# Assign unique ID to each set of identical somatic SVs
somatic_sv_recurrence <- mark_recurrence(sv_bedpe %>% filter(predicted_class == "SOMATIC"), 0)

# Count the number of SVs within `slop` distance of each SV, using the more stringent distance calculation described above
all_sv_recurrence <- calc_sv_frequency(mark_recurrence(sv_bedpe, slop), slop = slop, n_cores = n_cores)

sv_freq <- somatic_sv_recurrence %>%
  select(name, SV_DUP_ID) %>%
  left_join(all_sv_recurrence, by = "name") %>%
  group_by(SV_DUP_ID) %>%
  summarize(somatic_frequency = n(), all_frequency = unique(frequency))


pct_recurrent <- sum(sv_freq$all_frequency > 1) / nrow(sv_freq)
message(paste("Before filtering, ", round(pct_recurrent * 100, 2), '% of SVs are recurrent'))

if (pct_recurrent <= target_recurrence) {
  message("Data already meets target recurrence. Saving original data frame as deduplicated version.")
  saveRDS(sv_bedpe, file.path(paste0("classified_bedpe_dedup_slop_", slop, "_target_recurrence_", target_recurrence, ".rds")))
} else {
  num_unique_somatic <- nrow(sv_freq)
  dup_stats <- sv_freq %>%
    filter(all_frequency > 1) %>%
    arrange(desc(all_frequency)) %>%
    mutate(pct_recurrent = (n() - row_number()) / num_unique_somatic)

  # Start removing most recurrent SV groups until target recurrence is met
  cutoff_idx <- min(which(dup_stats$pct_recurrent <= target_recurrence))
  sv_groups_to_remove <- dup_stats$SV_DUP_ID[1:cutoff_idx]
  svs_to_remove <- somatic_sv_recurrence %>%
    filter(SV_DUP_ID %in% sv_groups_to_remove) %>%
    pull(name)
  pct_recurrent <- dup_stats$pct_recurrent[cutoff_idx]

  message(paste0("Re-labeling ", length(svs_to_remove), " SVs (", length(sv_groups_to_remove), " unique SVs)"))

  # Mark selected SVs as GERMLINE in predicted_class column
  sv_bedpe_dedup <- sv_bedpe %>%
    mutate(predicted_class = if_else(name %in% svs_to_remove, "GERMLINE", predicted_class))

  message(paste0("After filtering, ", round(pct_recurrent * 100, 2), '% of SVs are recurrent'))

  saveRDS(sv_bedpe_dedup, file.path(paste0("classified_bedpe_dedup_slop_", slop, "_target_recurrence_", target_recurrence, ".rds")))
}

# Algorithm details:
# - Percent-recurrence is calculated as (# of unique somatic SVs with at least one other SV within `slop` distance) / (total # of somatic unique SVs)
# 1) Unique somatic SVs are identified by grouping together SVs with exactly-matching breakpoints (same chrom1, chrom2, start1, and start2)
#    - This accounts for SVs that have their breakpoints listed in opposite order
# 2) Number of SVs within `slop` distance of each SV is calculated using a two step process:
#    a) First, SVs are grouped into components where each SV is within `slop` distance of at least one other SV in the component.
#       - Here, each breakpoint is allowed to be offset by `slop` bp
#    b) Then, for each SV, the number of SVs within `slop` distance is calculated by only comparing to other SVs in the same component,
#       and using the more stringent distance calculation where the sum of the offset distances of the two breakpoints of an SV must be less than or equal to `slop`
#   - The first step is done to speed up the second step, since SVs in different components are guaranteed to be more than `slop` distance away from each other
# 3) Unique somatic SVs are arranged in descending order of the number of SVs within `slop` distance, and the most recurrent SVs are removed until the target recurrence rate is met
