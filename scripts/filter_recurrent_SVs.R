library(dplyr)
library(GenomicRanges)
library(igraph)

## This script identifies SVs where both breakends overlap (within a specified slop)
## It then arranges SVs based on how many times they occur and removes the most recurrent ones
## until the target recurrence rate is met.

# Load your SV BEDPE data frame here
# sv_bedpe <- readRDS("../data//CCLE356_classified.rds")

# Window within which an SV is considered 'matching' to another SV
slop <- 10
target_recurrence <- 0.02

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
o11 <- findOverlaps(gr1, gr1, maxgap = slop, type = "any")
o22 <- findOverlaps(gr2, gr2, maxgap = slop, type = "any")
df11 <- as.data.frame(o11)[, c("queryHits","subjectHits")]
colnames(df11) <- c("i","j")
df22 <- as.data.frame(o22)[, c("queryHits","subjectHits")]
colnames(df22) <- c("i","j")

## Inner join to only keep pairs that match on both ends
## Remove self-matches (e.g. (1,1))
## Remove double-counting of pairs (e.g. (1,2) and (2,1))
direct_pairs <- inner_join(df11, df22, by = c("i","j")) %>%
  filter(i < j)

# Overlaps where both ends match in swapped orientation (1-2 and 2-1)
o12 <- findOverlaps(gr1, gr2, maxgap = slop, type = "any")
o21 <- findOverlaps(gr2, gr1, maxgap = slop, type = "any")
df12 <- as.data.frame(o12)[, c("queryHits","subjectHits")]
colnames(df12) <- c("i","j")
df21 <- as.data.frame(o21)[, c("queryHits","subjectHits")]
colnames(df21) <- c("i","j")

swapped_pairs <- inner_join(df12, df21, by = c("i","j")) %>%
  dplyr::filter(i < j)

# Combine edges
edges <- bind_rows(direct_pairs, swapped_pairs) %>%
  dplyr::distinct()

if (nrow(edges) == 0) {
  message("No overlapping SVs found at slop ", slop,'. Saving original data frame as deduplicated version.')
  saveRDS(sv_bedpe, file.path(paste0("classified_bedpe_dedup_slop_", slop, "target_recurrence", target_recurrence, ".rds")))
  break
}

g_dup <- graph_from_data_frame(
  edges,
  directed = FALSE,
  vertices = data.frame(sv_id = sv_bedpe$sv_id)
)

comp <- components(g_dup)
sv_bedpe <- sv_bedpe %>%
    dplyr::mutate(SV_DUP_ID = comp$membership) %>%
    group_by(SV_DUP_ID) %>%
    dplyr::mutate(GROUP_SIZE = n()) %>%
    ungroup()

num_somatic <- nrow(sv_bedpe %>% filter(predicted_class == "SOMATIC"))
num_somatic_dup <- nrow(sv_bedpe %>% filter(predicted_class == "SOMATIC", GROUP_SIZE >1))
pct_recurrent <- num_somatic_dup / num_somatic

# Save and exit if already below target recurrence
message(paste("Before filtering, ", round(pct_recurrent * 100, 2), '% of SVs are recurrent'))
if (pct_recurrent <= target_recurrence) {
  sv_bedpe <- sv_bedpe %>%
    select(-sv_id, SV_DUP_ID, GROUP_SIZE)
  saveRDS(sv_bedpe, file.path(paste0("classified_bedpe_dedup_slop_", slop, "target_recurrence", target_recurrence, ".rds")))
  break
}

dup_stats <- sv_bedpe %>%
  # Filter down to somatic SVs that are marked as duplicates
  filter(GROUP_SIZE > 1, predicted_class == "SOMATIC") %>%
  # Count number of occurences in somatic SVs
  group_by(SV_DUP_ID) %>%
  summarise(N_OCCURRENCES  = n()) %>%
  arrange(desc(N_OCCURRENCES)) %>%
  mutate(N_REMOVE = cumsum(N_OCCURRENCES)) %>%
  mutate(PCT_RECURRENT = (num_somatic_dup - N_REMOVE) / (num_somatic - N_REMOVE))

cutoff_idx <- min(which(dup_stats$PCT_RECURRENT <= target_recurrence))
svs_groups_to_remove <- dup_stats$SV_DUP_ID[1:cutoff_idx]
pct_recurrent <- dup_stats$PCT_RECURRENT[cutoff_idx]

sv_bedpe_dedup <- sv_bedpe %>%
  mutate(predicted_class = if_else(SV_DUP_ID %in% svs_groups_to_remove, "GERMLINE", predicted_class)) %>%
  select(-sv_id, SV_DUP_ID, GROUP_SIZE)

message(paste("After filtering, ", round(pct_recurrent * 100, 2), '% of SVs are recurrent'))

saveRDS(sv_bedpe_dedup, file.path(paste0("classified_bedpe_dedup_slop_", slop, "target_recurrence", target_recurrence, ".rds")))