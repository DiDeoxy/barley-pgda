import::from(dbscan, "hdbscan")
import::from(magrittr, "%>%")
import::from(readr, "read_rds", "write_rds")
import::from(SeqArray, "seqClose", "seqSetFilter", "seqGetData", "seqOpen")
import::from(SNPRelate, "snpgdsDiss")

# snp_set <- (read_rds("/workspace/barley/results/snp_sets.rds") %>% do.call(rbind, .))$id
# barley <- seqOpen("/workspace/barley/gds/imputed.gds")
# diss <- snpgdsDiss(barley, snp.id = snp_set, num.thread = 7)
# seqClose(barley)
# write_rds(diss, "/workspace/barley/results/diss.rds")

diss <- read_rds("/workspace/barley/results/diss.rds")

hdbscan <- hdbscan(diss$diss, minPts = 100)

write_rds(hdbscan, "/workspace/barley/results/hdbscan.rds")

################################################################################
# grm

grm_diss <- read_rds("/workspace/barley/results/grm_diss.rds")

grm_dist <- as.dist(grm_diss)

write_rds(grm_dist, "/workspace/barley/results/grm_dist.rds")

hdbscan <- hdbscan(grm_dist, minPts = 10)

write_rds(hdbscan, "/workspace/barley/results/grm_hdbscan_10.rds")
