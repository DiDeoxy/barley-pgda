library(parallel)
library(SeqArray)
library(SNPRelate)
library(magrittr, "%>%")
import::from(plyr, "llply")
import::from(readr, "read_tsv", "read_rds", "write_rds")

wheat <- snpgdsOpen("/workspace/data/intermediate/gds/sample_pruned_phys.gds")

chroms <- outer(as.character(1:7), c("A", "B", "D"), paste, sep = "") %>% 
t() %>% as.vector()

chrom <- as.character(read.gdsn(index.gdsn(wheat, "snp.chromosome")))

for (i in 1:21) {
    chrom[which(chrom == i)] <- chroms[i]
}

var_sets <- as.character(read.gdsn(index.gdsn(wheat, "snp.id"))) %>%
    split(chrom)
phys_pos <- as.numeric(read.gdsn(index.gdsn(wheat, "snp.position"))) %>%
    split(chrom)

size <- 10

ld_mats <- lapply(var_sets, function (var_set) {
    snpgdsLDMat(wheat, slide = size, snp.id = var_set, num.thread = 8)$LD
})

mean_lds <- lapply(ld_mats, function (ld_mat) {
    intervals <- cbind(
        c(rep(1, size), 2:(ncol(ld_mat) - size)), 1:(ncol(ld_mat) - 1)
    ) %>% split(., seq(nrow(.)))
    c(0, mclapply(intervals, function (interval) {
        cols <- interval[1]:interval[2]
        rows <- (interval[2] - cols) + 1
        lds <- abs(ld_mat[cbind(rows, cols)])
        lds[is.nan(lds)] <- NA
        mean_ld <- 0.1 / (lds %>% mean(na.rm = TRUE))
        if (is.nan(mean_ld)) {
            return(0)
        } else {
            return(mean_ld)
        }
    }, mc.cores = 1) %>% unlist() %>% as.numeric() %>% cumsum())
})

png(
   "/workspace/barley/results/wheat_phys_vs_recip_mean_LD.png",
  family = "Times New Roman", width = 660, height = 1400, pointsize = 5,
  units = "mm", res = 192
)
par(mfrow = c(7, 3))
for (i in 1:21) {
    plot(phys_pos[[i]], mean_lds[[i]])
}
dev.off()
