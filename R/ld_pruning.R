library(parallel)
library(SeqArray)
library(SNPRelate)
library(magrittr, "%>%")
import::from(plyr, "llply")
import::from(readr, "read_tsv", "read_rds", "write_rds")

barley <- seqOpen("/workspace/barley/gds/imputed.gds")

chrom <- as.character(read.gdsn(index.gdsn(barley, "chromosome"))) %>%
    as.factor()
var_sets <- as.character(read.gdsn(index.gdsn(barley, "variant.id"))) %>%
    split(chrom)
phys_pos <- as.numeric(read.gdsn(index.gdsn(barley, "position"))) %>%
    split(chrom)

size <- 20

ld_mats <- lapply(var_sets, function (var_set) {
    snpgdsLDMat(barley, slide = size, snp.id = var_set, num.thread = 8)$LD
})

gen_pos <- lapply(ld_mats, function (ld_mat) {
    intervals <- cbind(
        c(rep(1, size), 2:(ncol(ld_mat) - (size - 1))), 1:(ncol(ld_mat))
    ) %>% split(., seq(nrow(.)))
    mclapply(intervals, function (interval) {
        cols <- interval[1]:interval[2]
        rows <- (interval[2] - cols) + 1
        lds <- abs(ld_mat[cbind(rows, cols)])
        lds[is.nan(lds)] <- NA
        mean_ld <- 0.001 / (lds %>% mean(na.rm = TRUE))
        if (is.nan(mean_ld)) {
            return(1)
        } else {
            return(mean_ld)
        }
    }, mc.cores = 8) %>% unlist() %>% as.numeric() %>% cumsum()
})

cbind(
    c(rep(1, size), 2:(col - (size - 1))), 1:(col)
)

png(
   "/workspace/barley/results/phys_vs_recip_mean_LD.png",
  family = "Times New Roman", width = 220, height = 1400, pointsize = 5,
  units = "mm", res = 192
)
par(mfrow = c(7, 1))
for (i in 1:7) {
    plot(phys_pos[[i]], gen_pos[[i]])
}
dev.off()

make_rows <- function(cols, start, size) {
    cols <- c(start:cols[1], cols)
    if (length(cols) > size) {
        rows <- c((length(cols) - size):1, 1:size)
    } else {
        rows <- 1:size
    }
    cbind(rows, cols)
    
}

densities <- lapply(ld_mats, function (ld_mat) {
    intervals <- cbind(
        c(rep(1, size + 1), 2:(ncol(ld_mat) - (size))),
        c((size):(ncol(ld_mat)  - 1), rep(ncol(ld_mat) , size))
    ) %>% split(., seq(nrow(.)))
    mclapply(intervals, function (interval) {
        rng <- interval[1]:interval[2]
        if (length(rng) < (size * 2) && interval[1] == 1) {
            cols <- rep((interval[2] - size) + 1, size - 1)
            indices <- make_rows(cols, interval[1], size)
        } else if (length(rng) <= (size * 2) && interval[2] == ncol(ld_mat)) {
            cols <- rep((interval[2] - size), size - 1)
            indices <- make_rows(cols, interval[1], size)
        } else {
            cols <- c(
                interval[1]:(interval[2] - size),
                rep(interval[2] - (size - 1), size)
            )
            rows <- c(size:1, 1:size)
            indices <- cbind(rows, cols)
        }
        lds <- abs(ld_mat[indices])
        lds[is.nan(lds)] <- NA
        density <- mean(lds, na.rm = TRUE)
        if (is.nan(density)) {
            return(1)
        } else {
            return(density)
        }
    }, mc.cores = 8) %>% unlist() %>% as.numeric()
})

library(igraph)

prop = 0.5
window <- 10
threshold <- 0.99
snp_set <- mclapply(1:length(densities), function (i) {
    chrom <- data.frame(id = var_sets[[i]], density = densities[[i]], pos = gen_pos[[i]])
    n <- (nrow(chrom) * prop) %>% floor()
    while (nrow(chrom) > n) {
        densest <- chrom$id[(which.max(chrom$density) - window):(which.max(chrom$density) + window)] %>% as.character()
        comps <- combn(densest, 2)
        el <- comps[, (abs(snpgdsLDMat(barley, snp.id = densest, num.thread = 8)$LD) > threshold) %>% as.vector() %>% na.omit()] %>% t()
        largest_cliques(graph_from_edgelist(el, directed = FALSE)) %>% print()
    }
}, mc.cores = 1)

combn(1:3, 2)