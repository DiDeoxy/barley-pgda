library(igraph)
library(parallel)
library(SeqArray)
library(SNPRelate)
library(tibble)
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

window <- 10
size <- window * 2 + 1

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

################################################################################

# make_rows <- function(cols, start, size) {
#     cols <- c(start:cols[1], cols)
#     if (length(cols) > size) {
#         rows <- c((length(cols) - size):1, 1:size)
#     } else {
#         rows <- 1:size
#     }
#     cbind(rows, cols)
    
# }

# calc_densities <- function (ld_mat, cores) {
#     size <- (nrow(ld_mat) - 1) / 2
#     intervals <- cbind(
#         c(rep(1, size + 1), 2:(ncol(ld_mat) - (size))),
#         c(size:(ncol(ld_mat) - 1), rep(ncol(ld_mat), size))
#     ) %>% split(., seq(nrow(.)))
#     # print(intervals)
#     mclapply(intervals, function (interval) {
#         rng <- interval[1]:interval[2]
#         if (length(rng) < (size * 2) && interval[1] == 1) {
#             cols <- rep((interval[2] - size) + 1, size - 1)
#             indices <- make_rows(cols, interval[1], size)
#         } else if (length(rng) <= (size * 2) && interval[2] == ncol(ld_mat)) {
#             cols <- rep((interval[2] - size), size - 1)
#             indices <- make_rows(cols, interval[1], size)
#         } else {
#             cols <- c(
#                 interval[1]:(interval[2] - size),
#                 rep(interval[2] - (size - 1), size)
#             )
#             rows <- c(size:1, 1:size)
#             indices <- cbind(rows, cols)
#         }
#         lds <- abs(ld_mat[indices])
#         lds[is.nan(lds)] <- NA
#         density <- mean(lds, na.rm = TRUE)
#         if (is.nan(density)) {
#             return(1)
#         } else {
#             return(density)
#         }
#     }, mc.cores = cores) %>% unlist() %>% as.numeric()
# }

calc_densities <- function (ld_mat, cores) {
    int_size <- min(size, ncol(ld_mat))
    mclapply(1:ncol(ld_mat), function (i) {
        if (i - window > 0 && i + window < ncol(ld_mat)) {
            # print(1)
            cols <- c((i - window):i, rep(i, window - 1))
            rows <- c(window:1, 1:(window))
        } else if (i - window > 0) {
            # print(2)
            cols <- (ncol(ld_mat) - (int_size - 1)):(ncol(ld_mat) - 1)
            cols[which(cols > i)] <- i
            if (i == ncol(ld_mat)) {
                rows <- c(sum(cols < i):1)
            } else {
                rows <- c(sum(cols < i):1, 1:sum(cols == i))
            }
        } else {
            # print(3)
            cols <- 1:(int_size - 1)
            cols[which(cols > i)] <- i
            if (i == 1) {
                rows <- c(1:sum(cols == i))
            } else {
                rows <- c(sum(cols < i):1, 1:sum(cols == i))
            }
        }
        indices <- cbind(rows, cols)
        lds <- abs(ld_mat[indices])
        lds[is.nan(lds)] <- NA
        density <- mean(lds, na.rm = TRUE)
        if (is.nan(density)) {
            return(1)
        } else {
            return(density)
        }
    }, mc.cores = cores) %>% unlist() %>% as.numeric()
}

window <- 50
size <- window * 2 + 1

ld_mats <- lapply(var_sets, function (var_set) {
    snpgdsLDMat(barley, slide = size, snp.id = var_set, num.thread = 8)$LD
})

densities <- lapply(ld_mats, calc_densities, 8)

ld_prune <- function (i) {
    chrom <- tibble(id = var_sets[[i]], density = densities[[i]])
    n <- (nrow(chrom) * prop) %>% floor()
    while (nrow(chrom) > n) {
        # find the current densest marker
        densest <- which.max(chrom$density)
        # create a window around the densest marker
        markers <- (densest - window):(densest + window)
        # select those snps that actually exist in chrom and extract their ids
        markers <- chrom$id[markers[which((markers >= 1) & (markers <= nrow(chrom)))]]
        # find all combinations of size 2 of markers in this window
        comps <- combn(markers, 2)
        # get the LD mat of these markers, find which ones are in LD greater
        # than the threshold value, transform it into a vector columnwise and
        # omit NAs from the vector. This creates a vector with indices that
        # match the combos at the same indices in comps. Identify which indices
        # of the vector are true and subset the columns of comps by this and
        # transpose to create the edge list (which variants are in LD above the
        # threshold to each other)
        el <- comps[, which((abs(snpgdsLDMat(barley, snp.id = markers, num.thread = 1)$LD) > threshold) %>% as.vector() %>% na.omit())] %>% t()
        # use the edgelist to idenify the largest clique in the window
        clique <- largest_cliques(graph_from_edgelist(el, directed = FALSE))[[1]]$name
        # identify the left and rightmost markers in the clique, auto transforms
        # into numeric
        clique_range <- range(clique)
        # identify all possible markers in this range
        covered <- clique_range[1]:clique_range[2] %>% as.character()
        # identify those markers in this range that are still present in the
        # chrom data frame
        present <- which(chrom$id %in% covered)
        # identify the densest marker in the clique since this may be different
        # from the densest marker in the window
        clique_densest <- chrom$id[present][which.max(chrom$density[present])]
        # remove all but the densest marker in the clique range from chrom
        chrom <- chrom[-which(chrom$id %in% covered[-which(covered %in% clique_densest)]), ]
        # identify which row the densest marker now occupies in chrom
        densest_row <- which(chrom$id %in% clique_densest)
        # identify the markers whose densities need to be recalcualted due to
        # the removal of the markers in the clique range
        minor_range_rows <- (densest_row - window):(densest_row + window)
        # it is possible for the range to extend beyond the rows of chrom
        # therfore i remove these from the range
        minor_range_rows <- minor_range_rows[which((minor_range_rows >= 1) & (minor_range_rows <= nrow(chrom)))]
        # to calculate densities we must consider markers beyond the range of
        # those whose density need to be recalculated as the density of a marker
        # includes those both upstream and downstream up to window away
        major_range_rows <- (densest_row - window * 2):(densest_row + window * 2)
        # same issue as with minor_range_rows
        major_range_rows <- major_range_rows[which((major_range_rows >= 1) & (major_range_rows <= nrow(chrom)))]
        # get the ids of these markers
        snpset_ids <- chrom$id[major_range_rows]
        # calculate the LD of the snpset markers to each other
        ld <- snpgdsLDMat(barley, slide = size, snp.id = snpset_ids, num.thread = 1, verbose = FALSE)$LD
        # identify the indices of the markers in snpset_id that correspond to
        # the ones that need to be updated
        to_change <- which(major_range_rows %in% minor_range_rows)
        # calculate the densities of these markers and update them
        chrom$density[major_range_rows[to_change]] <- calc_densities(ld, 1)[to_change]
        print((nrow(chrom) / length(var_sets[[i]])))
    }
    chrom
}

prop <- 0.05
threshold <- 0.7
snp_sets <- mclapply(1:length(densities), ld_prune, mc.cores = 8)

which(! var_sets[[1]] %in% snp_sets[[1]]$id)
################################################################################
png(
   "/workspace/barley/results/phys_vs_recip_mean_LD_pruned.png",
  family = "Times New Roman", width = 220, height = 1400, pointsize = 5,
  units = "mm", res = 192
)
par(mfrow = c(7, 1))
for (i in 1:7) {
    plot(phys_pos[[i]][var_sets[[i]] %in% snp_sets[[i]]$id], gen_pos[[i]][var_sets[[i]] %in% snp_sets[[i]]$id])
}
dev.off()

pca <- snpgdsPCA(barley, snp.id=(snp_sets %>% do.call(rbind, .))$id, num.thread = 8)