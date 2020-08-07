library(countrycode)
library(igraph)
library(parallel)
library(plotly)
library(SeqArray)
library(SNPRelate)
library(tibble)
library(magrittr, "%>%")
import::from(htmlwidgets, "saveWidget")
import::from(pgda, "snpgds_parse")
import::from(plyr, "llply")
import::from(readr, "read_tsv", "read_rds", "write_rds")
import::from(stringr, "str_c")
source("/workspace/repos/wheat-pgda/R/colours.R")

wheat <- snpgdsOpen("/workspace/data/intermediate/gds/sample_pruned_phys.gds", allow.fork = TRUE)

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

pane <- 10
window <- pane * 2 + 1

ld_mats <- lapply(var_sets, function (var_set) {
    snpgdsLDMat(wheat, slide = window, snp.id = var_set, num.thread = 8)$LD
})

gen_pos <- lapply(ld_mats, function (ld_mat) {
    intervals <- cbind(
        c(rep(1, window), 2:(ncol(ld_mat) - (window - 1))), 1:(ncol(ld_mat))
    ) %>% split(., seq(nrow(.)))
    mclapply(intervals, function (interval) {
        cols <- interval[1]:interval[2]
        rows <- (interval[2] - cols) + 1
        lds <- abs(ld_mat[cbind(rows, cols)])
        lds[is.nan(lds)] <- NA
        mean_ld <- 0.1 / (lds %>% mean(na.rm = TRUE))
        if (is.nan(mean_ld)) {
            return(1)
        } else {
            return(mean_ld)
        }
    }, mc.cores = 8) %>% unlist() %>% as.numeric() %>% cumsum()
})

png(
   "/workspace/barley/results/wheat_phys_vs_recip_mean_LD.png",
  family = "Times New Roman", width = 660, height = 1400, pointsize = 5,
  units = "mm", res = 192
)
par(mfrow = c(7, 3))
for (i in 1:21) {
    plot(phys_pos[[i]], gen_pos[[i]])
}
dev.off()

# ################################################################################

# calc_densities <- function (ld_mat, cores) {
#     int_window <- min(window, ncol(ld_mat))
#     mclapply(1:ncol(ld_mat), function (i) {
#         if (i - pane >= 0 && i + pane <= ncol(ld_mat)) {
#             # print("a")
#             cols <- c((i - pane):i, rep(i, pane - 1))
#             rows <- c(pane:1, 1:(pane))
#         } else if (i - pane < 0) {
#             # print("b")
#             cols <- 1:(int_window - 1)
#             # print(cols)
#             if (i == 1) {
#                 rows <- c(1:length(cols))
#             } else {
#                 rows <- c(sum(cols < i):1, 1:sum(cols >= i))
#             }
#             cols[which(cols > i)] <- i
#         } else {
#             # print("c")
#             cols <- (ncol(ld_mat) - (int_window - 1)):(ncol(ld_mat) - 1)
#             if (i == ncol(ld_mat)) {
#                 rows <- c(length(cols):1)
#             } else {
#                 rows <- c(sum(cols < i):1, 1:sum(cols >= i))
#             }
#             cols[which(cols > i)] <- i
#         }
#         indices <- cbind(rows, cols)
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

# densities <- lapply(ld_mats, calc_densities, 1)

################################################################################

ld_prune <- function (i) {
    chrom <- tibble(id = var_sets[[i]], density = densities[[i]], gen_pos = gen_pos[[i]])
    n <- (nrow(chrom) * prop) %>% floor()

    while (nrow(chrom) > n && sum(chrom$density > 0)) {
        # this section finds the largest clique near the densest marker and 
        # removes all but the densest marker in the clique

        # find the current densest marker
        densest <- which.max(chrom$density)
        # create a window around the densest marker
        if (densest - pane2 > 0 && densest + pane2 < nrow(chrom)) {
            # print(11)
            cols <- (densest - pane2):(densest + pane2)
        } else if (densest - pane2 < 0) {
            # print(22)
            cols <- 1:min(window2, nrow(chrom))
        } else {
            # print(33)
            cols <- max(1, (nrow(chrom) - (pane2 * 2))):nrow(chrom)
        }
        # print(cols)
        dists <- abs(chrom$gen_pos[cols] - chrom$gen_pos[densest])
        cols <- cols[which(dists <= max_dist)]
        markers <- chrom$id[cols]

        if (length(markers) <= 1) {
            # print(1)
            # due to the distance cutoff there may be only one marker
            # we therefore set its density to 0 and skip to the next iteration
            # of the loop. This ensures that the next iteration will move on to
            # the next densest marker
            chrom$density[densest] <- 0
            next
        } else {
            # find all combinations of window 2 of markers in this pane
            comps <- combn(markers, 2)
            # get the LD mat of these markers, find which ones are in LD greater
            # than the threshold value, transform it into a vector columnwise and
            # omit NAs from the vector. This creates a vector with indices that
            # match the combos at the same indices in comps. Identify which indices
            # of the vector are true and subset the columns of comps by this and
            # transpose to create the edge list (which variants are in LD above the
            # threshold to each other)
            el <- comps[, which((abs(snpgdsLDMat(wheat, snp.id = markers, num.thread = 1, verbose = FALSE)$LD) > threshold) %>% as.vector() %>% na.omit())] %>% t()
            # use the edgelist to idenify the largest clique in the pane
            cliques <- largest_cliques(graph_from_edgelist(el, directed = FALSE))
            # print(cliques)
            if (length(cliques) > 0) {
                # print(2)
                # identify the left and rightmost markers in the clique, auto transforms
                # into numeric
                present <- range(which(chrom$id %in% cliques[[1]]$name))
                # identify the densest marker in the clique since this may be different
                # from the densest marker in the window
                clique_densest <- which.max(chrom$density[present])
                clique_densest_id <- chrom$id[present][clique_densest]
                # remove all but the densest marker in the clique range from chrom
                chrom <- chrom[-(present)[-(clique_densest)], ]
            } else {
                # print(3)
                chrom$density[densest] <- 0
                next
            }
        }

        ########################################################################
        
        # this section calculates the new densities of the markers near, and
        # including the densest marker in the clique

        # identify which row the densest marker now occupies in chrom
        densest_row <- which(chrom$id %in% clique_densest_id)
        # identify the markers whose densities need to be recalcualted due to
        # the removal of the markers in the clique range
        minor_range_rows <- (densest_row - pane):(densest_row + pane)
        # it is possible for the range to extend beyond the rows of chrom
        # therfore i remove these from the range
        minor_range_rows <- minor_range_rows[which((minor_range_rows >= 1) & (minor_range_rows <= nrow(chrom)))]
        # to calculate densities we must consider markers beyond the range of
        # those whose density need to be recalculated as the density of a marker
        # includes those both upstream and downstream up to pane away
        major_range_rows <- (densest_row - pane * 2):(densest_row + pane * 2)
        # same issue as with minor_range_rows
        major_range_rows <- major_range_rows[which((major_range_rows >= 1) & (major_range_rows <= nrow(chrom)))]
        # get the ids of these markers
        snpset_ids <- chrom$id[major_range_rows]
        # calculate the LD of the snpset markers to each other
        ld <- snpgdsLDMat(wheat, slide = window, snp.id = snpset_ids, num.thread = 1, verbose = FALSE)$LD
        # identify the indices of the markers in snpset_id that correspond to
        # the ones that need to be updated
        to_change <- which(major_range_rows %in% minor_range_rows)
        # calculate the densities of these markers and update them
        chrom$density[major_range_rows[to_change]] <- calc_densities(ld, 1)[to_change]
        # print(1)
        print((nrow(chrom) / length(var_sets[[i]])))
    }
    chrom
}

pane2 <- 50
window2 <- pane2 * 2 + 1

prop <- 0.05
threshold <- 0.7
max_dist <- 10

snp_sets <- mclapply(1:length(var_sets), ld_prune, mc.cores = 8)

################################################################################
png(
   "/workspace/barley/results/wheat_phys_vs_recip_mean_LD_pruned.png",
  family = "Times New Roman", width = 660, height = 1400, pointsize = 5,
  units = "mm", res = 192
)
par(mfrow = c(7, 3))
for (i in 1:length(var_sets)) {
    plot(phys_pos[[i]][var_sets[[i]] %in% snp_sets[[i]]$id], gen_pos[[i]][var_sets[[i]] %in% snp_sets[[i]]$id])
}
dev.off()

pca <- snpgdsPCA(wheat, snp.id = (snp_sets %>% do.call(rbind, .))$id, num.thread = 8)

mtg <- read.gdsn(index.gdsn(wheat, "samp_annot"))$mtg
colours_mtg <- colour_set[c(1, 2, 7, 16, 12, 4, 9, 22)]

scatter <- plot_ly() %>%
  add_markers(
    x = ~pca$eigenvect[, 1], y = ~pca$eigenvect[, 2], z = ~pca$eigenvect[, 3],
    color = mtg, colors = colours_mtg, marker = list(window = 5)
  )

saveWidget(
  as_widget(scatter),
  "/workspace/barley/results/pca.html"
)

snpgdsClose(wheat)