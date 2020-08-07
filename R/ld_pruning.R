library(igraph)
library(parallel)
library(plotly)
library(SeqArray)
library(SNPRelate)
library(tibble)
import::from(magrittr, "%>%")
import::from(htmlwidgets, "saveWidget")
import::from(readr, "write_rds")
source("/workspace/repos/barley-pgda/R/colours.R")

barley <- seqOpen("/workspace/barley/gds/imputed.gds")

chrom <- as.character(read.gdsn(index.gdsn(barley, "chromosome"))) %>%
    as.factor()
var_sets <- as.character(read.gdsn(index.gdsn(barley, "variant.id"))) %>%
    split(chrom)
phys_pos <- as.numeric(read.gdsn(index.gdsn(barley, "position"))) %>%
    split(chrom)

pane <- 5
window <- pane * 2 + 1

ld_mats <- lapply(var_sets, function (var_set) {
    snpgdsLDMat(barley, slide = window, snp.id = var_set, num.thread = 8)$LD
})

################################################################################

calc_density <- function (ld_mat, cores) {
  int_window <- min(window, ncol(ld_mat))
  mclapply(1:ncol(ld_mat), function (i) {
    if (i - pane >= 0 && i + pane <= ncol(ld_mat)) {
      cols <- c((i - pane):i, rep(i, pane - 1))
      rows <- c(pane:1, 1:(pane))
    } else if (i - pane < 0) {
      cols <- 1:(int_window - 1)
      if (i == 1) {
        rows <- c(1:length(cols))
      } else {
        rows <- c(sum(cols < i):1, 1:sum(cols >= i))
      }
      cols[which(cols > i)] <- i
    } else {
      cols <- (
          ncol(ld_mat) - (int_window - 1)
      ):(
          ncol(ld_mat) - 1
      )
      if (i == ncol(ld_mat)) {
        rows <- c(length(cols):1)
      } else {
        rows <- c(sum(cols < i):1, 1:sum(cols >= i))
      }
      cols[which(cols > i)] <- i
    }
    indices <- cbind(rows, cols)
    lds <- abs(ld_mat[indices])
    lds[is.nan(lds)] <- NA
    density <- mean(lds, na.rm = TRUE)
    if (is.nan(density)) {
        return(0)
    } else {
        return(density)
    }
  }, mc.cores = cores) %>% unlist() %>% as.numeric()
}

################################################################################

density <- lapply(ld_mats, calc_density, 7)

png(
   "/workspace/barley/results/phys_vs_density.png",
  family = "Times New Roman", width = 220, height = 1400, pointsize = 5,
  units = "mm", res = 192
)
par(mfrow = c(7, 1))
for (i in 1:7) {
    plot(phys_pos[[i]], density[[i]])
}
dev.off()

sparsity <- lapply(density, function (chrom) { -log(chrom) })

png(
   "/workspace/barley/results/phys_vs_sparsity.png",
  family = "Times New Roman", width = 220, height = 1400, pointsize = 5,
  units = "mm", res = 192
)
par(mfrow = c(7, 1))
for (i in 1:7) {
    plot(phys_pos[[i]], sparsity[[i]])
}
dev.off()

gen_pos <- lapply(sparsity, cumsum)

png(
   "/workspace/barley/results/phys_vs_gen_pos.png",
  family = "Times New Roman", width = 220, height = 1400, pointsize = 5,
  units = "mm", res = 192
)
par(mfrow = c(7, 1))
for (i in 1:7) {
    plot(phys_pos[[i]], gen_pos[[i]])
}
dev.off()


################################################################################

ld_prune <- function (i) {
  chrom <- tibble(id = var_sets[[i]], density = density[[i]], gen_pos = gen_pos[[i]])
  n <- (nrow(chrom) * prop) %>% floor()

  while (nrow(chrom) > n) {
    # this section finds the largest clique near the densest marker and 
    # removes all but the densest marker in the clique

    # find the current densest marker
    densest <- which.max(chrom$density)

    if (chrom$density[densest] == 0) {
      # print(0)
      return(chrom)
    } else {
      # create a window around the densest marker
      if (densest - pane2 > 0 && densest + pane2 < nrow(chrom)) {
        # print(11)
        cols <- (densest - pane2):(densest + pane2)
      } else if (densest - pane2 < 0) {
        # print(22)
        cols <- 1:window2
      } else {
        # print(33)
        cols <- (nrow(chrom) - (pane2 * 2)):nrow(chrom)
      }
      dists <- abs(chrom$gen_pos[cols] - chrom$gen_pos[densest])
      cols <- cols[which(dists <= max_dist)]
      markers <- chrom$id[cols]

      if (length(markers) <= 1) {
        # print(1)
        # due to the distance cutoff there may be only one marker
        # we therefore set its density to Inf and skip to the next iteration
        # of the loop. This ensures that the next iteration will move on to
        # the next densest marker
        chrom$density[which(chrom$id == markers)] <- 0
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
        el <- comps[, which((abs(snpgdsLDMat(barley, snp.id = markers, num.thread = 1)$LD) > threshold) %>% as.vector() %>% na.omit())] %>% t()
        # use the edgelist to idenify the largest clique in the pane
        cliques <- largest_cliques(graph_from_edgelist(el, directed = FALSE))
        if (length(cliques) > 0) {
          # print(2)
          # identify the left and rightmost markers in the clique, auto transforms
          # into numeric
          clique_range <- range(cliques[[1]]$name)
          # identify all possible markers in this range
          covered <- clique_range[1]:clique_range[2] %>% as.character()
          # identify those markers in this range that are still present in the
          # chrom data frame
          present <- which(chrom$id %in% covered)
          print(length(present))
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
      ld <- snpgdsLDMat(barley, slide = window, snp.id = snpset_ids, num.thread = 1, verbose = FALSE)$LD
      # identify the indices of the markers in snpset_id that correspond to
      # the ones that need to be updated
      to_change <- which(major_range_rows %in% minor_range_rows)
      # calculate the densities of these markers and update them
      chrom$density[major_range_rows[to_change]] <- calc_density(ld, 1)[to_change]
      # print(1)
      print((nrow(chrom) / length(var_sets[[i]])))
    }
  }
  chrom
}

################################################################################

pane2 <- 50
window2 <- pane2 * 2 + 1

prop <- 0.05
threshold <- 0.7
max_dist <- 1000

snp_sets <- mclapply(1:length(density), ld_prune, mc.cores = 7)

seqClose(barley)

write_rds(snp_sets, "/workspace/barley/results/snp_sets.rds")

png(
   "/workspace/barley/results/phys_vs_gen_pos_pruned.png",
  family = "Times New Roman", width = 220, height = 1400, pointsize = 5,
  units = "mm", res = 192
)
par(mfrow = c(7, 1))
for (i in 1:7) {
  plot(
    phys_pos[[i]][var_sets[[i]] %in% snp_sets[[i]]$id],
    gen_pos[[i]][var_sets[[i]] %in% snp_sets[[i]]$id]
  )
}
dev.off()

png(
   "/workspace/barley/results/phys_vs_sparsity_pruned.png",
  family = "Times New Roman", width = 220, height = 1400, pointsize = 5,
  units = "mm", res = 192
)
par(mfrow = c(7, 1))
for (i in 1:7) {
  plot(
    phys_pos[[i]][var_sets[[i]] %in% snp_sets[[i]]$id],
    sparsity[[i]][var_sets[[i]] %in% snp_sets[[i]]$id]
  )
}
dev.off()