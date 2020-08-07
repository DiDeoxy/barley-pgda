import::from(countrycode, "countrycode")
import::from(htmlwidgets, "saveWidget")
import::from(gdsfmt, "index.gdsn", "read.gdsn")
import::from(magrittr, "%>%")
import::from(plotly, "add_markers", "as_widget", "plot_ly")
import::from(readr, "read_rds", "write_rds", "read_tsv")
import::from(SeqArray, "seqClose", "seqOpen")
import::from(SNPRelate, "snpgdsGRM")
source("/workspace/repos/barley-pgda/R/colours.R")

snp_set <- (read_rds("/workspace/barley/results/snp_sets.rds") %>% do.call(rbind, .))$id
hdbscan <- read_rds("/workspace/barley/results/hdbscan.rds")
grm_hdbscan <- read_rds("/workspace/barley/results/grm_hdbscan_10.rds")
passp <- read_tsv("/workspace/barley/barley_GBS_SNP_data/mascher@IPK-GATERSLEBEN.DE/variant_matrices/barley_GBS_SNP_data/180906_sample_information.tsv")

barley <- seqOpen("/workspace/barley/gds/imputed.gds")
barley_cults <- match(as.character(read.gdsn(index.gdsn(barley, "sample.id"))), passp$ENA_accession_id)
regions <- countrycode(passp$country_of_origin[barley_cults], "iso3c", "region") %>% as.factor()

# barley_ib_grm <- snpgdsGRM(
#     barley, snp.id = snp_set, method = "IndivBeta", num.thread = 7
# )

seqClose(barley)

write_rds(barley_ib_grm, "/workspace/barley/results/ib_grm.rds")

barley_grm_dist <- as.dist(1 - (barley_ib_grm$grm / max(barley_ib_grm$grm)))

write_rds(barley_grm_dist, "/workspace/barley/results/grm_dist.rds")

barley_ib_grm_eig <- eigen(barley_ib_grm$grm)

# write_rds(barley_ib_grm_eig, "/workspace/barley/results/eigen.rds")

barley_ib_grm_eig <- read_rds("/workspace/barley/results/eigen.rds")

scatter <- plot_ly() %>%
  add_markers(
    x = ~barley_ib_grm_eig$vectors[, 1], y = ~barley_ib_grm_eig$vectors[, 2], z = ~barley_ib_grm_eig$vectors[, 3],
    color = regions, colors = clrs, marker = list(size = 5)
  )

saveWidget(
  as_widget(scatter),
  "/workspace/barley/results/grm_eig_region.html"
)

scatter <- plot_ly() %>%
  add_markers(
    x = ~barley_ib_grm_eig$vectors[, 1], y = ~barley_ib_grm_eig$vectors[, 2], z = ~barley_ib_grm_eig$vectors[, 3],
    color = passp$annual_growth_habit[barley_cults], colors = clrs, marker = list(size = 5)
  )

saveWidget(
  as_widget(scatter),
  "/workspace/barley/results/grm_eig_growth_habit.html"
)

scatter <- plot_ly() %>%
  add_markers(
    x = ~barley_ib_grm_eig$vectors[, 1], y = ~barley_ib_grm_eig$vectors[, 2], z = ~barley_ib_grm_eig$vectors[, 3],
    color = passp$row_type[barley_cults], colors = clrs, marker = list(size = 5)
  )

saveWidget(
  as_widget(scatter),
  "/workspace/barley/results/grm_eig_row_type.html"
)

scatter <- plot_ly() %>%
  add_markers(
    x = ~barley_ib_grm_eig$vectors[, 1], y = ~barley_ib_grm_eig$vectors[, 2], z = ~barley_ib_grm_eig$vectors[, 3],
    color = hdbscan$cluster, colors = clrs, marker = list(size = 5)
  )

saveWidget(
  as_widget(scatter),
  "/workspace/barley/results/grm_eig_hdbscan.html"
)

scatter <- plot_ly() %>%
  add_markers(
    x = ~barley_ib_grm_eig$vectors[, 1], y = ~barley_ib_grm_eig$vectors[, 2], z = ~barley_ib_grm_eig$vectors[, 3],
    color = grm_hdbscan$cluster, colors = clrs, marker = list(size = 5)
  )

saveWidget(
  as_widget(scatter),
  "/workspace/barley/results/grm_eig_grm_hdbscan_10.html"
)