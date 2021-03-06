import::from(countrycode, "countrycode")
import::from(htmlwidgets, "saveWidget")
import::from(gdsfmt, "index.gdsn", "read.gdsn")
import::from(magrittr, "%>%")
import::from(plotly, "add_markers", "as_widget", "plot_ly")
import::from(readr, "read_rds", "write_rds", "read_tsv")
import::from(SeqArray, "seqClose", "seqOpen")
import::from(SNPRelate,"snpgdsPCA")
source("/workspace/repos/barley-pgda/R/colours.R")

snp_set <- (read_rds("/workspace/barley/results/snp_sets.rds") %>% do.call(rbind, .))$id
hdbscan <- read_rds("/workspace/barley/results/hdbscan.rds")
passp <- read_tsv("/workspace/barley/barley_GBS_SNP_data/mascher@IPK-GATERSLEBEN.DE/variant_matrices/barley_GBS_SNP_data/180906_sample_information.tsv")

barley <- seqOpen("/workspace/barley/gds/imputed.gds")
barley_cults <- match(as.character(read.gdsn(index.gdsn(barley, "sample.id"))), passp$ENA_accession_id)
regions <- countrycode(passp$country_of_origin[barley_cults], "iso3c", "region") %>% as.factor()

# pca <- snpgdsPCA(barley, snp.id = snp_set, num.thread = 8)

seqClose(barley)

# write_rds(pca, file = "/workspace/barley/results/ld_prunded_pca.rds")

pca <- read_rds("/workspace/barley/results/ld_prunded_pca.rds")

png()
plot(cumsum(pca$varprop))
head(pca$varprop)

################################################################################

scatter <- plot_ly() %>%
  add_markers(
    x = ~pca$eigenvect[, 1], y = ~pca$eigenvect[, 2], z = ~pca$eigenvect[, 3],
    color = regions, colors = clrs, marker = list(size = 5)
  )

saveWidget(
  as_widget(scatter),
  "/workspace/barley/results/pca_region.html"
)

scatter <- plot_ly() %>%
  add_markers(
    x = ~pca$eigenvect[, 1], y = ~pca$eigenvect[, 2], z = ~pca$eigenvect[, 3],
    color = passp$annual_growth_habit[barley_cults], colors = clrs, marker = list(size = 5)
  )

saveWidget(
  as_widget(scatter),
  "/workspace/barley/results/pca_growth_habit.html"
)

scatter <- plot_ly() %>%
  add_markers(
    x = ~pca$eigenvect[, 1], y = ~pca$eigenvect[, 2], z = ~pca$eigenvect[, 3],
    color = passp$row_type[barley_cults], colors = clrs, marker = list(size = 5)
  )

saveWidget(
  as_widget(scatter),
  "/workspace/barley/results/pca_row_type.html"
)

scatter <- plot_ly() %>%
  add_markers(
    x = ~pca$eigenvect[, 1], y = ~pca$eigenvect[, 2], z = ~pca$eigenvect[, 3],
    color = hdbscan$cluster, colors = clrs, marker = list(size = 5)
  )

saveWidget(
  as_widget(scatter),
  "/workspace/barley/results/pca_hdbscan.html"
)