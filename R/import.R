library(countrycode)
library(SeqArray)
library(SNPRelate)
library(magrittr, "%>%")
import::from(htmlwidgets, "saveWidget")
import::from(readr, "read_tsv", "read_rds", "write_rds")
import::from(stringr, "str_c")
source("/workspace/repos/barley-pgda/R/colours.R")

seqVCF2GDS(
    str_c(
        "/workspace/barley/barley_GBS_SNP_data/mascher@IPK-GATERSLEBEN.DE/",
        "variant_matrices/barley_GBS_SNP_data/",
        "180606_GBS_imputed_matrix_for_GWAS_20458_samples_306049_SNPs.vcf.gz"
    ),
    "/workspace/barley/gds/imputed.gds"
)

barley <- seqOpen("/workspace/barley/gds/imputed.gds")

set.seed(1000)

# may try different LD thresholds for sensitivity analysis
snpset <- snpgdsLDpruning(barley, maf = 0.05, missing.rate = 0.10, 
  slide.max.bp = 1e7, ld.threshold = 0.7)
snpset.id <- unlist(snpset)

pca <- snpgdsPCA(barley, snp.id=snpset.id, num.thread = 8)

write_rds(pca, path = "/workspace/barley/results/pca.rds")
pca <- read_rds("/workspace/barley/results/pca.rds")

passp <- read_tsv("/workspace/barley/barley_GBS_SNP_data/mascher@IPK-GATERSLEBEN.DE/variant_matrices/barley_GBS_SNP_data/180906_sample_information.tsv")

barley_cults <- match(as.character(read.gdsn(index.gdsn(barley, "sample.id"))), passp$ENA_accession_id)

regions <- countrycode(passp$country_of_origin[barley_cults], "iso3c", "region") %>% as.factor()
png(
   "/workspace/barley/results/pca.png",
  family = "Times New Roman", width = 500, height = 500, pointsize = 5,
  units = "mm", res = 192
)
plot(pca, eig=1:4, pch=20, col = clrs[regions])
dev.off()
dim(pca$eigenvect)

scatter <- plot_ly() %>%
  add_markers(
    x = ~pca$eigenvect[, 1], y = ~pca$eigenvect[, 2], z = ~pca$eigenvect[, 3],
    color = regions, colors = clrs, marker = list(size = 5)
  )

saveWidget(
  as_widget(scatter),
  "/workspace/barley/results/pca.html"
)

barley_diss <- snpgdsDiss(barley, snp.id = snpset.id, num.thread = 8)
write_rds(barley_diss, path = "/workspace/barley/results/barley_diss.rds")
barley_hc <- snpgdsHCluster(barley_diss, sample.id = NULL, need.mat = FALSE, hang = 0.25)
write_rds(barley_hc, path = "/workspace/barley/results/barley_hc.rds")

png(
   "/workspace/barley/results/pca.png",
  family = "Times New Roman", width = 500, height = 500, pointsize = 5,
  units = "mm", res = 192
)
plot(pca, eig=1:4, pch=20, col = clrs[regions])
dev.off()
dim(pca$eigenvect)


knn_cl <- kmeans(barley_diss$diss, 1000, iter.max=20)
names()

dim(knn_cl$centers)