import::from(countrycode, "countrycode")
import::from(gdsfmt, "index.gdsn", "read.gdsn")
import::from(magrittr, "%>%")
import::from(readr, "read_rds", "read_tsv", "write_rds")
import::from(SeqArray, "seqClose", "seqOpen")
import::from(SNPRelate, "snpgdsAdmixPlot", "snpgdsAdmixProp")
source("/workspace/repos/barley-pgda/R/colours.R")

pca <- read_rds("/workspace/barley/results/ld_prunded_pca.rds")
hdbscan <- read_rds("/workspace/barley/results/hdbscan.rds")
passp <- read_tsv("/workspace/barley/barley_GBS_SNP_data/mascher@IPK-GATERSLEBEN.DE/variant_matrices/barley_GBS_SNP_data/180906_sample_information.tsv")

barley <- seqOpen("/workspace/barley/gds/imputed.gds")

barley_cults <- match(as.character(read.gdsn(index.gdsn(barley, "sample.id"))), passp$ENA_accession_id)
regions <- countrycode(passp$country_of_origin[barley_cults], "iso3c", "region") %>% as.factor()

sample_ids <- as.character(read.gdsn(index.gdsn(barley, "sample.id")))
samples_by_cluster <- split(sample_ids, hdbscan$cluster)

admix_prop <- snpgdsAdmixProp(pca, samples_by_cluster)

seqClose(barley)

png(
   "/workspace/barley/results/admix_prop.png",
  family = "Times New Roman", width = 1000, height = 1000, pointsize = 5,
  units = "mm", res = 192
)
snpgdsAdmixPlot(admix_prop, group = hdbscan$cluster, col = clrs[regions])
dev.off()