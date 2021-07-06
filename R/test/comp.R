pheno <- read.csv("/workspace/barley/phenotyping/GRINcore_2012_Agron_Aberdeen.csv", stringsAsFactors = FALSE)
all_grin <- read.csv("/workspace/barley/phenotyping/all_grin.csv", header = FALSE, stringsAsFactors = FALSE)
gbs <- read.csv("/workspace/barley/barley_GBS_SNP_data/mascher@IPK-GATERSLEBEN.DE/variant_matrices/barley_GBS_SNP_data/180906_sample_information.tsv", sep = "\t", stringsAsFactors = FALSE)
usda <- read.csv("/workspace/barley/phenotyping/journal.pone.0094688.s005.csv", stringsAsFactors = FALSE)
spring <- read.delim("/workspace/barley/phenotyping/spring_barley.tsv", sep = " ", stringsAsFactors = FALSE, header = FALSE)


names(pheno)
names(gbs)
names(usda)
usda$Accession

pheno$accession_name %in% all_grin$V1

pheno$accession_name[pheno$accession_name %in% gbs$name]
all_grin$V1[all_grin$V1 %in% gbs$name]
spring$V2[spring$V2 %in% gbs$name]

approx <- lapply(spring$V2[which(spring$V2 != "")], function (name) {
    matches <- agrep(name, gbs$name)
    if (length(matches) > 0 && length(matches) <= 5) {
        paste(gbs$name[matches], collapse = "/")
    }
})

approx <- lapply(all_grin$V1[which(all_grin$V1 != "")], function (name) {
    matches <- agrep(name, gbs$name)
    if (length(matches) > 0 && length(matches) <= 5) {
        paste(gbs$name[matches], collapse = "/")
    }
})

matched <- all_grin$V1[which(unlist(lapply(approx, function(matches) { if (length(matches)) { TRUE } else { FALSE } })))]

matches <- unlist(Filter(length, approx))

cbind(matched, matches)

write.table(cbind(matched, matches), "/workspace/barley/phenotyping/grin_fuzzy.tsv", sep = "\t", row.names = FALSE, col.names = FALSE)

currated_grin_gbs <- read.delim("/workspace/barley/phenotyping/grin_fuzzy copy.tsv", stringsAsFactors = FALSE)
write.table(currated_grin_gbs$Grin, "/workspace/barley/phenotyping/grin_only_fuzzy.tsv", quote = FALSE, col.names = FALSE, row.names = FALSE)

sum(usda$Accession %in% gbs$name)

gbs_names <- as.character(gbs$name[which(gbs$country_of_origin == "USA")])
gbs_names <- gbs_names[which(gbs_names != "")]

write.csv(gbs_names, "/workspace/barley/gbs_accessions.csv", quote = FALSE, row.names = FALSE)