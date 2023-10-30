library(data.table)
library(Biostrings)

pyke <- fread("pyke_filtered.csv")
alleles <- unique(pyke[, allele])

multiallelic <- fread("multiallelic_data.csv")
multiallelic_alleles <- unique(multiallelic[, allele])

human_proteome <- readAAStringSet("./uniprot-compressed_true_download_true_format_fasta_query__28_28prote-2022.11.20-21.11.23.09.fasta")
protein_sequences <- as.data.frame(human_proteome)$x
sequence_length <- c(8, 9, 10, 11)
i <- 0 
negatives_list <- lapply(protein_sequences, function(sequence) {
  len <- sample(sequence_length, 1)
  peptide_list <- strsplit(sequence, paste0("(?<=.{", len , "})"), perl = TRUE)[[1]]
  peptide_list <- peptide_list[!grepl(c("U|X"), peptide_list)]
  peptide_list <- peptide_list[nchar(peptide_list)==max(nchar(peptide_list))]
  i <<- i + 1
  message(paste(i, "/", length(protein_sequences)),"\r",appendLF=FALSE)
  peptide_list
})
negatives_list <- unique(unlist(negatives_list))
negatives_list <- negatives_list[nchar(negatives_list) > 7 & nchar(negatives_list) < 12]
netagives_list_8 <- unique(sample(negatives_list[nchar(negatives_list) == 8], 500000))
netagives_list_9 <- unique(sample(negatives_list[nchar(negatives_list) == 9], 500000))
netagives_list_10 <- unique(sample(negatives_list[nchar(negatives_list) == 10], 500000))
netagives_list_11 <- unique(sample(negatives_list[nchar(negatives_list) == 11], 500000))
netagives_list_8 <- rep(netagives_list_8, length(alleles))
netagives_list_9 <- rep(netagives_list_9, length(alleles))
netagives_list_10 <- rep(netagives_list_10, length(alleles))
netagives_list_11 <- rep(netagives_list_11, length(alleles))

alleles <- rep(alleles, each=500000)
motif_peptides_8 <- data.table(peptide = netagives_list_8, allele = alleles)
motif_peptides_9 <- data.table(peptide = netagives_list_9, allele = alleles)
motif_peptides_10 <- data.table(peptide = netagives_list_10, allele = alleles)
motif_peptides_11 <- data.table(peptide = netagives_list_11, allele = alleles)
fwrite(motif_peptides_8, "./motif_data8.csv")
fwrite(motif_peptides_9, "./motif_data9.csv")
fwrite(motif_peptides_10, "./motif_data10.csv")
fwrite(motif_peptides_11, "./motif_data11.csv")