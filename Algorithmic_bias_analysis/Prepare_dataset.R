library(data.table)
library(Biostrings)

pyke <- fread("./Pyke.csv")
pyke <- pyke[modification == "na"]
pyke[, label := 1]
fwrite(pyke, "./pyke_filtered.csv")

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
negatives_list <- sample(negatives_list, dim(pyke)[1]*50)
negative_alleles <- pyke[, allele]
negative_alleles <- rep(negative_alleles, each=50)
negatives <- data.table(peptide = negatives_list, allele = negative_alleles)
negatives[, label := 0]
testing_data <- rbindlist(list(pyke[, .(peptide, allele, label)], negatives))
testing_data <- testing_data[nchar(peptide) > 7 & nchar(peptide) < 12]
fwrite(testing_data, "./testing_data.csv")