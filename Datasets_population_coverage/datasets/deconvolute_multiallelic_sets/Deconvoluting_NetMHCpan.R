library(data.table)

cell_lines <- list.files("./Deconvoluted_cell_lines/")
resss <- lapply(cell_lines, function(cell_line) {
  alleles <- list.files(paste0("./Deconvoluted_cell_lines/", cell_line, "/netmhcpan_res/"))
  DT_list <- lapply(alleles, function(allele) {
    res <- fread(paste0("./Deconvoluted_cell_lines/", cell_line, "/netmhcpan_res/", allele))[, .(V2, V4)]
    allele <- strsplit(allele, "_")[[1]][1]
    print(allele)
    names(res) <- c("peptide", allele)
    res
  })
  MergedDT <- Reduce(function(...) merge(..., by = c("peptide")), DT_list)
  MergedDT <- melt.data.table(MergedDT, id.vars = c("peptide"), variable.name = "allele")
  MergedDT <- MergedDT[MergedDT[, .I[value == max(value)], by=.(peptide)]$V1][, .(peptide, allele)]
  MergedDT
})
resss <- rbindlist(resss)
fwrite(resss, "./Deconvoluted_NetMHCpan_multiallelic_training_dataset.csv")