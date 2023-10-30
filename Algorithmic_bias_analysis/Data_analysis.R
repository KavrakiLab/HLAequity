library(data.table)
library(yardstick)
library(magrittr)
library(tidyverse)
library(cutpointr)
library(caret)
library(rnaturalearthdata)
library(RColorBrewer)
library(ggplot2)
library(viridis)
theme_set(theme_bw())

# MHCFlurry analysis
mhcflurry_preds <- fread("mhcflurry_preds.csv")
positives <- dim(mhcflurry_preds[label == 1])[1]
mhcflurry_preds[, label := as.factor(label)]
alleles <- unique(mhcflurry_preds[, allele])

# Rank analysis (FOOP)
mhcflurry_motifs <- fread("mhcflurry_motifs.csv")
mhcflurry_motifs <- mhcflurry_motifs[, .(peptide, allele, mhcflurry_presentation_score)]
mhcflurry_motifs <- mhcflurry_motifs[order(mhcflurry_presentation_score), .(peptide, allele, mhcflurry_presentation_score), by = allele]
mhcflurry_motifs <- split(mhcflurry_motifs,mhcflurry_motifs$allele)

mhcflurry_preds <- fread("mhcflurry_preds.csv")[label == 1]
mhcflurry_preds <- mhcflurry_preds[, .(peptide, allele, mhcflurry_presentation_score)]

i <- 0
rank <- by(mhcflurry_preds, seq_len(nrow(mhcflurry_preds)), function(row) {
  i <<- i + 1
  message(paste(i, "/", nrow(mhcflurry_preds)),"\r",appendLF=FALSE)
  (500000 - findInterval(row$mhcflurry_presentation_score, mhcflurry_motifs[[row$allele]]$mhcflurry_presentation_score))/5000
})
mhcflurry_preds[, rank := as.numeric(rank)]

mhcflurry_preds_list <- split(mhcflurry_preds, mhcflurry_preds$allele)
allele_order <- names(mhcflurry_preds_list)
fraction <- sapply(mhcflurry_preds_list, function(DT) {
  size <- dim(DT)[1]
  observed_peptides <- dim(DT[rank < 0.1])[1]
  observed_peptides/size
})
fraction_of_observed_peptides <- data.table(allele = allele_order, fraction_of_observed_peptides = fraction)

mhcflurry_multiallelic <- fread("mhcflurry_multiallelic.csv")
mhcflurry_multiallelic <- mhcflurry_multiallelic[, .(peptide, allele, mhcflurry_presentation_score)]
i <- 0
rank <- by(mhcflurry_multiallelic, seq_len(nrow(mhcflurry_multiallelic)), function(row) {
  i <<- i + 1
  message(paste(i, "/", nrow(mhcflurry_multiallelic)),"\r",appendLF=FALSE)
  (500000 - findInterval(row$mhcflurry_presentation_score, mhcflurry_motifs[[row$allele]]$mhcflurry_presentation_score))/5000
})
mhcflurry_multiallelic[, rank := as.numeric(rank)]
mhcflurry_multiallelic_filtered <- mhcflurry_multiallelic[rank <= 0.5]
mhcflurry_multiallelic_filtered <- mhcflurry_multiallelic_filtered[mhcflurry_multiallelic_filtered[, .I[which.min(rank)], by = peptide]$V1]

mhcflurry_whole <- rbindlist(list(mhcflurry_preds, mhcflurry_multiallelic_filtered))
mhcflurry_whole_list <- split(mhcflurry_whole, mhcflurry_whole$allele)
allele_order <- names(mhcflurry_whole_list)
fraction <- sapply(mhcflurry_whole_list, function(DT) {
  size <- dim(DT)[1]
  observed_peptides <- dim(DT[rank < 0.1])[1]
  observed_peptides/size
})

fraction_of_observed_peptides_multiallelic <- data.table(allele = allele_order, fraction_of_observed_peptides = fraction)
fraction_of_observed_peptides_multiallelic[, allele := substr(allele, 5, nchar(allele))]
fraction_of_observed_peptides_multiallelic[, Allele := paste0(substr(allele, 1, 1), '*', substr(allele, 2, nchar(allele)))]
fraction_of_observed_peptides_multiallelic[, allele := NULL]


# PPV Analysis
mhcflurry_preds <- fread("mhcflurry_preds.csv")[label == 1]
mhcflurry_multiallelic_filtered[, label := 1]
mhcflurry_motifs <- fread("mhcflurry_motifs.csv")
mhcflurry_motifs <- mhcflurry_motifs[, .(peptide, allele, mhcflurry_presentation_score)]
mhcflurry_motifs[, label := 0]
mhcflurry_whole <- rbindlist(list(mhcflurry_preds[, .(peptide, allele, mhcflurry_presentation_score, label)],
                                  mhcflurry_multiallelic_filtered[, .(peptide, allele, mhcflurry_presentation_score, label)],
                                  mhcflurry_motifs))
alleles <- unique(mhcflurry_whole[, allele])
ppv_list_presentation_all <- sapply(alleles, function(x) {
  DT <- mhcflurry_whole[allele == x]
  positives <- dim(DT[label == 1])[1]
  ppv <- ((DT[order(-mhcflurry_presentation_score)] %>% head(positives))[, label] %>% as.character %>% as.integer %>% sum)/positives
  ppv
})
ppv_list_presentation_all <- data.table(allele = alleles, ppv = ppv_list_presentation_all)
ppv_list_presentation_all[, allele := substr(allele, 5, nchar(allele))]
ppv_list_presentation_all[, Allele := paste0(substr(allele, 1, 1), '*', substr(allele, 2, nchar(allele)))]
ppv_list_presentation_all[, allele := NULL]

numerical_results <- merge(ppv_list_presentation_all, fraction_of_observed_peptides_multiallelic, on = 'Allele')
numerical_results <- melt.data.table(numerical_results, id.vars = "Allele")
numerical_results[, locus := substr(Allele, 1, 1)]

rm(mhcflurry_motifs)
rm(mhcflurry_preds)
rm(mhcflurry_multiallelic)
rm(mhcflurry_whole)
rm(mhcflurry_multiallelic_filtered)
rm(mhcflurry_preds_list)
rm(mhcflurry_whole_list)
rm(rank)
gc()
gc()

ggplot(numerical_results, aes(x=Allele, y=value, fill=locus)) +
  geom_bar(stat = "identity")+
  facet_grid(rows = vars(variable), cols = vars(locus), scales = "free") +
  scale_fill_viridis_d() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12, vjust = 0.5), plot.title = element_text(hjust = 0.5, size = 20, face = "bold"), 
        strip.text.x = element_text(size = 16, face = "italic"), strip.text.y = element_text(size = 16, face = "italic"),
        axis.title.y = element_text(size = rel(1.3)), axis.title.x = element_text(size = rel(1.3)),
        axis.text.y = element_text(angle = 0, hjust = 1, size = 12)) +
  xlab("")

# AUPRC -PPV analysis
mhcflurry_preds <- fread("mhcflurry_preds.csv")[, .(peptide, allele, label)]
mhcflurry_preds[, allele := paste0(substr(allele, 1, 5), '*', substr(allele, 6, 10))]
mhcflurry_preds[, label := as.factor(label)]
netmhcpan_preds <- fread("netmhcpan_preds.csv")
netmhcpan_preds <- na.omit(netmhcpan_preds[mhcflurry_preds, on = c("peptide", "allele"), nomatch = 0])
positives <- dim(netmhcpan_preds[label == 1])[1]
alleles <- unique(netmhcpan_preds[, allele])

# Rank analysis
netmhcpan_motifs <- fread("netmhcpan_motifs.csv")
netmhcpan_motifs_BA <- netmhcpan_motifs[order(BA), .(peptide, allele, BA), by = allele]
netmhcpan_motifs_BA <- split(netmhcpan_motifs_BA, netmhcpan_motifs_BA$allele)
netmhcpan_motifs_EL <- netmhcpan_motifs[order(EL), .(peptide, allele, EL), by = allele]
netmhcpan_motifs_EL <- split(netmhcpan_motifs_EL, netmhcpan_motifs_EL$allele)

netmhcpan_preds <- fread("netmhcpan_preds.csv")
netmhcpan_preds <- na.omit(netmhcpan_preds[mhcflurry_preds, on = c("peptide", "allele"), nomatch = 0])
netmhcpan_preds <- netmhcpan_preds[label == 1]
netmhcpan_preds <- netmhcpan_preds[, .(peptide, allele, BA, EL)]

i <- 0
rank_BA <- by(netmhcpan_preds, seq_len(nrow(netmhcpan_preds)), function(row) {
  i <<- i + 1
  message(paste(i, "/", nrow(netmhcpan_preds)),"\r",appendLF=FALSE)
  (500000 - findInterval(row$BA, netmhcpan_motifs_BA[[row$allele]]$BA))/5000
})
i <- 0
rank_EL <- by(netmhcpan_preds, seq_len(nrow(netmhcpan_preds)), function(row) {
  i <<- i + 1
  message(paste(i, "/", nrow(netmhcpan_preds)),"\r",appendLF=FALSE)
  (500000 - findInterval(row$EL, netmhcpan_motifs_EL[[row$allele]]$EL))/5000
})
netmhcpan_preds[, rank_BA := as.numeric(rank_BA)]
netmhcpan_preds[, rank_EL := as.numeric(rank_EL)]

netmhcpan_preds_list <- split(netmhcpan_preds, netmhcpan_preds$allele)
allele_order <- names(netmhcpan_preds_list)
fraction_BA <- sapply(netmhcpan_preds_list, function(DT) {
  size <- dim(DT)[1]
  observed_peptides <- dim(DT[rank_BA < 0.1])[1]
  observed_peptides/size
})
fraction_EL <- sapply(netmhcpan_preds_list, function(DT) {
  size <- dim(DT)[1]
  observed_peptides <- dim(DT[rank_EL < 0.1])[1]
  observed_peptides/size
})
fraction_of_observed_peptides <- data.table(allele = allele_order, fraction_of_observed_peptides_BA = fraction_BA, fraction_of_observed_peptides_EL = fraction_EL)

netmhcpan_multiallelic <- fread("netmhcpan_multiallelic.csv")
i <- 0
rank_BA <- by(netmhcpan_multiallelic, seq_len(nrow(netmhcpan_multiallelic)), function(row) {
  i <<- i + 1
  message(paste(i, "/", nrow(netmhcpan_multiallelic)),"\r",appendLF=FALSE)
  (500000 - findInterval(row$BA, netmhcpan_motifs_BA[[row$allele]]$BA))/5000
})
i <- 0
rank_EL <- by(netmhcpan_multiallelic, seq_len(nrow(netmhcpan_multiallelic)), function(row) {
  i <<- i + 1
  message(paste(i, "/", nrow(netmhcpan_multiallelic)),"\r",appendLF=FALSE)
  (500000 - findInterval(row$EL, netmhcpan_motifs_EL[[row$allele]]$EL))/5000
})
netmhcpan_multiallelic[, rank_BA := as.numeric(rank_BA)]
netmhcpan_multiallelic[, rank_EL := as.numeric(rank_EL)]
netmhcpan_multiallelic_filtered_BA <- netmhcpan_multiallelic[rank_BA <= 0.5]
netmhcpan_multiallelic_filtered_EL <- netmhcpan_multiallelic[rank_EL <= 0.5]
netmhcpan_multiallelic_filtered_BA <- netmhcpan_multiallelic_filtered_BA[netmhcpan_multiallelic_filtered_BA[, .I[which.min(rank_BA)], by = peptide]$V1]
netmhcpan_multiallelic_filtered_EL <- netmhcpan_multiallelic_filtered_EL[netmhcpan_multiallelic_filtered_EL[, .I[which.min(rank_EL)], by = peptide]$V1]

netmhcpan_whole_BA <- rbindlist(list(netmhcpan_preds, netmhcpan_multiallelic_filtered_BA), use.names=TRUE)
netmhcpan_whole_list_BA <- split(netmhcpan_whole_BA, netmhcpan_whole_BA$allele)
allele_order <- names(netmhcpan_whole_list_BA)
fraction_BA <- sapply(netmhcpan_whole_list_BA, function(DT) {
  size <- dim(DT)[1]
  observed_peptides <- dim(DT[rank_BA < 0.1])[1]
  observed_peptides/size
})
fraction_of_observed_peptides_multiallelic_BA <- data.table(Allele = allele_order, fraction_of_observed_peptides = fraction_BA)
fraction_of_observed_peptides_multiallelic_BA[, Allele := substr(Allele, 5, nchar(Allele))]

rm(netmhcpan_motifs)
rm(netmhcpan_preds_list)
rm(netmhcpan_motifs_BA)
rm(netmhcpan_motifs_EL)
rm(netmhcpan_multiallelic)
rm(netmhcpan_whole_BA)
rm(netmhcpan_whole_list_BA)
rm(rank_BA)
rm(rank_EL)
gc()
gc()

netmhcpan_preds_ppv <- fread("netmhcpan_preds.csv")
netmhcpan_preds_ppv <- na.omit(netmhcpan_preds_ppv[mhcflurry_preds, on = c("peptide", "allele"), nomatch = 0])[label == 1]
netmhcpan_multiallelic_filtered_BA[, label := 1]
netmhcpan_motifs_BA <- fread("netmhcpan_motifs.csv")
netmhcpan_motifs_BA <- netmhcpan_motifs_BA[, .(peptide, allele, BA)]
netmhcpan_motifs_BA[, label := 0]
netmhcpan_whole_BA <- rbindlist(list(netmhcpan_preds_ppv[, .(peptide, allele, BA, label)],
                                     netmhcpan_multiallelic_filtered_BA[, .(peptide, allele, BA, label)],
                                     netmhcpan_motifs_BA))
alleles <- unique(netmhcpan_whole_BA[, allele])
rm(netmhcpan_multiallelic_filtered_BA)
rm(netmhcpan_motifs_BA)
rm(netmhcpan_preds_ppv)
gc()
gc()
ppv_list_presentation_all_BA <- sapply(alleles, function(x) {
  DT <- netmhcpan_whole_BA[allele == x]
  positives <- dim(DT[label == 1])[1]
  ppv <- ((DT[order(-BA)] %>% head(positives))[, label] %>% as.character %>% as.integer %>% sum)/positives
  ppv
})
ppv_list_presentation_all_BA <- data.table(Allele = alleles, ppv = ppv_list_presentation_all_BA)
ppv_list_presentation_all_BA[, Allele := substr(Allele, 5, nchar(Allele))]

netmhcpan_whole_EL <- rbindlist(list(netmhcpan_preds, netmhcpan_multiallelic_filtered_EL), use.names=TRUE)
netmhcpan_whole_list_EL <- split(netmhcpan_whole_EL, netmhcpan_whole_EL$allele)
allele_order <- names(netmhcpan_whole_list_EL)
fraction_EL <- sapply(netmhcpan_whole_list_EL, function(DT) {
  size <- dim(DT)[1]
  observed_peptides <- dim(DT[rank_EL < 0.1])[1]
  observed_peptides/size
})
fraction_of_observed_peptides_multiallelic_EL <- data.table(Allele = allele_order, fraction_of_observed_peptides = fraction_EL)
fraction_of_observed_peptides_multiallelic_EL[, Allele := substr(Allele, 5, nchar(Allele))]

netmhcpan_preds_ppv <- fread("netmhcpan_preds.csv")
netmhcpan_preds_ppv <- na.omit(netmhcpan_preds_ppv[mhcflurry_preds, on = c("peptide", "allele"), nomatch = 0])[label == 1]
netmhcpan_multiallelic_filtered_EL[, label := 1]
netmhcpan_motifs <- fread("netmhcpan_motifs.csv")
netmhcpan_motifs_EL <- netmhcpan_motifs[, .(peptide, allele, EL)]
netmhcpan_motifs_EL[, label := 0]
netmhcpan_whole_EL <- rbindlist(list(netmhcpan_preds_ppv[, .(peptide, allele, EL, label)],
                                     netmhcpan_multiallelic_filtered_EL[, .(peptide, allele, EL, label)],
                                     netmhcpan_motifs_EL))
alleles <- unique(netmhcpan_whole_EL[, allele])
rm(netmhcpan_multiallelic_filtered_EL)
rm(netmhcpan_motifs_EL)
rm(netmhcpan_preds_ppv)
gc()
gc()
ppv_list_presentation_all_EL <- sapply(alleles, function(x) {
  DT <- netmhcpan_whole_EL[allele == x]
  positives <- dim(DT[label == 1])[1]
  ppv <- ((DT[order(-EL)] %>% head(positives))[, label] %>% as.character %>% as.integer %>% sum)/positives
  ppv
})
ppv_list_presentation_all_EL <- data.table(Allele = alleles, ppv = ppv_list_presentation_all_EL)
ppv_list_presentation_all_EL[, Allele := substr(Allele, 5, nchar(Allele))]

numerical_results_BA <- merge(ppv_list_presentation_all_BA, fraction_of_observed_peptides_multiallelic_BA, on = 'Allele')
numerical_results_BA <- melt.data.table(numerical_results_BA, id.vars = "Allele")
numerical_results_BA[, locus := substr(Allele, 1, 1)]
ggplot(numerical_results_BA, aes(x=Allele, y=value, fill=locus)) +
  geom_bar(stat = "identity")+
  facet_grid(rows = vars(variable), cols = vars(locus), scales = "free") +
  scale_fill_manual(values = c("#554F66", "#C2D1A4", "#FA3456")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12, vjust = 0.5), plot.title = element_text(hjust = 0.5, size = 20, face = "bold"), 
        strip.text.x = element_text(size = 16, face = "italic"), strip.text.y = element_text(size = 16, face = "italic"),
        axis.title.y = element_text(size = rel(1.3)), axis.title.x = element_text(size = rel(1.3)),
        axis.text.y = element_text(angle = 0, hjust = 1, size = 12)) +
  xlab("")

numerical_results_EL <- merge(ppv_list_presentation_all_EL, fraction_of_observed_peptides_multiallelic_BA, on = 'Allele')
numerical_results_EL <- melt.data.table(numerical_results_EL, id.vars = "Allele")
numerical_results_EL[, locus := substr(Allele, 1, 1)]
ggplot(numerical_results_EL, aes(x=Allele, y=value, fill=locus)) +
  geom_bar(stat = "identity")+
  facet_grid(rows = vars(variable), cols = vars(locus), scales = "free") +
  scale_fill_viridis_d() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12, vjust = 0.5), plot.title = element_text(hjust = 0.5, size = 20, face = "bold"), 
        strip.text.x = element_text(size = 16, face = "italic"), strip.text.y = element_text(size = 16, face = "italic"),
        axis.title.y = element_text(size = rel(1.3)), axis.title.x = element_text(size = rel(1.3)),
        axis.text.y = element_text(angle = 0, hjust = 1, size = 12)) +
  xlab("")

numerical_results[, method := "MHCFlurry2.0"]
numerical_results[, paired := seq(1, 118)]
numerical_results_EL[, method := "NetMHCpan4.1"]
numerical_results_EL[, paired := seq(1, 118)]
numerical_results_whole <- rbindlist(list(numerical_results, numerical_results_EL))
numerical_results_whole[, variable := ifelse(variable == "fraction_of_observed_peptides", "FOOP", "PPV")]
numerical_results_whole[, variable := factor(variable, levels = c("PPV", "FOOP"))]
ggplot(numerical_results_whole, aes(x=Allele, y=value, shape=locus, colour=method)) +
  geom_point(size = 4) +
  geom_line(aes(group = paired), colour = 'black') +
  geom_vline(xintercept = 17.5, linetype="dashed", color = "blue", linewidth=1.5) +
  geom_vline(xintercept = 45.5, linetype="dashed", color = "blue", linewidth=1.5) +
  facet_wrap(vars(variable), nrow = 2) +
  scale_fill_viridis_d() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12, vjust = 0.5), plot.title = element_text(hjust = 0.5, size = 20, face = "bold"), 
        strip.text.x = element_text(size = 16, face = "italic"), strip.text.y = element_text(size = 16, face = "italic"),
        axis.title.y = element_text(size = rel(1.3)), axis.title.x = element_text(size = rel(1.3)),
        axis.text.y = element_text(angle = 0, hjust = 1, size = 12), legend.position="none") +
  xlab("") +
  ylab("")

head(numerical_results[, lapply(.SD,sum), by=.(Allele), .SDcols = c("value")][order(value)], 10)
head(numerical_results_EL[, lapply(.SD,sum), by=.(Allele), .SDcols = c("value")][order(value)], 10)

tail(numerical_results[, lapply(.SD,sum), by=.(Allele), .SDcols = c("value")][order(value)], 10)
tail(numerical_results_EL[, lapply(.SD,sum), by=.(Allele), .SDcols = c("value")][order(value)], 10)
#devtools::install_github("schliebs/ggoxford")
library(ggoxford)
library(ggtext)
A3001 <- data.table(
  Population = c("India West Coast Parsi", " Guinea-Bissau Bijago", "Switzerland Lucerne population", "Italy"),
  Country = c("India", "Guinea-Bissau", "Switzerland", "Italy"),
  Percentage = c(16.0, 14.1, 4.38, 3.40)
)
A3001[, iso3 := countrycode::countrycode(Country, origin = "country.name.en", destination = "iso3c")]
A3001[, iso3 := factor(iso3, levels = c("ITA", "CHE", "GNB", "IND"))]
A3001 <- tibble(A3001)
ggplot(data = A3001 ,
       aes(x = Percentage, y = iso3)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 2.5, linetype="dashed", color = "blue", linewidth=1.5) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, vjust = 0.5), plot.title = element_text(hjust = 0.5, size = 20, face = "bold"), 
        strip.text.x = element_text(size = 16, face = "italic"), strip.text.y = element_text(size = 16, face = "italic"),
        axis.title.y = element_text(size = rel(1.3)), axis.title.x = element_text(size = rel(1.3)),
        axis.text.y = element_text(angle = 0, hjust = 1, size = 12), plot.background = element_rect(colour = "black", fill=NA, linewidth=3)) +
  geom_axis_flags(breaks = A3001$iso3,
                  labels = A3001$Population,
                  country_icons = A3001$iso3,
                  width = 30,
                  axis = "y",
                  lineheight = 2,
                  fontface = "bold"
  ) +
  ylab("Population") +
  xlab("Percentage (%)")

A0252 <- data.table(
  Population = c("Iran Kurds", "Brazil Rio State Caucasians", "Spanish regions"),
  Country = c("Iran", "Brazil", "Spain"),
  Percentage = c(6.7, 0.12, 0.02)
)
A0252[, iso3 := countrycode::countrycode(Country, origin = "country.name.en", destination = "iso3c")]
A0252[, iso3 := factor(iso3, levels = c("ESP", "BRA", "IRN"))]
A0252 <- tibble(A0252)
ggplot(data = A0252 ,
       aes(x = Percentage, y = iso3)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 2.5, linetype="dashed", color = "blue", linewidth=1.5) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, vjust = 0.5), plot.title = element_text(hjust = 0.5, size = 20, face = "bold"), 
        strip.text.x = element_text(size = 16, face = "italic"), strip.text.y = element_text(size = 16, face = "italic"),
        axis.title.y = element_text(size = rel(1.3)), axis.title.x = element_text(size = rel(1.3)),
        axis.text.y = element_text(angle = 0, hjust = 1, size = 12), plot.background = element_rect(colour = "black", fill=NA, linewidth=3)) +
  geom_axis_flags(breaks = A0252$iso3,
                  labels = A0252$Population,
                  country_icons = A0252$iso3,
                  width = 30,
                  axis = "y",
                  lineheight = 2,
                  fontface = "bold"
  ) +
  ylab("Population") +
  xlab("Percentage (%)")

A2402 <- data.table(
  Population = c("Ecuador Cayapas", "New Caledonia", "Japan", "Romania"),
  Country = c("Ecuador", "New Caledonia", "Japan", "Romania"),
  Percentage = c(61.4, 60.7, 36.4, 23.8)
)
A2402[, iso3 := countrycode::countrycode(Country, origin = "country.name.en", destination = "iso3c")]
A2402[, iso3 := factor(iso3, levels = c("ROU", "JPN", "NCL", "ECU"))]
A2402 <- tibble(A2402)
ggplot(data = A2402 ,
       aes(x = Percentage, y = iso3)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 2.5, linetype="dashed", color = "blue", linewidth=1.5) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, vjust = 0.5), plot.title = element_text(hjust = 0.5, size = 20, face = "bold"), 
        strip.text.x = element_text(size = 16, face = "italic"), strip.text.y = element_text(size = 16, face = "italic"),
        axis.title.y = element_text(size = rel(1.3)), axis.title.x = element_text(size = rel(1.3)),
        axis.text.y = element_text(angle = 0, hjust = 1, size = 12), plot.background = element_rect(colour = "black", fill=NA, linewidth=3)) +
  geom_axis_flags(breaks = A2402$iso3,
                  labels = A2402$Population,
                  country_icons = A2402$iso3,
                  width = 30,
                  axis = "y",
                  lineheight = 2,
                  fontface = "bold"
  ) +
  ylab("Population") +
  xlab("Percentage (%)")

A0101 <- data.table(
  Population = c("India Tamil Nadu", "Kenya Nandi", "Ireland South", "USA Caucasian Population"),
  Country = c("India", "Kenya", "Ireland", "USA"),
  Percentage = c(15.93, 11.8, 33.33, 25)
)
A0101[, iso3 := countrycode::countrycode(Country, origin = "country.name.en", destination = "iso3c")]
A0101[, iso3 := factor(iso3, levels = c("USA", "IRL", "KEN", "IND"))]
A0101 <- tibble(A0101)
ggplot(data = A0101 ,
       aes(x = Percentage, y = iso3)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 2.5, linetype="dashed", color = "blue", linewidth=1.5) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, vjust = 0.5), plot.title = element_text(hjust = 0.5, size = 20, face = "bold"), 
        strip.text.x = element_text(size = 16, face = "italic"), strip.text.y = element_text(size = 16, face = "italic"),
        axis.title.y = element_text(size = rel(1.3)), axis.title.x = element_text(size = rel(1.3)),
        axis.text.y = element_text(angle = 0, hjust = 1, size = 12), plot.background = element_rect(colour = "black", fill=NA, linewidth=3)) +
  geom_axis_flags(breaks = A0101$iso3,
                  labels = A0101$Population,
                  country_icons = A0101$iso3,
                  width = 30,
                  axis = "y",
                  lineheight = 2,
                  fontface = "bold"
  ) +
  ylab("Population") +
  xlab("Percentage (%)")

C0401 <- data.table(
  Population = c("Mexico Oaxaca Mixe", "Australia Yuendumu Aborigine", "Greece", "Italy South"),
  Country = c("Mexico", "Australia", "Greece", "Italy"),
  Percentage = c(37.8, 26.8, 19.28, 18.88)
)
C0401[, iso3 := countrycode::countrycode(Country, origin = "country.name.en", destination = "iso3c")]
C0401[, iso3 := factor(iso3, levels = c("ITA", "GRC", "AUS", "MEX"))]
C0401 <- tibble(C0401)
ggplot(data = C0401 ,
       aes(x = Percentage, y = iso3)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 2.5, linetype="dashed", color = "blue", linewidth=1.5) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, vjust = 0.5), plot.title = element_text(hjust = 0.5, size = 20, face = "bold"), 
        strip.text.x = element_text(size = 16, face = "italic"), strip.text.y = element_text(size = 16, face = "italic"),
        axis.title.y = element_text(size = rel(1.3)), axis.title.x = element_text(size = rel(1.3)),
        axis.text.y = element_text(angle = 0, hjust = 1, size = 12), plot.background = element_rect(colour = "black", fill=NA, linewidth=3)) +
  geom_axis_flags(breaks = C0401$iso3,
                  labels = C0401$Population,
                  country_icons = C0401$iso3,
                  width = 30,
                  axis = "y",
                  lineheight = 2,
                  fontface = "bold"
  ) +
  ylab("Population") +
  xlab("Percentage (%)")

C0722 <- data.table(
  Population = c("Tanzania Maasai", "South Africa Caucasians", "Switzerland Geneva population"),
  Country = c("Tanzania", "South Africa", "Switzerland"),
  Percentage = c(0.78, 0.5, 0.31)
)
C0722[, iso3 := countrycode::countrycode(Country, origin = "country.name.en", destination = "iso3c")]
C0722[, iso3 := factor(iso3, levels = c("CHE", "ZAF", "TZA"))]
C0722 <- tibble(C0722)
ggplot(data = C0722,
       aes(x = Percentage, y = iso3)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 2.5, linetype="dashed", color = "blue", linewidth=1.5) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, vjust = 0.5), plot.title = element_text(hjust = 0.5, size = 20, face = "bold"), 
        strip.text.x = element_text(size = 16, face = "italic"), strip.text.y = element_text(size = 16, face = "italic"),
        axis.title.y = element_text(size = rel(1.3)), axis.title.x = element_text(size = rel(1.3)),
        axis.text.y = element_text(angle = 0, hjust = 1, size = 12), plot.background = element_rect(colour = "black", fill=NA, linewidth=3)) +
  geom_axis_flags(breaks = C0722$iso3,
                  labels = C0722$Population,
                  country_icons = C0722$iso3,
                  width = 30,
                  axis = "y",
                  lineheight = 2,
                  fontface = "bold"
  ) +
  ylab("Population") +
  xlab("Percentage (%)")

B3801 <- data.table(
  Population = c("Tunisia", "Thailand", "Israel Askhenazi Jews", "Italy North population"),
  Country = c("Tunisia", "Thailand", "Israel", "Italy"),
  Percentage = c(3.0, 2.8, 16.7, 6.7)
)
B3801[, iso3 := countrycode::countrycode(Country, origin = "country.name.en", destination = "iso3c")]
B3801[, iso3 := factor(iso3, levels = c("ITA", "ISR", "THA", "TUN"))]
B3801 <- tibble(B3801)
ggplot(data = B3801,
       aes(x = Percentage, y = iso3)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 2.5, linetype="dashed", color = "blue", linewidth=1.5) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, vjust = 0.5), plot.title = element_text(hjust = 0.5, size = 20, face = "bold"), 
        strip.text.x = element_text(size = 16, face = "italic"), strip.text.y = element_text(size = 16, face = "italic"),
        axis.title.y = element_text(size = rel(1.3)), axis.title.x = element_text(size = rel(1.3)),
        axis.text.y = element_text(angle = 0, hjust = 1, size = 12), plot.background = element_rect(colour = "black", fill=NA, linewidth=3)) +
  geom_axis_flags(breaks = B3801$iso3,
                  labels = B3801$Population,
                  country_icons = B3801$iso3,
                  width = 30,
                  axis = "y",
                  lineheight = 2,
                  fontface = "bold"
  ) +
  ylab("Population") +
  xlab("Percentage (%)")

B1513 <- data.table(
  Population = c("Malaysia Mandailing", "Indonesia Java population", "Netherlands Leiden Population"),
  Country = c("Malaysia", "Indonesia", "Netherlands"),
  Percentage = c(27.8, 12.5, 0.01)
)
B1513[, iso3 := countrycode::countrycode(Country, origin = "country.name.en", destination = "iso3c")]
B1513[, iso3 := factor(iso3, levels = c("NLD", "IDN", "MYS"))]
B1513 <- tibble(B1513)
ggplot(data = B1513,
       aes(x = Percentage, y = iso3)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 1.5, linetype="dashed", color = "blue", linewidth=1.5) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, vjust = 0.5), plot.title = element_text(hjust = 0.5, size = 20, face = "bold"), 
        strip.text.x = element_text(size = 16, face = "italic"), strip.text.y = element_text(size = 16, face = "italic"),
        axis.title.y = element_text(size = rel(1.3)), axis.title.x = element_text(size = rel(1.3)),
        axis.text.y = element_text(angle = 0, hjust = 1, size = 12), plot.background = element_rect(colour = "black", fill=NA, linewidth=3)) +
  geom_axis_flags(breaks = B1513$iso3,
                  labels = B1513$Population,
                  country_icons = B1513$iso3,
                  width = 30,
                  axis = "y",
                  lineheight = 2,
                  fontface = "bold"
  ) +
  ylab("Population") +
  xlab("Percentage (%)")

B0801 <- data.table(
  Population = c("Ireland South", "USA San Francisco Caucasian", "Oman", "India New Delhi population"),
  Country = c("Ireland", "USA", "Oman", "India"),
  Percentage = c(18.2, 15.1, 11.0, 8.3)
)
B0801[, iso3 := countrycode::countrycode(Country, origin = "country.name.en", destination = "iso3c")]
B0801[, iso3 := factor(iso3, levels = c("USA", "IRL", "IND", "OMN"))]
B0801 <- tibble(B0801)
ggplot(data = B0801,
       aes(x = Percentage, y = iso3)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 2.5, linetype="dashed", color = "blue", linewidth=1.5) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, vjust = 0.5), plot.title = element_text(hjust = 0.5, size = 20, face = "bold"), 
        strip.text.x = element_text(size = 16, face = "italic"), strip.text.y = element_text(size = 16, face = "italic"),
        axis.title.y = element_text(size = rel(1.3)), axis.title.x = element_text(size = rel(1.3)),
        axis.text.y = element_text(angle = 0, hjust = 1, size = 12), plot.background = element_rect(colour = "black", fill=NA, linewidth=3)) +
  geom_axis_flags(breaks = B0801$iso3,
                  labels = B0801$Population,
                  country_icons = B0801$iso3,
                  width = 30,
                  axis = "y",
                  lineheight = 2,
                  fontface = "bold"
  ) +
  ylab("Population") +
  xlab("Percentage (%)")