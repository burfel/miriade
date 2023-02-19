##################################################
## Project: MIRIADE
## Script purpose: Integrate different biomarker sources in MIRIADE and annotate them with stable IDs
## Date: 07.01.2022, 22.09.2022
## Author: Marek Ostaszewski, Felicia Burtscher
##################################################

options(stringsAsFactors = F)
library(here)

source(here("functions", "data_processing_functions.R"))

datasets_root_directory <- define_datasets_root()

### Read the VUMC Olink biomarkers
# vumc_ol <- read.table("_notgit/cleaned_VUMC_olink.tsv", sep = "\t", header = T,
#                    quote = "") ### 'quote' need due to a special sign in 'Name'
vumc_ol <- read.table(file.path(datasets_root_directory, "Olink/20201216_OLINK_Data/cleaned_VUMC_olink.tsv"), sep = "\t", header = T,
                      quote = "") ### 'quote' need due to a special sign in 'Name'

### Read the VUMC Mass Spec biomarkers
# vumc_ms <- read.table("_notgit/cleaned_VUMC_mspec.tsv", sep = "\t", header = T)
vumc_ms <- read.table(file.path(datasets_root_directory, "VUMC/cleaned_VUMC_mspec.tsv"), sep = "\t", header = T)

### Read the KTH biomarkers
# kth <- read.table("_notgit/cleaned_kth_entries.tsv", sep = "\t", header = T)
kth <- read.table(file.path(datasets_root_directory, "KTH/cleaned_KTH_entries.tsv"), sep = "\t", header = T)

### Read the ADNI biomarkers (Tijms et al, 2020, Brain)
# adni <- read.table("_notgit/cleaned_adni_entries.tsv", sep = "\t", header = T)
adni <- read.table(file.path(datasets_root_directory, "Brain 2020 AD/cleaned_ADNI_entries.tsv"), sep = "\t", header = T)

### Read the VUMC biomarkers (Tijms et al, 2020, Brain)
# emif <- read.table("_notgit/cleaned_emif_entries.tsv", sep = "\t", header = T)
emif <- read.table(file.path(datasets_root_directory, "Brain 2020 AD/cleaned_EMIF_entries.tsv"), sep = "\t", header = T)

### First, combine p values across all datasets, using UniProt as the index column
common_cols <- c("UniProt", "p.val", "disease")

### Pool all the entries together, add source information
combined <- do.call(rbind, list(cbind(vumc_ol[,common_cols], source = "VUMC Olink"),
                                cbind(vumc_ms[,common_cols], source = "VUMC MSpec"),
                                cbind(kth[,common_cols], source = "KTH"),
                                cbind(adni[,common_cols], source = "ADNI"),
                                cbind(emif[,common_cols], source = "EMIF")))

### Combine p values, thanks to Armin Rauschenberger
### for explaining and implementing the Fisher and Simes method

fisher_method <- function(pval) {
  1-stats::pchisq(q  = sum(-2*log(pval[!is.na(pval)])),
                  df = 2*sum(!is.na(pval)))
}

# # package "mppa" no longer available (25/09/2022)
# library(dplyr)
# install.packages("mppa")
# library("mppa")

# install.packages("/Users/felicia.burtscher/Downloads/mppa",
#                  repos = NULL, type = "source")

library(installr)
# install.packages("https://cran.r-project.org/src/contrib/Archive/mppa/mppa_1.0.tar.gz",verbose=TRUE)
library(remotes)
install_version("mppa","1.0")
library("mppa")

# Simes <- function(x){
#   p.vals <- apply(x, 2, function(z) t.test(z)$p.value) # Compute variable-wise pvalues
#   p <- ncol(x)
#   p.Simes <- p * min(sort(p.vals)/seq_along(p.vals)) # Compute the Simes statistic
#   return(c(pvalue=p.Simes))
# }

### Calculate for both Fisher and Simes, create adjusted p values
combined <- dplyr::group_by(combined, disease, UniProt) %>%
  dplyr::summarise(p.val_fisher = fisher_method(p.val),
                   p.val_sims = mppa::simes.test(p.val),
                   sources = paste(sort(source), collapse = ",")) %>%
  dplyr::ungroup() %>% dplyr::group_by(disease) %>%
  dplyr::mutate(adj.p.val_fisher = p.adjust(p.val_fisher, method = "BY"), 
                adj.p.val_sims = p.adjust(p.val_sims, method = "BY"))

### Tables of disease-specific entries 
ALZ_df <- dplyr::select(combined, UniProt, adj.p.val_fisher, disease, sources) %>%
  dplyr::filter(disease == "ALZ") %>%
  dplyr::rename(adj.p = adj.p.val_fisher) %>%
  dplyr::filter(adj.p <= 0.1) %>% distinct()
  
FTD_df <- dplyr::select(combined, UniProt, adj.p.val_sims, disease, sources) %>%
  dplyr::filter(disease == "FTD") %>%
  dplyr::rename(adj.p = adj.p.val_sims) %>%
  dplyr::filter(adj.p <= 0.1) %>% distinct()

DLB_df <- dplyr::select(combined, UniProt, adj.p.val_sims, disease, sources) %>%
  dplyr::filter(disease == "DLB") %>%
  dplyr::rename(adj.p = adj.p.val_sims) %>%
  dplyr::filter(adj.p <= 0.1) %>% distinct()

### Create a table of proteins specific for each of the tables
specific_BMs <- rbind(
  dplyr::filter(ALZ_df, UniProt %in% setdiff(ALZ_df$UniProt, c(FTD_df$UniProt, DLB_df$UniProt))),
  dplyr::filter(FTD_df, UniProt %in% setdiff(FTD_df$UniProt, c(ALZ_df$UniProt, DLB_df$UniProt))),
  dplyr::filter(DLB_df, UniProt %in% setdiff(DLB_df$UniProt, c(FTD_df$UniProt, ALZ_df$UniProt)))
  )

### Create a table of proteins common between the tables
common_BMs <- rbind(
  dplyr::filter(ALZ_df, UniProt %in% union(FTD_df$UniProt, DLB_df$UniProt)),
  dplyr::filter(FTD_df, UniProt %in% union(ALZ_df$UniProt, DLB_df$UniProt)),
  dplyr::filter(DLB_df, UniProt %in% union(FTD_df$UniProt, ALZ_df$UniProt))
) 
### Cleanup the common table by merging same UniProt entries
common_BMs <- dplyr::group_by(common_BMs, UniProt) %>%
  dplyr::summarise(adj.p = min(adj.p),
                   sources = paste(paste0(disease,":",sources), collapse = ";"),
                   disease = paste(disease, collapse = ","))

### Create a common table of candidate proteins to be used downstream
candidate_BMs <- rbind(specific_BMs, common_BMs)

### Gene names and Entrez ids, for downstream data analysis/enrichment
### Retrieve HGNCs and NCBI gene ids from AnnotationDbi, using UniProt ids
library(AnnotationDbi)
# BiocManager::install("org.Hs.eg.db")
BiocManager::install("org.Hs.eg.db", force = TRUE)
id_mapping <- AnnotationDbi::select(org.Hs.eg.db, candidate_BMs$UniProt,
                                    columns = c("SYMBOL", "ENTREZID"), keytype = "UNIPROT")
### Manual tweaks to the mapping, there are some inconsistencies
id_mapping[id_mapping$UNIPROT == "Q8WWJ7", 2:3] <- c("EPHB6", "2051")
id_mapping[id_mapping$UNIPROT == "Q8N423", 2:3] <- c("LILRB2", "10288")
id_mapping[id_mapping$UNIPROT == "Q8WWJ7", 2:3] <- c("ACAN", "176")
id_mapping <- dplyr::filter(id_mapping, !(SYMBOL %in% c("TNFSF12-TNFSF13", "C4B_2", "CALM2", "CALM3")))

### Create the annotated set of biomarkers
annotated_BMs <- merge(candidate_BMs, id_mapping, by.x = "UniProt", by.y = "UNIPROT")
# write.table(annotated_BMs, file = "_notgit/annotated_candidate_biomarkers.tsv",
#             sep = "\t", quote = F, row.names = F, col.names = T)
write.table(annotated_BMs, file = file.path(datasets_root_directory, "annotated_candidate_biomarkers.tsv"),
            sep = "\t", quote = F, row.names = F, col.names = T)
