##################################################
## Project: MIRIADE
## Script purpose: Clean the data table for Mass Spec entries
## Date: 07.01.2022, 28.09.2022
## Author: Marek Ostaszewski, Felicia Burtscher
##################################################

options(stringsAsFactors = F)

library(dplyr)
library(readxl)
library(here)

source(here("functions", "data_processing_functions.R"))

datasets_root_directory <- define_datasets_root()

### Read the results Excel file for FTD
# ftd_mspec <- readxl::read_xlsx("_notgit/VUMC/WBI_CSF_Candidates_FC1.5_P0.05_DP0.5_FTD_DLB.xlsx", sheet = 1) %>%
ftd_mspec <- readxl::read_xlsx(file.path(datasets_root_directory, "VUMC/WBI_CSF_Candidates_FC1.5_P0.05_DP0.5_FTD_DLB.xlsx"), sheet = 1) %>%
  dplyr::rename(UniProt = Accession, Name = Gene, p.val = Pval) %>%
  dplyr::select(UniProt, Name, p.val)

### Read the results Excel file for DLB
# dlb_mspec <- readxl::read_xlsx("_notgit/VUMC/WBI_CSF_Candidates_FC1.5_P0.05_DP0.5_FTD_DLB.xlsx", sheet = 2) %>%
dlb_mspec <- readxl::read_xlsx(file.path(datasets_root_directory, "VUMC/WBI_CSF_Candidates_FC1.5_P0.05_DP0.5_FTD_DLB.xlsx"), sheet = 2) %>%
  dplyr::rename(UniProt = Accession, Name = Gene, p.val = Pval) %>%
  dplyr::select(UniProt, Name, p.val)

### Combine into a single data table
vumc_mspec <- rbind(data.frame(ftd_mspec, disease = "FTD"),
                    data.frame(dlb_mspec, disease = "DLB"))

# write.table(vumc_mspec, file = "_notgit/cleaned_VUMC_mspec.tsv", sep = "\t",
write.table(vumc_mspec, file = file.path(datasets_root_directory, "VUMC/cleaned_VUMC_mspec.tsv"), sep = "\t",
            row.names = F, col.names = T, quote = F)

### separate file according to disease type -- to feed into MetaCore etc.
vumc_mspec_FTD <- dplyr::filter(vumc_mspec, disease %in% c("FTD"))
vumc_mspec_DLB <- dplyr::filter(vumc_mspec, disease %in% c("DLB"))

write.table(vumc_mspec_FTD, file = file.path(datasets_root_directory, "VUMC/cleaned_VUMC_mspec_FTD.tsv"), sep = "\t",
            row.names = F, col.names = T, quote = F)
write.table(vumc_mspec_DLB, file = file.path(datasets_root_directory, "VUMC/cleaned_VUMC_mspec_DLB.tsv"), sep = "\t",
            row.names = F, col.names = T, quote = F)

