##################################################
## Project: MIRIADE
## Script purpose: Clean the data table for Olink entries
## Date: 07.01.2022, 15.09.2022
## Author: Marek Ostaszewski, Felicia Burtscher
##################################################

options(stringsAsFactors = F)

library(dplyr)
library(readxl)
library(here)

source(here("functions", "data_processing_functions.R"))

datasets_root_directory <- define_datasets_root()

### Read the olink mapping file
# olp <- readxl::read_xlsx("_notgit/VUMC/Olink/Olink_proteins_list.xlsx") %>%
olp <- readxl::read_xlsx(file.path(datasets_root_directory, "Olink/20201216_OLINK_Data/Olink_proteins_list_Marek.xlsx")) %>%
  dplyr::select(Assay, UniProt) %>% 
  dplyr::distinct() 

replace_add <- function(table, multiup, up_reps) {
  row <- table[table$UniProt == multiup,]
  for(r in 1:length(up_reps)) {
    row$UniProt <- up_reps[r]
    table <- rbind(table, row)
  }
  table[table$UniProt != multiup,]
}

### Manual cleanup of some incorrect UniProt entries
olp$UniProt[olp$UniProt == "O43521-2"] <- "O43521"
olp$UniProt[olp$UniProt == "P49763-3"] <- "P49763"
olp <- replace_add(olp, "P29460/P29459", up_reps = c("P29460", "P29459"))
olp <- replace_add(olp, "Q8NEV9/Q14213", up_reps = c("Q8NEV9", "Q14213"))
olp <- replace_add(olp, "Q29983-Q29980", up_reps = c("Q29983", "Q29980"))
olp <- replace_add(olp, "Q11128/P21217", up_reps = c("Q11128", "P21217"))

### Read the datasets (sorted by increasing adj. p-values), add a position to each
### the Alzheimers dataset
# alzp <- read.csv("_notgit/VUMC/Olink/FVOL_AD_CON.csv", sep = ",",
#                  col.names = c("Name", "Effect", "p.val", "adj.p")) %>%
#   dplyr::mutate(disease = "ALZ")
alzp <- read.csv(file.path(datasets_root_directory, "Olink/20201216_OLINK_Data/FVOL_AD_CON.csv"), sep = ",",
                 col.names = c("Name", "Effect", "p.val", "adj.p")) %>%
  dplyr::mutate(disease = "ALZ")


### the Frontotemporal Dementia dataset
# ftdp <- read.csv("_notgit/VUMC/Olink/FVOL_FTD_CON.csv", sep = ",",
#                  col.names = c("Name", "Effect", "p.val", "adj.p")) %>%
#   dplyr::mutate(disease = "FTD")
ftdp <- read.csv(file.path(datasets_root_directory, "Olink/20201216_OLINK_Data/FVOL_FTD_CON.csv"), sep = ",",
                 col.names = c("Name", "Effect", "p.val", "adj.p")) %>%
  dplyr::mutate(disease = "FTD")

### the Dementia with Lewy Bodies dataset
# dlbp <- read.csv("_notgit/VUMC/Olink/FVOL_DLB_CON.csv", sep = ",",
#                  col.names = c("Name", "Effect", "p.val", "adj.p")) %>%
#   dplyr::mutate(disease = "DLB")
dlbp <- read.csv(file.path(datasets_root_directory, "Olink/20201216_OLINK_Data/FVOL_DLB_CON.csv"), sep = ",",
                 col.names = c("Name", "Effect", "p.val", "adj.p")) %>%
  dplyr::mutate(disease = "DLB")

all <- rbind(rbind(alzp, ftdp), dlbp)

### Read cutoffs file, get UniProts that have less than 20% LOD percentage
# cutoffs <- readxl::read_xlsx("_notgit/VUMC/Olink/Perc_below_LOD_per_assay.xlsx")
cutoffs <- readxl::read_xlsx(file.path(datasets_root_directory, "Olink/20201216_OLINK_Data/Perc_below_LOD_per_assay.xlsx"))
cutoff_ups <- unique(cutoffs[cutoffs$percentage_above_LOD < 20, "UniProt"]) %>% pull()

### Add UniProt ids based on assay and filter our NA UniProts,
### exclude low LOD UniProt ids from the dataset
all_up <- merge(x = all, y = unique(olp[,c("Assay", "UniProt")]),
                by.x = "Name", by.y = "Assay") %>%
  dplyr::filter(!is.na(UniProt) & UniProt != "NA") %>%
  dplyr::filter(!(UniProt %in% cutoff_ups)) %>%
  dplyr::group_by(UniProt, disease) %>%
  dplyr::summarise(p.val = min(p.val))

write.table(all_up, file = file.path(datasets_root_directory, "Olink/20201216_OLINK_Data/cleaned_VUMC_olink.tsv"), sep = "\t",
            row.names = F, col.names = T, quote = F)

### separate file according to disease type -- to feed into MetaCore etc.
all_up_ALZ <- dplyr::filter(all_up, disease %in% c("ALZ"))
all_up_DLB <- dplyr::filter(all_up, disease %in% c("DLB"))
all_up_FTD <- dplyr::filter(all_up, disease %in% c("FTD"))

write.table(all_up_ALZ, file = file.path(datasets_root_directory, "Olink/20201216_OLINK_Data/cleaned_VUMC_olink_AD.tsv"), sep = "\t",
            row.names = F, col.names = T, quote = F)
write.table(all_up_DLB, file = file.path(datasets_root_directory, "Olink/20201216_OLINK_Data/cleaned_VUMC_olink_DLB.tsv"), sep = "\t",
            row.names = F, col.names = T, quote = F)
write.table(all_up_FTD, file = file.path(datasets_root_directory, "Olink/20201216_OLINK_Data/cleaned_VUMC_olink_FTD.tsv"), sep = "\t",
            row.names = F, col.names = T, quote = F)

