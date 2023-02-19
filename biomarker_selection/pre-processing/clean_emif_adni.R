##################################################
## Project: MIRIADE
## Script purpose: Clean the data table for ADNI entries (Tijms et al, Brain, 2020)
## Date: 10.01.2022, 15.09.2022
## Author: Marek Ostaszewski, Felicia Burtscher
##################################################

options(stringsAsFactors = F)

library(readxl)
library(dplyr)
library(here)

source(here("functions", "data_processing_functions.R"))

datasets_root_directory <- define_datasets_root()
### Read the Supplementary Materials Excel sheet, Supplementary table 2
### Contains the combined list of preselected ADNI and EMIF entries
### We need to skip the first two lines and select columns 1:6 (metadata) and 12:13 (AD vs Controls)
### Focus only on those selected for clustering

supp <- readxl::read_xlsx(file.path(datasets_root_directory, "Brain 2020 AD/awaa325-suppl_data/brain-2020-00570-File008.xlsx"),
                          skip = 2,
                          sheet = "Supplementary table 2")[, c(1:6, 12:13)] %>%
  dplyr::filter(`Selected for clustering` == "Yes") %>%
  dplyr::rename(Name = Gene, UniProt = Uniprot, p.val = `p value...13`) %>%
  dplyr::select(Name, UniProt, p.val, Cohort)

### Convert the p.val column to numeric, "<.001 values become 0.001"
final <- dplyr::mutate(supp, p.val = case_when(p.val == "<.001" ~ "0.001",
                                                   TRUE ~ p.val)) %>%
  dplyr::mutate(p.val = as.numeric(p.val))

### Handle multientries
to_split <- grep("; ", final$UniProt)
for(t in to_split) {
  mentry <- final[t,]
  for(s in strsplit(mentry$UniProt, split = "; ")[[1]]) {
    mentry$UniProt <- s
    final <- rbind(final, mentry)
  }
}
final <- final[-to_split,]

### Adjust the p values, separately for both cohorts
final <- cbind(final, adj.p = final$p.val, disease = "ALZ")
final$adj.p[final$Cohort == "ADNI"] <- p.adjust(final$p.val[final$Cohort == "ADNI"])
final$adj.p[final$Cohort == "EMIF-AD MBD"] <- p.adjust(final$p.val[final$Cohort == "EMIF-AD MBD"])

### Write down separately EMIF and ADNI data
# write.table(dplyr::filter(final, Cohort == "ADNI"), file = "_notgit/cleaned_ADNI_entries.tsv",
#             sep = "\t", quote = F, row.names = F, col.names = T)
# write.table(dplyr::filter(final, Cohort == "EMIF-AD MBD"), file = "_notgit/cleaned_EMIF_entries.tsv",
#             sep = "\t", quote = F, row.names = F, col.names = T)

write.table(dplyr::filter(final, Cohort == "ADNI"), file = file.path(datasets_root_directory, "Brain 2020 AD/cleaned_ADNI_entries.tsv"),
            sep = "\t", quote = F, row.names = F, col.names = T)
write.table(dplyr::filter(final, Cohort == "EMIF-AD MBD"), file = file.path(datasets_root_directory, "Brain 2020 AD/cleaned_EMIF_entries.tsv"),
            sep = "\t", quote = F, row.names = F, col.names = T)

