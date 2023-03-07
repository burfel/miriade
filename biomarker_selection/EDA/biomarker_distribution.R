#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
## Project: MIRIADE
## Script purpose: Plot out biomarker distributions according to different
##                 sortings of the patients
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
library(here)
library(dplyr)

source(here("biomarker_selection", "EDA", "meld_all_fractured_datasets.R"))
source(here("functions", "mapping_functions.R"))

### Read the KTH table
kth <- readxl::read_xlsx(file.path(datasets_root_directory,
                                   "KTH/KTH AD dataset for MIRIADE biomarker selection sample info and results.xlsx"))

### Change the names into HGNC symbols (prefixes)
colnames(kth) <- sapply(colnames(kth), take_only_pre_underscore_substring)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
####                            Extract proteins                            ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

olink_proteins <- data.frame(Name = colnames(olink)[5:ncol(olink)]) %>%
  map_hgnc_to_uniprot_and_cbind("UniProt") %>%
  na.omit()

kth_proteins <- data.frame(Name = colnames(kth)[5:ncol(kth)]) %>%
  map_hgnc_to_uniprot_and_cbind("UniProt") %>%
  na.omit()

emif_proteins <- data.frame(UniProt = colnames(emif)[6:ncol(emif)]) %>%
  map_uniprot_to_hgnc_and_cbind("Name") %>%
  na.omit() %>%
  dplyr::select(UniProt, Name)

proteins_in_all_three_datasets <-
  dplyr::inner_join(olink_proteins, kth_proteins,
                    by = "UniProt", suffix = c("", "2")) %>%
  dplyr::select(-c("Name2")) %>%
  dplyr::inner_join(emif_proteins,
                    by = "UniProt", suffix = c("", "2")) %>%
  dplyr::select(-c("Name2"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
####     Filter out non overlapping proteins + additional processing        ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# use sapply to create a named vector of the hgnc and their new uniprots
protein_renaming_vector <- sapply(proteins_in_all_three_datasets$UniProt, function(x) {
  proteins_in_all_three_datasets$Name[which(proteins_in_all_three_datasets$UniProt == x)]
})

# Create a mapping rule as a named vector for the diagnosis to have a unified format
diagnosis_mapping <- c("CN" = "Control",
                       "NL" = "Control",
                       "AD dementia" = "AD")

filtered_olink <- olink %>%
  dplyr::select(SampleId, dx, age_gr, sex, all_of(proteins_in_all_three_datasets$Name)) %>%
  dplyr::arrange(age_gr) %>%
  dplyr::rename(Age = "age_gr", Diagnosis = "dx", Gender = "sex", Subject = "SampleId") %>%
  dplyr::rename(!!!protein_renaming_vector) %>%
  dplyr::mutate(Diagnosis = recode(Diagnosis, !!!diagnosis_mapping)) %>%
  dplyr::filter(Diagnosis == "Control" | Diagnosis == "AD") %>%
  dplyr::mutate(Id = row_number())

filtered_kth <- kth %>%
  dplyr::select(Diagnosis, Age, Gender, all_of(proteins_in_all_three_datasets$Name)) %>%
  dplyr::arrange(Age) %>%
  dplyr::rename(!!!protein_renaming_vector) %>%
  dplyr::filter(Diagnosis == "Control" | Diagnosis == "AD") %>%
  dplyr::mutate(Id = row_number())

filtered_emif <- emif %>%
  dplyr::select(SubjectId, Diagnosis, Age, Gender, all_of(proteins_in_all_three_datasets$UniProt)) %>%
  dplyr::arrange(Age) %>%
  dplyr::rename(Subject = "SubjectId") %>%
  dplyr::mutate(Diagnosis = recode(Diagnosis, !!!diagnosis_mapping)) %>%
  dplyr::filter(Diagnosis == "Control" | Diagnosis == "AD") %>%
  dplyr::mutate(Id = row_number())

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
####               Bind the three datasets and then pivot long              ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# We need to make sure the "Subject" column is of the same type
filtered_olink$Subject <- as.character(filtered_olink$Subject)

df_all <- bind_rows(filtered_olink %>% mutate(cohort = "olink"),
                    filtered_kth %>% mutate(cohort = "kth"),
                    filtered_emif %>% mutate(cohort = "emif"))

long_all <- df_all %>%
  tidyr::pivot_longer(cols = 5:12, names_to = "UniProt")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
####                             Normalize                                  ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

means <- long_all %>%
  group_by(cohort, UniProt) %>%
  summarize(mean_value = mean(value, na.rm = TRUE)) %>%
  ungroup()

long_all <- long_all %>%
  left_join(means, by = c("cohort", "UniProt")) %>%
  mutate(normalized_value = value / mean_value) %>%
  dplyr::select(-mean_value)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
####                             Plotting                                   ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

plot_biomarker_distribution <- function(long_all, uniprot)
{
  uniprot_data <- long_all %>%
    dplyr::filter(UniProt == uniprot)
  
  ggplot(uniprot_data, aes(x = Id, y = normalized_value, color = Diagnosis)) +
    geom_line() +
    labs(title = paste("Distribution for", uniprot),
         subtitle = "where the subjects were sorted by age and the Id is the row number") +
    facet_wrap(~ cohort, ncol = 3, scales = "free_x")
}

# Per unique UniProt create a plot
biomarker_distribution_plots <-
  lapply(unique(long_all$UniProt),
         function(uniprot) plot_biomarker_distribution(long_all, uniprot))
