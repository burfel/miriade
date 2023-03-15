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

emif_proteins <- data.frame(Name = colnames(emif)[6:ncol(emif)]) %>%
  map_hgnc_to_uniprot_and_cbind("UniProt") %>%
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

protein_renaming_to_hgnc_vector <- sapply(proteins_in_all_three_datasets$Name, function(x) {
  proteins_in_all_three_datasets$UniProt[which(proteins_in_all_three_datasets$Name == x)]
})

# Create a mapping rule as a named vector for the diagnosis to have a unified format
diagnosis_mapping <- c("CN" = "Control",
                       "NL" = "Control",
                       "AD dementia" = "AD")

spread_ids <- function(df)
{
  # Calculate number of rows for each group
  n_ad <- sum(df$Diagnosis == "AD")
  n_control <- sum(df$Diagnosis == "Control")
  
  # Calculate maximum range of integers to use
  max_range <- max(n_ad, n_control)
  
  # Add new column with sequence of integers
  df <- df %>%
    group_by(Diagnosis) %>%
    mutate(Id = ifelse(Diagnosis == "AD", 
                             seq(1, max_range, length.out = n_ad),
                             seq(1, max_range, length.out = n_control))) %>%
    ungroup()
  return(df)
}

filtered_olink <- olink %>%
  dplyr::select(SampleId, Diagnosis, Age, Gender, all_of(proteins_in_all_three_datasets$Name)) %>%
  dplyr::arrange(Age) %>%
  dplyr::filter(Diagnosis == "Control" | Diagnosis == "AD") %>%
  arrange(Diagnosis != "AD") %>%
  spread_ids()

filtered_kth <- kth %>%
  dplyr::select(Diagnosis, Age, Gender, all_of(proteins_in_all_three_datasets$Name)) %>%
  dplyr::arrange(Age) %>%
  dplyr::filter(Diagnosis == "Control" | Diagnosis == "AD") %>%
  arrange(Diagnosis != "AD") %>%
  spread_ids()

filtered_emif <- emif %>%
  dplyr::select(SubjectId, Diagnosis, Age, Gender, all_of(proteins_in_all_three_datasets$Name)) %>%
  dplyr::arrange(Age) %>%
  dplyr::rename(Subject = "SubjectId") %>%
  #dplyr::rename(!!!protein_renaming_to_hgnc_vector) %>%
  dplyr::filter(Diagnosis == "Control" | Diagnosis == "AD") %>%
  arrange(Diagnosis != "AD") %>%
  spread_ids()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
####               Bind the three datasets and then pivot long              ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# We need to make sure the "Subject" column is of the same type
filtered_olink$Subject <- as.character(filtered_olink$Subject)

df_all <- bind_rows(filtered_olink %>% mutate(cohort = "olink"),
                    filtered_kth %>% mutate(cohort = "kth"),
                    filtered_emif %>% mutate(cohort = "emif"))

long_all <- df_all %>%
  #tidyr::pivot_longer(cols = 5:12, names_to = "UniProt")
  tidyr::pivot_longer(cols = 5:12, names_to = "HGNC")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
####                             Normalize                                  ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

means <- long_all %>%
  #group_by(cohort, UniProt) %>%
  group_by(cohort, HGNC) %>%
  summarize(mean_value = mean(value, na.rm = TRUE)) %>%
  ungroup()

long_all <- long_all %>%
  #left_join(means, by = c("cohort", "UniProt")) %>%
  left_join(means, by = c("cohort", "HGNC")) %>%
  mutate(normalized_value = value / mean_value) %>%
  dplyr::select(-mean_value)

source(here("functions", "auxiliary_statistical_functions.R"))
breaks_and_labels <- prepare_breaks_and_labels_by_quantiles(long_all$Age, 3)

long_all <- long_all %>%
  convert_age_column_to_age_group_column(sym("Age"), sym("Age_Group"),
                                         cutting_breaks = breaks_and_labels$Breaks,
                                         cutting_labels = breaks_and_labels$Labels) %>%
  dplyr::select(Diagnosis, Age_Group, everything())

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
####                             Plotting                                   ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

plot_biomarker_distribution <- function(long_all, protein)
{
  require(ggplot2)
  uniprot_data <- long_all %>%
    #dplyr::filter(UniProt == uniprot)
    dplyr::filter(HGNC == protein)
  plots <- vector("list")
  
  plots[["Scatter by cohort"]] <-
    ggplot(uniprot_data, aes(x = Id, y = normalized_value, color = Diagnosis)) +
    geom_point() +
    labs(title = paste("Distribution for", protein),
         subtitle = "where the subjects were sorted by age and the Id was built so the diagnoses are approx. similarly spread") +
    facet_wrap(~ cohort , ncol = 3, scales = "free_x")
  
  plots[["Scatter by cohort and age group"]] <-
    ggplot(uniprot_data, aes(x = Id, y = normalized_value, color = Diagnosis)) +
    geom_point() +
    labs(title = paste("Distribution for", protein),
         subtitle = "where the subjects were sorted by age and the Id was built so the diagnoses are approx. similarly spread") +
    facet_wrap(~ cohort + Age_Group, ncol = 3, scales = "free_x")
  
  plots[["Boxplot by cohort"]] <-
    ggplot(uniprot_data, aes(x = Diagnosis, y = normalized_value, color = Diagnosis)) +
    geom_boxplot() +
    labs(title = paste("Distribution for", protein),
         subtitle = "where the subjects were sorted by age") +
    facet_wrap(~ cohort, ncol = 3, scales = "free_x")
  
  plots[["Violin by cohort"]] <-
    ggplot(uniprot_data, aes(x = Diagnosis, y = normalized_value, color = Diagnosis)) +
    geom_violin() +
    labs(title = paste("Distribution for", protein),
         subtitle = "where the subjects were sorted by age") +
    facet_wrap(~ cohort, ncol = 3, scales = "free_x")
  
  return(plots)
}

# Per unique UniProt create a plot
biomarker_distribution_plots <-
  #lapply(unique(long_all$UniProt),
  #       function(uniprot) plot_biomarker_distribution(long_all, uniprot))
  lapply(unique(long_all$HGNC),
         function(hgnc) plot_biomarker_distribution(long_all, hgnc))
names(biomarker_distribution_plots) <- unique(long_all$HGNC)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
## Try for specfifc proteins EXPERIMENTAL ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#\

std_olink <- olink %>%
  dplyr::filter(Diagnosis == "Control" | Diagnosis == "AD") %>%
  arrange(Diagnosis != "AD") %>%
  spread_ids() %>%
  dplyr::select(Id, everything())

std_kth <- kth %>%
  dplyr::select(Diagnosis, Age, Gender, everything()) %>%
  dplyr::arrange(Age) %>%
  dplyr::filter(Diagnosis == "Control" | Diagnosis == "AD") %>%
  arrange(Diagnosis != "AD") %>%
  spread_ids() %>%
  dplyr::select(Id, everything())

std_emif <- emif %>%
  dplyr::select(SubjectId, Diagnosis, Age, Gender, everything()) %>%
  dplyr::arrange(Age) %>%
  dplyr::rename(Subject = "SubjectId") %>%
  dplyr::filter(Diagnosis == "Control" | Diagnosis == "AD") %>%
  arrange(Diagnosis != "AD") %>%
  spread_ids() %>%
  dplyr::select(Id, everything())

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
####               Bind the three datasets and then pivot long              ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# We need to make sure the "Subject" column is of the same type
std_olink$Subject <- as.character(std_olink$Subject)

std_all <- bind_rows(std_olink %>% mutate(cohort = "olink"),
                    std_kth %>% mutate(cohort = "kth"),
                    std_emif %>% mutate(cohort = "emif")) %>%
  dplyr::select(Id, cohort, everything(), -Subject, -SampleId)

long_std <- std_all %>%
  tidyr::pivot_longer(cols = 7:ncol(std_all), names_to = "HGNC") %>%
  na.omit()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
####                             Normalize                                  ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

means <- long_std %>%
  group_by(cohort, HGNC) %>%
  summarize(mean_value = mean(value, na.rm = TRUE)) %>%
  ungroup()

long_std <- long_std %>%
  left_join(means, by = c("cohort", "HGNC")) %>%
  mutate(normalized_value = value / mean_value) %>%
  dplyr::select(-mean_value)

proteins_to_check <- c("DDAH1", "NPTX2", "SPON1", "PEBP1", "ENO2", "ENO1", "NPTXR",
                       "TREM1", "YWHAE", "YWHAZ", "GAP43", "IL1B", "EDN1", "VGF",
                       "NEFM", "GAD1", "THOP1", "ABL1", "MANBA",
                       "APP", "SOD1", "MMP_1", "MMP_10")

specific_biomarkers_distribution <-
  lapply(proteins_to_check, function(hgnc) plot_biomarker_distribution(long_std, hgnc))
names(specific_biomarkers_distribution) <- proteins_to_check
# Recommend to remove when done with due to size
remove(specific_biomarkers_distribution)
