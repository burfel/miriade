library(here)
library(dplyr)

# This loads the melded emif and olink dataframes to be used later in the file
# It also calls to set the dataset root directory
source(here("biomarker_selection", "EDA", "meld_all_fractured_datasets.R"))
source(here("functions", "auxiliary_statistical_functions.R"))
source(here("functions", "boxplot_functions.R"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
## Processing the datasets ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
### KTH ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
### Read the table
kth <- readxl::read_xlsx(file.path(datasets_root_directory,
                                   "KTH/KTH AD dataset for MIRIADE biomarker selection sample info and results.xlsx"))

### Change the names into HGNC symbols (prefixes)
colnames(kth) <- sapply(colnames(kth), take_only_pre_underscore_substring)

breaks_and_labels <- prepare_breaks_and_labels_by_quantiles(kth$Age, 3)
mutated_kth <- kth %>%
  convert_age_column_to_age_group_column(sym("Age"), sym("Age_Group"),
                                         cutting_breaks = breaks_and_labels$Breaks,
                                         cutting_labels = breaks_and_labels$Labels) %>%
  dplyr::filter(Diagnosis == "AD" | Diagnosis == "Control") %>%
  dplyr::mutate(Diagnosis_Age_Group = paste(Diagnosis, Age_Group)) %>%
  dplyr::select(Diagnosis_Age_Group, Age_Group, everything())

print("The number of patients per age group (in KTH) is:")
mutated_kth %>% dplyr::group_by(Age_Group) %>% dplyr::summarise(n())


### Melt the dataframe, so there's only one readout variable
melted_kth <-
  reshape2::melt(kth, id = 1:9,
                 variable.name = "HGNC_Symbol", value.name = "Value") %>%
  convert_age_column_to_age_group_column(sym("Age"), sym("Age_Group"),
                                         cutting_breaks = breaks_and_labels$Breaks,
                                         cutting_labels = breaks_and_labels$Labels) %>%
  dplyr::select(Diagnosis, Age_Group, Gender, HGNC_Symbol, Value)

#### Adjusting to have Diagnosis + age group combined column                ####

supercharged_kth_melt <-
  melted_kth %>%
  dplyr::filter(Diagnosis == "AD" | Diagnosis == "Control") %>%
  dplyr::mutate(Diagnosis_Age_Group = paste(Diagnosis, Age_Group)) %>%
  dplyr::select(Diagnosis_Age_Group, Gender, HGNC_Symbol, Value)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
### Emif ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

source(here("functions", "mapping_functions.R"))
emif_proteins <- data.frame(UniProt = colnames(emif)[6:ncol(emif)]) %>%
  map_uniprot_to_hgnc_and_cbind("Name") %>%
  na.omit() %>%
  dplyr::select(UniProt, Name)
# use sapply to create a named vector of the hgnc and their new uniprots
protein_renaming_vector <- sapply(emif_proteins$Name, function(x) {
  emif_proteins$UniProt[which(emif_proteins$Name == x)]
})

emif$Age <- round(emif$Age) # Must round, otherwise the age groups are classified incorrectly

breaks_and_labels <- prepare_breaks_and_labels_by_quantiles(emif$Age, 3)

mutated_emif <- emif %>%
  convert_age_column_to_age_group_column(sym("Age"), sym("Age_Group"),
                                         cutting_breaks = breaks_and_labels$Breaks,
                                         cutting_labels = breaks_and_labels$Labels) %>%
  dplyr::filter(Diagnosis == "AD" | Diagnosis == "NL") %>%
  dplyr::mutate(Diagnosis = case_when(Diagnosis == "NL" ~ "Control", TRUE ~ Diagnosis)) %>%
  dplyr::mutate(Diagnosis_Age_Group = paste(Diagnosis, Age_Group)) %>%
  dplyr::select(Diagnosis_Age_Group, Age_Group, everything()) %>%
  dplyr::rename(!!!protein_renaming_vector)

print("The number of patients per age group (in Emif) is:")
mutated_emif %>% dplyr::group_by(Age_Group) %>% dplyr::summarise(n())

### Melt the emif dataframe, so there's only one readout variable
melted_emif <- emif %>%
  dplyr::rename(!!!protein_renaming_vector) %>%
  reshape2::melt(id = 1:5,
                 variable.name = "HGNC_Symbol", value.name = "Value") %>%
  convert_age_column_to_age_group_column(sym("Age"), sym("Age_Group"),
                                         cutting_breaks = breaks_and_labels$Breaks,
                                         cutting_labels = breaks_and_labels$Labels) %>%
  dplyr::select(Diagnosis, Age_Group, Gender, HGNC_Symbol, Value)

supercharged_emif_melt <-
  melted_emif %>%
  dplyr::mutate(Diagnosis = case_when(Diagnosis == "NL" ~ "Control", TRUE ~ Diagnosis)) %>%
  dplyr::filter(Diagnosis == "AD" | Diagnosis == "Control") %>%
  dplyr::mutate(Diagnosis_Age_Group = paste(Diagnosis, Age_Group)) %>%
  dplyr::select(Diagnosis_Age_Group, Gender, HGNC_Symbol, Value)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
### Olink ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

breaks_and_labels <- prepare_breaks_and_labels_by_quantiles(olink$age_gr, 3)
mutated_olink <- olink %>%
  convert_age_column_to_age_group_column(sym("age_gr"), sym("Age_Group"),
                                         cutting_breaks = breaks_and_labels$Breaks,
                                         cutting_labels = breaks_and_labels$Labels) %>%
  dplyr::filter(dx == "AD dementia" | dx == "CN") %>%
  dplyr::mutate(Diagnosis = case_when(dx == "CN" ~ "Control", dx == "AD dementia" ~ "AD", TRUE ~ dx)) %>%
  dplyr::mutate(Diagnosis_Age_Group = paste(Diagnosis, Age_Group)) %>%
  dplyr::select(Diagnosis_Age_Group, Diagnosis, Age_Group, everything(), -dx)

print("The number of patients per age group (in Olink) is:")
mutated_olink %>% dplyr::group_by(Age_Group) %>% dplyr::summarise(n())

### Melt the olink dataframe, so there's only one readout variable
melted_olink <-
  reshape2::melt(olink, id = 1:4,
                 variable.name = "HGNC_Symbol", value.name = "Value") %>%
  convert_age_column_to_age_group_column(sym("age_gr"), sym("Age_Group"),
                                         cutting_breaks = breaks_and_labels$Breaks,
                                         cutting_labels = breaks_and_labels$Labels) %>%
  dplyr::select(dx, Age_Group, sex, HGNC_Symbol, Value)

supercharged_olink_melt <-
  melted_olink %>%
  dplyr::mutate(Diagnosis = case_when(dx == "CN" ~ "Control", dx == "AD dementia" ~ "AD", TRUE ~ dx)) %>%
  dplyr::filter(Diagnosis == "AD" | Diagnosis == "Control") %>%
  dplyr::mutate(Diagnosis_Age_Group = paste(Diagnosis, Age_Group)) %>%
  dplyr::select(Diagnosis_Age_Group, sex, HGNC_Symbol, Value)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
### Combine the datasets ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#### First, we bring all 3 datasets to a standard form ####
kth <- kth %>%
  dplyr::select(-Class) %>%
  dplyr::mutate(Gender = case_when(Gender == 'M' ~ 'm', Gender == 'F' ~ 'f', TRUE ~ Gender)) %>%
  dplyr::filter(Diagnosis == "AD" | Diagnosis == "Control")

# Create a mapping rule as a named vector for the diagnosis to have a unified format
diagnosis_mapping <- c("CN" = "Control",
                       "NL" = "Control",
                       "AD dementia" = "AD")

emif <- emif %>%
  dplyr::select(-SubjectId, -Assay.ID) %>%
  dplyr::mutate(Diagnosis = recode(Diagnosis, !!!diagnosis_mapping)) %>%
  dplyr::filter(Diagnosis == "Control" | Diagnosis == "AD") %>%
  dplyr::rename(!!!protein_renaming_vector)

olink <- olink %>%
  dplyr::select(-SampleId) %>%
  dplyr::rename(Age = "age_gr", Diagnosis = "dx", Gender = "sex") %>%
  dplyr::mutate(Diagnosis = recode(Diagnosis, !!!diagnosis_mapping)) %>%
  dplyr::filter(Diagnosis == "Control" | Diagnosis == "AD")

#### Now, we bind_rows ####
combined_df <- bind_rows(kth, emif, olink)

#### Now, mutate and melt

breaks_and_labels <- prepare_breaks_and_labels_by_quantiles(combined_df$Age, 3)
mutated_combined_df <- combined_df %>%
  convert_age_column_to_age_group_column(sym("Age"), sym("Age_Group"),
                                         cutting_breaks = breaks_and_labels$Breaks,
                                         cutting_labels = breaks_and_labels$Labels) %>%
  dplyr::mutate(Diagnosis_Age_Group = paste(Diagnosis, Age_Group)) %>%
  dplyr::select(Diagnosis_Age_Group, Age_Group, everything())

print("The number of patients per age group (in KTH) is:")
mutated_combined_df %>% dplyr::group_by(Age_Group) %>% dplyr::summarise(n())


### Melt the dataframe, so there's only one readout variable
melted_combined_df <-
  reshape2::melt(combined_df, id = 1:3,
                 variable.name = "HGNC_Symbol", value.name = "Value") %>%
  convert_age_column_to_age_group_column(sym("Age"), sym("Age_Group"),
                                         cutting_breaks = breaks_and_labels$Breaks,
                                         cutting_labels = breaks_and_labels$Labels) %>%
  dplyr::select(Diagnosis, Age_Group, Gender, HGNC_Symbol, Value) #%>%
  #na.omit()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
## KW and wilcoxon ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
### KTH ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
### Run Kruskal-Wallis test for all biomarkers grouped by gender
kth_gender_grouped_pvals <- perform_kw_wilcoxon_according_to_grouping(
  melted_df = melted_kth,
  differentiating_feature_symbol = sym("HGNC_Symbol"),
  grouping_feature_symbol = sym("Gender"),
  measurement_symbol = sym("Value"),
  significance_limit = 0.05)

### Age_Group
kth_age_group_grouped_pvals <- perform_kw_wilcoxon_according_to_grouping(
  melted_df = melted_kth,
  differentiating_feature_symbol = sym("HGNC_Symbol"),
  grouping_feature_symbol = sym("Age_Group"),
  measurement_symbol = sym("Value"),
  significance_limit = 0.05)
### Diagnosis
kth_diagnosis_grouped_pvals <- perform_kw_wilcoxon_according_to_grouping(
  melted_df = melted_kth,
  differentiating_feature_symbol = sym("HGNC_Symbol"),
  grouping_feature_symbol = sym("Diagnosis"),
  measurement_symbol = sym("Value"),
  significance_limit = 0.05)

kth_summary <- melted_kth %>%
dplyr::summarise(mean = mean(Value),
                 sd = sd(Value),
                 median = median(Value),
                 .by = c(Diagnosis, Age_Group, Gender, HGNC_Symbol))

### Diagnosis + age group
kth_diagnosis_age_group_grouped_pvals <- perform_kw_wilcoxon_according_to_grouping(
  melted_df = supercharged_kth_melt,
  differentiating_feature_symbol = sym("HGNC_Symbol"),
  grouping_feature_symbol = sym("Diagnosis_Age_Group"),
  measurement_symbol = sym("Value"),
  significance_limit = 0.05)

significant_kth_diagnosis_age_group_proteins <-
  kth_diagnosis_age_group_grouped_pvals[[1]]$HGNC_Symbol

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
###                                EMIF                                    ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

### Gender
emif_gender_grouped_pvals <- perform_kw_wilcoxon_according_to_grouping(
  melted_df = melted_emif,
  differentiating_feature_symbol = sym("UniProt"),
  grouping_feature_symbol = sym("Gender"),
  measurement_symbol = sym("Value"),
  significance_limit = 0.05)
### Age group
emif_age_group_grouped_pvals <- perform_kw_wilcoxon_according_to_grouping(
  melted_df = melted_emif,
  differentiating_feature_symbol = sym("UniProt"),
  grouping_feature_symbol = sym("Age_Group"),
  measurement_symbol = sym("Value"),
  significance_limit = 0.05)
### Diagnosis
emif_diagnosis_grouped_pvals <- perform_kw_wilcoxon_according_to_grouping(
  melted_df = melted_emif,
  differentiating_feature_symbol = sym("UniProt"),
  grouping_feature_symbol = sym("Diagnosis"),
  measurement_symbol = sym("Value"),
  significance_limit = 0.05)
### Diagnosis + age group
emif_diagnosis_age_group_grouped_pvals <- perform_kw_wilcoxon_according_to_grouping(
  melted_df = supercharged_emif_melt,
  differentiating_feature_symbol = sym("HGNC_Symbol"),
  grouping_feature_symbol = sym("Diagnosis_Age_Group"),
  measurement_symbol = sym("Value"),
  significance_limit = 0.05)

significant_emif_diagnosis_age_group_proteins <-
  emif_diagnosis_age_group_grouped_pvals[[1]]$HGNC_Symbol

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
###                                 OLINK                                  ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
### Gender
olink_gender_grouped_pvals <- perform_kw_wilcoxon_according_to_grouping(
  melted_df = melted_olink,
  differentiating_feature_symbol = sym("Gene_name"),
  grouping_feature_symbol = sym("sex"),
  measurement_symbol = sym("Value"),
  significance_limit = 0.05)
### Age group
olink_age_group_grouped_pvals <- perform_kw_wilcoxon_according_to_grouping(
  melted_df = melted_olink,
  differentiating_feature_symbol = sym("Gene_name"),
  grouping_feature_symbol = sym("Age_Group"),
  measurement_symbol = sym("Value"),
  significance_limit = 0.05)
### Diagnosis
olink_diagnosis_grouped_pvals <- perform_kw_wilcoxon_according_to_grouping(
  melted_df = melted_olink,
  differentiating_feature_symbol = sym("Gene_name"),
  grouping_feature_symbol = sym("dx"),
  measurement_symbol = sym("Value"),
  significance_limit = 0.05)
### Diagnosis + age group
olink_diagnosis_age_group_grouped_pvals <- perform_kw_wilcoxon_according_to_grouping(
  melted_df = supercharged_olink_melt,
  differentiating_feature_symbol = sym("HGNC_Symbol"),
  grouping_feature_symbol = sym("Diagnosis_Age_Group"),
  measurement_symbol = sym("Value"),
  significance_limit = 0.05)

significant_olink_diagnosis_age_group_proteins <-
  olink_diagnosis_age_group_grouped_pvals[[1]]$HGNC_Symbol

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
### Combined datasets ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
###Gender
combined_gender_grouped_pvals <- perform_kw_wilcoxon_according_to_grouping(
  melted_df = melted_combined_df,
  differentiating_feature_symbol = sym("HGNC_Symbol"),
  grouping_feature_symbol = sym("Gender"),
  measurement_symbol = sym("Value"),
  significance_limit = 0.05)

### Age_Group
combined_age_group_grouped_pvals <- perform_kw_wilcoxon_according_to_grouping(
  melted_df = melted_combined_df,
  differentiating_feature_symbol = sym("HGNC_Symbol"),
  grouping_feature_symbol = sym("Age_Group"),
  measurement_symbol = sym("Value"),
  significance_limit = 0.05)
### Diagnosis
combined_diagnosis_grouped_pvals <- perform_kw_wilcoxon_according_to_grouping(
  melted_df = melted_combined_df,
  differentiating_feature_symbol = sym("HGNC_Symbol"),
  grouping_feature_symbol = sym("Diagnosis"),
  measurement_symbol = sym("Value"),
  significance_limit = 0.05)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
## Box plots for diagnosis + age group combinations ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# KTH

kth_boxplots <-
  create_box_plots_of_all_columns_starting_at_column_number(
    mutated_kth %>%
      dplyr::select(Diagnosis_Age_Group, all_of(significant_kth_diagnosis_age_group_proteins)),
    "KTH", "Diagnosis_Age_Group", 2)

# EMIF

emif_boxplots <-
  create_box_plots_of_all_columns_starting_at_column_number(
    mutated_emif %>%
      dplyr::select(Diagnosis_Age_Group, all_of(significant_emif_diagnosis_age_group_proteins)),
    "Emif", "Diagnosis_Age_Group", 2)

# Olink
olink_boxplots <-
  create_box_plots_of_all_columns_starting_at_column_number(
    mutated_olink %>%
      dplyr::select(Diagnosis_Age_Group, all_of(significant_olink_diagnosis_age_group_proteins)),
    "Olink", "Diagnosis_Age_Group", 2)
# Recommended to remove once done because of size
remove(olink_boxplots)
