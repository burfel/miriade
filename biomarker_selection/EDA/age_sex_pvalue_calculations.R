library(here)
library(dplyr)

# This calls meld_all_fractured_datasets, sets the dataset root directory,
# and prepares kth (/emif/olink), mutated_kth(/emif/olink),
# melted_kth(/emif/olink) and supercharged_kth(/emif/olink)_melt to be used
# in this file 
source(here("biomarker_selection", "EDA", "prep_datasets_for_kw_wilcoxon.R"))
source(here("functions", "boxplot_functions.R"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
## Processing the datasets ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
### Combine the datasets ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#### First, we bring all 3 datasets to a standard form ####
kth <- kth %>%
  dplyr::filter(Diagnosis == "AD" | Diagnosis == "Control")

emif <- emif %>%
  dplyr::select(-SubjectId, -Assay.ID) %>%
  dplyr::filter(Diagnosis == "Control" | Diagnosis == "AD") %>%
  dplyr::rename(!!!protein_renaming_vector)

olink <- olink %>%
  dplyr::select(-SampleId) %>%
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
