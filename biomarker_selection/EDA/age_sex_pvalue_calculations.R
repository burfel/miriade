library(here)
library(dplyr)

# This loads the melded emif and olink dataframes to be used later in the file
# It also calls to set the dataset root directory
source(here("biomarker_selection", "EDA", "meld_all_fractured_datasets.R"))
source(here("functions", "auxiliary_statistical_functions.R"))
source(here("functions", "boxplot_functions.R"))

### Read the table
kth <- readxl::read_xlsx(file.path(datasets_root_directory,
                                   "KTH/KTH AD dataset for MIRIADE biomarker selection sample info and results.xlsx"))

### Change the names into HGNC symbols (prefixes)
colnames(kth) <- sapply(colnames(kth), take_only_pre_underscore_substring)

### Melt the dataframe, so there's only one readout variable
melted_kth <-
  reshape2::melt(kth, id = 1:9,
                 variable.name = "HGNC_Symbol", value.name = "median_ab_readout") %>%
  #dplyr::filter(Diagnosis == "AD") %>%
  convert_age_column_to_age_group_column(sym("Age"), sym("Age_Group")) %>%
  dplyr::select(Diagnosis, Age_Group, Gender, HGNC_Symbol, median_ab_readout)


### Run Kruskal-Wallis test for all biomarkers grouped by gender
kth_gender_grouped_pvals <- perform_kw_wilcoxon_according_to_grouping(
  melted_df = melted_kth,
  differentiating_feature_symbol = sym("HGNC_Symbol"),
  grouping_feature_symbol = sym("Gender"),
  measurement_symbol = sym("median_ab_readout"),
  significance_limit = 0.05)

### Age_Group
kth_age_group_grouped_pvals <- perform_kw_wilcoxon_according_to_grouping(
  melted_df = melted_kth,
  differentiating_feature_symbol = sym("HGNC_Symbol"),
  grouping_feature_symbol = sym("Age_Group"),
  measurement_symbol = sym("median_ab_readout"),
  significance_limit = 0.05)
### Diagnosis
kth_diagnosis_grouped_pvals <- perform_kw_wilcoxon_according_to_grouping(
  melted_df = melted_kth,
  differentiating_feature_symbol = sym("HGNC_Symbol"),
  grouping_feature_symbol = sym("Diagnosis"),
  measurement_symbol = sym("median_ab_readout"),
  significance_limit = 0.05)

kth_summary <- melted_kth %>%
dplyr::summarise(mean = mean(median_ab_readout),
                 sd = sd(median_ab_readout),
                 median = median(median_ab_readout),
                 .by = c(Diagnosis, Age_Group, Gender, HGNC_Symbol))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### Adjusting to have Diagnosis + age group combined column                ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

supercharged_melt <-
  melted_kth %>%
  #dplyr::filter(Diagnosis == "AD" | Diagnosis == "Control") %>%
  dplyr::mutate(Diagnosis_Age_Group = paste(Diagnosis, Age_Group)) %>%
  dplyr::select(Diagnosis_Age_Group, Gender, HGNC_Symbol, median_ab_readout)

### Diagnosis + age group
kth_diagnosisage_group_grouped_pvals <- perform_kw_wilcoxon_according_to_grouping(
  melted_df = supercharged_melt,
  differentiating_feature_symbol = sym("HGNC_Symbol"),
  grouping_feature_symbol = sym("Diagnosis_Age_Group"),
  measurement_symbol = sym("median_ab_readout"),
  significance_limit = 0.05)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
####                                EMIF                                    ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

### Melt the emif dataframe, so there's only one readout variable
melted_emif <-
  reshape2::melt(emif, id = 1:5,
                 variable.name = "UniProt", value.name = "Value") %>%
  convert_age_column_to_age_group_column(sym("Age"), sym("Age_Group")) %>%
  dplyr::select(Diagnosis, Age_Group, Gender, UniProt, Value)

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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
####                                 OLINK                                  ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

### Melt the olink dataframe, so there's only one readout variable
melted_olink <-
  reshape2::melt(olink, id = 1:4,
                 variable.name = "Gene_name", value.name = "Value") %>%
  convert_age_column_to_age_group_column(sym("age_gr"), sym("Age_Group")) %>%
  dplyr::select(dx, Age_Group, sex, Gene_name, Value)

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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### Box plots for diagnosis + age group combinations ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# KTH
mutated_kth <- kth %>%
  convert_age_column_to_age_group_column(sym("Age"), sym("Age_Group")) %>%
  dplyr::filter(Diagnosis == "AD" | Diagnosis == "Control") %>%
  dplyr::mutate(Diagnosis_Age_Group = paste(Diagnosis, Age_Group)) %>%
  dplyr::select(Diagnosis_Age_Group, Age_Group, everything())

kth_boxplots <-
  create_box_plots_of_all_columns_starting_at_column_number(mutated_kth, "Diagnosis_Age_Group", 7)

# EMIF
source(here("functions", "mapping_functions.R"))
emif_proteins <- data.frame(UniProt = colnames(emif)[6:ncol(emif)]) %>%
  map_uniprot_to_hgnc_and_cbind("Name") %>%
  na.omit() %>%
  dplyr::select(UniProt, Name)
# use sapply to create a named vector of the hgnc and their new uniprots
protein_renaming_vector <- sapply(emif_proteins$Name, function(x) {
  emif_proteins$UniProt[which(emif_proteins$Name == x)]
})
mutated_emif <- emif %>%
  convert_age_column_to_age_group_column(sym("Age"), sym("Age_Group")) %>%
  dplyr::filter(Diagnosis == "AD" | Diagnosis == "NL") %>%
  dplyr::mutate(Diagnosis = case_when(Diagnosis == "NL" ~ "Control", TRUE ~ Diagnosis)) %>%
  dplyr::mutate(Diagnosis_Age_Group = paste(Diagnosis, Age_Group)) %>%
  dplyr::select(Diagnosis_Age_Group, Age_Group, everything()) %>%
  dplyr::rename(!!!protein_renaming_vector)

emif_boxplots <-
  create_box_plots_of_all_columns_starting_at_column_number(mutated_emif, "Diagnosis_Age_Group", 8)
# Recommended to remove once done because of size
remove(emif_boxplots)
# Olink
mutated_olink <- olink %>%
  convert_age_column_to_age_group_column(sym("age_gr"), sym("Age_Group")) %>%
  dplyr::filter(dx == "AD dementia" | dx == "CN") %>%
  dplyr::mutate(Diagnosis = case_when(dx == "CN" ~ "Control", dx == "AD dementia" ~ "AD", TRUE ~ dx)) %>%
  dplyr::mutate(Diagnosis_Age_Group = paste(Diagnosis, Age_Group)) %>%
  dplyr::select(Diagnosis_Age_Group, Diagnosis, Age_Group, everything(), -dx)

olink_boxplots <-
  create_box_plots_of_all_columns_starting_at_column_number(mutated_olink, "Diagnosis_Age_Group", 7)
# Recommended to remove once done because of size
remove(olink_boxplots)
