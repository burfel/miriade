library(here)
library(corrr)
library(dplyr)
library(ggplot2)

source(here("functions", "data_processing_functions.R"))
source(here("functions", "heatmap_functions.R"))

datasets_root_directory <- define_datasets_root()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
####                                Olink                                   ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

olink <- read_df_from_file(file.path(datasets_root_directory,
                                     "Olink/OLINK basic data files",
                                     "protein_data_LODmax_10percent_outliers_rm_missing_imputed_16112020.tsv"))

olink_corr <- corrr::correlate(olink[,4:ncol(olink)]) %>%
  stretch()

plot_correlation_heatmap(olink_corr, "Olink")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
####                                KTH                                     ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
kth <- readxl::read_xlsx(file.path(datasets_root_directory,
                                   "KTH/KTH AD dataset for MIRIADE biomarker selection sample info and results.xlsx")) %>%
  dplyr::select(-c("Class", "Diagnosis", "Age", "Gender", "ApoE4", "MMSE_score"))
### Change the names into HGNC symbols (prefixes)
colnames(kth) <- sapply(colnames(kth), take_only_pre_underscore_substring)

kth_corr <- corrr::correlate(kth) %>% stretch()

plot_correlation_heatmap(kth_corr, "KTH")
