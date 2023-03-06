library(here)
library(corrr)
library(dplyr)
library(ggplot2)

source(here("functions", "data_processing_functions.R"))
source(here("functions", "auxiliary_statistical_functions.R"))

datasets_root_directory <- define_datasets_root()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
####                                Olink                                   ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

olink <- read_df_from_file(file.path(datasets_root_directory,
                                     "Olink/OLINK basic data files",
                                     "protein_data_LODmax_10percent_outliers_rm_missing_imputed_16112020.tsv"))

olink_corr <- corrr::correlate(olink[,4:ncol(olink)]) %>%
  stretch()

olink_corr %>%
  ggplot(aes(x, y, fill = r)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
####                                KTH                                     ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
kth <- readxl::read_xlsx(file.path(datasets_root_directory,
                                   "KTH/KTH AD dataset for MIRIADE biomarker selection sample info and results.xlsx")) %>%
  dplyr::select(-c("Class", "Diagnosis", "Age", "Gender", "ApoE4", "MMSE_score"))
### Change the names into HGNC symbols (prefixes)
colnames(kth) <- sapply(colnames(kth), take_only_pre_underscore_substring)

kth_corr <- corrr::correlate(kth) %>% stretch()

kth_corr %>%
  ggplot(aes(x, y, fill = r)) +
  geom_tile()  +
  labs(x = NULL, y = NULL, fill = "Pearson's\nCorrelation", title="Correlations in KTH") +
  # map a red, white and blue color scale to correspond to -1:1 sequential gradient
  scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", limits=c(-1,1)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  # remove excess space on x and y axes
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0))
