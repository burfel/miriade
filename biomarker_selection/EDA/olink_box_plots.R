# 15.11.2022
# visualising the original Olink data as box plots
library(here)
source(here("functions", "data_processing_functions.R"))
source(here("functions", "boxplot_functions.R"))

datasets_root_directory <- define_datasets_root()

olink_basic <- read.table(file.path(datasets_root_directory, "Olink/OLINK basic data files/protein_data_LODmax_10percent_outliers_rm_missing_imputed_16112020.txt"),
                          sep = "\t", header = T)

# This ends up being ~4.7 GB
olink_basic_box_plots <- create_box_plots_of_all_columns_starting_at_column_number(olink_basic, 6)
# Use an index between 1 and 807 (for proteins)
## TODO: perhaps implementing a dictionary to call uniprot-ids (once converted to uniprot ids) or standard gene names
olink_basic_box_plots[[67]]
# Recommend to remove it when done using it
remove(olink_basic_box_plots)
