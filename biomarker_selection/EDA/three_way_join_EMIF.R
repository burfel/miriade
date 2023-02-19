##################################################
## Project: MIRIADE
## Script purpose: Perform data analysis on
## Date: 16.11.2022
## Author: Felicia Burtscher
##################################################
options(stringsAsFactors = F)

library(dplyr)
library(readxl)
library(ggplot2)
library(here)
source(here("functions", "data_processing_functions.R"))

datasets_root_directory <- define_datasets_root()

#three_way_join_append
### read the 3 files
emif_protein_data <- read.table(file.path(datasets_root_directory, "EMIF-AD MBD Study/data_protein.tsv"),
                                sep = "\t", header = T)
emif_samples_data <- read.table(file.path(datasets_root_directory, "EMIF-AD MBD Study/samples.tsv"),
                                          sep = "\t", header = T)
emif_clinical_data <- read.table(file.path(datasets_root_directory, "EMIF-AD MBD Study/Clinical_EMIFMBD_MIRIADE.csv"),
                                 sep = ";", header = T)

three_way_join_append <- function(target_df, intermediate_df, source_df, 
                                  target_to_intermediate_matchers, 
                                  intermediate_to_source_matchers, appended_columns) {
  temp_df <- dplyr::inner_join(intermediate_df, source_df, by = intermediate_to_source_matchers)
  new_df <- dplyr::inner_join(target_df, temp_df, by = target_to_intermediate_matchers) %>%
    dplyr::select(names(target_df), all_of(appended_columns))
  return(new_df)
}

# '.' replaces ' ' for the column names
# We do lose (558843 - 513389) rows out of 558843 in the process
emif_protein_data_with_gender <- three_way_join_append(target_df = emif_protein_data, intermediate_df = emif_samples_data,
                                source_df = emif_clinical_data, 
                                target_to_intermediate_matchers =  c("Assay.ID"),
                                intermediate_to_source_matchers =  c("Subject.ID" = "SubjectId"),
                                appended_column = "Gender") %>%
  dplyr::mutate(Gender = case_when(Gender == 0 ~ 'm', Gender == 1 ~ 'f'))
