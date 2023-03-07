##################################################
## Project: MIRIADE
## Script purpose: Explore the Olink dataset
## Date: 03.10.2022
## Author: Felicia Burtscher
##################################################
library(here)
library(dplyr)
source(here("biomarker_selection", "EDA", "significance_olink_prep.R"))
source(here("functions", "mapping_functions.R"))
source(here("functions", "comparison_functions.R"))

# Read the KTH
kth <- read.table(file.path(datasets_root_directory, "KTH/cleaned_KTH_entries.tsv"),
                  sep = "\t", header = T, quote = "")
adni <- read.table(file.path(datasets_root_directory, "Brain 2020 AD/cleaned_ADNI_entries.tsv"),
                   sep = "\t", header = T, quote = "")
emif <- read.table(file.path(datasets_root_directory, "Brain 2020 AD/cleaned_EMIF_entries.tsv"),
                   sep = "\t", header = T, quote = "") # Has two duplicate uniprots
mspec <- read.table(file.path(datasets_root_directory, "VUMC/cleaned_VUMC_mspec.tsv"),
                    sep = "\t", header = T, quote = "")

vumc_ol <- map_uniprot_to_hgnc_and_cbind(vumc_ol, "Name")

# Compare olink and kth for AZ
alz_olink_and_kth <- compare_datasets(
  datasets = list(vumc_ol, kth),
  dataset_names = c('olink', 'kth'),
  should_filter_by_disease = c(TRUE, FALSE),
  additional_columns = c("Name"),
  disease_to_filter_by = "ALZ")

# Now do the same thing we did with olink vs. kth with the rest
alz_olink_and_adni <- compare_datasets(
  datasets = list(vumc_ol, adni),
  dataset_names = c('olink', 'adni'),
  should_filter_by_disease = c(TRUE, FALSE),
  additional_columns = c("Name"),
  disease_to_filter_by = "ALZ")

alz_olink_and_emif <- compare_datasets(
  datasets = list(vumc_ol, emif),
  dataset_names = c('olink', 'emif'),
  should_filter_by_disease = c(TRUE, FALSE),
  additional_columns = c("Name"),
  disease_to_filter_by = "ALZ")

alz_adni_and_emif <- compare_datasets(
  datasets = list(adni, emif),
  dataset_names = c('adni', 'emif'),
  should_filter_by_disease = c(FALSE, FALSE),
  additional_columns = c("Name"))

# All 4 ALZ datasets
alz_olink_adni_emif_kth <- compare_datasets(
  datasets = list(vumc_ol, adni, emif, kth),
  dataset_names = c('olink', 'adni', 'emif', 'kth'),
  should_filter_by_disease = c(TRUE, FALSE, FALSE, FALSE),
  additional_columns = c("Name"),
  disease_to_filter_by = "ALZ")

# Without KTH
alz_olink_adni_emif <- compare_datasets(
  datasets = list(vumc_ol, adni, emif),
  dataset_names = c('olink', 'adni', 'emif'),
  should_filter_by_disease = c(TRUE, FALSE, FALSE),
  additional_columns = c("Name"),
  disease_to_filter_by = "ALZ")

# Now Olink and mspec
dlb_olink_mspec <- compare_datasets(
  datasets = list(vumc_ol, mspec),
  dataset_names = c('olink', 'mspec'),
  should_filter_by_disease = c(TRUE, TRUE),
  additional_columns = c("Name"),
  disease_to_filter_by = "DLB")

ftd_olink_mspec <- compare_datasets(
  datasets = list(vumc_ol, mspec),
  dataset_names = c('olink', 'mspec'),
  should_filter_by_disease = c(TRUE, TRUE),
  additional_columns = c("Name"),
  disease_to_filter_by = "FTD")

# Now adni and kth
alz_adni_kth <- compare_datasets(
  datasets = list(adni, kth),
  dataset_names = c('adni', 'kth'),
  should_filter_by_disease = c(FALSE, FALSE),
  additional_columns = c("Name"))

# Now emif and kth
alz_emif_kth <- compare_datasets(
  datasets = list(emif, kth),
  dataset_names = c('emif', 'kth'),
  should_filter_by_disease = c(FALSE, FALSE),
  additional_columns = c("Name"))
