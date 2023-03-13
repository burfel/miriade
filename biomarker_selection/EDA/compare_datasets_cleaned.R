##################################################
## Project: MIRIADE
## Script purpose: Explore the Olink dataset
## Date: 03.10.2022
## Author: Felicia Burtscher
##################################################
# DISCLAIMER: This script is using "cleaned" datasets, i.e. datasets that were
# already processed, underwent statistical tests, and then only significant
# proteins accordig to those tests were preserved. Therefore, it does not contain
# the whole data from the datasets.
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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
####                             Venn diagrams                              ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
library(ggvenn)

generate_venn_data <- function(dataset_comparison_list)
{
  dataset_names <- dataset_comparison_list$dataset_names
  result <-
  sapply(dataset_names,
         function(dataset_name)
           dplyr::filter(dataset_comparison_list$raw, Source == dataset_name) %>%
           dplyr::select(UniProt))
  names(result) <- dataset_names

  return(result)
}

# For the colors, we want: Yellow for olink, red for kth, blue for adni,
# green for emif, and purple for mspec
alz_olink_and_kth$venn_diagram <-
  ggvenn(generate_venn_data(alz_olink_and_kth),
         fill_color = c("yellow", "red")) +
  labs(title = "Olink and KTH overlap", subtitle = "Alzheimer")
alz_olink_and_adni$venn_diagram <-
  ggvenn(generate_venn_data(alz_olink_and_adni),
         fill_color = c("yellow", "blue")) +
  labs(title = "Olink and Adni overlap", subtitle = "Alzheimer")
alz_olink_and_emif$venn_diagram <-
  ggvenn(generate_venn_data(alz_olink_and_emif),
         fill_color = c("yellow", "green")) +
  labs(title = "Olink and Emif overlap", subtitle = "Alzheimer")
alz_adni_and_emif$venn_diagram <-
  ggvenn(generate_venn_data(alz_adni_and_emif),
         fill_color = c("blue", "green")) +
  labs(title = "Adni and Emif overlap", subtitle = "Alzheimer")
alz_olink_adni_emif_kth$venn_diagram <-
  ggvenn(generate_venn_data(alz_olink_adni_emif_kth),
         fill_color = c("yellow", "blue", "green", "red")) +
  labs(title = "Olink, Adni, Emif and KTH overlap", subtitle = "Alzheimer")
alz_olink_adni_emif$venn_diagram <-
  ggvenn(generate_venn_data(alz_olink_adni_emif),
         fill_color = c("yellow", "blue", "green")) +
  labs(title = "Olink, Adni and Emif overlap", subtitle = "Alzheimer")
dlb_olink_mspec$venn_diagram <-
  ggvenn(generate_venn_data(dlb_olink_mspec),
         fill_color = c("yellow", "purple")) +
  labs(title = "Olink and Mspec overlap", subtitle = "DLB")
ftd_olink_mspec$venn_diagram <-
  ggvenn(generate_venn_data(ftd_olink_mspec),
         fill_color = c("yellow", "purple")) +
  labs(title = "Olink and Mspec overlap", subtitle = "FTD")
alz_adni_kth$venn_diagram <-
  ggvenn(generate_venn_data(alz_adni_kth),
         fill_color = c("blue", "red")) +
  labs(title = "Adni and KTH overlap", subtitle = "Alzheimer")
alz_emif_kth$venn_diagram <-
  ggvenn(generate_venn_data(alz_emif_kth),
         fill_color = c("green", "red")) +
  labs(title = "Emif and KTH overlap", subtitle = "Alzheimer")
# showing the plots
alz_olink_and_kth$venn_diagram
alz_olink_and_adni$venn_diagram
alz_olink_and_emif$venn_diagram
alz_adni_and_emif$venn_diagram
alz_olink_adni_emif_kth$venn_diagram
alz_olink_adni_emif$venn_diagram
dlb_olink_mspec$venn_diagram
ftd_olink_mspec$venn_diagram
alz_adni_kth$venn_diagram
alz_emif_kth$venn_diagram
