##################################################
## Project: MIRIADE
## Script purpose: Explore the Olink dataset
## Date: 03.10.2022
## Author: Felicia Burtscher
##################################################
options(stringsAsFactors = F)

library(dplyr)
library(readxl)
library(ggplot2)
library(here)
source(here("biomarker_selection", "EDA", "significance_olink_prep.R"))

# Read the KTH
kth <- read.table("/Users/felicia.burtscher/Documents/UL/DATASETS/KTH/cleaned_KTH_entries.tsv",
                  sep = "\t", header = T, quote = "")
adni <- read.table("/Users/felicia.burtscher/Documents/UL/DATASETS/Brain 2020 AD/cleaned_ADNI_entries.tsv",
                   sep = "\t", header = T, quote = "")
emif <- read.table("/Users/felicia.burtscher/Documents/UL/DATASETS/Brain 2020 AD/cleaned_EMIF_entries.tsv",
                   sep = "\t", header = T, quote = "")
mspec <- read.table("/Users/felicia.burtscher/Documents/UL/DATASETS/VUMC/cleaned_VUMC_mspec.tsv",
                    sep = "\t", header = T, quote = "")

# Preprocess datasets to have the unified format of UniProt and Pvalue, possibly after filtering
# some of them according to disease.
# datasets and should_filter_by_disease must be lists/vectors of the same length
# datasets HAS TO BE A LIST OF DATAFRAMES
# Any dataframe for which should_filter_by_disease is true needs to have a disease column
# If any of should_filter_by_disease is true, disease_to_filter_by should be a string that appears in the disease column
preprocess_datasets <- function(datasets, should_filter_by_disease, disease_to_filter_by = NULL) {
  datasets_length <- length(datasets)
  if (datasets_length != length(should_filter_by_disease)) {
    stop("datasets and should_filter_by_disease must be lists/vectors of the same length")
  }
  processed_datasets <- vector(mode = "list", length = datasets_length)
  ind <- 1
  for (dataset in datasets) {
    if(should_filter_by_disease[ind]) {
      dataset <- dataset %>% dplyr::filter(disease == disease_to_filter_by)
    }
    dataset <- dataset %>% dplyr::select(UniProt, p.val)
    processed_datasets[[ind]] <- dataset
    ind <- ind + 1
  }
  return(processed_datasets)
}

# Binds datasets together and produces a dataframe with UniProt p_values and a column for the source of that row
# TODO: fold the mutate into this function once we realize how to use case_when with expressions built
# on the fly with indices unknown in advance
bind_datasets_together <- function(processed_datasets, should_filter_by_disease, disease_to_filter_by = NULL) {
  # Unfortunately, failed to do the renaming of source here... might have to repeat it
  return(dplyr::bind_rows(processed_datasets, .id = "Source")
  )
}

# Compare olink and kth for AZ
alz_olink_and_kth <- bind_datasets_together(
  preprocess_datasets(
  datasets = list(vumc_ol,kth),
  should_filter_by_disease = c(TRUE,FALSE),
  disease_to_filter_by = "ALZ")) %>%
  dplyr::mutate(Source = case_when(Source == 1 ~ 'olink', Source == 2 ~ 'kth'))

# This function takes a dataframe with uniprots and p values from more than one source
# filters and returns according to the list of limiting dataframes
filter_according_to_UniProt_masks <- function(data_frame, masking_frames_vector) {
  filtered_data_frame <- data_frame
  for (mask in masking_frames_vector) {
    filtered_data_frame <- filtered_data_frame %>%
      dplyr::filter(UniProt %in% mask$UniProt)
  }
  return(filtered_data_frame)
}

alz_olink_and_kth_intersection_list <- filter_according_to_UniProt_masks(alz_olink_and_kth, list(kth, vumc_ol))

# Takes a datafrane with a UniProt column, and extracts the distinct UniProts
extract_unique_uniprots <- function(dataframe) {
  return(dataframe %>%
           dplyr::select(UniProt) %>% dplyr::distinct(UniProt))
}

alz_olink_and_kth_intersection_uniprots <- extract_unique_uniprots(alz_olink_and_kth_intersection_list)

preprocessed_olink_and_adni <- preprocess_datasets(
  datasets = list(vumc_ol,adni),
  should_filter_by_disease = c(TRUE,FALSE),
  disease_to_filter_by = "ALZ")
# Now do the same thing we did with olink vs. kth with the rest
alz_olink_and_adni <- bind_datasets_together(preprocessed_olink_and_adni) %>%
  dplyr::mutate(Source = case_when(Source == 1 ~ 'olink', Source == 2 ~ 'adni'))

alz_olink_and_adni_intersection_list <- filter_according_to_UniProt_masks(alz_olink_and_adni,
                                                                          preprocessed_olink_and_adni)
alz_olink_and_adni_intersection_uniprots <- extract_unique_uniprots(alz_olink_and_adni_intersection_list)
alz_olink_and_adni_intersection_significant_in_both_datasets <- 
  dplyr::inner_join(preprocessed_olink_and_adni[[1]], preprocessed_olink_and_adni[[2]], by = "UniProt",
                   suffix = c("_olink", "_adni")) %>%
  dplyr::filter(p.val_olink <= 0.05 & p.val_adni <=0.05)
remove(preprocessed_olink_and_adni)

preprocessed_olink_and_emif <- preprocess_datasets(
  datasets = list(vumc_ol,emif),
  should_filter_by_disease = c(TRUE,FALSE),
  disease_to_filter_by = "ALZ")

alz_olink_and_emif <- bind_datasets_together(preprocessed_olink_and_emif) %>%
  dplyr::mutate(Source = case_when(Source == 1 ~ 'olink', Source == 2 ~ 'emif'))

alz_olink_and_emif_intersection_list <- filter_according_to_UniProt_masks(alz_olink_and_emif,
                                                                          preprocessed_olink_and_emif)
alz_olink_and_emif_intersection_uniprots <- extract_unique_uniprots(alz_olink_and_emif_intersection_list)
alz_olink_and_emif_intersection_significant_in_both_datasets <- 
  dplyr::inner_join(preprocessed_olink_and_emif[[1]], preprocessed_olink_and_emif[[2]], by = "UniProt",
                    suffix = c("_olink", "_emif")) %>%
  dplyr::filter(p.val_olink <= 0.05 & p.val_emif <=0.05)
remove(preprocessed_olink_and_emif)

preprocessed_adni_and_emif <- preprocess_datasets(
  datasets = list(adni,emif),
  should_filter_by_disease = c(FALSE,FALSE))

alz_adni_and_emif <- bind_datasets_together(preprocessed_adni_and_emif) %>%
  dplyr::mutate(Source = case_when(Source == 1 ~ 'adni', Source == 2 ~ 'emif'))

alz_adni_and_emif_intersection_list <- filter_according_to_UniProt_masks(alz_adni_and_emif,
                                                                         preprocessed_adni_and_emif)
alz_adni_and_emif_intersection_uniprots <- extract_unique_uniprots(alz_adni_and_emif_intersection_list)
alz_adni_and_emif_intersection_significant_in_both_datasets <- 
  dplyr::inner_join(preprocessed_adni_and_emif[[1]], preprocessed_adni_and_emif[[2]], by = "UniProt",
                    suffix = c("_adni", "_emif")) %>%
  dplyr::filter(p.val_adni <= 0.05 & p.val_emif <=0.05)
remove(preprocessed_adni_and_emif)

# All 4 ALZ datasets
preprocessed_olink_adni_emif_kth <- preprocess_datasets(
  datasets = list(vumc_ol,adni,emif,kth),
  should_filter_by_disease = c(TRUE,FALSE,FALSE,FALSE),
  disease_to_filter_by = "ALZ")
alz_olink_adni_emif_kth <- bind_datasets_together(preprocessed_olink_adni_emif_kth) %>%
  dplyr::mutate(Source = case_when(Source == 1 ~ 'olink', Source == 2 ~ 'adni', Source == 3 ~ 'emif', Source == 4 ~ 'kth'))

alz_olink_adni_emif_kth_intersection_list <- filter_according_to_UniProt_masks(alz_olink_adni_emif_kth,
                                                                               preprocessed_olink_adni_emif_kth)
remove(preprocessed_olink_adni_emif_kth)

# Without KTH
preprocessed_olink_adni_emif <- preprocess_datasets(
  datasets = list(vumc_ol,adni,emif),
  should_filter_by_disease = c(TRUE,FALSE,FALSE),
  disease_to_filter_by = "ALZ")
alz_olink_adni_emif <- bind_datasets_together(preprocessed_olink_adni_emif) %>%
  dplyr::mutate(Source = case_when(Source == 1 ~ 'olink', Source == 2 ~ 'adni', Source == 3 ~ 'emif'))

alz_olink_adni_emif_intersection_list <- filter_according_to_UniProt_masks(alz_olink_adni_emif,
                                                                           preprocessed_olink_adni_emif)
alz_olink_adni_emif_intersection_uniprots <- extract_unique_uniprots(alz_olink_adni_emif_intersection_list)
remove(preprocessed_olink_adni_emif)

# Now Olink and mspec
preprocessed_olink_mspec_dlb <- preprocess_datasets(
  datasets = list(vumc_ol,mspec),
  should_filter_by_disease = c(TRUE,TRUE),
  disease_to_filter_by = "DLB")
dlb_olink_mspec <- bind_datasets_together(preprocessed_olink_mspec_dlb) %>%
  dplyr::mutate(Source = case_when(Source == 1 ~ 'olink', Source == 2 ~ 'mspec'))
dlb_olink_mspec_intersection_list <- filter_according_to_UniProt_masks(dlb_olink_mspec, preprocessed_olink_mspec_dlb)
dlb_olink_mspec_intersection_uniprots <- extract_unique_uniprots(dlb_olink_mspec_intersection_list)
dlb_olink_mspec_intersection_significant_in_both_datasets <- 
  dplyr::inner_join(preprocessed_olink_mspec_dlb[[1]], preprocessed_olink_mspec_dlb[[2]], by = "UniProt",
                    suffix = c("_olink", "_mspec")) %>%
  dplyr::filter(p.val_olink <= 0.05 & p.val_mspec <=0.05)
remove(preprocessed_olink_mspec_dlb)

preprocessed_olink_mspec_ftd <- preprocess_datasets(
  datasets = list(vumc_ol,mspec),
  should_filter_by_disease = c(TRUE,TRUE),
  disease_to_filter_by = "FTD")
ftd_olink_mspec <- bind_datasets_together(preprocessed_olink_mspec_ftd) %>%
  dplyr::mutate(Source = case_when(Source == 1 ~ 'olink', Source == 2 ~ 'mspec'))
ftd_olink_mspec_intersection_list <- filter_according_to_UniProt_masks(ftd_olink_mspec, preprocessed_olink_mspec_ftd)
ftd_olink_mspec_intersection_uniprots <- extract_unique_uniprots(ftd_olink_mspec_intersection_list)
ftd_olink_mspec_intersection_significant_in_both_datasets <- 
  dplyr::inner_join(preprocessed_olink_mspec_ftd[[1]], preprocessed_olink_mspec_ftd[[2]], by = "UniProt",
                    suffix = c("_olink", "_mspec")) %>%
  dplyr::filter(p.val_olink <= 0.05 & p.val_mspec <=0.05)
remove(preprocessed_olink_mspec_ftd)
