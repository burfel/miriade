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
alz_olink_and_kth <- bind_datasets_together(
  preprocess_datasets(
  datasets = list(vumc_ol, kth),
  should_filter_by_disease = c(TRUE, FALSE),
  additional_columns = c("Name"),
  disease_to_filter_by = "ALZ")) %>%
  dplyr::mutate(Source = case_when(Source == 1 ~ 'olink', Source == 2 ~ 'kth'))

alz_olink_and_kth_intersection_list <- filter_according_to_UniProt_masks(alz_olink_and_kth, list(kth, vumc_ol))

alz_olink_and_kth_intersection_uniprots <- extract_unique_uniprots(alz_olink_and_kth_intersection_list, c("Name"))

preprocessed_olink_and_adni <- preprocess_datasets(
  datasets = list(vumc_ol, adni),
  should_filter_by_disease = c(TRUE, FALSE),
  additional_columns = c("Name"),
  disease_to_filter_by = "ALZ")
# Now do the same thing we did with olink vs. kth with the rest
alz_olink_and_adni <- bind_datasets_together(preprocessed_olink_and_adni) %>%
  dplyr::mutate(Source = case_when(Source == 1 ~ 'olink', Source == 2 ~ 'adni'))

alz_olink_and_adni_intersection_list <- filter_according_to_UniProt_masks(alz_olink_and_adni,
                                                                          preprocessed_olink_and_adni)
alz_olink_and_adni_intersection_uniprots <- extract_unique_uniprots(alz_olink_and_adni_intersection_list, c("Name"))
alz_olink_and_adni_intersection_significant_in_both_datasets <- 
  dplyr::inner_join(preprocessed_olink_and_adni[[1]], preprocessed_olink_and_adni[[2]], by = "UniProt",
                   suffix = c("_olink", "_adni")) %>%
  dplyr::filter(p.val_olink <= 0.05 & p.val_adni <=0.05)
remove(preprocessed_olink_and_adni)

preprocessed_olink_and_emif <- preprocess_datasets(
  datasets = list(vumc_ol, emif),
  should_filter_by_disease = c(TRUE, FALSE),
  additional_columns = c("Name"),
  disease_to_filter_by = "ALZ")

alz_olink_and_emif <- bind_datasets_together(preprocessed_olink_and_emif) %>%
  dplyr::mutate(Source = case_when(Source == 1 ~ 'olink', Source == 2 ~ 'emif'))

alz_olink_and_emif_intersection_list <- filter_according_to_UniProt_masks(alz_olink_and_emif,
                                                                          preprocessed_olink_and_emif)
alz_olink_and_emif_intersection_uniprots <- extract_unique_uniprots(alz_olink_and_emif_intersection_list, c("Name"))
alz_olink_and_emif_intersection_significant_in_both_datasets <- 
  dplyr::inner_join(preprocessed_olink_and_emif[[1]], preprocessed_olink_and_emif[[2]], by = "UniProt",
                    suffix = c("_olink", "_emif")) %>%
  dplyr::filter(p.val_olink <= 0.05 & p.val_emif <=0.05)
remove(preprocessed_olink_and_emif)

preprocessed_adni_and_emif <- preprocess_datasets(
  datasets = list(adni, emif),
  should_filter_by_disease = c(FALSE, FALSE),
  additional_columns = c("Name"))

alz_adni_and_emif <- bind_datasets_together(preprocessed_adni_and_emif) %>%
  dplyr::mutate(Source = case_when(Source == 1 ~ 'adni', Source == 2 ~ 'emif'))

alz_adni_and_emif_intersection_list <- filter_according_to_UniProt_masks(alz_adni_and_emif,
                                                                         preprocessed_adni_and_emif)
alz_adni_and_emif_intersection_uniprots <- extract_unique_uniprots(alz_adni_and_emif_intersection_list, c("Name"))
alz_adni_and_emif_intersection_significant_in_both_datasets <- 
  dplyr::inner_join(preprocessed_adni_and_emif[[1]], preprocessed_adni_and_emif[[2]], by = "UniProt",
                    suffix = c("_adni", "_emif")) %>%
  dplyr::filter(p.val_adni <= 0.05 & p.val_emif <=0.05)
remove(preprocessed_adni_and_emif)

# All 4 ALZ datasets
preprocessed_olink_adni_emif_kth <- preprocess_datasets(
  datasets = list(vumc_ol, adni, emif, kth),
  should_filter_by_disease = c(TRUE, FALSE, FALSE, FALSE),
  additional_columns = c("Name"),
  disease_to_filter_by = "ALZ")
alz_olink_adni_emif_kth <- bind_datasets_together(preprocessed_olink_adni_emif_kth) %>%
  dplyr::mutate(Source = case_when(Source == 1 ~ 'olink', Source == 2 ~ 'adni', Source == 3 ~ 'emif', Source == 4 ~ 'kth'))

alz_olink_adni_emif_kth_intersection_list <- filter_according_to_UniProt_masks(alz_olink_adni_emif_kth,
                                                                               preprocessed_olink_adni_emif_kth)
remove(preprocessed_olink_adni_emif_kth)

# Without KTH
preprocessed_olink_adni_emif <- preprocess_datasets(
  datasets = list(vumc_ol, adni, emif),
  should_filter_by_disease = c(TRUE, FALSE, FALSE),
  additional_columns = c("Name"),
  disease_to_filter_by = "ALZ")
alz_olink_adni_emif <- bind_datasets_together(preprocessed_olink_adni_emif) %>%
  dplyr::mutate(Source = case_when(Source == 1 ~ 'olink', Source == 2 ~ 'adni', Source == 3 ~ 'emif'))

alz_olink_adni_emif_intersection_list <- filter_according_to_UniProt_masks(alz_olink_adni_emif,
                                                                           preprocessed_olink_adni_emif)
alz_olink_adni_emif_intersection_uniprots <- extract_unique_uniprots(alz_olink_adni_emif_intersection_list, c("Name"))
remove(preprocessed_olink_adni_emif)

# Now Olink and mspec
preprocessed_olink_mspec_dlb <- preprocess_datasets(
  datasets = list(vumc_ol, mspec),
  should_filter_by_disease = c(TRUE, TRUE),
  additional_columns = c("Name"),
  disease_to_filter_by = "DLB")
dlb_olink_mspec <- bind_datasets_together(preprocessed_olink_mspec_dlb) %>%
  dplyr::mutate(Source = case_when(Source == 1 ~ 'olink', Source == 2 ~ 'mspec'))
dlb_olink_mspec_intersection_list <- filter_according_to_UniProt_masks(dlb_olink_mspec, preprocessed_olink_mspec_dlb)
dlb_olink_mspec_intersection_uniprots <- extract_unique_uniprots(dlb_olink_mspec_intersection_list, c("Name"))
dlb_olink_mspec_intersection_significant_in_both_datasets <- 
  dplyr::inner_join(preprocessed_olink_mspec_dlb[[1]], preprocessed_olink_mspec_dlb[[2]], by = "UniProt",
                    suffix = c("_olink", "_mspec")) %>%
  dplyr::filter(p.val_olink <= 0.05 & p.val_mspec <=0.05)
remove(preprocessed_olink_mspec_dlb)

preprocessed_olink_mspec_ftd <- preprocess_datasets(
  datasets = list(vumc_ol, mspec),
  should_filter_by_disease = c(TRUE, TRUE),
  additional_columns = c("Name"),
  disease_to_filter_by = "FTD")
ftd_olink_mspec <- bind_datasets_together(preprocessed_olink_mspec_ftd) %>%
  dplyr::mutate(Source = case_when(Source == 1 ~ 'olink', Source == 2 ~ 'mspec'))
ftd_olink_mspec_intersection_list <- filter_according_to_UniProt_masks(ftd_olink_mspec, preprocessed_olink_mspec_ftd)
ftd_olink_mspec_intersection_uniprots <- extract_unique_uniprots(ftd_olink_mspec_intersection_list, c("Name"))
ftd_olink_mspec_intersection_significant_in_both_datasets <- 
  dplyr::inner_join(preprocessed_olink_mspec_ftd[[1]], preprocessed_olink_mspec_ftd[[2]], by = "UniProt",
                    suffix = c("_olink", "_mspec")) %>%
  dplyr::filter(p.val_olink <= 0.05 & p.val_mspec <=0.05)
remove(preprocessed_olink_mspec_ftd)

# Now adni and kth
preprocessed_adni_and_kth <- preprocess_datasets(
  datasets = list(adni, kth),
  should_filter_by_disease = c(FALSE, FALSE),
  additional_columns = c("Name"))
alz_adni_kth <- bind_datasets_together(preprocessed_adni_and_kth) %>%
  dplyr::mutate(Source = case_when(Source == 1 ~ 'adni', Source == 2 ~ 'kth'))
alz_adni_kth_intersection_list <-
  filter_according_to_UniProt_masks(alz_adni_kth, preprocessed_adni_and_kth)
alz_adni_kth_intersection_uniprots <- extract_unique_uniprots(alz_adni_kth_intersection_list, c("Name"))
remove(preprocessed_adni_and_kth)

# Now emif and kth
preprocessed_emif_and_kth <- preprocess_datasets(
  datasets = list(emif, kth),
  should_filter_by_disease = c(FALSE, FALSE),
  additional_columns = c("Name"))
alz_emif_kth <- bind_datasets_together(preprocessed_emif_and_kth) %>%
  dplyr::mutate(Source = case_when(Source == 1 ~ 'emif', Source == 2 ~ 'kth'))
alz_emif_kth_intersection_list <-
  filter_according_to_UniProt_masks(alz_emif_kth, preprocessed_emif_and_kth)
alz_emif_kth_intersection_uniprots <- extract_unique_uniprots(alz_emif_kth_intersection_list, c("Name"))
remove(preprocessed_emif_and_kth)
