#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
## Project: MIRIADE
## Script purpose: Compare the non cleaned datasets
## Date: 13.03.2023
## Author: Felicia Burtscher
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
library(here)
library(dplyr)
source(here("functions", "mapping_functions.R"))
source(here("functions", "comparison_functions.R"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
## Redux comparison ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

compare_datasets_redux <- function(datasets, dataset_names)
{
  datasets_length <- length(datasets)
  if (datasets_length != length(dataset_names)) {
    stop("datasets and dataset_names must be the same length")
  }

  bound_datasets <- bind_datasets_together(datasets) %>%
    dplyr::mutate(Source = case_source(Source, dataset_names))

  unique_proteins <- bound_datasets %>%
    dplyr::distinct(HGNC_Symbol, .keep_all = TRUE) %>%
    dplyr::select(HGNC_Symbol)

  return(list(
    dataset_names = dataset_names,
    raw = bound_datasets,
    unique_proteins = unique_proteins))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
## Perform comparisons ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Compare olink and kth for AZ
alz_olink_and_kth <- compare_datasets_redux(
  datasets = list(olink_control_vs_ad_pvals, kth_control_vs_ad_pvals),
  dataset_names = c('olink', 'kth'))

# Now do the same thing we did with olink vs. kth with the rest
alz_olink_and_adni <- compare_datasets_redux(
  datasets = list(olink_control_vs_ad_pvals, adni_control_vs_ad),
  dataset_names = c('olink', 'adni'))

alz_olink_and_emif <- compare_datasets_redux(
  datasets = list(olink_control_vs_ad_pvals, emif_control_vs_ad_pvals),
  dataset_names = c('olink', 'emif'))

alz_adni_and_emif <- compare_datasets_redux(
  datasets = list(adni_control_vs_ad, emif_control_vs_ad_pvals),
  dataset_names = c('adni', 'emif'))

# All 4 ALZ datasets
alz_olink_adni_emif_kth <- compare_datasets_redux(
  datasets = list(olink_control_vs_ad_pvals, adni_control_vs_ad, emif_control_vs_ad_pvals, kth_control_vs_ad_pvals),
  dataset_names = c('olink', 'adni', 'emif', 'kth'))

# Without KTH
alz_olink_adni_emif <- compare_datasets_redux(
  datasets = list(olink_control_vs_ad_pvals, adni_control_vs_ad, emif_control_vs_ad_pvals),
  dataset_names = c('olink', 'adni', 'emif'))

# Now Olink and mspec
# dlb_olink_mspec <- compare_datasets_redux(
#   datasets = list(vumc_ol, mspec),
#   dataset_names = c('olink', 'mspec'),
#   should_filter_by_disease = c(TRUE, TRUE),
#   additional_columns = c("Name"),
#   disease_to_filter_by = "DLB")
# 
# ftd_olink_mspec <- compare_datasets_redux(
#   datasets = list(vumc_ol, mspec),
#   dataset_names = c('olink', 'mspec'),
#   should_filter_by_disease = c(TRUE, TRUE),
#   additional_columns = c("Name"),
#   disease_to_filter_by = "FTD")

# Now adni and kth
alz_adni_kth <- compare_datasets_redux(
  datasets = list(adni_control_vs_ad, kth_control_vs_ad_pvals),
  dataset_names = c('adni', 'kth'))

# Now emif and kth
alz_emif_kth <- compare_datasets_redux(
  datasets = list(emif_control_vs_ad_pvals, kth_control_vs_ad_pvals),
  dataset_names = c('emif', 'kth'))

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
           dplyr::select(HGNC_Symbol))
  names(result) <- dataset_names

  return(result)
}

# For the colors, we want: Yellow for olink, red for kth, blue for adni,
# green for emif, and purple for mspec
alz_olink_and_kth$venn_diagram <-
  ggvenn(generate_venn_data(alz_olink_and_kth),
         fill_color = c("yellow", "red")) +
  labs(title = "Olink and KTH overlap", subtitle = "AD")
alz_olink_and_adni$venn_diagram <-
  ggvenn(generate_venn_data(alz_olink_and_adni),
         fill_color = c("yellow", "blue")) +
  labs(title = "Olink and Adni overlap", subtitle = "AD")
alz_olink_and_emif$venn_diagram <-
  ggvenn(generate_venn_data(alz_olink_and_emif),
         fill_color = c("yellow", "green")) +
  labs(title = "Olink and Emif overlap", subtitle = "AD")
alz_adni_and_emif$venn_diagram <-
  ggvenn(generate_venn_data(alz_adni_and_emif),
         fill_color = c("blue", "green")) +
  labs(title = "Adni and Emif overlap", subtitle = "AD")
alz_olink_adni_emif_kth$venn_diagram <-
  ggvenn(generate_venn_data(alz_olink_adni_emif_kth),
         fill_color = c("yellow", "blue", "green", "red")) +
  labs(title = "Olink, Adni, Emif and KTH overlap", subtitle = "AD")
alz_olink_adni_emif$venn_diagram <-
  ggvenn(generate_venn_data(alz_olink_adni_emif),
         fill_color = c("yellow", "blue", "green")) +
  labs(title = "Olink, Adni and Emif overlap", subtitle = "AD")
# dlb_olink_mspec$venn_diagram <-
#   ggvenn(generate_venn_data(dlb_olink_mspec),
#          fill_color = c("yellow", "purple")) +
#   labs(title = "Olink and Mspec overlap", subtitle = "DLB")
# ftd_olink_mspec$venn_diagram <-
#   ggvenn(generate_venn_data(ftd_olink_mspec),
#          fill_color = c("yellow", "purple")) +
#   labs(title = "Olink and Mspec overlap", subtitle = "FTD")
alz_adni_kth$venn_diagram <-
  ggvenn(generate_venn_data(alz_adni_kth),
         fill_color = c("blue", "red")) +
  labs(title = "Adni and KTH overlap", subtitle = "AD")
alz_emif_kth$venn_diagram <-
  ggvenn(generate_venn_data(alz_emif_kth),
         fill_color = c("green", "red")) +
  labs(title = "Emif and KTH overlap", subtitle = "AD")
# showing the plots
alz_olink_and_kth$venn_diagram
alz_olink_and_adni$venn_diagram
alz_olink_and_emif$venn_diagram
alz_adni_and_emif$venn_diagram
alz_olink_adni_emif_kth$venn_diagram
alz_olink_adni_emif$venn_diagram
# dlb_olink_mspec$venn_diagram
# ftd_olink_mspec$venn_diagram
alz_adni_kth$venn_diagram
alz_emif_kth$venn_diagram
