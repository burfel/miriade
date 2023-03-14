#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
## Project: MIRIADE
## Script purpose: Perform meta analysis on all datasets
## Date: 14.03.2023
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
library(here)
library(dplyr)
# The following are sourced by sourcing `prep_datasets_for_comparisons.R`
# `prep_datasets_for_kw_wilcoxon.R`, `auxiliary_statistical_functions.R`,
# `meld_all_fractured_datasets.R`, `data_processing_functions.R`
# and `dataset_root_directory` is defined
source(here("biomarker_selection", "EDA", "prep_datasets_for_comparisons.R"))

### First, combine p values across all AD datasets, using HGNC_Symbol as the index column
common_cols <- c("HGNC_Symbol", "p.val", "adj.p")

### Pool all the entries together, add source information
combined_alz <-
  do.call(rbind,
          list(cbind(olink_control_vs_ad_pvals[,common_cols], source = "Olink"),
               cbind(kth_control_vs_ad_pvals[,common_cols], source = "KTH"),
               cbind(adni_control_vs_ad[,common_cols], source = "ADNI"),
               cbind(emif_control_vs_ad_pvals[,common_cols], source = "EMIF")))

if(!require("mppa")) {
  library(remotes)
  install_version("mppa","1.0")
}
library("mppa")
library("poolr")

meta_alz <- dplyr::group_by(combined_alz, HGNC_Symbol) %>%
  dplyr::summarise(p.val_fisher = poolr::fisher(p.val)$p,
                   p.val_simes = mppa::simes.test(p.val),
                   p.val_stouffer = poolr::stouffer(p.val)$p,
                   p.val_invchisq = poolr::invchisq(p.val)$p,
                   p.val_binomtest = poolr::binomtest(p.val)$p,
                   p.val_bonferroni = poolr::bonferroni(p.val)$p,
                   p.val_tippett = poolr::tippett(p.val)$p,
                   sources = paste(sort(source), collapse = ",")) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(adj.p.val_fisher = p.adjust(p.val_fisher, method = "BY"), 
                adj.p.val_simes = p.adjust(p.val_simes, method = "BY"),
                adj.p.val_stouffer = p.adjust(p.val_stouffer, method = "BY"),
                adj.p.val_invchisq = p.adjust(p.val_invchisq, method = "BY"),
                adj.p.val_binomtest = p.adjust(p.val_binomtest, method = "BY"),
                adj.p.val_bonferroni = p.adjust(p.val_bonferroni, method = "BY"),
                adj.p.val_tippett = p.adjust(p.val_tippett, method = "BY")) %>%
  dplyr::select(-sources, everything(), sources)
