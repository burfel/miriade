#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
## Project: MIRIADE
## Script purpose: Refine the preparation from `prep_datasets_for_kw_wilcoxon`
##                 further by performing kw wilcoxon with a significance_limit = 1
##                 to get p values according to control vs disease
## Date: 13.03.2023
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
library(here)
library(dplyr)

# This calls meld_all_fractured_datasets, sets the dataset root directory,
# and prepares kth (/emif/olink), mutated_kth(/emif/olink),
# melted_kth(/emif/olink) and supercharged_kth(/emif/olink)_melt to be used
# in this file 
source(here("biomarker_selection", "EDA", "prep_datasets_for_kw_wilcoxon.R"))

# DISCLAIMER: We will use `significance_limit = 1` to retain all proteins

kth_control_vs_ad_pvals <- perform_kw_wilcoxon_according_to_grouping(
  melted_df = melted_kth,
  differentiating_feature_symbol = sym("HGNC_Symbol"),
  grouping_feature_symbol = sym("Diagnosis"),
  measurement_symbol = sym("Value"),
  significance_limit = 1)[["Control_vs_AD"]]

emif_control_vs_ad_pvals <- perform_kw_wilcoxon_according_to_grouping(
  melted_df = melted_emif,
  differentiating_feature_symbol = sym("HGNC_Symbol"),
  grouping_feature_symbol = sym("Diagnosis"),
  measurement_symbol = sym("Value"),
  significance_limit = 1)[["AD_vs_Control"]]

olink_diagnosis_pvals <- perform_kw_wilcoxon_according_to_grouping(
  melted_df = melted_olink,
  differentiating_feature_symbol = sym("HGNC_Symbol"),
  grouping_feature_symbol = sym("Diagnosis"),
  measurement_symbol = sym("Value"),
  significance_limit = 1)

olink_control_vs_ad_pvals <- olink_diagnosis_pvals[["AD_vs_Control"]]
olink_control_vs_dlb_pvals <- olink_diagnosis_pvals[["DLB_vs_Control"]]
olink_control_vs_ftd_pvals <- olink_diagnosis_pvals[["FTD_vs_Control"]]

### Brain 2020 AD Adni and Emif thing

### Read the Supplementary Materials Excel sheet, Supplementary table 2
### Contains the combined list of preselected ADNI and EMIF entries
### We need to skip the first two lines and select columns 1:6 (metadata) and 13
### (AD vs Controls p values).
### Focus only on those selected for clustering
supp_control_vs_ad <- readxl::read_xlsx(file.path(datasets_root_directory, "Brain 2020 AD/awaa325-suppl_data/brain-2020-00570-File008.xlsx"),
                          skip = 2,
                          sheet = "Supplementary table 2")[, c(1:6, 13)] %>%
  dplyr::filter(`Selected for clustering` == "Yes") %>%
  dplyr::rename(HGNC_Symbol = Gene, UniProt = Uniprot, p.val = `p value...13`) %>%
  dplyr::select(HGNC_Symbol, UniProt, p.val, Cohort) %>%
  ### Convert the p.val column to numeric, "<.001 values become 0.001"
  dplyr::mutate(p.val = case_when(p.val == "<.001" ~ "0.001",
                                        TRUE ~ p.val)) %>%
  dplyr::mutate(p.val = as.numeric(p.val))
# Remove 2 rows that had two uniprots separated by ';'
supp_control_vs_ad <- supp_control_vs_ad[-grep("; ", supp_control_vs_ad$UniProt),]
### Adjust the p values, separately for both cohorts
supp_control_vs_ad <- cbind(supp_control_vs_ad, adj.p = supp_control_vs_ad$p.val)
# TODO: FIX THESE TO ACTUALLY BE MEANINGFUL
supp_control_vs_ad$adj.p[supp_control_vs_ad$Cohort == "ADNI"] <- p.adjust(supp_control_vs_ad$p.val[supp_control_vs_ad$Cohort == "ADNI"])
supp_control_vs_ad$adj.p[supp_control_vs_ad$Cohort == "EMIF-AD MBD"] <- p.adjust(supp_control_vs_ad$p.val[supp_control_vs_ad$Cohort == "EMIF-AD MBD"])

adni_control_vs_ad <- dplyr::filter(supp_control_vs_ad, Cohort == "ADNI")
emif_brain2020_control_vs_ad <- dplyr::filter(supp_control_vs_ad, Cohort == "EMIF-AD MBD")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
### mspec ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
### Read the results Excel file for FTD
ftd_mspec <- readxl::read_xlsx(file.path(datasets_root_directory, "VUMC/WBI_CSF_Candidates_FC1.5_P0.05_DP0.5_FTD_DLB.xlsx"), sheet = 1) %>%
  dplyr::rename(UniProt = Accession, HGNC_Symbol = Gene, p.val = Pval) %>%
  dplyr::select(UniProt, HGNC_Symbol, p.val) %>%
  dplyr::mutate(adj.p = p.adjust(p.val, method = "BH"))

### Read the results Excel file for DLB
dlb_mspec <- readxl::read_xlsx(file.path(datasets_root_directory, "VUMC/WBI_CSF_Candidates_FC1.5_P0.05_DP0.5_FTD_DLB.xlsx"), sheet = 2) %>%
  dplyr::rename(UniProt = Accession, HGNC_Symbol = Gene, p.val = Pval) %>%
  dplyr::select(UniProt, HGNC_Symbol, p.val) %>%
  dplyr::mutate(adj.p = p.adjust(p.val, method = "BH"))
