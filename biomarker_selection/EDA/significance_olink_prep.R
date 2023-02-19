library(here)
source(here("functions", "data_processing_functions.R"))

datasets_root_directory <- define_datasets_root()

### Read the VUMC Olink biomarkers
vumc_ol <- read.table(file.path(datasets_root_directory, "Olink/20201216_OLINK_Data/cleaned_VUMC_olink.tsv"),
                      sep = "\t", header = T,
                      quote = "") # 'quote' needed due to a special sign in 'Name'

# 3D p-value space of the 3 diseases
# points corresponding to proteins should be labeled
# missing values lead to skipping that protein
# Now we want to pivot so that each unique protein will have 3 values for p values of each disease
# we would also like to have a column that will help color the plot so that proteins significant for
# all diseases and proteins significant to only one will be easily identified
significance_limit <- 0.05 # below ~0.00002322709 significance limit we lose the significant for all category
combined_3d_p_value_per_protein <- tidyr::pivot_wider(vumc_ol, names_from = "disease", values_from = "p.val") %>%
  dplyr::mutate(Significance = case_when(
    (ALZ <= significance_limit & DLB <= significance_limit & FTD <= significance_limit) ~ "Significant for all",
    (ALZ <= significance_limit & DLB > significance_limit & FTD > significance_limit) ~ "Significant for ALZ",
    (ALZ > significance_limit & DLB <= significance_limit & FTD > significance_limit) ~ "Significant for DLB",
    (ALZ > significance_limit & DLB > significance_limit & FTD <= significance_limit) ~ "Significant for FTD",
    TRUE ~ "Not interesting"
  ))

## TODO: Possibly add enum implementation from https://stackoverflow.com/questions/33838392/enum-like-arguments-in-r for significance_type
# Filters only the rows with the Significance matching significance_type and leaves only the UniProt column
extract_uniprots_according_to_significance_type <- function(protein_pvalues, significance_type) {
  protein_pvalues %>%
    dplyr::filter(Significance == significance_type) %>%
    dplyr::select(UniProt)
}

# A list of UniProts significant for all 3 diseases
biomarkers_for_all_dementias <-
  extract_uniprots_according_to_significance_type(combined_3d_p_value_per_protein,
                                                  "Significant for all")
# A list of UniProts significant for ALZ
biomarkers_for_ALZ <-
  extract_uniprots_according_to_significance_type(combined_3d_p_value_per_protein,
                                                  "Significant for ALZ")
# A list of UniProts significant for DLB
biomarkers_for_DLB <-
  extract_uniprots_according_to_significance_type(combined_3d_p_value_per_protein,
                                                  "Significant for DLB")

# A list of UniProts significant for FTD
biomarkers_for_FTD <-
  extract_uniprots_according_to_significance_type(combined_3d_p_value_per_protein,
                                                  "Significant for FTD")
