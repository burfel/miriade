### Read the VUMC Olink biomarkers
vumc_ol <- read.table("/Users/felicia.burtscher/Documents/UL/DATASETS/Olink/20201216_OLINK_Data/cleaned_VUMC_olink.tsv",
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

if(!require("plotly")) {
  install.packages("plotly")
}
library("plotly")
plot_ly(combined_3d_p_value_per_protein, x = ~ALZ, y = ~DLB, z = ~FTD, text = ~UniProt) %>%
  add_markers(color = ~Significance,
              hovertemplate = paste('<b>UniProt</b>: %{text}',
                                    '<br><b>ALZ p value</b>: %{x}',
                                    '<br><b>DLB p value</b>: %{y}',
                                    '<br><b>FTD p value</b>: %{z}'),) %>%
  layout(title = 'P value of proteins for 3 diseases',
         scene = list(
           xaxis = list(title = "ALZ p value"),
           yaxis = list(title = "DLB p value"),
           zaxis = list(title = "FTD p value")
         )
  )

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

# TODO: Extract the ggplot call etc. to a function
# Plot 4: 1-d curve of the 3 diseases
# we give the proteins an ordering, a panel of proteins according to the cleaned_VUMC_olink.tsv dataset
# vumc_ol needs to filter out the non important ones
# vumc_significants <- vumc_ol %>%
#   dplyr::filter(UniProt %in% biomarkers_for_all_dementias$UniProt |
#                   UniProt %in% biomarkers_for_ALZ$UniProt |
#                   UniProt %in% biomarkers_for_DLB$UniProt |
#                   UniProt %in% biomarkers_for_FTD$UniProt)
# pval_vs_uniprot <- ggplot(vumc_significants, aes(x = UniProt, y = p.val, group = disease, color = disease)) +
#                             geom_line() + geom_point()
# pval_vs_uniprot
pval_vs_uniprot_significant <- ggplot(vumc_ol %>%
                                        dplyr::filter(UniProt %in% biomarkers_for_all_dementias$UniProt)
                                      , aes(x = UniProt, y = p.val, group = disease, color = disease)) +
  geom_line() + geom_point() +
  labs(title = "P values for proteins significant for all 3 diseases")
pval_vs_uniprot_significant

pval_vs_uniprot_ALZ_significant <- ggplot(vumc_ol %>%
                                            dplyr::filter(UniProt %in% biomarkers_for_ALZ$UniProt &
                                                            disease == "ALZ")
                                          , aes(x = UniProt, y = p.val, group = disease)) +
  geom_line() + geom_point() +
  labs(title = "P values for ALZ for proteins significant for ALZ")
pval_vs_uniprot_ALZ_significant

pval_vs_uniprot_DLB_significant <- ggplot(vumc_ol %>%
                                            dplyr::filter(UniProt %in% biomarkers_for_DLB$UniProt &
                                                            disease == "DLB")
                                          , aes(x = UniProt, y = p.val, group = disease)) +
  geom_line() + geom_point() +
  labs(title = "P values for DLB for proteins significant for DLB")
pval_vs_uniprot_DLB_significant

pval_vs_uniprot_FTD_significant <- ggplot(vumc_ol %>%
                                            dplyr::filter(UniProt %in% biomarkers_for_FTD$UniProt &
                                                            disease == "FTD")
                                          , aes(x = UniProt, y = p.val, group = disease)) +
  geom_line() + geom_point() +
  labs(title = "P values for FTD for proteins significant for FTD")
pval_vs_uniprot_FTD_significant
