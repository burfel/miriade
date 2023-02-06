library(here)
source(here("biomarker_selection", "EDA", "significance_olink_prep.R"))

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
