library(here)
source(here("biomarker_selection", "EDA", "compare_datasets.R"))

# Do the ALZ plot with KTH
pval_vs_uniprot_ALZ_kth <- ggplot(alz_olink_and_kth
                                  , aes(x = UniProt, y = p.val, color = Source)) +
  geom_line() + geom_point() +
  labs(title = "P values for ALZ including kth")
pval_vs_uniprot_ALZ_kth

pval_vs_uniprot_ALZ_kth_significant <- ggplot(alz_olink_and_kth %>%
                                                dplyr::filter(p.val <= 0.05)
                                              , aes(x = UniProt, y = p.val, color = Source)) +
  geom_vline(xintercept = kth$UniProt[which(kth$UniProt %in% vumc_ol$UniProt)]) + geom_point() +
  labs(title = "P values for ALZ for proteins significant for ALZ in olink and kth")
pval_vs_uniprot_ALZ_kth_significant

pval_vs_uniprot_olink_kth_intersection <- ggplot(alz_olink_and_kth_intersection_list
                                                 , aes(x = UniProt, y = p.val, color = Source)) + geom_point() +
  labs(title = "P values for ALZ for proteins for ALZ both in olink and kth")
pval_vs_uniprot_olink_kth_intersection

# This function takes a dataframe with uniprots and p values from more than one source
# and filters according to the list of limiting dataframes. Also, takes a title
create_plot_of_significant_according_to_limiting_lists <- function(data_frame_to_be_plotted, masking_frames_vector, title) {
  return(ggplot(filter_according_to_UniProt_masks(data_frame_to_be_plotted, masking_frames_vector),
                aes(x = UniProt, y = p.val, color = Source)) + geom_point() +
           labs(title = title))
}

pval_vs_uniprot_olink_kth_intersection_ALZ_significant <- create_plot_of_significant_according_to_limiting_lists(
  alz_olink_and_kth_intersection_list, list(biomarkers_for_ALZ), "P values for ALZ for significnat proteins for ALZ both in olink and kth"
)
pval_vs_uniprot_olink_kth_intersection_ALZ_significant

pval_vs_uniprot_olink_adni_intersection_ALZ_significant <- create_plot_of_significant_according_to_limiting_lists(
  alz_olink_and_adni_intersection_list, list(biomarkers_for_ALZ, alz_olink_and_adni_intersection_significant_in_both_datasets),
  "P values for ALZ for significant proteins for ALZ both in olink and adni"
)
pval_vs_uniprot_olink_adni_intersection_ALZ_significant

pval_vs_uniprot_olink_emif_intersection_ALZ_significant <- create_plot_of_significant_according_to_limiting_lists(
  alz_olink_and_emif_intersection_list, list(biomarkers_for_ALZ, alz_olink_and_emif_intersection_significant_in_both_datasets),
  "P values for ALZ for significant proteins for ALZ both in olink and emif"
)
pval_vs_uniprot_olink_emif_intersection_ALZ_significant

pval_vs_uniprot_adni_emif_intersection_ALZ_significant <- create_plot_of_significant_according_to_limiting_lists(
  alz_adni_and_emif_intersection_list, list(biomarkers_for_ALZ, alz_adni_and_emif_intersection_significant_in_both_datasets),
  "P values for ALZ for significant proteins for ALZ both in adni and emif"
)
pval_vs_uniprot_adni_emif_intersection_ALZ_significant
