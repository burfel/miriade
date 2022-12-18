library(dplyr)
biomarker_network_functions_path <- "/Users/felicia.burtscher/Documents/UL/GITHUB/miriade/disease_map/biomarker-network_functions.R"
source(biomarker_network_functions_path)
from_metacore_path <- "/Users/felicia.burtscher/Documents/UL/GITHUB/miriade/biomarker_selection/networks/from-metacore"
miriade_genes <- readxl::read_xls(paste(from_metacore_path, 
                                        "all-MIRIADE-biomarker-candidates_analyse-networks_genes.xls", 
                                        sep="/"), skip = 2) %>%
  dplyr::select(`SwissProt IDs`, `Network Object Name`) %>%
  # Remove useless "No ID" rows
  dplyr::filter(`SwissProt IDs`!="No ID") %>%
  dplyr::distinct(`Network Object Name`, .keep_all = TRUE)

miriade_interactions <- readxl::read_xls(paste(from_metacore_path, 
                                        "all-MIRIADE-biomarker-candidates_analyse-networks_interactions.xls", 
                                        sep="/"), skip = 2) %>%
  dplyr::select(`Network Object \"FROM\"`, `Network Object \"TO\"`)

miriade_uniprot_interactions <- dplyr::inner_join(miriade_interactions, miriade_genes, 
                  by = c("Network Object \"FROM\"" = "Network Object Name")) %>%
  dplyr::rename(UniProt_FROM = `SwissProt IDs`) %>%
  dplyr::inner_join(miriade_genes, 
                   by = c("Network Object \"TO\"" = "Network Object Name")) %>%
  dplyr::rename(UniProt_TO = `SwissProt IDs`) %>%
  dplyr::select(UniProt_FROM, UniProt_TO) %>%
  # Separate the ';' delimited lists of UniProts into separate rows TWICE:
  # once for the FROM column, and once for the TO column
  tidyr::separate_rows(UniProt_FROM,sep=";") %>% 
  tidyr::separate_rows(UniProt_TO,sep=";") %>%
  dplyr::distinct()
# Produce a list of all UniProts in the file
uniprots_miriade <- union(
  dplyr::distinct(miriade_uniprot_interactions, UniProt_TO)[["UniProt_TO"]],
  dplyr::distinct(miriade_uniprot_interactions, UniProt_FROM)[["UniProt_FROM"]]
  )
# DANGER! - Do not run locally. Too much for the poor machine
# generate the igraph finally
#ig <- graph_from_edge_df_filtered_by_genes(miriade_uniprot_interactions,
#                                           source_col_index = 1,
#                                           target_col_index = 2)
#draw_igraph(ig, list("miriade" = uniprots_miriade), list("miriade"))

# dqrng for fast sampling
library(dqrng)
# Generate the proper combinations of genes by which type(s) they belong to
vertices_by_types <- generate_all_vertex_type_combinations(list("miriade" = uniprots_miriade))
color_seq <- generate_coloring_scheme_for_graph(vertices_by_types)
ig <- graph_from_edge_df_filtered_by_genes(
  miriade_uniprot_interactions[dqsample(nrow(miriade_uniprot_interactions), 20),],
                                           source_col_index = 1,
                                           target_col_index = 2) %>%
  prepare_graph_attributes_according_to_vertices(vertices_by_types, color_seq)

draw_igraph(ig, extract_type_list_from_vertices_by_types(vertices_by_types), color_seq)

# 09DEC2022 - write the graph to file
write_graph(ig, "/Users/felicia.burtscher/Documents/UL/GITHUB/miriade/disease_map/test-metacore-network.graphml", format = c("graphml"))

#############
## TEST UNIONING GRAPHS
#############
test_interactions <- miriade_uniprot_interactions[dqsample(nrow(miriade_uniprot_interactions), 20),]
test_interactions2 <- miriade_uniprot_interactions[dqsample(nrow(miriade_uniprot_interactions), 20),]
test_source_vector <- c("test_interactions", "test_interactions2")

interactions <- list(test_interactions, test_interactions2)
names(interactions) <- test_source_vector

test_v1 <- unlist(list(test_interactions[[1]], test_interactions[[2]]))
test_v2 <- unlist(list(test_interactions2[[1]], test_interactions2[[2]]))

vertices_by_test <- generate_all_vertex_type_combinations(
  list(
    "test1" = test_v1,
    "test2" = test_v2
    ))
v_color_seq <- generate_color_sequence(length(vertices_by_test) + 1)
graph_list = list()
for (i in 1:length(interactions)) {
  graph_list[[i]] <- graph_from_edge_df_filtered_by_genes(
    interactions[[i]],
    source_col_index = 1,
    target_col_index = 2) %>%
    adjust_vertices_attributes_according_to_type(vertices_by_test, v_color_seq)
}
