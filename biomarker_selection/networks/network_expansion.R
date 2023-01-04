library(dplyr)
biomarker_network_functions_path <- "/Users/felicia.burtscher/Documents/UL/GITHUB/miriade/functions/biomarker-network_functions.R"
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
  dplyr::select(`Network Object \"FROM\"`, `Network Object \"TO\"`) %>%
  dplyr::distinct()

miriade_uniprot_interactions <- translate_interactions_using_dictionary(
        interactions_df = miriade_interactions, 
        interaction_columns = c("Network Object \"FROM\"", "Network Object \"TO\""),
        dictionary_df = miriade_genes,
        dictionary_columns = c("Network Object Name", "SwissProt IDs"),
        new_interaction_names = c("UniProt_FROM", "UniProt_TO"))
# Produce a list of all UniProts in the file
uniprots_miriade <- extract_all_distinct_objects(miriade_uniprot_interactions,
                                     c("UniProt_TO", "UniProt_FROM"))
# Produce a list of all network objects
network_objects_miriade <- extract_all_distinct_objects(miriade_interactions,
                                    c("Network Object \"FROM\"", "Network Object \"TO\""))
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
