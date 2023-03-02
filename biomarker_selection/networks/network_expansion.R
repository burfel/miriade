library(dplyr)
library(here)
source(here("functions", "biomarker-network_functions.R"))
source(here("functions", "data_processing_functions.R"))
from_metacore_path <- here("biomarker_selection", "networks", "from-metacore")

miriade_interactions <- get_all_interactions_that_have_backing_uniprots(
  dataset_directory = from_metacore_path,
  interactions_file = "all-MIRIADE-biomarker-candidates_analyse-networks_interactions.xls",
  dictionary_file = "all-MIRIADE-biomarker-candidates_analyse-networks_genes.xls",
  return_as_uniprots = T
)

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
# vertices_by_types <- generate_all_vertex_type_combinations(list("miriade" = uniprots_miriade))
# color_seq <- generate_coloring_scheme_for_graph(vertices_by_types)
# ig <- graph_from_edge_df_filtered_by_genes(
#   miriade_uniprot_interactions[dqsample(nrow(miriade_uniprot_interactions), 20),],
#                                            source_col_index = 1,
#                                            target_col_index = 2) %>%
#   prepare_graph_attributes_according_to_vertices(vertices_by_types, color_seq)
#
# draw_igraph(ig, extract_type_list_from_vertices_by_types(vertices_by_types), color_seq)

# 09DEC2022 - write the graph to file
# write_graph(ig, "/Users/felicia.burtscher/Documents/UL/GITHUB/miriade/disease_map/test-metacore-network.graphml", format = c("graphml"))

###########################################################################
# Taking the csf metacore interactions and combining with the miriade ones
###########################################################################
csf_interactions <- readxl::read_xls(here(from_metacore_path,
                                           "pooled-CSF-proteins_interactions.xls"),
                                           skip = 2) %>%
  dplyr::select(`Network Object \"FROM\"`, `Network Object \"TO\"`) %>%
  dplyr::distinct() # 21 NAs
# Produce a list of all network objects for csf
network_objects_csf <- extract_all_distinct_objects(csf_interactions,
                                                        c("Network Object \"FROM\"", "Network Object \"TO\""))
sources <- c("miriade", "csf")
# Prepare the list of dataframes
interaction_dfs <- list(miriade_interactions, csf_interactions)
# prepare the simple named list of vertices by sources
vertices <- list(network_objects_miriade, network_objects_csf)
names(vertices) <- sources
names(interaction_dfs) <- sources # Not necessary, but helpful for inspection and debug
# prepare the extended list of vertices by source
vertices_by_source <- generate_all_vertex_type_combinations(vertices)
color_sequence <- generate_color_sequence(length(vertices_by_source) + 1)

graph_list = list()
for (i in 1:2) {
  graph_list[[i]] <- graph_from_edge_df_filtered_by_genes(
    interaction_dfs[[i]],
    source_col_index = 1,
    target_col_index = 2) %>%
    add_source_to_edges(sources[[i]])
}

extended_source_list <- extract_type_list_from_vertices_by_types(vertices_by_source)

unified_graph <- unify_graphs(graph_list, vertices_by_source,
                              extended_source_list[-1], color_sequence)

write_graph(unified_graph,
            here("biomarker_selection", "networks", "/miriade_csf.graphml"),
            format = "graphml")

data_framed_graph <- igraph::as_data_frame(unified_graph)