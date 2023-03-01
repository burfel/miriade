library(dplyr)
library(here)
source(here("functions", "biomarker-network_functions.R"))
source(here("functions", "data_processing_functions.R"))
from_metacore_path <- here("biomarker_selection", "networks", "from-metacore")

# Whether to work with uniprots or network names
work_with_uniprots <- FALSE

miriade_interactions <- get_all_interactions_that_have_backing_uniprots(
  dataset_directory = from_metacore_path,
  interactions_file = "all-MIRIADE-biomarker-candidates_analyse-networks_interactions.xls",
  dictionary_file = "all-MIRIADE-biomarker-candidates_analyse-networks_genes.xls",
  return_as_uniprots = work_with_uniprots
) # 5 NAs


###########################################################################
# Taking the csf metacore interactions and combining with the miriade ones
###########################################################################
csf_interactions <- get_all_interactions_that_have_backing_uniprots(
  dataset_directory = from_metacore_path,
  interactions_file = "pooled-CSF-proteins_network_new_interactions.xls",
  dictionary_file = "pooled-CSF-proteins_network_new_genes.xls",
  return_as_uniprots = work_with_uniprots
) # 21 NAs

sources <- c("miriade", "csf")

# Prepare the list of dataframes
interaction_dfs <- named_list(names_vector = sources,
                              miriade_interactions, csf_interactions)

# prepare the simple named list of vertices by sources
vertices <- extract_vertices_for_all_graphs(interaction_dfs,
                                            return_as_uniprots = work_with_uniprots)
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

draw_igraph(unified_graph, extended_source_list, color_sequence)
write_graph(unified_graph,
            here("biomarker_selection", "networks", "/miriade_csf.graphml"),
            format = "graphml")
# data_framed_graph <- igraph::as_data_frame(unified_graph)

########################
# Community detection  #
######################### Random walk does only 4 step by default, resulting in 662 communities.
# With 100 steps, 371 communities are detected
# random_walk_communities <- igraph::walktrap.community(unified_graph, steps = 100)
# edge_betweenness_communities <- edge.betweenness.community(unified_graph)
# fast_greedy_communities <- fastgreedy.community(unified_graph) # ONLY FOR Undirected Graphs
# DO NOT USE - Has run for several dozen minutes with no result
# label_propagation_communities <- label.propagation.community(unified_graph)
# DO NOT USE - The optimal method is causing some memory leak
# optimal_communities <- igraph::optimal.community(unified_graph)
# infomap_communities <- igraph:: infomap.community(unified_graph)
# Now we take the induced subgraph of only the intersection vertices
vertex_intersection_graph <- igraph::induced_subgraph(unified_graph, V(unified_graph)[V(unified_graph)$type == "miriade+csf"])
v_random_walk_communities <- igraph::walktrap.community(vertex_intersection_graph, steps = 100)
# v_edge_betweenness_communities <- edge.betweenness.community(vertex_intersection_graph)
# v_infomap_communities <- igraph:: infomap.community(vertex_intersection_graph)

enriched_communities <-
  Filter(is_df_populated,
         sapply(communities(v_random_walk_communities), community_enrichment))
# Aggregate according to Gene so that all of its terms are concatenated with ';'
aggregated_enriched_communities <-
  lapply(enriched_communities, function(df) df %>% 
  group_by(Gene) %>% 
  summarise(Term = paste(Term, collapse = ";")))
