###########
# This auxiliary function splits an interaction dataframe into randomly
# sampled parts, annotates each part as a source and then performs a union
# on the parts and possibly draws it
#
# @param source_df A dataframe convertible to a graph via 
#        graph_from_edge_df_filtered_by_genes
# @param source_and_target_columns the indices and names of the 2 columns
#        representing the source and target of the edges
# @param number_of_parts the number of parts sampled from the whole df
#        2 by default
# @param size_of_parts the size of parts randomly sampled from the whole df
#        20 by default
# @param should_draw_graph whether or not it should plot the graph with
#        draw_igraph. FALSE by default
#
# @return the unified_graph
###########
split_dataframe_and_rejoin_resulting_graphs <-
  function(source_df, source_and_target_columns,
           number_of_parts = 2, size_of_parts = 20,
           should_draw_graph = FALSE) {
    require(dqrng)
    test_dfs <- list()
    test_sources <- vector(mode = "character")
    test_vertices <- list()
    for (i in 1:number_of_parts) {
      test_dfs[[i]] <- source_df[dqsample(nrow(source_df), size_of_parts),]
      test_sources[[i]] <- paste("test_df_", i, sep="")
      test_vertices[[i]] <- extract_all_distinct_objects(test_dfs[[i]],
                                      source_and_target_columns$name)
    }

    names(test_dfs) <- test_sources
    names(test_vertices) <- test_sources
    
    vertices_by_test <- generate_all_vertex_type_combinations(test_vertices)
    v_color_seq <- generate_color_sequence(length(vertices_by_test) + 1)

    graph_list = list()
    for (i in 1:number_of_parts) {
      graph_list[[i]] <- graph_from_edge_df_filtered_by_genes(
        test_dfs[[i]],
        source_col_index = source_and_target_columns$index[[1]],
        target_col_index = source_and_target_columns$index[[2]]) %>%
        add_source_to_edges(test_sources[[i]])
    }

    extended_source_list <- extract_type_list_from_vertices_by_types(vertices_by_test)

    unified_graph <- unify_graphs(graph_list, vertices_by_test,
                                  extended_source_list[-1], v_color_seq)
    if (should_draw_graph) {
      draw_igraph(unified_graph, extended_source_list, v_color_seq)
    }
    return(unified_graph)
  }

###########
# This test tests whether the graph functions from 'biomarker-network_functions'
# do their job correctly when performing a union.
# It optionally draws the graph
# It returns a dataframed version of the graph for easier inspection
###########
should_produce_properly_colored_and_typed_graph_when_doing_union <-
  function(number_of_parts = 2, size_of_parts = 20,
           should_draw_graph = FALSE) {
from_metacore_path <- "/Users/felicia.burtscher/Documents/UL/GITHUB/miriade/biomarker_selection/networks/from-metacore"
miriade_interactions <- readxl::read_xls(paste(from_metacore_path, 
                                               "all-MIRIADE-biomarker-candidates_analyse-networks_interactions.xls", 
                                               sep="/"), skip = 2) %>%
  dplyr::select(`Network Object \"FROM\"`, `Network Object \"TO\"`) %>%
  dplyr::distinct()
source_and_target_columns <- data.frame(index = c(1,2), name = c("Network Object \"FROM\"", "Network Object \"TO\""))
test_unified_graph <- split_dataframe_and_rejoin_resulting_graphs(miriade_interactions,
                                                                 source_and_target_columns,
                                                                 number_of_parts,
                                                                 size_of_parts,
                                                                 should_draw_graph)
return(as_data_frame(test_unified_graph))
}

###########################################
# Default driver code for graph union test#
###########################################
library(dplyr)
biomarker_network_functions_path <- "/Users/felicia.burtscher/Documents/UL/GITHUB/miriade/functions/biomarker-network_functions.R"
source(biomarker_network_functions_path)
data_processing_functions_path <- "/Users/felicia.burtscher/Documents/UL/GITHUB/miriade/functions/data_processing_functions.R"
source(data_processing_functions_path)
data_framed_graph <- should_produce_properly_colored_and_typed_graph_when_doing_union()
