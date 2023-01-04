###########
# Changes the colors of the edges of the graph
#
# color_sequence should have at least 1 more element than vertices_by_types
#
# @return the colored graph
###########
color_edges_based_on_vertices <- function(graph, vertices_by_types, color_sequence) {
  length <- length(vertices_by_types)
  logical_indexing_lists <- generate_logical_indexing_lists(graph, vertices_by_types)
  
  for (i in 1:(length + 1)) {
    source_vertex_group <- V(graph)[logical_indexing_lists[[i]]]
    for (j in i:(length + 1)) {
      target_vertex_group <- V(graph)[logical_indexing_lists[[j]]]
      E(graph)[source_vertex_group %--% target_vertex_group]$color <- color_sequence[[(i - 1) * (length + 1) + j]]
    }
  }
  return(graph)
}

#######################################################
# AUXILIARY FUNCTION FOR color_edges_based_on_vertices#
#######################################################
generate_logical_indexing_lists <- function(graph, vertices_by_types) {
  length <- length(vertices_by_types)
  logical_indexing <- list()
  # Additionally create a new grouping of vertices that are in non of the types (possibly empty)
  no_type_logical_indexing <- TRUE # To initialize it so the first loop iteration will return !(V(graph)$name %in% vertices_by_types[[1]])
  for (i in 1:length) {
    no_type_logical_indexing <- !(V(graph)$name %in% vertices_by_types[[i]]) & no_type_logical_indexing
    current_type_indexing <- V(graph)$name %in% vertices_by_types[[i]]
    logical_indexing[[i]] <- current_type_indexing
  }
  append(logical_indexing, list(no_type_logical_indexing), after = 0)
}

###########
# Changes the colors of the edges of the graph based on their source
#
# @param graph
# @param sources the list of sources. In order to properly color edges that come
#.       from more than one source use the expanded source list based on
# @param color_sequence should have at least 1 more element than sources
#
# @return the colored graph
###########
color_edges_based_on_sources <- function(graph, sources, color_sequence) {
  length <- length(sources)
  
  for (i in 1:length) {
    E(graph)[E(graph)$source == sources[[i]]]$color <- color_sequence[[i + 1]]
  }
  return(graph)
}

###########
# Adds the source attribute to the graph's edged.
# This function is a thin wrapper to allow using the '%>%' operator
# to perform this operation
#
# @return the graph with the source added to the edges
###########
add_source_to_edges <- function(graph, source) {
  E(graph)$source <- source
  return(graph)
}

###########
# This function takes a union of two graphs with a source attribute in their
# edges, and fixes the unified graph to have only one (correct) source
# attribute in the edges.
# This function also creates a union source if edges appear in both sources
# under the assumption that all sources are unique
#
# @return the graph with the edge sources converged into
#         one source attribute and the source1&source2 attributes removed
###########
converge_edge_sources <- function(graph) {
  exists_in_source1 <- !is.na(E(graph)$source_1)
  exists_in_source2 <- !is.na(E(graph)$source_2)
  only_in_source_1 <- exists_in_source1 & !exists_in_source2
  only_in_source_2 <- exists_in_source2 & !exists_in_source1
  in_both <- exists_in_source1 & exists_in_source2
  E(graph)[only_in_source_1]$source <-
    E(graph)[only_in_source_1]$source_1
  E(graph)[only_in_source_2]$source <-
    E(graph)[only_in_source_2]$source_2
  E(graph)[in_both]$source <-
    paste(E(graph)[in_both]$source_1, E(graph)[in_both]$source_2, sep = '+')
  return(delete_edge_attr(graph, "source_1") %>%
           delete_edge_attr("source_2"))
}