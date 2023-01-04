###########
# Adjusts the attributes (color, size and type) of the vertices according to their type
#
# color_sequence should have at least 1 more element than vertices_by_types
#
# @return the colored and resized graph
###########
adjust_vertices_attributes_according_to_type <- function(graph, vertices_by_types, color_sequence) {
  length <- length(vertices_by_types)
  # First the color of vertices that are somehow not associated with anything
  V(graph)$color <- color_sequence[[1]]
  V(graph)$size <- 15
  V(graph)$type <- "None"
  for (i in 1:length) {
    V(graph)[V(graph)$name %in% vertices_by_types[[i]]]$color <- color_sequence[[i + 1]]
    V(graph)[V(graph)$name %in% vertices_by_types[[i]]]$size <- 22
    V(graph)[V(graph)$name %in% vertices_by_types[[i]]]$type <- names(vertices_by_types)[[i]]
  }
  return(graph)
}

##############################################################
# AUX function for the biomarker network graph prepatation
# generating a list of sets of combinations of original types
##############################################################
generate_all_vertex_type_combinations <- function(vertices_by_types) {
  if(!require('sets')) {
    install.packages('sets')
    library('sets')
  }
  length <- length(vertices_by_types)
  # Generate the permutations of which sets to take from
  permutations <- expand.grid(rep(list(0:1), length(vertices_by_types)))
  number_of_permutations <- nrow(permutations)
  vertex_sets <- list()
  all_used_vertices <- set()
  for (perm in 2:number_of_permutations) {
    # Start from the permutation containing all sets to be able to preemptively trim duplicates
    permutation <- number_of_permutations - perm + 2
    temp_set <- set()
    key <- ""
    for (i in 1:length) {
      if(permutations[[i]][[permutation]]) {
        key <- expand_key(key, vertices_by_types, i)
        temp_set <- expand_vertex_set(temp_set, vertices_by_types, i)
      }
    }
    # Remove from temp_set vertices that are already included in one of the sets
    temp_set <- temp_set - (temp_set & all_used_vertices)
    # Record in all_used_vertices the genes that are now in this set so they will not be used again in another set
    all_used_vertices <- all_used_vertices | temp_set
    vertex_sets[key] <- set(temp_set)
  }
  # In order to get the magrittr pipe back
  library(magrittr)
  return(vertex_sets)
}

### AUX function
expand_key <- function(key, vertices_by_types, index) {
  if(nchar(key)) {
    key <- paste(key, "+", sep="")
  }
  return(paste(key, names(vertices_by_types)[[index]], sep=""))
}

### AUX function
expand_vertex_set <- function(temp_set, vertices_by_types, index) {
  if(set_is_empty(temp_set)) {
    temp_set <- as.set(vertices_by_types[[index]])
  }
  else {
    # The intersection of what's in the set and the next to be included vertex list
    temp_set <- temp_set & as.set(vertices_by_types[[index]])
  }
  return(temp_set)
}