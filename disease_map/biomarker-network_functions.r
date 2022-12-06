################################################################################
# Takes an edge dataframe and extracts from it a graph with edges where at
# least one node is in the vertex_list
#
# @param edge_df is the dataframe with the edge data
# @param vertex_list NULL if no filtering desired or a list of vertices otherwise
# @param source_col_index column number for the source
# @param target_col_index column number for the target
#
# @return an igraph graph with edges containing the vertices
################################################################################
graph_from_edge_df_filtered_by_genes <- function(edge_df, vertex_list = NULL, 
                                                 source_col_index,
                                                 target_col_index) {
  require(igraph)
  if(is.null(vertex_list)) {
    filtered_edge_df <- edge_df
  }
  else {
    filtered_edge_df <- edge_df[edge_df[,source_col_index] %in% vertex_list 
                                | edge_df[,target_col_index] %in% vertex_list,]
  }
  graph <- graph_from_edgelist(as.matrix(
      filtered_edge_df[,c(source_col_index,target_col_index)]))
  return(graph)
}

#
# Generates the color sequence for the graph with those vertices
# DO NOT USE FOR NOW
# REASON -> no I generate the vertices by combination types before calling
# this function, so for the vertices only I would have way too many colors...
generate_coloring_scheme_for_graph <- function(vertices_by_types) {
  type_length <- length(vertices_by_types)
  # Number of combinations for the edges (i.e. combos of type of source and type of target)
  number_of_edge_combinations <- (type_length + 1) * (type_length + 1)
  color_seq <- generate_color_sequence(number_of_edge_combinations)
  return(color_seq)
}

################################################################################
# Generate a sequence of colors to use in plots.
# 
# Thanks https://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r
################################################################################
generate_color_sequence <- function(desired_number_of_colors) {
  set.seed(1493)
  # The general group (or first group in general) should be in a gray color to
  # easily distinguish between the important bits from the non important bits
  color_list <- c("#D3D3D3") # #d3d3d3 is the hex for lightgray
  if (desired_number_of_colors <= 75) {
    # Using RColorBrewer
    # The method I found using RColorBrewer can generate up to 74 distinct colors
    require(RColorBrewer)
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    brewer_colors <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    color_list <- append(color_list, 
                         sample(brewer_colors, desired_number_of_colors - 1))
  }
  else if (desired_number_of_colors <= 434) {
    # Using grDevices::colors()
    # The method I found using RColorBrewer can generate up to 433 distinct colors
    r_named_colors <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
    color_list <- append(color_list, 
                         sample(r_named_colors, desired_number_of_colors - 1))
  }
  else {
    # Using randomcoloR
    require(randomcoloR)
    color_list <- append(color_list, 
                         distinctColorPalette(desired_number_of_colors - 1))
  }
  return(color_list)
}

################################################################################
# Takes a graph and a list of sets of vertices according to all type combinations
# and uses that list to assign colors to the vertices and edges according to type.
#
# @param graph 
# @param vertices_by_types The list of sets of vertices of all combinations
#        of basic types
# @param color_seq the color sequence to be used by this graph
#
# @return an igraph graph with edges containing the vertices and proper coloring
################################################################################
prepare_graph_coloring_and_sizing_according_to_vertices <- function(graph, vertices_by_types, color_seq) {
  graph <- color_and_size_vertices(graph, vertices_by_types, color_seq)
  graph <- color_and_size_edges_based_on_vertices(graph, vertices_by_types, color_seq)
  return(graph)
}

###########
# Changes the colors and sizes of the vertices of the graph
#
# color_sequence should have at least 1 more element than vertices_by_types
#
# @return the colored and resized graph
###########
color_and_size_vertices <- function(graph, vertices_by_types, color_sequence) {
  length <- length(vertices_by_types)
  # First the color of vertices that are somehow not associated with anything
  V(graph)$color <- color_sequence[[1]]
  V(graph)$size <- 15
  for (i in 1:length) {
    V(graph)[V(graph)$name %in% vertices_by_types[[i]]]$color <- color_sequence[[i + 1]]
    V(graph)[V(graph)$name %in% vertices_by_types[[i]]]$size <- 22
  }
  return(graph)
}

###########
# Changes the colors and sizes of the edges of the graph
#
# color_sequence should have at least 1 more element than vertices_by_types
#
# @return the colored and resized graph
########### TODO: Figure out how to address all combinations of types
color_and_size_edges_based_on_vertices <- function(graph, vertices_by_types, color_sequence) {
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
# Extract a list of all combinations types               #
#
# @param vertices_by_types The list of sets of vertices of all combinations
#        of basic types
#######################################################
extract_type_list_from_vertices_by_types <- function(vertices_by_types) {
  # Use the list of named sets to have an updated type_list
  type_list <- as.list(names(vertices_by_types))
  # In order to properly align the colors with the non type and the types
  type_list <- append(type_list, "None", after = 0)
  return(type_list)
}

#######################################################
# Creating graph layout and plotting it               #
#
# @param graph
# @param type_list including the None type
# @param color_seq the color sequence to be used by this graph
#######################################################
draw_igraph <- function(graph, type_list, color_seq) {
  # this ensures the starting random position is the same
  # for the layouts that use a random starting position
  set.seed(1492) 
  l <- layout_with_dh(graph, coords = NULL, maxiter = 10, 
                      fineiter = max(10, log2(vcount(graph))), 
                      cool.fact = 0.75, weight.node.dist = 19,
                      weight.border = 0, 
                      #weight.edge.lengths = edge_density(graph)/17,
                      weight.edge.lengths = edge_density(graph)/22,
                      weight.edge.crossings = 100 - sqrt(edge_density(graph)),
                      weight.node.edge.dist = 0.2* (1 - edge_density(graph)))
  
  plot.igraph(graph, layout=l, 
              edge.width=1, edge.arrow.size = 0.05, #edge.length = 10, 
              vertex.label.cex=0.6, 
              vertex.label.family="Helvetica",
              vertex.label.font=1,
              vertex.shape="circle",)
  #edge.arrow.mode=2,)
  #edge.curved=TRUE,)
  
  legend("topleft",bty = "n",
         legend=type_list,
         fill=color_seq, border=NA,
         title = "Vertices")
}

### AUX function for the draw_igraph generating a list of sets of combinations of original types
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

##############################################
# AUXILIARY FUNCTION FOR color_and_size_edges#
##############################################
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