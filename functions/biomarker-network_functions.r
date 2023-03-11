library(here)
source(here("functions", "graph_vertices_functions.R"))
source(here("functions", "graph_edges_functions.R"))
source(here("functions", "data_processing_functions.R"))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Takes an edge dataframe and extracts from it a graph with edges where at
# least one node is in the vertex_list
#
# @param edge_df is the dataframe with the edge data
# @param vertex_list NULL if no filtering desired or a list of vertices otherwise
# @param source_col_index column number for the source
# @param target_col_index column number for the target
#
# @return an igraph graph with edges containing the vertices
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
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
  # There was a bug where the presence of NA caused an error.
  # Will therefore remove any edge that has NA as a vertex
  filtered_edge_df <- na.omit(filtered_edge_df)
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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Generate a sequence of colors to use in plots.
# 
# Thanks https://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Takes a graph and a list of sets of vertices according to all type combinations
# and uses that list to assign attributes (e.g. colors) to the vertices and edges 
# according to type.
#
# @param graph 
# @param vertices_by_types The list of sets of vertices of all combinations
#        of basic types
# @param color_seq the color sequence to be used by this graph
#
# @return an igraph graph with edges containing the vertices and proper coloring
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
prepare_graph_attributes_according_to_vertices <- function(graph, vertices_by_types, color_seq) {
  graph <- adjust_vertices_attributes_according_to_type(graph, vertices_by_types, color_seq)
  graph <- color_edges_based_on_vertices(graph, vertices_by_types, color_seq)
  return(graph)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Extract a list of all combinations types               #
#
# @param vertices_by_types The list of sets of vertices of all combinations
#        of basic types
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
extract_type_list_from_vertices_by_types <- function(vertices_by_types) {
  # Use the list of named sets to have an updated type_list
  type_list <- as.list(names(vertices_by_types))
  # In order to properly align the colors with the non type and the types
  type_list <- append(type_list, "None", after = 0)
  return(type_list)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Creating graph layout and plotting it               #
#
# @param graph
# @param type_list including the None type
# @param color_seq the color sequence to be used by this graph
# @param comm A list derived from a communities object or an empty list
#             (default) if no communities are to be shown
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
draw_igraph <- function(graph, type_list, color_seq, comm = list()) {
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
              vertex.shape="circle",
              mark.groups = comm)
  #edge.arrow.mode=2,)
  #edge.curved=TRUE,)
  
  legend("topleft",bty = "n",
         legend=type_list,
         fill=color_seq, border=NA,
         title = "Vertices")
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# This function takes a list of graphs and generates a unified graph
# with the correct sources for the edges
# and the correct coloring of the vertices and edges
#
# @param graph_list
# @param vertices_by_type named list of vertices by their extended type
# @param extended_source_list
# @param color_sequence at least as big as
#        max(length(vertices_by_type), length(extended_source_list))
#
# @return the graph with the edge sources converged into
#         one source attribute and the source1&source2 attributes removed
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
unify_graphs <- function(graph_list, vertices_by_type, extended_source_list,
                         color_sequence) {
  number_of_graphs <- length(graph_list)
  unified_graph <- graph_list[[1]]
  if (number_of_graphs > 1) {
    for (i in 2:number_of_graphs) {
      unified_graph <- union(unified_graph,graph_list[[i]]) %>%
        converge_edge_sources()
    }
  }

  unified_graph <- unified_graph %>%
    adjust_vertices_attributes_according_to_type(vertices_by_type, color_sequence) %>%
    color_edges_based_on_sources(extended_source_list, color_sequence)
  return(unified_graph)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#' Getting all network object interactions that have valid uniprots for the
#' network objects
#'
#' @param dataset_directory The directory where the interactions dataset and
#'        and the dictionary (network object to uniprot) dataset sit
#' @param interactions_file The file where the interactions between the
#'        network objects are
#' @param dictionary_file The file translating network objects to uniprots
#' @param interaction_columns The column names of the source and target of
#'        the interactions
#' @param dictionary_columns The column names of the network object to uniprot
#'        translation
#' @param uniprot_interaction_names The name to give to the uniprot based
#'        interaction columns
#' @param return_as_uniprots whether to return the interaction df with uniprots
#'
#' @return A df that contains the interactions either as target and source
#'         network objects or target and source uniprots.
#'         The rows are unique.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
get_all_interactions_that_have_backing_uniprots <- function(
    dataset_directory,
    interactions_file,
    dictionary_file,
    interaction_columns = c("Network Object \"FROM\"", "Network Object \"TO\""),
    dictionary_columns = c("Network Object Name", "Input IDs"),
    uniprot_interaction_names = c("UniProt_FROM", "UniProt_TO"),
    return_as_uniprots = FALSE)
{
  require(dplyr)
  require(here)
  dictionary <- readxl::read_xls(here(dataset_directory, dictionary_file), skip = 2) %>%
    dplyr::select(all_of(dictionary_columns)) %>%
    # Remove useless "No ID" rows
    dplyr::filter(!!sym(dictionary_columns[[2]])!="No ID") %>%
    dplyr::distinct(!!sym(dictionary_columns[[1]]), .keep_all = TRUE) %>%
    dplyr::filter(!grepl(';', !!sym(dictionary_columns[[2]])))

  interactions <- readxl::read_xls(here(from_metacore_path, interactions_file),
                                   skip = 2) %>%
    dplyr::select(all_of(interaction_columns)) %>%
    dplyr::distinct()

  uniprot_interactions <- translate_interactions_using_dictionary(
    interactions_df = interactions,
    interaction_columns = interaction_columns,
    dictionary_df = dictionary,
    dictionary_columns = dictionary_columns,
    new_interaction_names = uniprot_interaction_names,
    only_interaction_columns = F
  )

  # After cleaning with the uniprots, let's continue with the network objects
  # or uniprot, depending on return_as_uniprots
  if (return_as_uniprots) {
    returning_columns <- uniprot_interaction_names
  }
  else {
    returning_columns <- interaction_columns
  }

  return(uniprot_interactions %>%
    dplyr::select(all_of(returning_columns)) %>%
    dplyr::distinct())
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#' Extract vertices for all graphs
#' 
#' This function takes the interaction dfs representing edge lists of graphs
#' and extracts all unique vertices from them
#'
#' @param interaction_dfs Named list of interaction dfs
#' @param interaction_columns Names of source and target of interactions
#' @param uniprot_interaction_columns Names of uniprot source and target of interactions
#' @param return_as_uniprots Whether we should return uniprots or not
#'
#' @return A named list with the vertices of each interaction df. The names
#'         from the interaction_dfs are inherited correctly
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
extract_vertices_for_all_graphs <- function(
    interaction_dfs,
    interaction_columns = c("Network Object \"FROM\"", "Network Object \"TO\""),
    uniprot_interaction_columns = c("UniProt_FROM", "UniProt_TO"),
    return_as_uniprots = FALSE
    )
{
  if (return_as_uniprots) {
    return(lapply(interaction_dfs, extract_all_distinct_objects, uniprot_interaction_columns))
  } else {
    return(lapply(interaction_dfs, extract_all_distinct_objects, interaction_columns))
  }
}

# Auxiliary function to check if dataframe is populated
is_df_populated <- function(df)
{
  return(nrow(df) > 0)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# This 
#
# @param community
# @param adj.p.cutoff 
# @param term_reject which terms to reject
# @param pathway_datasets char vector of pathway databases
#
# @return 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
community_enrichment <- function(community, adj.p.cutoff = 0.1,
                                      term_reject = "infection|cancer|carcinoma|virus|[v|V]iral|trypanosomiasis|Malaria|Tuberculosis|Melanoma|leukemia|Influenza|disease|Amoebiasis",
                                      pathway_datasets = c("KEGG_2021_Human", "MSigDB_Hallmark_2020", "Reactome_2022", "HDSigDB_Human_2021", "WikiPathway_2021_Human"))
{
  require(enrichR)
  require(dplyr)
  ### KEGG enrichment
  enrc <- enrichR::enrichr(community, databases = pathway_datasets)
  # Due to a bug where sometimes empty members are of type logical causing
  # `bind_rows` to fail, we filter empty ones
  filtered_enrc <-
    lapply(enrc, function(x) dplyr::filter(x, Adjusted.P.value < adj.p.cutoff &
                                             !grepl(term_reject,Term)))
  enrichments <-
    bind_rows(Filter(is_df_populated, filtered_enrc))
  ### Turn pathway-genes (one-many) association table into a flat pathway-gene (one-one) table
  if(nrow(enrichments) > 0)
  {
    return((do.call(rbind, apply(enrichments, 1, function(x) data.frame(Term = x["Term"],
                                                                Gene = unlist(strsplit(x["Genes"], split = ";")),
                                                                row.names = NULL))) %>%
             dplyr::group_by(Gene) %>%
             dplyr::arrange(.by_group = T))[,c(2,1)])
  }
  else
  {
    return(data.frame())
  }
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#' Apply enrichment to all communities with a fail safe for EnirchR timeouts
#'
#' @param communities_list List of communities
#' @param i The starting index to start appednding enrichments. Either 1 or
#'          wherever it got stuck if you accidentally stopped it
#'
#' @return The (non empty) enriched communities
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
apply_enrichment_to_communities <- function(communities_list, i, enriched_communities) {
  num_of_communities <- length(communities_list)
  while (i <= num_of_communities) {
    cat("Index", i, "\n")
    tryCatch({enriched_communities[[i]] <-
      community_enrichment(communities_list[[i]])
    i <- i + 1},
    error = function(e) {
      cat("Got stuck at index", i, "\n")
    }
    )
  }
  names(enriched_communities) <- 1:num_of_communities
  enriched_communities <-
    Filter(is_df_populated, enriched_communities)
}
