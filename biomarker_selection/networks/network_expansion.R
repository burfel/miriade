library(dplyr)
library(here)
library(igraph)
library(future)
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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Taking the csf metacore interactions and combining with the miriade ones
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### Community detection  ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
random_walk_communities <- igraph::walktrap.community(unified_graph, steps = 100)
# Number of communities with more than two nodes
print(length(Filter(function(x) length(x)>2, communities(random_walk_communities))))
# Random walk does only 4 step by default, resulting in 662 communities.
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
# Number of communities with more than two nodes
print(length(Filter(function(x) length(x)>2, communities(v_random_walk_communities))))
# v_edge_betweenness_communities <- edge.betweenness.community(vertex_intersection_graph)
# v_infomap_communities <- igraph:: infomap.community(vertex_intersection_graph)

# Vertex intersection graph community enrichment
datasets_root_directory <<- define_datasets_root()
vec_filepath <- file.path(datasets_root_directory,
                          "Enriched_pathways/VertexIntersectionRandomWalkEnrichedCommunities.xlsx")
if(!file.exists(vec_filepath)) {
  v_enriched_communities <-
    Filter(is_df_populated,
           sapply(communities(v_random_walk_communities), community_enrichment))
  write_df_list_to_xlsx(v_enriched_communities, vec_filepath)
} else {
  require(XLConnect)
  # Read the Excel file into a workbook object
  wb <- loadWorkbook(vec_filepath)

  # Get the names of all sheets in the workbook
  sheet_names <- getSheets(wb)

  # Loop through the sheets and read each one into a data frame
  v_enriched_communities <- lapply(sheet_names, function(sheet) readWorksheet(wb, sheet = sheet))
  names(v_enriched_communities) <- sheet_names

  remove(wb, sheet_names)
}

# Unified graph community enrichment
uec_filepath <- file.path(datasets_root_directory,
                          "Enriched_pathways/UnifiedRandomWalkEnrichedCommunities.xlsx")
if(!file.exists(uec_filepath)) {
  tryCatch({
    i <<- length(enriched_communities) + 1
  }, error = function(e) {
    enriched_communities <<- vector("list")
    i <<- 1
  })
  enriched_communities <-
    apply_enrichment_to_communities(communities(random_walk_communities),
                                    i, enriched_communities)
  write_df_list_to_xlsx(enriched_communities, uec_filepath)
} else {
  require(XLConnect)
  # Read the Excel file into a workbook object
  wb <- loadWorkbook(uec_filepath)

  # Get the names of all sheets in the workbook
  sheet_names <- getSheets(wb)

  # Loop through the sheets and read each one into a data frame
  enriched_communities <- lapply(sheet_names, function(sheet) readWorksheet(wb, sheet = sheet))
  names(enriched_communities) <- sheet_names

  remove(wb, sheet_names)
}
# Aggregate according to Gene so that all of its terms are concatenated with ';'
aggregated_enriched_communities <-
  lapply(enriched_communities, function(df) df %>% 
           group_by(Gene) %>% 
           summarise(Term = paste(Term, collapse = ";")))
# Aggregate according to Gene so that all of its terms are concatenated with ';'
aggregated_v_enriched_communities <-
  lapply(v_enriched_communities, function(df) df %>% 
  group_by(Gene) %>% 
  summarise(Term = paste(Term, collapse = ";")))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### Removing lonely nodes ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

all_lonely_nodes <- which(degree(vertex_intersection_graph)==0)
vertex_intersection_no_lonely_graph <- delete.vertices(vertex_intersection_graph, all_lonely_nodes)
vnl_random_walk_communities <- igraph::walktrap.community(vertex_intersection_no_lonely_graph, steps = 100)
# Unified graph community enrichment
vnlec_filepath <- file.path(datasets_root_directory,
                          "Enriched_pathways/VertexIntersectionNoLonelyRandomWalkEnrichedCommunities.xlsx")
if(!file.exists(vnlec_filepath)) {
  tryCatch({
    i <<- length(vnl_enriched_communities) + 1
  }, error = function(e) {
    vnl_enriched_communities <<- vector("list")
    i <<- 1
  })
  vnl_enriched_communities <-
    apply_enrichment_to_communities(communities(vnl_random_walk_communities),
                                    i, vnl_enriched_communities)
  write_df_list_to_xlsx(vnl_enriched_communities, vnlec_filepath)
} else {
  require(XLConnect)
  # Read the Excel file into a workbook object
  wb <- loadWorkbook(vnlec_filepath)

  # Get the names of all sheets in the workbook
  sheet_names <- getSheets(wb)

  # Loop through the sheets and read each one into a data frame
  vnl_enriched_communities <- lapply(sheet_names, function(sheet) readWorksheet(wb, sheet = sheet))
  names(vnl_enriched_communities) <- sheet_names

  remove(wb, sheet_names)
}
aggregated_vnl_enriched_communities <-
  lapply(vnl_enriched_communities, function(df) df %>%
           group_by(Gene) %>%
           summarise(Term = paste(Term, collapse = ";")))

pathways<-unique(unlist(lapply(vnl_enriched_communities, function(x) x$Term)))

pathway_related_subgraph <- function(graph, pathway, enriched_vertices)
{
  vertices <- lapply(enriched_vertices, function(x) x[x$Term == pathway]$UniProt)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### enrichment comparison ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
enriched_unified <- community_enrichment(V(unified_graph)$name)
aggregated_enriched_unified <- enriched_unified %>%
  group_by(Gene) %>%
  summarise(Term = paste(Term, collapse = ";"))
aggregated_enriched_communities <- bind_rows(aggregated_enriched_communities)
unified_pathways <- unique(enriched_unified$Term)
communities_pathways <- unique(bind_rows(enriched_communities)$Term)
pathway_df <- as.data.frame(rbind(cbind(unified_pathways, "Unified graph"), cbind(communities_pathways, "Communities")))
colnames(pathway_df) <- c("Pathway", "Source")
aggregated_pathway_df <- pathway_df %>%
  dplyr::group_by(Pathway) %>%
  dplyr::summarise(Source = paste(Source, collapse = ";"))
pathways_in_both <- aggregated_pathway_df %>% dplyr::filter(Source == "Unified graph;Communities")
num_of_pathways_unique_to_communities <-
  length(communities_pathways) - nrow(pathways_in_both)
pathways_only_in_communities <-
  dplyr::filter(aggregated_pathway_df, Source == "Communities")$Pathway

color_source_nodes <- function(pathway_overlap_graph,
                               sources = c("Unified graph", "Communities")) {
  V(pathway_overlap_graph)$color <- "orange"
    V(pathway_overlap_graph)[V(pathway_overlap_graph)$name %in% sources]$color <- "green"
      return(pathway_overlap_graph)
}

pathway_overlap_graph <- pathway_df %>%
  dplyr::filter(!(Pathway %in% pathways_only_in_communities)) %>%
  rbind(c(paste("... and", num_of_pathways_unique_to_communities, "additional pathways"),
          "Communities")) %>%
  graph_from_edge_df_filtered_by_genes(source_col_index = 1, target_col_index = 2) %>%
  color_source_nodes()
draw_igraph(pathway_overlap_graph)
