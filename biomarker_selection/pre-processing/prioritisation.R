##################################################
## Project: MIRIADE
## Script purpose: Interpret MIRIADE biomarker lists according to known mechanisms
## in pathways, databases etc.
## Date: 21.01.2022, 03.10.2022
## Author: Marek Ostaszewski, Felicia Burtscher
##################################################

options(stringsAsFactors = F)
library(here)

source(here("functions", "data_processing_functions.R"))

datasets_root_directory <- define_datasets_root()

### Read the harmonised, annotated biomarkers
# annotated_BMs <- read.table("_notgit/annotated_candidate_biomarkers.tsv", sep = "\t", header = T)
annotated_BMs <- read.table(file.path(datasets_root_directory, "annotated_candidate_biomarkers.tsv"),
                            sep = "\t", header = T)

library(dplyr)

### DisGeNet
### Read the DisGeNET file, can be replaced with an API call for better reproducibility (as any changes are not captured in .tsv file)
### Use only the C10 mapping, "Neurological disorders" and score > 0.1
# dgn_c10 <- read.table("_notgit/knowledge_bases/all_gene_disease_associations.tsv",
dgn_c10 <- read.table(file.path(datasets_root_directory, "knowledge_bases/curated_gene_disease_associations.tsv"),
                      sep = "\t", quote = "", comment.char = "", header = T) %>%
  dplyr::filter(grepl("(^C10;)|(^C10$)|(;C10$)|(;C10;)", diseaseClass) & score > 0.3)
### Integrate with the main table, count how many hits each of the BMs has
annotated_BMs <- dplyr::mutate(annotated_BMs,
                                DGN_hits = sapply(SYMBOL,
                                                  function(x) sum(dgn_c10$geneSymbol == x)))

### Disease Maps
### Access selected Disease Maps and cound how often a given BM is present in these maps
### A wrapper function: for a given diseae map url, get the default map, 
### get elements of all diagrams and extract UniProt ids

### Load facilitating functions from minerva_access.R file
source("https://gitlab.lcsb.uni.lu/minerva/minervar/-/raw/master/Rscripts/access_functions.R")

### A convenience function to return UniProt ids per Disease Map
wrap_map <- function(dmap_url, project = NULL) {
  ### Get MINERVA elements
  components <- get_map_components(dmap_url, project_id = project,
                                   e_columns = "id,references", r_columns = "id")
  ### return UNIPROT references in MINERVA diagrams
  return(get_components_annotations(components, "UNIPROT"))
}

### Selected Disease Maps: PD, Ageing and AD
pd_ups <- wrap_map("https://pdmap.uni.lu/minerva/api/")
ag_ups <- wrap_map("https://progeria.uni.lu/minerva/api/")
ad_ups <- wrap_map("https://minerva-dev.lcsb.uni.lu/minerva/api/","alzpath_8APR")

### A combined look-up table
dm_ups <- table(unlist(c(pd_ups, ag_ups, ad_ups)))

### Integrate with the main table
annotated_BMs <- dplyr::mutate(annotated_BMs, 
                               DMaps_hits = sapply(UniProt, function(x) dm_ups[x])) %>%
  dplyr::mutate(DMaps_hits = ifelse(is.na(DMaps_hits), 0, DMaps_hits)) ### Clean up the NAs

sort_table <- function(bm_table, dis) {
  dplyr::filter(bm_table, disease %in% dis) %>%
    arrange(adj.p, desc(DGN_hits), desc(DMaps_hits))
}

alz_bms <- sort_table(annotated_BMs, "ALZ")
ftd_bms <- sort_table(annotated_BMs, "FTD")
dlb_bms <- sort_table(annotated_BMs, "DLB")
com_bms <- sort_table(annotated_BMs, c("ALZ,FTD,DLB","ALZ,DLB","ALZ,FTD","FTD,DLB"))

###
### Tissue location based on the Human Protein Atlas
###

# hpa_brain_gtex <- read.table("_notgit/rna_brain_gtex.tsv", sep = "\t", header = T) %>%
hpa_brain_gtex <- read.table(file.path(datasets_root_directory, "HPA/rna_brain_gtex.tsv"),
                             sep = "\t", header = T) %>%
  dplyr::filter(!(Brain.region %in% c("amygdala", "pituitary gland", "retina"))). ## TODO: READ UP ON THE BIOLOGY
# hpa_brain_fant <- read.table("_notgit/rna_brain_fantom.tsv", sep = "\t", header = T) %>%
hpa_brain_fant <- read.table(file.path(datasets_root_directory, "HPA/rna_brain_fantom.tsv"),
                             sep = "\t", header = T) %>%
  dplyr::filter(!(Brain.region %in% c("amygdala", "pituitary gland", "retina",
                                      "medulla oblongata", "pons")))

### Get HPA-GTEX entries that have values >= 3, but not those where "spinal cord" has the maximum expression
hpa_brain_gtex_genes <- dplyr::group_by(hpa_brain_gtex, Gene, Gene.name) %>%
  dplyr::filter(max(nTPM) >= 3) %>%
  dplyr::select(Gene.name) %>%
  dplyr::distinct() %>% pull()

library(OmnipathR)
op_narrow <- OmnipathR::import_all_interactions() %>%
  dplyr::filter(consensus_direction == 1 & n_references > 1 & (consensus_inhibition + consensus_stimulation) > 0) %>%
  dplyr::select(-starts_with("consensus_"), -starts_with("is_"))

library(igraph)
op_graph <- function(omnipath_net, bms) {
  fop <- dplyr::select(omnipath_net, source_genesymbol, target_genesymbol) %>%
    dplyr::filter(source_genesymbol %in% bms$SYMBOL | target_genesymbol %in% bms$SYMBOL)
  g <- igraph::graph.edgelist(as.matrix(fop), directed = TRUE)
  bmvs <- igraph::V(g)[name %in% bms$SYMBOL]
  return(list(graph = g, 
              inout = data.frame(name = bmvs$name, 
                                 OmniPath_out = igraph::degree(g, v = bmvs, mode = "out"),
                                 OmniPath_in = igraph::degree(g, v = bmvs, mode = "in"))))
}

###
### Combine enrichment analysis of the biomarkers with previous information about connectivity and brain tissue expression
###

library(enrichR)
ageing <- c("Aging_Perturbations_from_GEO_down", "Aging_Perturbations_from_GEO_up")
pathways <- c("KEGG_2021_Human", "MSigDB_Hallmark_2020")

library(ReactomePA)
reac_reject <- ""

### A wrapper function to run enrichment and OmniPath lookup, and combine the analysis
combined_analysis <- function(bms, adj.p.cutoff = 0.1,
                    kegg_reject = "infection|cancer|carcinoma|virus|[v|V]iral|trypanosomiasis|Malaria|Tuberculosis|Melanoma|leukemia|Influenza|disease|Amoebiasis") {
  ### KEGG enrichment
  enrc <- enrichR::enrichr(bms$SYMBOL, databases = pathways)
  kegg <- dplyr::filter(enrc$KEGG_2021_Human, Adjusted.P.value < adj.p.cutoff & !grepl(kegg_reject,Term))
  ### Turn pathway-genes (one-many) association table into a flat pathway-gene (one-one) table
  kegg_flat <- do.call(rbind, apply(kegg, 1, function(x) data.frame(Term = x["Term"],
                                                                    Genes = unlist(strsplit(x["Genes"], split = ";")),
                                                                    row.names = NULL)))
  ### Use the flat table to compute pathway hits and pathway names per gene
  bms <- dplyr::mutate(bms, KEGG_hits = sapply(SYMBOL, function(x) sum(kegg_flat$Genes == x))) %>%
    dplyr::mutate(bms, KEGG_pathways = sapply(SYMBOL, function(x) paste(with(kegg_flat, Term[Genes == x]), collapse = ",")))

  ### Reactome enrichment
  reac <- ReactomePA::enrichPathway(bms$ENTREZID)
  reac <- dplyr::filter(reac@result, p.adjust < adj.p.cutoff)
  ### As above, for Reactome results
  reac_flat <- do.call(rbind, apply(reac, 1, function(x) data.frame(Description = x["Description"], 
                                                                    geneID = unlist(strsplit(x["geneID"], split = "/")),
                                                                    row.names = NULL)))
  ### Use the flat table to compute pathway hits and pathway names per gene
  bms <- dplyr::mutate(bms, Reactome_hits = sapply(ENTREZID, function(x) sum(reac_flat$geneID == x))) %>%
    dplyr::mutate(bms, Reactome_pathways = sapply(ENTREZID, function(x) paste(with(reac_flat, Description[geneID == x]), collapse = ",")))

  ### OmniPath connectivity
  bms_g <- op_graph(op_narrow, bms)

  # ### Use the graph degree table to enrich the table, preferred over "merge"
  # bms <- dplyr::mutate(bms, OmniPath_out = sapply(SYMBOL, function(x) bms_g$inout[x,2])) %>%
  #   dplyr::mutate(bms, OmniPath_in = sapply(SYMBOL, function(x) bms_g$inout[x,3]))

  ### Wrap things up
  return(list(kegg = kegg, reac = reac, opco = bms_g, enriched_bms = bms))
}

## TODO: UNDERSTAND (BIOLOGICAL) REASONING HERE
alz_sum <- combined_analysis(alz_bms)
alz_select <- dplyr::filter(alz_sum[["enriched_bms"]], SYMBOL %in% c("OLR1", "GLO1", "SNAP29", "KYAT1"))
ftd_sum <- combined_analysis(ftd_bms)
ftd_select <- dplyr::filter(ftd_sum[["enriched_bms"]], SYMBOL %in% c("RTN4R", "ROBO2", "SMPD1", "CTSF"))
dlb_sum <- combined_analysis(dlb_bms)
dlb_select <- dplyr::filter(dlb_sum[["enriched_bms"]], SYMBOL %in% c("PRCP", "COL4A1", "IDUA"))
com_sum <- combined_analysis(com_bms)
com_select <- dplyr::filter(com_sum[["enriched_bms"]], SYMBOL %in% c("SDC4", "PAM", "FKBP4", "PRDX1"))

write.table(rbind(alz_select, rbind(ftd_select, rbind(dlb_select,com_select))),
            file = "selected_new_biomarkers.tsv", sep = "\t", quote = F, row.names = F)

write.table(rbind(alz_sum$enriched_bms, rbind(ftd_sum$enriched_bms, rbind(dlb_sum$enriched_bms,com_sum$enriched_bms))),
            file = "all_new_biomarkers.tsv", sep = "\t", quote = F, row.names = F)

### Visualisation of the pathway enrichment overlaps

postfilter <- function(summary) {
  dplyr::filter(summary[["enriched_bms"]], (DGN_hits + DMaps_hits > 0) & (KEGG_hits + Reactome_hits) > 0)
}

ftd_narrow <- postfilter(ftd_sum)

kegg_net <- rbind(
  cbind(alz_sum$kegg$Term, "ALZ"), cbind(ftd_sum$kegg$Term, "FTD"), cbind(dlb_sum$kegg$Term, "DLB")
)

g <- igraph::graph.edgelist(kegg_net, directed = FALSE)
plot(g)

reac_net <- rbind(
  cbind(alz_sum$reac$Description, "ALZ"), cbind(ftd_sum$reac$Description, "FTD"), cbind(dlb_sum$reac$Description, "DLB")
)

reac_g <- igraph::graph.edgelist(reac_net, directed = FALSE)
plot(reac_g)



##  further in-depth study of the disease biomarker candidates

#ftds <- read.table("_notgit/ftd candidates")
ftds <- read.table(file.path(datasets_root_directory, "Olink/20201216_OLINK_Data/cleaned_VUMC_olink_FTD.tsv"),
                   header = TRUE)

dplyr::filter(ftd_bms, UniProt %in% ftds[,1])

curl_head <- "curl -X GET --header 'Accept: application/json' 'https://www.ebi.ac.uk/europepmc/annotations_api/annotationsByRelationship?firstEntity="
curl_tail <- "&secondEntity=dementia'"

eupmc <- lapply(ftd_sum$enriched_bms$SYMBOL, function(g) fromJSON(system(paste0(curl_head,g,curl_tail), intern = T)))
# Assuming gene order is preserved
ftd_prep_bms <- cbind(ftd_sum$enriched_bms, TM_hits = sapply(eupmc, function(x) length(x$articles$source)))

hitcols <- c("DGN_hits", "DMaps_hits", "KEGG_hits", "Reactome_hits", "TM_hits")


ftd_supported <- dplyr::filter(ftd_prep_bms, rowSums(ftd_prep_bms[,hitcols]) > 0)
ftd_low_evidence <- dplyr::filter(ftd_prep_bms, rowSums(ftd_prep_bms[,hitcols]) == 0)


proteindb_head <- "https://www.proteomicsdb.org/proteomicsdb/logic/api/proteinexpression.xsodata/InputParams(PROTEINFILTER='"
proteindb_tail <- "',MS_LEVEL=1,TISSUE_ID_SELECTION='',TISSUE_CATEGORY_SELECTION='tissue;fluid',SCOPE_SELECTION=1,GROUP_BY_TISSUE=1,CALCULATION_METHOD=0,EXP_ID=-1)/Results?$select=UNIQUE_IDENTIFIER,TISSUE_ID,TISSUE_NAME,TISSUE_SAP_SYNONYM,SAMPLE_ID,SAMPLE_NAME,AFFINITY_PURIFICATION,EXPERIMENT_ID,EXPERIMENT_NAME,EXPERIMENT_SCOPE,EXPERIMENT_SCOPE_NAME,PROJECT_ID,PROJECT_NAME,PROJECT_STATUS,UNNORMALIZED_INTENSITY,NORMALIZED_INTENSITY,MIN_NORMALIZED_INTENSITY,MAX_NORMALIZED_INTENSITY,SAMPLES&$format=json"

ftd_tissues <- lapply(ftd_prep_bms$UniProt, function(x) fromJSON(paste0(proteindb_head, x, proteindb_tail)))

### Source: Human Protein Atlas, https://www.proteinatlas.org/about/download
# hpa_ihch_tissues <- read.table("_notgit/normal_tissue.tsv", sep = "\t", header = T)
hpa_ihch_tissues <- read.table(file.path(datasets_root_directory, "HPA/normal_tissue.tsv"),
                               sep = "\t", header = T)
# hpa_rna_tissues <- read.table("_notgit/rna_tissue_consensus.tsv", sep = "\t", header = T)
hpa_rna_tissues <- read.table(file.path(datasets_root_directory, "HPA/rna_tissue_consensus.tsv"),
                              sep = "\t", header = T)

table(hpa_ihch_tissues$Reliability)

hpa_ihch_tissues_slim <- dplyr::filter(hpa_ihch_tissues, Reliability != "Uncertain" & Level != "Not detected") %>%
  dplyr::group_by(Gene.name) %>% dplyr::summarise(Tissues = paste(unique(Tissue), collapse = ","))

unique(hpa_ihch_tissues$Gene.name)

# dea_preds <- read.table("_notgit/biomarkers_annotated_EV_features.csv", sep = ",", header = T) %>%
#   dplyr::select(-dplyr::matches("^[A-Z]$")) %>%
#   dplyr::select(-dplyr::ends_with(c("_exposed", "_netsurfp2", "_UP", "_all", "_MSD")))

