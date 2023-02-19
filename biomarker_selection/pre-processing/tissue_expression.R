##################################################
## Project: MIRIADE
## Script purpose: Combine protein tissue levels for selected MIRIADE biomarkers
## from public sources
## Date: 02.06.2022, 26.09.2022
## Author: Marek Ostaszewski, Felicia Burtscher
##################################################
library(here)

source(here("functions", "data_processing_functions.R"))

datasets_root_directory <- define_datasets_root()

### Previously pooled biomarkers
# annotated_BMs <- read.table(file = "_notgit/annotated_candidate_biomarkers.tsv", sep = "\t", header = T)
annotated_BMs <- read.table(file = file.path(datasets_root_directory, "annotated_candidate_biomarkers.tsv"), sep = "\t", header = T)


### Human Protein Atlas sources, loaded locally, obtained at: https://www.proteinatlas.org/about/download
### HPA dataset load
### Protein expression (ICH): https://www.proteinatlas.org/download/normal_tissue.tsv.zip
# normal_tissue_ich <- read.table(file = "_notgit/normal_tissue.tsv", sep = "\t", header = T, quote = "")
normal_tissue_ich <- read.table(file = file.path(datasets_root_directory, "HPA/normal_tissue.tsv"), sep = "\t", header = T, quote = "")
### Gene expression (RNASeq): https://www.proteinatlas.org/download/rna_tissue_consensus.tsv.zip
# normal_tissue_RNAseq <- read.table(file = "_notgit/rna_tissue_consensus.tsv", sep = "\t", header = T, quote = "")
normal_tissue_RNAseq <- read.table(file = file.path(datasets_root_directory, "HPA/rna_tissue_consensus.tsv"), sep = "\t", header = T, quote = "")
### CNS tissues, for aggregation
hpa_cns_tissues <- c("amygdala", "basal ganglia", "caudate", "cerebellum", "cerebral cortex", 
                     "dorsal raphe", "hippocampus", "hippocampal formation", "hypothalamus", "medulla oblongata", "midbrain",
                     "olfactory bulb", "pituitary gland", "pons", "substantia nigra","thalamus", "white matter",
                     "spinal cord")

library(dplyr)

### Filter the HPA ICH data per gene symbol and aggregate CNS-related data
aggregated_ich <- normal_tissue_ich %>%
  ### Remove non-biomarker entries
  dplyr::filter(Gene.name %in% annotated_BMs$SYMBOL) %>%
  ### Keep only values 'High' and 'Medium'
  dplyr::filter(Reliability != "Uncertain" & Level %in%  c("High", "Medium")) %>%
  ### For better comparison, turn categprical 'High' and 'Medium' into numeric '2' and '1'
  dplyr::mutate(Level = dplyr::case_when(Level == "High" ~ 2,
                                         TRUE ~ 1)) %>%
  ### Aggregate CNS tissues (see the grouping above)
  dplyr::mutate(Tissue = dplyr::case_when(Tissue %in% hpa_cns_tissues ~ "CNS",
                                          TRUE ~ Tissue)) %>% 
  ### Take the max value for the group
  dplyr::group_by(Gene.name, Tissue) %>% dplyr::summarise(Level = max(Level)) %>%
  ### Summarise on the gene level
  dplyr::summarise(cns_level_ich = Level[Tissue == "CNS"],
                   n_entries_hi_mid_ich = n(),
                   entries_hi_mid_ich = paste(sort(Tissue), collapse = ","))

### Filter the HPA RNAseq data per gene symbol and aggregate CNS-related data
aggregated_rna <- normal_tissue_RNAseq %>%
  ### Remove non-biomarker entries
  dplyr::filter(Gene.name %in% annotated_BMs$SYMBOL) %>%
  ### Aggregate CNS tissues (see the grouping above)
  dplyr::mutate(Tissue = dplyr::case_when(Tissue %in% hpa_cns_tissues ~ "CNS",
                                          TRUE ~ Tissue)) %>%
  ### Take the max value for the group
  dplyr::group_by(Gene.name, Tissue) %>% dplyr::summarise(norm_exp = max(nTPM)) %>%
  ### Percentage within gene group
  dplyr::mutate(perc_max = norm_exp/max(norm_exp)) %>%
  ### Filter by percentage within gene group
  dplyr::filter(perc_max >= 0.5) %>%
  ### Summarise on the gene level
  dplyr::summarise(cns_perc_max_rna = perc_max[Tissue == "CNS"],
                   n_entries_50perc_max_rna = n(),
                   entries_50perc_max_rna = paste(sort(Tissue), collapse = ","))

### Load data from Protein DB, if the table does not exist - regenerate
### Data acquired from "https://www.proteomicsdb.org"
# if(!file.exists("_notgit/proteindb_tissues_annotated_BM.tsv")) {
if(!file.exists(file.path(datasets_root_directory, "proteindb_tissues_annotated_BM.tsv"))) {
  ### Head and tail of a curl query for a single UniProt
  proteindb_head <- "https://www.proteomicsdb.org/proteomicsdb/logic/api/proteinexpression.xsodata/InputParams(PROTEINFILTER='"
  proteindb_tail <- "',MS_LEVEL=1,TISSUE_ID_SELECTION='',TISSUE_CATEGORY_SELECTION='tissue;fluid',SCOPE_SELECTION=1,GROUP_BY_TISSUE=1,CALCULATION_METHOD=0,EXP_ID=-1)/Results?$select=UNIQUE_IDENTIFIER,TISSUE_ID,TISSUE_NAME,TISSUE_SAP_SYNONYM,UNNORMALIZED_INTENSITY,NORMALIZED_INTENSITY,MIN_NORMALIZED_INTENSITY,MAX_NORMALIZED_INTENSITY&$format=json"
  library(jsonlite)
  ### Run API request for each UniProt of the annotated biomarkers
  proteindb_tissues <- lapply(annotated_BMs$UniProt, function(x) fromJSON(paste0(proteindb_head, x, proteindb_tail), flatten = F))
  names(proteindb_tissues) <- annotated_BMs$UniProt
  proteindb_tab <- data.frame()
  ### Selected names in the table; dplyr chokes on a nested data.frame, didn't want to use tidyverse just for this
  sel_names <- c("UNIQUE_IDENTIFIER", "TISSUE_ID", "TISSUE_NAME", "TISSUE_SAP_SYNONYM", "NORMALIZED_INTENSITY")
  ### Remove empty entries, it's much easier to work this way
  proteindb_tissues <- proteindb_tissues[sapply(proteindb_tissues, function(x) length(x$d$results)) > 0]
  ### Bind all the tables for each entry, harmonise for only selected column names
  proteindb_tab <- do.call(rbind, lapply(proteindb_tissues, function(x) x$d$results %>% dplyr::select(sel_names)))
  ### Write the source JSON and the tidied table
  write.table(proteindb_tab, file = file.path(datasets_root_directory, "proteindb_tissues_annotated_BM.tsv"),
              sep = "\t", col.names = T, row.names = F, quote = F)
  # write_json(proteindb_tissues, "_notgit/proteindb_tissues")
  write_json(proteindb_tissues, file.path(datasets_root_directory, "proteindb_tissues.json"))
} else {
  # proteindb_tab <- read.table("_notgit/proteindb_tissues_annotated_BM.tsv", sep = "\t", header = T)
  proteindb_tab <- read.table(file.path(datasets_root_directory, "proteindb_tissues_annotated_BM.tsv"),
                              sep = "\t", header = T)
}

### Filter the HPA RNAseq data per gene symbol and aggregate CNS-related data
aggregated_proteindb <- proteindb_tab %>% mutate(NORMALIZED_INTENSITY = as.numeric(NORMALIZED_INTENSITY)) %>%
  ### Remove non-biomarker entries
  dplyr::filter(UNIQUE_IDENTIFIER %in% annotated_BMs$UniProt) %>%
  ### Aggregate CNS tissues (see the grouping above)
  dplyr::mutate(TISSUE_SAP_SYNONYM = dplyr::case_when(TISSUE_SAP_SYNONYM %in% c("brain", "cerebrospinal_fluid") ~ "CNS",
                                                      TRUE ~ TISSUE_SAP_SYNONYM)) %>%
  ### Take the max value for the group
  dplyr::group_by(UNIQUE_IDENTIFIER, TISSUE_SAP_SYNONYM) %>% dplyr::summarise(norm_int = max(NORMALIZED_INTENSITY)) %>%
  ### Percentage within gene group
  dplyr::mutate(perc_max = norm_int/max(norm_int)) %>%
  ### Filter by percentage within gene group
  dplyr::filter(perc_max >= 0.5) %>%
  ### Summarise on the gene level
  dplyr::summarise(cns_perc_max_protdb = perc_max[TISSUE_SAP_SYNONYM == "CNS"],
                   n_entries_50perc_max_protdb = n(),
                   entries_50perc_max_protdb = paste(sort(TISSUE_SAP_SYNONYM), collapse = ","))
### Chain merge to combine the three tables
summary_table_tissue_BMs <- merge(annotated_BMs, aggregated_ich, by.x = "SYMBOL", by.y = "Gene.name", all.x = TRUE)
summary_table_tissue_BMs <- merge(summary_table_tissue_BMs, aggregated_rna, by.x = "SYMBOL", by.y = "Gene.name", all.x = TRUE)
summary_table_tissue_BMs <- merge(summary_table_tissue_BMs, aggregated_proteindb, by.x = "UniProt", by.y = "UNIQUE_IDENTIFIER", all.x = TRUE)

### Write the table down
# write.table(summary_table_tissue_BMs, file = "_notgit/summary_tissues_annotated_BM.tsv",
write.table(summary_table_tissue_BMs, file = file.path(datasets_root_directory, "summary_tissues_annotated_BM.tsv"),
            sep = "\t", col.names = T, row.names = F, quote = F)

tmp <- cbind(is.na(summary_table_tissue_BMs$cns_level_ich),
             is.na(summary_table_tissue_BMs$cns_perc_max_rna),
             is.na(summary_table_tissue_BMs$cns_perc_max_protdb))

table(rowSums(tmp))

