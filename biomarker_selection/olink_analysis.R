##################################################
## Project: MIRIADE
## Script purpose: Combine expression levels from Olink tables, integrate prior knowledge
## Date: 21.12.2020
## Author: Marek Ostaszewski, Felicia Burtscher
##################################################

setwd("~/Documents/_PhD/UNILUX/CODING/miriade/biomarker_selection")

library(reshape2)
library("scales")

options(stringsAsFactors = F)

### Read the olink mapping file
#olp <- read.csv("Data/Olink_proteins_list_Marek.csv", sep = ";")
olp <- read.csv("Data/Olink_proteins_list_Marek.csv")

### Read the datasets (sorted by increasing adj. p-values), add a position to each
### the Alzheimers dataset
alzp <- read.csv("Data/FVOL_AD_CON.csv", sep = ",", 
                 col.names = c("Name", "Eff.ALZ", "p.val.ALZ", "adj.p.ALZ"))
alzp <- cbind(alzp, pos.ALZ = 1:nrow(alzp)) # fixed from pos.AlZ to pos.ALZ 17.12.2020 23:29

### the Frontotemporal Dementia dataset
ftdp <- read.csv("Data/FVOL_FTD_CON.csv", sep = ",",
                 col.names = c("Name", "Eff.FTD", "p.val.FTD", "adj.p.FTD"))
ftdp <- cbind(ftdp, pos.FTD = 1:nrow(ftdp))

### the Dementia with Lewy Bodies dataset
dlbp <- read.csv("Data/FVOL_DLB_CON.csv", sep = ",",
                 col.names = c("Name", "Eff.DLB", "p.val.DLB", "adj.p.DLB"))
dlbp <- cbind(dlbp, pos.DLB = 1:nrow(dlbp))

all <- merge(merge(alzp, ftdp, by = "Name"), dlbp, by = "Name")
# Which of the names in all are not in the protein list?
to_change <- which(!(all$Name %in% olp$Assay))
# None  -- all good.

### A convenience function to help split multiple ids in the dataset
split_multi_id <- function(ftab, multi_id, split_ids) {
  row <- which(ftab$Name == multi_id)
  ftab$Name[row] <- split_ids[1]
  if(length(split_ids) > 1) {
    for(i in 2:length(split_ids)) {
      rbind(ftab, ftab[row,])
      ftab$Name[row] <- split_ids[i]
    }
  }
  return(ftab)
}

### ids to split are: "FUT3/FUT5" "IL-27"     "IL12"      "MIC-A/B"
### here handled in-code, for bigger instances it should be a mapping file
all <- split_multi_id(all, "FUT3/FUT5", c("FUT3", "FUT5"))
all <- split_multi_id(all, "IL-27", c("IL-27A", "IL-27B")) ## different isoforms??
all <- split_multi_id(all, "IL12", c("IL12A", "IL12B"))
all <- split_multi_id(all, "MIC-A/B", c("MIC-A", "MIC-B"))

### 0. Add UniProt ids
all_up <- merge(x = all, y = unique(olp[,c("Assay", "UniProt")]), by.x = "Name", by.y = "Assay")

### 1. Prune NA UniProts
all_up <- all_up[!is.na(all_up$UniProt),]

### 2. Cleaning
### Read cutoffs file, get UniProts that have less than 20% LOD percentage
cutoffs <- read.csv("Data/Perc_below_LOD_per_assay.csv", sep = ";")
cutoff_ups <- unique(cutoffs[cutoffs$percentage_above_LOD < 20, "UniProt"])
### Exclude low LOD UniProt ids from the dataset
all_up <- all_up[!(all_up$UniProt %in% cutoff_ups),]

### Trim by adjusted pvalue, keep records that have at least one adj.p < 0.1
pass <- rowSums(all_up[,c("adj.p.ALZ", "adj.p.FTD", "adj.p.DLB")] < 0.1) > 0
all_up <- all_up[pass,]

### 3. Merging with prior knowledge
###
### 3.1 DisGeNet
message("DisGeNET integration...")
### Add DisGeNet gene mapping to the main table
up2dgn <- read.table("Data/mapa_geneid_4_uniprot_crossref.tsv", header = T)
all_up_dgn <- merge(all_up, up2dgn, by.x = "UniProt", by.y = "UniProtKB", all.x = F, all.y = F) ## 2 entries get lost??
### Read the DisGeNET file, can be replaced with an API call for better reproducibility (as any changes are not captured in .tsv file)
dgn <- read.table("Data/curated_gene_disease_associations.tsv", sep = "\t", quote = "", comment.char = "", header = T)
### Use only the C10 mapping, "Neurological disorders"
dgn_c10 <- dgn[grep("(^C10;)|(^C10$)|(;C10$)|(;C10;)", dgn$diseaseClass),]
### Integrate with the main table
all_up_dgn_c10 <- cbind(all_up_dgn, DGN_hits = sapply(all_up_dgn$GENEID, function(x) sum(dgn_c10$geneId == x)))

### 3.2 PathwayStudio
message("PathwayStudio integration...")
library(xml2)
### RNEF is the native XML format of the PathwayStudio
message("Loading RNEF")
miriade_rnef <- read_xml("Data/PathwayStudio_protein_disease_associations.rnef")
### Reducing the size of the file
invisible(sapply(xml_find_all(miriade_rnef, "//attr[@index]"), xml_remove))

message("Processing RNEF...")
### A convenience function to find specific UniProts and get their urns
find_urns <- function(funiprot, fxml) {
  xml_attr(xml_parent(xml_find_first(fxml, paste0("//node/attr[@name='Swiss-Prot Accession' and @value='",funiprot,"']"))), "urn")
}
### Get urns for all UniProts
miriade_rnef_urns <- sapply(unique(all_up_dgn$UniProt), find_urns, miriade_rnef)

### Convenience function to find how many nodes have a given UniProt urn
find_nids <- function(furn, fxml) {
  if(is.na(furn)) { return(0) }
  sum(sapply(sapply(xml_find_all(fxml, paste0("//node[@urn='",furn,"']")), xml_attr, "local_id"),
         function(x) length(xml_find_all(fxml, paste0("//link[@ref='",x,"']")))))
}
### Find in how many protein-disease interactions a given UniProt is found
miriade_rnef_links <- sapply(miriade_rnef_urns, find_nids, miriade_rnef)

message("Integrating RNEF")
### Integrate with the main table
all_up_dgn_c10_ps <- cbind(all_up_dgn_c10, 
                           PS_hits = sapply(all_up_dgn$UniProt, function(x) miriade_rnef_links[x]))


### 3.3 Disease Maps
message("Disease maps integration...")
library(jsonlite)
### A convenience function to handle API queries
ask_GET <- function(furl, fask) {
  resp <- httr::GET(url = paste0(furl, fask),
                    httr::add_headers('Content-Type' = "application/x-www-form-urlencoded"))
  if(httr::status_code(resp) == 200) {
    return(httr::content(resp, as = "text"))
  }
  return(NULL)
}

### A convenience function to extract a given annotation type from references
extract_anno <- function(refs, type) {
  refs <- refs[sapply(refs, length) > 0]
  if(length(refs) == 0) { return(NULL) }
  ret <- sapply(refs, function(x) x[x$type == type, "resource"])
  ret <- ret[sapply(ret, function(x) ifelse(is.character(x) & length(x) > 0, TRUE, FALSE))]
}

### A convenience function to parse the annotations
extract_annos <- function(model, type) {
  ### For erroneous response, return NULL
  if(is.null(model)) { return(NULL) }
  ### Only elements that have annotations
  these_refs <- model$references[sapply(model$references, length) > 0]
  ### For empty list, return NULL
  if(length(these_refs) == 0) { return(NULL) }
  ### Get HGNC symbols, for elements that have annotations
  res <- extract_anno(these_refs, type)
  return(unlist(res))
}
### A wrapper function: for a given diseae map url, get the default map, 
### get elements of all diagrams and extract UniProt ids
wrap_map <- function(furl, fbase = NULL) {
  ### Get MINERVA elements
  if(is.null(fbase)) {
    ### Get configuration to obtain the latest (default) version
    cfg <- fromJSON(ask_GET(furl, "configuration/"))
    project_id <- cfg$options[cfg$options$type == "DEFAULT_MAP","value"]
    ### The address of the latest (default) build 
    mnv_base <- paste0(furl,"projects/",project_id,"/")
  } else {
    mnv_base <- paste0(furl,"projects/",fbase,"/")
  }
  message(paste0("Asking for diagrams in: ", mnv_base, "models/"))
  ### Get diagrams
  models <- fromJSON(ask_GET(mnv_base, "models/"), flatten = F)
  ### Get elements of diagrams
  model_elements <- lapply(models$idObject, 
                           function(x) fromJSON(ask_GET(paste0(mnv_base,"models/",x,"/"), 
                                                        "bioEntities/elements/?columns=id,name,type,references,complexId"), 
                                                flatten = F))
  names(model_elements) <- models$name
  return(as.vector(sapply(model_elements, extract_annos, "UNIPROT")))
}
### Ask for UniProt ids in each map
pd_ups <- table(unlist(wrap_map("https://pdmap.uni.lu/minerva/api/")))
ag_ups <- table(unlist(wrap_map("https://progeria.uni.lu/minerva/api/")))
ad_ups <- table(unlist(wrap_map("https://minerva-dev.lcsb.uni.lu/minerva/api/","alzpath_8APR")))
### Summarize hits
across_dmaps <- sapply(pd_ups[all_up_dgn_c10_ps$UniProt], function(x) ifelse(is.na(x), 0, x)) + 
                sapply(ag_ups[all_up_dgn_c10_ps$UniProt], function(x) ifelse(is.na(x), 0, x)) + 
                sapply(ad_ups[all_up_dgn_c10_ps$UniProt], function(x) ifelse(is.na(x), 0, x))
### Integrate with the main map
all_up_dgn_c10_ps_dmaps <- cbind(all_up_dgn_c10_ps, DMaps_hits = across_dmaps)

### Integration ends here, afterwards it's sorting and positioning of proteins based on their pvals/hits and disease association

### 4. Rank proteins based on their p-values and sum of prior knowledge hits
###
### A convenience function to simplify pvals, pooling them into fixed ranges;
### for rank-based sorts
floor_pvals <- function(fvec) {
  fvec <- as.numeric(fvec)
  fvec[fvec <= 0.0001] <- 0.0001
  fvec[fvec > 0.0001 & fvec <= 0.001] <- 0.001
  fvec[fvec > 0.001 & fvec <= 0.01] <- 0.01
  fvec[fvec > 0.01 & fvec <= 0.1] <- 0.1
  fvec[fvec > 0.1] <- 1
  fvec
}

### A convenience function to simplify adj pvals, either based on fixed ranges (above) or by their significant value
### for rank-based sorts
treat_pvals <- function(ftab, scale) {
# for scale=0, we use the bins above (coarse sorting) 
  if(scale == 0) {
    ## ADDED "RETURN" HERE (v.3)
    return(
      data.frame(ALZ = floor_pvals(ftab$adj.p.ALZ),
               FTD = floor_pvals(ftab$adj.p.FTD),
               DLB = floor_pvals(ftab$adj.p.DLB))
    )
  } else {
    return(
      data.frame(ALZ = signif(ftab$adj.p.ALZ, scale),
                 FTD = signif(ftab$adj.p.FTD, scale),
                 DLB = signif(ftab$adj.p.DLB, scale))
    )
  }
}

### For better display, columns that can be omitted
## TODO: CREATE A FUNCTION INSTEAD (taking column names and returning indices)
masked_columns <- c(3,4,7,8,11,12,15)

### Position-based sort, across diseases
simple_sort <- function(ftab) {
 # ftab[order(ftab$pos.ALZ + ftab$pos.FTD + ftab$pos.DLB,
#             -(type.convert(ftab$DGN_hits + ftab$PS_hits + ftab$DMaps_hits))),]
  ## added the temp variables 17.12.2020 23:09
  position_vector <- ftab$pos.ALZ + ftab$pos.FTD + ftab$pos.DLB
  hits_vector <- -(type.convert(ftab$DGN_hits + ftab$PS_hits + ftab$DMaps_hits))
  ftab[order(position_vector, hits_vector),]
}

### Adj.p-based sort, ALZ specific
# alz_p_sort <- function(ftab, scale = 1) {
#  sort <- treat_pvals(ftab, scale)
#  ftab[order(sort$ALZ, -sort$FTD*sort$DLB, -(ftab$DGN_hits + ftab$PS_hits + ftab$DMaps_hits)),]
# }
### Searching for more stable sorting:
# alz_p_sort <- function(ftab, scale = 1) {
#   sort <- treat_pvals(ftab, scale)
#   ftab[order(sort$ALZ, -(sort$FTD+sort$DLB), -(ftab$DGN_hits + ftab$PS_hits + ftab$DMaps_hits)),]
# }
alz_p_sort <- function(ftab, scale = 1) {
  sort <- treat_pvals(ftab, scale)
  ftab[order(sort$ALZ / sort$FTD + sort$ALZ / sort$DLB, sort$ALZ, -(ftab$DGN_hits + ftab$PS_hits + ftab$DMaps_hits)),]
}
# improve discrimination between diseases even further
# alz_p_sort <- function(ftab, scale = 1) {
#   sort <- treat_pvals(ftab, scale)
#   ftab[order((sort$FTD / sort$DLB + sort$DLB / sort$FTD)*(sort$ALZ / sort$FTD + sort$ALZ / sort$DLB), sort$ALZ, -(ftab$DGN_hits + ftab$PS_hits + ftab$DMaps_hits)),]
# }

### Adj.p-based sort, FTD specific
# ftd_p_sort <- function(ftab, scale = 1) {
#  sort <- treat_pvals(ftab, scale)
#  ftab[order(sort$FTD, -sort$ALZ*sort$DLB, -(ftab$DGN_hits + ftab$PS_hits + ftab$DMaps_hits)),]
# }
### Searching for more stable sorting:
# ftd_p_sort <- function(ftab, scale = 1) {
#   sort <- treat_pvals(ftab, scale)
#   ftab[order(sort$FTD, -(sort$ALZ+sort$DLB), -(ftab$DGN_hits + ftab$PS_hits + ftab$DMaps_hits)),]
# }
ftd_p_sort <- function(ftab, scale = 1) {
  sort <- treat_pvals(ftab, scale)
  ftab[order(sort$FTD / sort$ALZ + sort$FTD / sort$DLB, sort$FTD, -(ftab$DGN_hits + ftab$PS_hits + ftab$DMaps_hits)),]
}
# improve discrimination between diseases even further
# ftd_p_sort <- function(ftab, scale = 1) {
#   sort <- treat_pvals(ftab, scale)
#   ftab[order((sort$ALZ / sort$DLB + sort$DLB / sort$ALZ)*(sort$FTD / sort$ALZ + sort$FTD / sort$DLB), sort$FTD, -(ftab$DGN_hits + ftab$PS_hits + ftab$DMaps_hits)),]
# }

### Adj.p-based sort, DLB specific
# dlb_p_sort <- function(ftab, scale = 1) {
#  sort <- treat_pvals(ftab, scale)
#  ftab[order(sort$DLB, -sort$FTD*sort$ALZ, -(ftab$DGN_hits + ftab$PS_hits + ftab$DMaps_hits)),]
# }
### Searching for more stable sorting:
# dlb_p_sort <- function(ftab, scale = 1) {
#   sort <- treat_pvals(ftab, scale)
#   ftab[order(sort$DLB, -(sort$FTD+sort$ALZ), -(ftab$DGN_hits + ftab$PS_hits + ftab$DMaps_hits)),]
# }
dlb_p_sort <- function(ftab, scale = 1) {
  sort <- treat_pvals(ftab, scale)
  ftab[order(sort$DLB / sort$ALZ + sort$DLB / sort$FTD, sort$DLB, -(ftab$DGN_hits + ftab$PS_hits + ftab$DMaps_hits)),]
}
# improve discrimination between diseases even further
# dlb_p_sort <- function(ftab, scale = 1) {
#   sort <- treat_pvals(ftab, scale)
#   ftab[order((sort$ALZ / sort$FTD + sort$FTD / sort$ALZ)*(sort$DLB / sort$ALZ + sort$DLB / sort$FTD), sort$DLB, -(ftab$DGN_hits + ftab$PS_hits + ftab$DMaps_hits)),]
# }

### Display results
### Simple sort
head(simple_sort(all_up_dgn_c10_ps_dmaps[,-masked_columns]), n = 15)
### Adj.p sorts
head(alz_p_sort(all_up_dgn_c10_ps_dmaps[,-masked_columns], scale = 1), n = 15)
head(ftd_p_sort(all_up_dgn_c10_ps_dmaps[,-masked_columns], scale = 1), n = 15)
head(dlb_p_sort(all_up_dgn_c10_ps_dmaps[,-masked_columns], scale = 1), n = 15)

### Write results to file

write.table(simple_sort(all_up_dgn_c10_ps_dmaps),
            file = "Results/across_diseases_sorted.txt", sep = "\t", quote = F, row.names = F, col.names = T)

write.table(alz_p_sort(all_up_dgn_c10_ps_dmaps, scale = 1),
            file = "Results/ALZ_sorted_strict.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(alz_p_sort(all_up_dgn_c10_ps_dmaps, scale = 0),
            file = "Results/ALZ_sorted_relaxed.txt", sep = "\t", quote = F, row.names = F, col.names = T)

write.table(ftd_p_sort(all_up_dgn_c10_ps_dmaps, scale = 1),
            file = "Results/FTD_sorted_strict.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(ftd_p_sort(all_up_dgn_c10_ps_dmaps, scale = 0),
            file = "Results/FTD_sorted_relaxed.txt", sep = "\t", quote = F, row.names = F, col.names = T)

write.table(dlb_p_sort(all_up_dgn_c10_ps_dmaps, scale = 1),
            file = "Results/DLB_sorted_strict.txt", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(dlb_p_sort(all_up_dgn_c10_ps_dmaps, scale = 0),
            file = "Results/DLB_sorted_relaxed.txt", sep = "\t", quote = F, row.names = F, col.names = T)

### Generate Excel files from the sorting results
library(openxlsx)
### Defaults
options("openxlsx.borderColour" = "#4F80BD") 
options("openxlsx.border" = "TopBottom")

### Style for header
headerStyle <- createStyle(fontSize = 14, fontColour = "#FFFFFF", halign = "center", fgFill = "#4F81BD")
### Style for body
bodyStyleHL <- createStyle(textDecoration = "underline", fontColour = "#1400FB")
bodyStyleNUM <- createStyle(numFmt = "0.000", halign = "center")
bodyStyleSCI <- createStyle(numFmt = "SCIENTIFIC",halign = "center")
bodyStylePOS <- createStyle(fgFill = "#DAEEF3")

### Nice arrangement of sorted tables, adding hyperlinks to UniProt entries
add_sorted_data <- function(fwb, ftab, fname) {
  ## Add a worksheets
  addWorksheet(fwb, fname, gridLines = FALSE)
  this_sheet <- length(names(fwb))
  ## write data to worksheet
  writeData(fwb, sheet = this_sheet, ftab, rowNames = FALSE)
  ## Uniprot hyperlinks
  links <- paste0("https://www.uniprot.org/uniprot/", ftab$UniProt)
  names(links) <- ftab$UniProt
  class(links) <- "hyperlink"
  writeData(fwb, sheet = this_sheet, x = links, startCol = 1, startRow = 2)
  
  ## Add a style to the column headers 
  addStyle(fwb, sheet = this_sheet, headerStyle, rows = 1, cols = 1:ncol(ftab), gridExpand = TRUE)
  
  ## Add a style to the table body
  ## Hyperlinks (UniProt)
  addStyle(fwb, sheet = this_sheet, bodyStyleHL, rows = 2:nrow(ftab), cols = 1, gridExpand = TRUE)
  ## Numeric (Effect size)
  addStyle(fwb, sheet = this_sheet, bodyStyleNUM, rows = 2:nrow(ftab), cols = which(startsWith(colnames(ftab), "Eff.")), gridExpand = TRUE) 
  ## Scientific (pvalues)
  addStyle(fwb, sheet = this_sheet, bodyStyleSCI, rows = 2:nrow(ftab), cols = which(startsWith(colnames(ftab), "p.val") | startsWith(colnames(ftab), "adj.p")), gridExpand = TRUE) 
  ## With background (position)
  addStyle(fwb, sheet = this_sheet, bodyStylePOS, rows = 2:nrow(ftab), cols = which(startsWith(colnames(ftab), "pos.")), gridExpand = TRUE)
  
  ## Column width, row height
  setColWidths(fwb, 1, cols = 1:ncol(ftab), widths = rep(12, ncol(ftab)))
  setRowHeights(fwb, 1, rows = 1:nrow(ftab), heights = rep(20, nrow(ftab)))
  
  ## Modified workbook
  return(fwb)
}

### Skip the last column
all_excel <- all_up_dgn_c10_ps_dmaps[,-15]

# ####-------------------------------------------------------
# ### Choice of biomarker
# # 1st criterium: low p-value
# # 2nd criterium: lots of hits/ known + in lots of pathways
# # 3rd criterium: dicriminating between diseases
# 
# # initialise normed all matrix
# all_norm <- all_up_dgn
# 
# # 1. Normalise adjusted p-values to range (0,1)
# all_norm$adj.p.ALZ <- sapply(data.frame(all_norm$adj.p.ALZ), rescale, to = c(0, 1))
# all_norm$adj.p.FTD <- sapply(data.frame(all_norm$adj.p.FTD), rescale, to = c(0, 1))
# all_norm$adj.p.DLB <- sapply(data.frame(all_norm$adj.p.DLB), rescale, to = c(0, 1))
# 
# ## normalize(all_norm$adj.p.ALZ, method = "standardize", range = c(0, 1), margin = 2)
# ##nf_ALZ <- sum(all_norm$adj.p.ALZ) 
# ##nf_FTD <- sum(all_norm$adj.p.FTD)
# ##nf_DLB <- sum(all_norm$adj.p.DLB)
# ##apply(all_norm$adj.p.ALZ, 2, function(x) {(x / nf_ALZ)})
# 
# 
# # 2. Calculate distinctiveness score using p-values
# # TODO: might want to use some other (exponential function) instead of mean
# all_norm$distinct <- abs(all_norm$adj.p.ALZ - all_norm$adj.p.FTD)/2  + abs(all_norm$adj.p.ALZ - all_norm$adj.p.DLB)/2 + abs(all_norm$adj.p.FTD - all_norm$adj.p.DLB)/2
# # normalise distinctiveness score
# all_norm$distinct_norm <- sapply(data.frame(all_norm$distinct), rescale, to = c(0, 1))
# 
# # 3. Calculate hits count as average normalised hit counts from all 3 databases
# # TODO: CREATE DATFRAME OUT OF TABLES, load it in workspace
# # Read a txt files
# hits_data <- read.delim("ALZ_sorted_strict.txt")
# #hits_data <- read.delim("across_diseases_sorted.txt")
# all_norm$DGN_norm <- sapply(data.frame(hits_data$DGN_hits), rescale, to = c(0, 1))
# all_norm$PS_norm <- sapply(data.frame(hits_data$PS_hits), rescale, to = c(0, 1))
# all_norm$DMaps_norm <- sapply(data.frame(hits_data$DMaps_hits), rescale, to = c(0, 1))
# # weight the hits accordingly
# hwDGN <- 1
# hwPS <- 1
# hwDMaps <- 1
# all_norm$hitstotal_norm <- (hwDGN * all_norm$DGN_norm + hwPS * all_norm$PS_norm + hwDMaps * all_norm$DMaps_norm) / (hwDGN + hwPS + hwDMaps)
# 
# ## calculate biomarker score bm_score: we want to minimise p-value, maximise the distinctiveness score and maximise the hits count
# # choose weights/ covariates for bm_score
# w_pv <- 1
# w_distinct <- 1
# w_hits <- 1
# 
# bm_score_ALZ <- data.frame(all_norm$UniProt)
# bm_score_ALZ$Name <- all_norm$Name
# bm_score_ALZ$bm <- w_pv * (1 - all_norm$adj.p.ALZ) + w_distinct * all_norm$distinct_norm + w_hits * all_norm$hitstotal_norm
# # bm_score_FTD <- w_pv * (1 - all_norm$adj.p.FTD) + w_distinct * all_norm$distinct_norm + w_hits * all_norm$hitstotal_norm
# # bm_score_DLB <- w_pv * (1 - all_norm$adj.p.DLB) + w_distinct * all_norm$distinct_norm + w_hits * all_norm$hitstotal_norm
# # NOW PICK BIOMARKERS WITH HIGHEST SCORE
# 
# ####---------------------------------------
##### TODO: ADJUST subfolders
## Create a new workbook
wb <- createWorkbook("Results/Olink MIRIADE biomarker priorisation")
modifyBaseFont(wb, fontSize = 13, fontColour = "black", fontName = "Calibri")
### Add sheets with different tables
wb <- add_sorted_data(wb, simple_sort(all_excel), "Cross-disease")
wb <- add_sorted_data(wb, alz_p_sort(all_excel, scale = 1), "ALZ-specific, strict")
wb <- add_sorted_data(wb, alz_p_sort(all_excel, scale = 0), "ALZ-specific, relaxed")
wb <- add_sorted_data(wb, ftd_p_sort(all_excel, scale = 1), "FTD-specific, strict")
wb <- add_sorted_data(wb, ftd_p_sort(all_excel, scale = 0), "FTD-specific, relaxed")
wb <- add_sorted_data(wb, dlb_p_sort(all_excel, scale = 1), "DLB-specific, strict")
# FIXED MISTAKE scale=0
wb <- add_sorted_data(wb, dlb_p_sort(all_excel, scale = 0), "DLB-specific, relaxed")

#saveWorkbook(wb, file = "Results/MIRIADE_Olink_sorted_biomarkers.xlsx", overwrite = TRUE)
saveWorkbook(wb, file = "Results/MIRIADE_Olink_sorted_biomarkers_v2.xlsx", overwrite = TRUE)

