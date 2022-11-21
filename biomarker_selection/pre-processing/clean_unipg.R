##################################################
## Project: MIRIADE
## Script purpose: Caclculate significant biomarkers from the UNIPG dataset
## Date: 11.10.2021, 15.09.2022
## Authors: Marek Ostaszewski, Felicia Burtscher
##################################################

library(readxl)
library(reshape)
library(dplyr)
library(stats)
library(purrr)

### Read the table
# tab <- readxl::read_xlsx("_notgit/UNIPG/Metabolomic data.xlsx") %>%
#   dplyr::rename(Sample = ...1) %>%
#   dplyr::mutate(Diagnosis = sapply(Sample, function(x) tail(strsplit(x, split = "_")[[1]],1)),
#                 .before = Center)
tab <- readxl::read_xlsx("/Users/felicia.burtscher/Documents/UL/DATASETS/UNIPG/Metabolomic data.xlsx") %>%
  dplyr::rename(Sample = ...1) %>%
  dplyr::mutate(Diagnosis = sapply(Sample, function(x) tail(strsplit(x, split = "_")[[1]],1)),
                .before = Center)

### Melt the dataframe, so there's only one readout variable
tab2 <- reshape2::melt(tab, id = 1:3, variable.name = "Metabolite", value.name = "readout")

sapply(tab$Sample, function(x) tail(strsplit(x, split = "_")[[1]],1))

####
#### Kruskal-Wallis test for significant biomarkers
####

### Thanks to Veronica Codoni and Armin Rauschenberger for the tip
### We begin with Kruskal-Wallis test to see if there is any dominance in the group
### No distribution assumption needed

### A helper function to calculate Kruskal-Wallis for all the Diagnoses in the table, given the hgnc symbol
kw_for_diagnoses <- function(metabolite, table) {
  ### This *magic* one-liner filters by 'metabolite', splits 'readout' by 'Diagnosis'
  kw_tibbles <- tibble::as_tibble(table) %>% 
    dplyr::filter(Metabolite == metabolite) %>%                          
    dplyr::select(Diagnosis, readout) %>% 
    dplyr::group_split(Diagnosis, .keep = FALSE)
  ### We need to unlist the tibbles from the 'group_split' to run the kruskal.test
  stats::kruskal.test(lapply(kw_tibbles, unlist))$p.value
}

### Run Kruskal-Wallis test for all metabolites in 'table'
run_kw <- function(table) {
  kw_per_metabolite <- sapply(unique(table$Metabolite), kw_for_diagnoses, table)
  ### Adjust for multiple hypothesis testing
  kw_biomarkers <- as.character(unique(table$Metabolite)[p.adjust(kw_per_metabolite, method = "BH") < 0.1])
  ### Keep only the significant entries
  return(dplyr::filter(table, Metabolite %in% kw_biomarkers))
}

### Run Kruskal-Wallis test for sample comparisons
run_kw2 <- function(table) {
  kw_per_metabolite <- sapply(unique(table$Metabolite), kw_for_diagnoses, table)
  ### Adjust for multiple hypothesis testing
  kw_biomarkers <- as.character(unique(table$Metabolite)[p.adjust(kw_per_metabolite, method = "BH") < 0.1])
  ### Keep only the significant entries
  return(dplyr::filter(table, Metabolite %in% kw_biomarkers))
}

per_bms <- run_kw(dplyr::filter(tab2, Center == "Perugia"))
ams_bms <- run_kw(dplyr::filter(tab2, Center == "Amsterdam"))

# library(ggplot2)
# ggplot(dplyr::filter(tab2, Center == "Perugia" & Metabolite %in% c("ascorbate")), 
#        aes(Metabolite, readout, colour = Diagnosis)) + 
#   geom_boxplot() + 
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

####
#### Wilcoxon pairwise comparisons for biomarkers passing the KW test
####

### A helper function to calculate pairwise Wilcoxon tests for the given list of
### group pairs. This function is specific for the variable/value names assigned in reshape2::melt, above
test_pairwise_wilcoxon <- function(df, groups) {
  return(sapply(groups, function(g) with(df, 
                                         stats::wilcox.test(readout[Diagnosis == g[1]],
                                                            readout[Diagnosis == g[2]], 
                                                            exact = F)$p.value)))
}

run_wilcoxon <- function(table) {
  ### Define the pairs automatically based on the "Diagnosis" values
  pairs <- unique(table$Diagnosis) %>%
    combn(2, simplify = F) %>%
    purrr::set_names(purrr::map_chr(., ~ paste(., collapse = "_vs_")))
  print(pairs)
  ### Calculate Wilcoxon pairwise tests
  sub_kw_wilcox <- sapply(unique(table$Metabolite), 
                          function(met) test_pairwise_wilcoxon(table[table$Metabolite == met,], pairs))
  # colnames(sub_kw_wilcox) <- unique(table$Metabolite)
  
  ### Adjust for measurements
  sub_kw_adj_wilcox <- matrix(p.adjust(sub_kw_wilcox, method = "BH"), nrow = length(pairs))
  rownames(sub_kw_adj_wilcox) <- rownames(sub_kw_wilcox)
  colnames(sub_kw_adj_wilcox) <- colnames(sub_kw_wilcox)
  return(list(t1 = sub_kw_wilcox, t2 = sub_kw_adj_wilcox))
}


run_wilcoxon(per_bms)



### Select only Control vs AD
### Construct the final table
sub_final <- data.frame(Name = colnames(sub_kw_wilcox), 
                        p.val.ALZ = sub_kw_wilcox["Control_vs_AD",],
                        adj.p.ALZ = sub_kw_adj_wilcox["Control_vs_AD",])
### A helper function to get UniProt id for a given HGNC
HGNC_to_UniProt <- function(hgnc) {
  require(RCurl)
  require(jsonlite)
  hgnc_entry <- RCurl::httpGET(paste0("http://rest.genenames.org/fetch/symbol/", hgnc), 
                               httpheader = c(Accept = "application/json"))
  ups <- jsonlite::fromJSON(hgnc_entry)$response$docs$uniprot_ids
  return(unlist(ups))
}
### Add UniProt ids based on the HGNCs in the Name
final <- cbind(UniProt = sapply(sub_final$Name, HGNC_to_UniProt), sub_final)

# write.table(final, file = "_notgit/cleaned_KTH_entries.tsv", 
#             sep = "\t", quote = F, row.names = F, col.names = T)

write.table(final, file = "/Users/felicia.burtscher/Documents/UL/DATASETS/UNIPG/cleaned_UNIPG_entries.tsv", 
            sep = "\t", quote = F, row.names = F, col.names = T)


