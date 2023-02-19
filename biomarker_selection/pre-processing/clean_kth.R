##################################################
## Project: MIRIADE
## Script purpose: Caclculate significant biomarkers from the KTH dataset
## Date: 11.10.2021, 15.09.2022
## Authors: Marek Ostaszewski, Felicia Burtscher
##################################################

library(readxl)
library(reshape)
library(dplyr)
library(stats)
library(purrr)
library(here)

source(here("functions", "data_processing_functions.R"))

datasets_root_directory <- define_datasets_root()

### Read the table
# tab <- readxl::read_xlsx("_notgit/KTH/KTH AD dataset for MIRIADE biomarker selection sample info and results.xlsx")
tab <- readxl::read_xlsx(file.path(datasets_root_directory, "KTH/KTH AD dataset for MIRIADE biomarker selection sample info and results.xlsx"))

### Change the names into HGNC symbols (prefixes)
colnames(tab) <- sapply(colnames(tab), function(n) strsplit(n, "_")[[1]][1])

### Melt the dataframe, so there's only one readout variable
tab2 <- reshape2::melt(tab, id = 1:9, variable.name = "HGNC_Symbol", value.name = "median_ab_readout")

####
#### Kruskal-Wallis test for significant biomarkers
####

### Thanks to Veronica Codoni and Armin Rauschenberger for the tip
### We begin with Kruskal-Wallis test to see if there is any dominance in the group
### No distribution assumption needed

### A helper function to calculate Kruskal-Wallis for all the Diagnoses in the table, given the hgnc symbol
kw_for_diagnoses <- function(hgnc, table) {
  ### This *magic* one-liner filters by 'hgnc', splits 'median_ab_readout' by 'Diagnosis'
  kw_tibbles <- tibble::as_tibble(table) %>%
    dplyr::filter(HGNC_Symbol == hgnc) %>%
    dplyr::select(Diagnosis, median_ab_readout) %>%
    dplyr::group_split(Diagnosis, .keep = FALSE)
  ### We need to unlist the tibbles from the 'group_split' to run the kruskal.test
  stats::kruskal.test(lapply(kw_tibbles, unlist))$p.value
}

### Run Kruskal-Wallis test for all biomarkers
kw_per_hgnc <- sapply(unique(tab2$HGNC_Symbol), kw_for_diagnoses, tab2)
### Adjust for multiple hypothesis testing
kw_biomarkers <- as.character(unique(tab2$HGNC_Symbol)[p.adjust(kw_per_hgnc, method = "BH") < 0.1])
### Keep only the significant entries
sub_kw <- dplyr::filter(tab2, HGNC_Symbol %in% kw_biomarkers)

####
#### Wilcoxon pairwise comparisons for biomarkers passing the KW test
####

### A helper function to calculate pairwise Wilcoxon tests for the given list of
### group pairs. This function is specific for the variable/value names assigned in reshape2::melt, above
test_pairwise_wilcoxon <- function(df, groups) {
  return(sapply(groups, function(g) with(df,
                                         stats::wilcox.test(median_ab_readout[Diagnosis == g[1]],
                                                            median_ab_readout[Diagnosis == g[2]],
                                                            exact = F)$p.value)))
}

### Define the pairs automatically based on the "Diagnosis" values
pairs <- unique(tab2$Diagnosis) %>%
  combn(2, simplify = F) %>%
  purrr::set_names(purrr::map_chr(., ~ paste(., collapse = "_vs_")))


### Calculate Wilcoxon pairwise tests
sub_kw_wilcox <- sapply(unique(sub_kw$HGNC_Symbol),
                        function(hgnc) test_pairwise_wilcoxon(sub_kw[sub_kw$HGNC_Symbol == hgnc,], pairs))
colnames(sub_kw_wilcox) <- unique(sub_kw$HGNC_Symbol)

### Adjust for measurements
sub_kw_adj_wilcox <- matrix(p.adjust(sub_kw_wilcox, method = "BH"), nrow = 6)
rownames(sub_kw_adj_wilcox) <- rownames(sub_kw_wilcox)
colnames(sub_kw_adj_wilcox) <- colnames(sub_kw_wilcox)

### Select only Control vs AD
### Construct the final table
final <- data.frame(Name = colnames(sub_kw_wilcox),
                    p.val = sub_kw_wilcox["Control_vs_AD",],
                    adj.p = sub_kw_adj_wilcox["Control_vs_AD",],
                    disease = "ALZ")
### Use biomaRt to provide UniProt ids, more stable than curl call to HGNC API
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install('biomaRt')
library('biomaRt')
if(!("ensembl_h" %in% ls())) {
  ensembl_h = useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
}

id_mapping <- getBM(attributes = c("hgnc_symbol", "uniprotswissprot"),
                    filters = "hgnc_symbol",
                    values = final$Name,
                    mart = ensembl_h)

id_mapping <- dplyr::filter(id_mapping, uniprotswissprot != "")

### Add UniProt ids based on the HGNCs in the Name
final <- cbind(UniProt = sapply(final$Name, function(x) with(id_mapping, uniprotswissprot[hgnc_symbol == x][1])),
               final)

# write.table(final, file = "_notgit/cleaned_KTH_entries.tsv",
#             sep = "\t", quote = F, row.names = F, col.names = T)
write.table(final, file = file.path(datasets_root_directory, "KTH/cleaned_KTH_entries.tsv"),
            sep = "\t", quote = F, row.names = F, col.names = T)

