library(tidyverse)
library(here)
library(devtools)
library(dplyr)
install_version("ggplot2", version = "3.3.6")

source(here("functions", "data_processing_functions.R"))

datasets_root_directory <- define_datasets_root()

keep_first_uniprot <- function(string) {
  if ("," %in% string) {
    uniprots <- strsplit(string, ",")[[1]]
    uniprot1 <- uniprots[1]
  } else {
    uniprot1 <- string
  }
  return(uniprot1)
}

###---- MACRON-2018 DATASET -----###
# load dataset
Macron2018 <- read.csv(file.path(datasets_root_directory, "background-CSF-proteins/2018-Macron.csv"), sep=";", skip=1, header=TRUE)
## rename column
#colnames(Macron2018)[colnames(Macron2018)=="Protein.Accession.Number"] <- "UniProt"
# rename and select UniProt column
Macron2018 <- Macron2018 %>% dplyr::select("Protein.Accession.Number") %>%
  dplyr::rename(UniProt = "Protein.Accession.Number")

# count number of peptides per UniProt-ID
Macron2018 <- Macron2018 %>% dplyr::count(UniProt)
# drop duplicate entries
Macron2018 <- unique(Macron2018, by = "UniProt")
# keep first entry for proteins with more than one associated Uniprot ID
# Macron2018$UniProt <- keep_first_uniprot(Macron2018$UniProt)
Macron2018$UniProt <- sapply(strsplit(as.character(Macron2018$UniProt), ","), function(x) x[1])

# Number of unique proteins with at least 1 peptide
nrow(subset(Macron2018, Macron2018$n > 0))

# Number of unique proteins with at least 2 peptides
nrow(subset(Macron2018, Macron2018$n > 1))


###------- MACRON-2020 DATASET -------###
#### load dataset
Macron2020 <- read.csv(file.path(datasets_root_directory, "background-CSF-proteins/2020-Macron.csv"), sep=";", skip=1, header=TRUE)
Macron2020 <- Macron2020 %>% dplyr::select("Protein.Accession.Number") %>%
  dplyr::rename(UniProt = "Protein.Accession.Number")

# count number of peptides per UniProt-ID
Macron2020 <- Macron2020 %>% dplyr::count(UniProt)
# drop duplicate entries
Macron2020 <- unique(Macron2020, by = "UniProt")
# omit NA entries
Macron2020 <- na.omit(Macron2020)
# keep first entry for proteins with more than one associated Uniprot ID
Macron2020$UniProt <- sapply(strsplit(as.character(Macron2020$UniProt), ","), function(x) x[1])

# Number of unique proteins with at least 1 peptide
nrow(subset(Macron2020, Macron2020$n > 0))

# Number of unique proteins with at least 2 peptides
nrow(subset(Macron2020, Macron2020$n > 1))


###---------- ZHANG-2015 DATASET -------###
# load dataset
Zhang2015 <- read.csv(file.path(datasets_root_directory, "background-CSF-proteins/2015-Zhang.csv"), sep=";", header=TRUE)

# select the columns of interest
Zhang2015 <- Zhang2015 %>%
  dplyr::select("Accession", "Flow.through.Proteins", "Original.Proteins", "Bound.Proteins")
# omit NA entries
Zhang2015 <- na.omit(Zhang2015)
# find the maximum value for each row and save it as a new column
Zhang2015 <- Zhang2015 %>%
  dplyr::mutate("n" = pmax(Flow.through.Proteins,Original.Proteins, Bound.Proteins))
# drop the unwanted columns and rename the remaining column
Zhang2015 <- Zhang2015 %>%
  dplyr::select(Accession, n) %>%
  dplyr::rename(UniProt = Accession)

# Number of unique proteins with at least 1 peptide
nrow(subset(Zhang2015, Zhang2015$n > 0))

# Number of unique proteins with at least 2 peptides
nrow(subset(Zhang2015, Zhang2015$n > 1))
