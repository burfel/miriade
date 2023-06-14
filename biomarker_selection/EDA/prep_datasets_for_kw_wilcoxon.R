library(here)
library(dplyr)

# This loads the melded emif and olink dataframes to be used later in the file
# It also calls to set the dataset root directory
source(here("biomarker_selection", "EDA", "meld_all_fractured_datasets.R"))
source(here("functions", "auxiliary_statistical_functions.R"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
## Processing the datasets ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
### KTH ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
kth <- readxl::read_xlsx(file.path(datasets_root_directory,
                                   "KTH/KTH AD dataset for MIRIADE biomarker selection sample info and results.xlsx"))
### Change the names into HGNC symbols (prefixes)
colnames(kth) <- sapply(colnames(kth), take_only_pre_underscore_substring)

breaks_and_labels <- prepare_breaks_and_labels_by_quantiles(kth$Age, 3)
### Read the table
kth <- kth %>%
  dplyr::select(-Class) %>%
  dplyr::mutate(Gender = case_when(Gender == 'M' ~ 'm', Gender == 'F' ~ 'f', TRUE ~ Gender)) %>%
  convert_age_column_to_age_group_column(sym("Age"), sym("Age_Group"),
                                         cutting_breaks = breaks_and_labels$Breaks,
                                         cutting_labels = breaks_and_labels$Labels) %>%
  dplyr::select(Diagnosis, Age_Group, everything())

mutated_kth <- kth %>%
  dplyr::filter(Diagnosis == "AD" | Diagnosis == "Control") %>%
  dplyr::mutate(Diagnosis_Age_Group = paste(Diagnosis, Age_Group)) %>%
  dplyr::select(Diagnosis_Age_Group, Age_Group, everything())

print("The number of patients per age group (in KTH) is:")
kth %>% dplyr::group_by(Age_Group) %>% dplyr::summarise(n())

write.table(kth, file = file.path(datasets_root_directory, "KTH/KTH-raw-data.csv"),
            sep = "\t", quote = F, row.names = F, col.names = T)

### Melt the dataframe, so there's only one readout variable
melted_kth <-
  reshape2::melt(kth, id = 1:9,
                 variable.name = "HGNC_Symbol", value.name = "Value") %>%
  dplyr::select(Diagnosis, Age_Group, Gender, HGNC_Symbol, Value)

#### Adjusting to have Diagnosis + age group combined column                ####

supercharged_kth_melt <-
  melted_kth %>%
  dplyr::filter(Diagnosis == "AD" | Diagnosis == "Control") %>%
  dplyr::mutate(Diagnosis_Age_Group = paste(Diagnosis, Age_Group)) %>%
  dplyr::select(Diagnosis_Age_Group, Gender, HGNC_Symbol, Value)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
### Emif ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

emif$Age <- round(emif$Age) # Must round, otherwise the age groups are classified incorrectly

breaks_and_labels <- prepare_breaks_and_labels_by_quantiles(emif$Age, 3)

# Uniform diagnosis names
emif <- emif %>%
  dplyr::mutate(Diagnosis = case_when(Diagnosis == "NL" ~ "Control", TRUE ~ Diagnosis)) %>%
  convert_age_column_to_age_group_column(sym("Age"), sym("Age_Group"),
                                         cutting_breaks = breaks_and_labels$Breaks,
                                         cutting_labels = breaks_and_labels$Labels) %>%
  dplyr::select(SubjectId, Age_Group, everything())
# Reomve two rows that have NA as Diagnosis
emif <- emif[complete.cases(emif$Diagnosis),]

write.table(emif, file = file.path(datasets_root_directory, "EMIF-AD MBD Study/EMIF-raw-data.csv"),
            sep = "\t", quote = F, row.names = F, col.names = T)


source(here("functions", "mapping_functions.R"))
emif_proteins <- data.frame(UniProt = colnames(emif)[6:ncol(emif)]) %>%
  map_uniprot_to_hgnc_and_cbind("Name") %>%
  na.omit() %>%
  dplyr::select(UniProt, Name)
# use sapply to create a named vector of the hgnc and their new uniprots
protein_renaming_vector <- sapply(emif_proteins$Name, function(x) {
  emif_proteins$UniProt[which(emif_proteins$Name == x)]
})

mutated_emif <- emif %>%
  dplyr::filter(Diagnosis == "AD" | Diagnosis == "Control") %>%
  dplyr::mutate(Diagnosis_Age_Group = paste(Diagnosis, Age_Group)) %>%
  dplyr::select(Diagnosis_Age_Group, Age_Group, everything()) %>%
  dplyr::rename(!!!protein_renaming_vector)

print("The number of patients per age group (in Emif) is:")
mutated_emif %>% dplyr::group_by(Age_Group) %>% dplyr::summarise(n())

### Melt the emif dataframe, so there's only one readout variable
melted_emif <- emif %>%
  dplyr::rename(!!!protein_renaming_vector) %>%
  reshape2::melt(id = 1:6,
                 variable.name = "HGNC_Symbol", value.name = "Value") %>%
  dplyr::select(Diagnosis, Age_Group, Gender, HGNC_Symbol, Value)

supercharged_emif_melt <-
  melted_emif %>%
  dplyr::filter(Diagnosis == "AD" | Diagnosis == "Control") %>%
  dplyr::mutate(Diagnosis_Age_Group = paste(Diagnosis, Age_Group)) %>%
  dplyr::select(Diagnosis_Age_Group, Gender, HGNC_Symbol, Value)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
### Olink ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

breaks_and_labels <- prepare_breaks_and_labels_by_quantiles(olink$age_gr, 3)
# Uniform diagnosis names
olink <- olink %>%
  dplyr::mutate(Diagnosis = case_when(dx == "CN" ~ "Control", 
                                      dx == "AD dementia" ~ "AD", 
                                      TRUE ~ dx)) %>%
  dplyr::rename(Age = "age_gr", Gender = "sex") %>%
  convert_age_column_to_age_group_column(sym("Age"), sym("Age_Group"),
                                         cutting_breaks = breaks_and_labels$Breaks,
                                         cutting_labels = breaks_and_labels$Labels) %>%
  dplyr::select(SampleId, Diagnosis, Age_Group, everything(), -dx)

write.table(olink, file = file.path(datasets_root_directory, "Olink/OLINK basic data files/Olink-raw-data.csv"),
            sep = "\t", quote = F, row.names = F, col.names = T)

mutated_olink <- olink %>%
  dplyr::filter(Diagnosis == "AD" | Diagnosis == "Control") %>%
  dplyr::mutate(Diagnosis_Age_Group = paste(Diagnosis, Age_Group)) %>%
  dplyr::select(Diagnosis_Age_Group, Diagnosis, Age_Group, everything())

print("The number of patients per age group (in Olink) is:")
mutated_olink %>% dplyr::group_by(Age_Group) %>% dplyr::summarise(n())

### Melt the olink dataframe, so there's only one readout variable
melted_olink <-
  reshape2::melt(olink, id = 1:5,
                 variable.name = "HGNC_Symbol", value.name = "Value") %>%
  dplyr::select(Diagnosis, Age_Group, Gender, HGNC_Symbol, Value)

supercharged_olink_melt <-
  melted_olink %>%
  dplyr::filter(Diagnosis == "AD" | Diagnosis == "Control") %>%
  dplyr::mutate(Diagnosis_Age_Group = paste(Diagnosis, Age_Group)) %>%
  dplyr::select(Diagnosis_Age_Group, Gender, HGNC_Symbol, Value)
