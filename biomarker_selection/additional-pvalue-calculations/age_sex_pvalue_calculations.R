library(here)
library(dplyr)

source(here("functions", "data_processing_functions.R"))
source(here("functions", "auxiliary_statistical_functions.R"))

datasets_root_directory <- define_datasets_root()

### Read the table
kth <- readxl::read_xlsx(file.path(datasets_root_directory,
                                   "KTH/KTH AD dataset for MIRIADE biomarker selection sample info and results.xlsx"))

### Change the names into HGNC symbols (prefixes)
colnames(kth) <- sapply(colnames(kth), take_only_pre_underscore_substring)

### Melt the dataframe, so there's only one readout variable
melted_kth <-
  reshape2::melt(kth, id = 1:9,
                 variable.name = "HGNC_Symbol", value.name = "median_ab_readout") %>%
  #dplyr::filter(Diagnosis == "AD") %>%
  convert_age_column_to_age_group_column(sym("Age"), sym("Age_Group")) %>%
  dplyr::select(Diagnosis, Age_Group, Gender, HGNC_Symbol, median_ab_readout)


### Run Kruskal-Wallis test for all biomarkers grouped by gender
kw_per_hgnc_gender_grouped <-
  sapply(
    unique(melted_kth$HGNC_Symbol),
    kruskal_wallis_test_for_sample_groups,
    melted_kth,
    sym("HGNC_Symbol"),
    sym("Gender"),
    sym("median_ab_readout"))

significant_kw_per_hgnc_gender_grouped <-
  keep_only_significant_entries(
    melted_kth,
    differentiating_feature_symbol = sym("HGNC_Symbol"),
    kw_per_differentiating_feature = kw_per_hgnc_gender_grouped,
    significance_limit = 0.1)

### Age_Group
kw_per_hgnc_age_group_grouped <-
  sapply(
    unique(melted_kth$HGNC_Symbol),
    kruskal_wallis_test_for_sample_groups,
    melted_kth,
    sym("HGNC_Symbol"),
    sym("Age_Group"),
    sym("median_ab_readout"))

significant_kw_per_hgnc_age_group_grouped <-
  keep_only_significant_entries(
    melted_kth,
    differentiating_feature_symbol = sym("HGNC_Symbol"),
    kw_per_differentiating_feature = kw_per_hgnc_age_group_grouped,
    significance_limit = 0.05)

significant_kw_wilcox_age_group_grouped <-
  test_pairwise_wilcoxon(
    melted_kth,
    sym("HGNC_Symbol"),
    significant_kw_per_hgnc_age_group_grouped,
    sym("Age_Group"),
    sym("median_ab_readout"))

kth_summary <- melted_kth %>%
dplyr::summarise(mean = mean(median_ab_readout),
                 sd = sd(median_ab_readout),
                 median = median(median_ab_readout),
                 .by = c(Diagnosis, Age_Group, Gender, HGNC_Symbol))

boxplot(EDN1 ~ Diagnosis, data = kth)

################################################################################
# EXPERIMENTAL
################################################################################

supercharged_melt <-
  melted_kth %>%
  dplyr::filter(Diagnosis == "AD" | Diagnosis == "Control") %>%
  dplyr::mutate(Diagnosis_Age_Group = paste(Diagnosis, Age_Group)) %>%
  dplyr::select(Diagnosis_Age_Group, Gender, HGNC_Symbol, median_ab_readout)

kw_per_hgnc_diagnosis_age_group <-
  sapply(
    unique(supercharged_melt$HGNC_Symbol),
    kruskal_wallis_test_for_sample_groups,
    supercharged_melt,
    sym("HGNC_Symbol"),
    sym("Diagnosis_Age_Group"),
    sym("median_ab_readout"))

signif_kw_per_hgnc_diagnosis_age_group <-
  keep_only_significant_entries(
    supercharged_melt,
    differentiating_feature_symbol = sym("HGNC_Symbol"),
    kw_per_differentiating_feature = kw_per_hgnc_diagnosis_age_group,
    significance_limit = 0.05)

signif_kw_wilcox_diagnosis_age_group <-
  test_pairwise_wilcoxon(
    supercharged_melt,
    sym("HGNC_Symbol"),
    signif_kw_per_hgnc_diagnosis_age_group,
    sym("Diagnosis_Age_Group"),
    sym("median_ab_readout"))

mutated_kth <- kth %>%
  convert_age_column_to_age_group_column(sym("Age"), sym("Age_Group")) %>%
  dplyr::filter(Diagnosis == "AD" | Diagnosis == "Control") %>%
  dplyr::mutate(Diagnosis_Age_Group = paste(Diagnosis, Age_Group))

create_box_plot <- function(dataframe, x_column, y_column) {
  x_column_symbol <- sym(x_column)
  y_column_symbol <- sym(y_column)
  return(
    ggplot(dataframe, aes(x = !!x_column_symbol, y = !!y_column_symbol)) +
      geom_boxplot() +
      scale_x_discrete() #+
      #labs(y = y_column)
  )
}

create_box_plots_of_all_columns_starting_at_column_number <-
  function(dataframe, x_column, y_starting_column) {
  number_of_columns <- length(names(dataframe))
  y_column_names <- names(dataframe)[y_starting_column:number_of_columns]
  new_length <- number_of_columns - y_starting_column
  box_plots <- vector("list", length = new_length)
  for (i in 1:new_length) { # TODO: Possibly replace with sapply
    box_plots[[i]] <- create_box_plot(dataframe, x_column, y_column_names[[i]])
  }
  return(box_plots)
}
create_box_plot(mutated_kth, "Diagnosis_Age_Group", "BDP1")
boxplot(BDP1 ~ Diagnosis_Age_Group, data = mutated_kth)

boxplot(CHGB ~ Diagnosis_Age_Group, data = mutated_kth)

boxplot(DDAH1 ~ Diagnosis_Age_Group, data = mutated_kth)

boxplot(EDN1 ~ Diagnosis_Age_Group, data = mutated_kth)

boxplot(IL1B ~ Diagnosis_Age_Group, data = mutated_kth)

boxplot(MDH1 ~ Diagnosis_Age_Group, data = mutated_kth)

boxplot(PEBP1 ~ Diagnosis_Age_Group, data = mutated_kth)

boxplot(S100B ~ Diagnosis_Age_Group, data = mutated_kth)

################################################################################
#####################################EMIF#######################################
################################################################################
if(!file.exists(file.path(datasets_root_directory,
                          "EMIF-AD MBD Study/dcast_proteins.tsv"))) {
  melted_proteins <-
    read_df_from_file(file.path(datasets_root_directory,
                                "EMIF-AD MBD Study/data_protein.tsv"))
  casted_proteins <-
    reshape2::dcast(melted_proteins,
            Assay.ID ~ PEPTIDE.SEQUENCE,
            value.var = "VALUE") %>%
    dplyr::select_if(~!any(is.na(.)))
  
  write.table(casted_proteins, file = file.path(datasets_root_directory, "EMIF-AD MBD Study/dcast_proteins.tsv"),
              sep = "\t", row.names = F)
}

emif_dir <- "EMIF-AD MBD Study"

excludes <- list(
  c("Fu","Study","Diag_num", "VisitTime", "Studyname", "Center", "Eduy",
    "Eduy_IMP", "BL_Diaggroups", "AMYLOIDstatus", "ABmethod",
    "MMSE", "DateMMSE", "DateBloodCollection", "DateCSFCollection", "DateMRI",
    "DatePET","FU_available", "LastFU", "LastFU_Date", "LastFU_Diagnosis",
    "MCI_Convert", "CTR_Convert", "Plasma_ID", "MRI_ID", "CSF_ID", "Local_AB42",
    "Local_AB42_Abnormal", "Local_AB42_Cutoff", "Local_PTAU", "Local_PTAU_Abnormal",
    "Local_PTAU_Cutoff", "Local_TTAU", "Local_TTAU_Abnormal", "Local_TTAU_Cutoff",
    "AmyPET_SUVR", "AmyPET_SUVRAbnormal", "AmyPET_SUVRCutoff", "AmyPET_Method",
    "AmyPET_Tracer", "AmyPET_Manufacturer", "AB_Zscore", "Ptau_ASSAY_Zscore",
    "Ttau_ASSAY_Zscore", "Central_CSF_AB40", "Central_CSF_AB42",
    "Central_CSF_AB4240ratio", "Central_CSF_ratiodich", "Hippo_sum_FS_adj",
    "LOCAL_CSFAssay", "LOCAL_CSFAssayCompany", "CDR_Global", "IntervalAB_Cog",
    "IntervalAB_MRI", "IntervalAB_Blood", "IntervalMRI_Blood", "Sel_intervalABCog",
    "Sel_intervalABMRI", "Sel_intervalABBlood", "Sel_intervalCogMRI",
    "Sel_intervalCogBlood", "Sel_intervalMRIBlood", "APOEdich", "APOEdich_IMP",
    "APOE_Allele1", "APOE_Allele2", "APOE_flag"),
  c("Sample.Type", "Time.Point", "Tissue.Type", "Platform.ID", "Sample.Code"),
  c())

emif <-
  meld_fractured_dataset(
    dataset_paths = c(
      file.path(datasets_root_directory, emif_dir,"Clinical_EMIFMBD_MIRIADE.xlsx"),
      file.path(datasets_root_directory, emif_dir,"samples.tsv"),
      file.path(datasets_root_directory, emif_dir,"dcast_proteins.tsv")),
    by_columns = setNames(c("Subject.ID","Assay.ID"), c("SubjectId","Assay.ID")),
    exclude_columns = excludes) %>%
  dplyr::mutate(Gender = case_when(Gender == 0 ~ 'm', Gender == 1 ~ 'f', TRUE ~ as.character(Gender)))

### Melt the dataframe, so there's only one readout variable
melted_emif <-
  reshape2::melt(emif, id = 1:5,
                 variable.name = "UniProt", value.name = "Value") %>%
  #dplyr::filter(Diagnosis == "AD") %>%
  convert_age_column_to_age_group_column(sym("Age"), sym("Age_Group")) %>%
  dplyr::select(Diagnosis, Age_Group, Gender, UniProt, Value)

### Run Kruskal-Wallis test for all biomarkers grouped by gender
kw_per_uniprot_gender_grouped <-
  sapply(
    unique(melted_emif$UniProt),
    kruskal_wallis_test_for_sample_groups,
    melted_emif,
    sym("UniProt"),
    sym("Gender"),
    sym("Value"))

significant_kw_per_uniprot_gender_grouped <-
  keep_only_significant_entries(
    melted_emif,
    differentiating_feature_symbol = sym("UniProt"),
    kw_per_differentiating_feature = kw_per_uniprot_gender_grouped,
    significance_limit = 0.000001)

emif_kw_wilcox_gender_grouped <-
  test_pairwise_wilcoxon(
    melted_emif,
    sym("UniProt"),
    significant_kw_per_uniprot_gender_grouped,
    sym("Gender"),
    sym("Value"))

### Adjust for measurements
emif_kw_wilcox_gender_grouped_adj <- matrix(p.adjust(emif_kw_wilcox_gender_grouped, method = "BH"), nrow = 1)
rownames(emif_kw_wilcox_gender_grouped_adj) <- rownames(emif_kw_wilcox_gender_grouped)
colnames(emif_kw_wilcox_gender_grouped_adj) <- colnames(emif_kw_wilcox_gender_grouped)

### Construct the final table
final <- data.frame(Name = colnames(emif_kw_wilcox_gender_grouped),
                    p.val = t(emif_kw_wilcox_gender_grouped)[,1],
