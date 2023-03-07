#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
## Project: MIRIADE
## Script purpose: Prepare all the fractured datasets by melding them together
##                 and have a df for them
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
library(here)
library(dplyr)

source(here("functions", "data_processing_functions.R"))

datasets_root_directory <- define_datasets_root()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
####                                Read Emif                               ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
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
  remove(melted_proteins, casted_proteins)
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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
####                                 Read Olink                             ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
olink_dir <- "Olink/OLINK basic data files"

excludes <- list(
  c("X","Olink_project", "Olink_plate"),
  c("sam_ind", "X")
)

olink <-
  meld_fractured_dataset(
    dataset_paths = c(
      file.path(datasets_root_directory, olink_dir,"clin_bl_basic_04022020.csv"),
      file.path(datasets_root_directory, olink_dir,"protein_data_LODmax_10percent_outliers_rm_missing_imputed_16112020.tsv")),
    by_columns = setNames(c("SampleId"), c("SampleId")),
    exclude_columns = excludes)  %>%
  dplyr::mutate(sex = case_when(sex == 1 ~ 'm', sex == 2 ~ 'f', TRUE ~ sex))

remove(emif_dir, olink_dir, excludes)
