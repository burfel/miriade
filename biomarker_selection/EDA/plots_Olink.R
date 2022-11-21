##################################################
## Project: MIRIADE
## Script purpose: Explore the Olink dataset
## Date: 03.10.2022
## Author: Felicia Burtscher
##################################################
options(stringsAsFactors = F)

library(dplyr)
library(readxl)
library(ggplot2)

### read the VUMC Olink cohort file
olink_cohort <- read.table("/Users/felicia.burtscher/Documents/UL/DATASETS/Olink/OLINK basic data files/clin_bl_basic_04022020.csv", sep = ",", header = T) 

# select the three columns and filter out anything not the 3 diseases (AD, DLB or FTD)
age_sex_per_disease <- dplyr::select(olink_cohort, dx, age_gr, sex) %>% 
  dplyr::filter(dx == "AD dementia" | dx == "DLB" | dx == "FTD") %>%
  dplyr::mutate(sex = case_when(sex == 1 ~ 'm', sex == 2 ~ 'f', TRUE ~ sex))

# AUXILIARY FUNCTION FOR make_histogram_plots. DO NOT USE DIRECTLY.
generate_ggplot_object <- function(data, column_names) {
  x_column_symbol <- sym(column_names[["x"]])
  if("fill" %in% names(column_names)) {
    histogram <- ggplot(data, aes(x = !!x_column_symbol, fill = !!sym(column_names[["fill"]])))
  } else {
    histogram <- ggplot(data, aes(x = !!x_column_symbol))
  }
  return(histogram)
}
# AUXILIARY FUNCTION FOR make_histogram_plots. DO NOT USE DIRECTLY.
append_geom_histogram <- function(histogram, data, column_names, type, binwidth) {
  if (type == "count") {
    aesthethic <- aes(y=..count..)
  } else if (type == "frequency") {
    aesthethic <- aes(y=..density..)
  }
  if(is.numeric(data[[column_names[["x"]]]])) {
    histogram <- histogram + geom_histogram(aesthethic, position = "dodge", binwidth = binwidth)
  }
  else {
    histogram <- histogram + geom_histogram(position = "dodge", stat = "count")
  }
  return(histogram)
}
# AUXILIARY FUNCTION FOR make_histogram_plots. DO NOT USE DIRECTLY.
append_geom_density <- function(histogram, number_of_rows, type, binwidth) {
  if (is.null(binwidth)) {
    binwidth <- 30 # Since it's the default that would have been used
  }
  if (type == "count") {
    histogram <- histogram + geom_density(aes(y = ..density.. * (number_of_rows * binwidth) / 3), alpha=0.3)
  } else if (type == "frequency") {
    histogram <- histogram + geom_density(alpha=0.3)
  }
  return(histogram)
}
# AUXILIARY FUNCTION FOR make_histogram_plots. DO NOT USE DIRECTLY.
append_geom_facet <- function(histogram, facet_name) {
  facet_name_symbol <- sym(facet_name)
  histogram <- histogram + facet_grid(eval(expr(!!facet_name_symbol ~ .)))
  return(histogram)
}
# AUXILIARY FUNCTION FOR make_histogram_plots. DO NOT USE DIRECTLY.
append_histogram_labels <- function(histogram, titles, plot_title_suffix, type, are_there_facets) {
  plot_title <- paste("Patient", tolower(titles[["x"]]))
  if (type == "count") {
    plot_title <- paste(plot_title, "counts")
    y_title <- "# of patients"
  } else if (type == "frequency") {
    plot_title <- paste(plot_title, "frequencies")
    y_title <- "Frequency of patients"
  }
  if ("fill" %in% names(titles)) {
    plot_title <- paste(plot_title, "per", tolower(titles[["fill"]]))
  }
  if (are_there_facets) {
    plot_title <- paste(plot_title, "divided by", tolower(titles[["facet"]]))
  }
  plot_title <- paste(plot_title, plot_title_suffix)
  histogram <- histogram + labs(title = plot_title, x = titles[["x"]], y = y_title, fill = titles[["fill"]])
  return(histogram)
}
# Function to creaate histogram plots
# data - dataframe to be plotted
# column_names - a list containing the column names from which the x axis, fill and facet will take the data.
#                The fill is optional in case you want to split the data into several groups.
#                The facet is optional in case you want to split into several histograms each for a specific value from that column.
#                Should be of the format `list(x = ?, fill = ?, facet =?)`
# titles - a list containing the text for the titles of the x axis and fill legend. Same format as column_names (wihtout facet).
# plot_title_suffix - the plot title will use the previous to populate a standard title format but will need a suffix to differentiate similar plots
# should_density - a boolean whether the density should be plotted. Default is FALSE. A count density requires binwidth to be defined..(?)
# type - count or frequency. Default is count.
# binwidth - the binwidth. The default is NULL, in which case this will default to whatever produces 30 bins.
make_histogram_plots <- function(data, column_names, titles, plot_title_suffix = NULL, should_density = FALSE, type = "count", binwidth = NULL) {
  are_there_facets <- "facet" %in% names(column_names)
  histogram <- generate_ggplot_object(data, column_names) %>%
    append_geom_histogram(data, column_names, type, binwidth) %>%
    {`if`(should_density, append_geom_density(., nrow(data), type, binwidth), .)} %>%
    {`if`(are_there_facets, append_geom_facet(., column_names[["facet"]]), .)} %>%
    append_histogram_labels(titles, plot_title_suffix, type, are_there_facets)
  return(histogram)
}

# Plot 1: distribution or histogram: age across diseases (x-axis: age, y-axis: #patients), per disease (diseases in different colour)
age_binwidth <- 5
# age_freq_hist <- ggplot(age_sex_per_disease, aes(x = age_gr, fill = dx)) + 
#   geom_histogram(aes(y=..density..), position = "dodge", binwidth = age_binwidth) +
#   geom_density(alpha=0.3)+
#   labs(title = "Patient age frequencies per disease", x = "Age", y = "Frequency of patients", fill = "Disease")
age_freq_hist <- make_histogram_plots(age_sex_per_disease, 
                                                 column_names = list(x="age_gr", fill="dx"),
                                                 titles = list(x="Age", fill="Disease"),
                                                 plot_title_suffix = "for Olink",
                                                 should_density = TRUE,
                                                 type = "frequency",
                                                 binwidth = age_binwidth)
age_freq_hist
# age_count_hist <- ggplot(age_sex_per_disease, aes(x = age_gr, fill = dx)) + 
#   geom_histogram(aes(y=..count..), position = "dodge", binwidth = age_binwidth) +
#   geom_density(aes(y = ..density.. * (nrow(age_sex_per_disease) * age_binwidth) / 2), alpha=0.3)+
#   labs(title = "Patient age counts per disease", x = "Age", y = "# of patients", fill = "Disease")
age_count_hist <- make_histogram_plots(age_sex_per_disease,
                                                  column_names = list(x="age_gr", fill="dx"),
                                                  titles = list(x="Age", fill="Disease"),
                                                  plot_title_suffix = "for Olink",
                                                  should_density = TRUE,
                                                  type = "count",
                                                  binwidth = age_binwidth)
age_count_hist

# Plot 2: same as plot 1, but dissected according to sex, so 2 curves per disease
# sex_hist <- ggplot(age_sex_per_disease, aes(x = sex, fill = dx)) + 
#   geom_histogram(position = "dodge", stat = "count") +
#   labs(title = "Patient sex counts per disease", x = "Sex", y = "# of patients", fill = "Disease")
sex_hist <- make_histogram_plots(age_sex_per_disease,
                                            column_names = list(x="sex", fill="dx"),
                                            titles = list(x="Sex", fill="Disease"),
                                            plot_title_suffix = "for Olink")
sex_hist

# Plot 2.5: according to age and sex
# age_and_sex_freq_hist <- ggplot(age_sex_per_disease, aes(x = age_gr, fill = dx)) +
#   geom_histogram(aes(y=..density..), position = "dodge", binwidth = age_binwidth) +
#   geom_density(alpha=0.3) +
#   facet_grid(sex ~ .) +
#   labs(title = "Patient age frequencies per disease divided by sex", x = "Age", y = "Frequency of patients", fill = "Disease")
age_and_sex_freq_hist <- make_histogram_plots(age_sex_per_disease,
                                                         column_names = list(x="age_gr", fill="dx", facet="sex"),
                                                         titles = list(x="Age", fill="Disease", facet="Sex"),
                                                         plot_title_suffix = "for Olink",
                                                         should_density = TRUE,
                                                         type = "frequency",
                                                         binwidth = age_binwidth)
age_and_sex_freq_hist
# age_and_sex_count_hist <- ggplot(age_sex_per_disease, aes(x = age_gr, fill = dx)) + 
#   geom_histogram(aes(y=..count..), position = "dodge", binwidth = age_binwidth) +
#   geom_density(aes(y = ..density.. * (nrow(age_sex_per_disease) * age_binwidth) / 5), alpha=0.3) +
#   facet_grid(sex ~ .) +
#   labs(title = "Patient age counts per disease divided by sex", x = "Age", y = "# of patients", fill = "Disease")
age_and_sex_count_hist <- make_histogram_plots(age_sex_per_disease,
                                                         column_names = list(x="age_gr", fill="dx", facet="sex"),
                                                         titles = list(x="Age", fill="Disease", facet="Sex"),
                                                         plot_title_suffix = "for Olink",
                                                         should_density = TRUE,
                                                         type = "count",
                                                         binwidth = age_binwidth)
age_and_sex_count_hist

# Samey for KTH
kth_age_sex_disease <- read.table("/Users/felicia.burtscher/Documents/UL/DATASETS/KTH/KTH AD dataset for MIRIADE biomarker selection sample info and results.csv"
                                  , sep = ";", fileEncoding = "UTF-8", header = T) %>%
  dplyr::select(Class, Diagnosis, Age, Gender, Tau, Ptau, Ab42, ApoE4, MMSE_score)
  
kth_age_count_hist <- make_histogram_plots(kth_age_sex_disease,
                                                      column_names = list(x="Age", fill="Diagnosis"),
                                                      titles = list(x="Age", fill="Diagnosis"),
                                                      plot_title_suffix = "for KTH",
                                                      should_density = TRUE,
                                                      binwidth = age_binwidth)
kth_age_count_hist

kth_sex_hist <- make_histogram_plots(kth_age_sex_disease,
                                            column_names = list(x="Gender", fill="Diagnosis"),
                                            titles = list(x="Sex", fill="Diagnosis"),
                                            plot_title_suffix = "for KTH")
kth_sex_hist

kth_age_and_sex_count_hist <- make_histogram_plots(kth_age_sex_disease,
                                                          column_names = list(x="Age", fill="Diagnosis", facet="Gender"),
                                                          titles = list(x="Age", fill="Diagnosis", facet="Sex"),
                                                          plot_title_suffix = "for KTH",
                                                          should_density = TRUE,
                                                          type = "count",
                                                          binwidth = age_binwidth)
kth_age_and_sex_count_hist

emif_age_sex_disease <- read.table("/Users/felicia.burtscher/Documents/UL/DATASETS/EMIF-AD MBD Study/Clinical_EMIFMBD_MIRIADE.csv"
                                   , sep = ";", fileEncoding = "UTF-8", header = T, dec = ",") %>% 
  dplyr::filter(Age != "" & Diagnosis !="") %>%
  dplyr::mutate(Gender = case_when(Gender == 0 ~ 'm', Gender == 1 ~ 'f', TRUE ~ as.character(Gender)))

emif_age_count_hist <- make_histogram_plots(emif_age_sex_disease ,
                                           column_names = list(x="Age", fill="Diagnosis"),
                                           titles = list(x="Age", fill="Diagnosis"),
                                           plot_title_suffix = "for Emif",
                                           should_density = TRUE,
                                           binwidth = age_binwidth)
emif_age_count_hist

emif_sex_hist <- make_histogram_plots(emif_age_sex_disease,
                                     column_names = list(x="Gender", fill="Diagnosis"),
                                     titles = list(x="Sex", fill="Diagnosis"),
                                     plot_title_suffix = "for Emif")
emif_sex_hist

emif_age_and_sex_count_hist <- make_histogram_plots(emif_age_sex_disease,
                                                   column_names = list(x="Age", fill="Diagnosis", facet="Gender"),
                                                   titles = list(x="Age", fill="Diagnosis", facet="Sex"),
                                                   plot_title_suffix = "for Emif",
                                                   should_density = TRUE,
                                                   type = "count",
                                                   binwidth = age_binwidth)
emif_age_and_sex_count_hist

### Read the VUMC Olink biomarkers
vumc_ol <- read.table("/Users/felicia.burtscher/Documents/UL/DATASETS/Olink/20201216_OLINK_Data/cleaned_VUMC_olink.tsv", sep = "\t", header = T,
                      quote = "") ### 'quote' needed due to a special sign in 'Name'

# Plot 3: 3D p-value space of the 3 diseases
# points corresponding to proteins should be labeled
# missing values lead to skipping that protein
library(tidyr)
# Now we want to pivot so that each unique protein will have 3 values for p values of each disease
# we would also like to have a column that will help color the plot so that proteins significant for
# all diseases and proteins significant to only one will be easily identified
significance_limit <- 0.05 # below ~0.00002322709 significance limit we lose the significant for all category
combined_3d_p_value_per_protein <- tidyr::pivot_wider(vumc_ol, names_from = "disease", values_from = "p.val") %>%
  dplyr::mutate(Significance = case_when(
    (ALZ <= significance_limit & DLB <= significance_limit & FTD <= significance_limit) ~ "Significant for all",
    (ALZ <= significance_limit & DLB > significance_limit & FTD > significance_limit) ~ "Significant for ALZ",
    (ALZ > significance_limit & DLB <= significance_limit & FTD > significance_limit) ~ "Significant for DLB",
    (ALZ > significance_limit & DLB > significance_limit & FTD <= significance_limit) ~ "Significant for FTD",
    TRUE ~ "Not interesting"
  ))
if(!require("plotly")) {
  install.packages("plotly")
}
library("plotly")
plot_ly(combined_3d_p_value_per_protein, x = ~ALZ, y = ~DLB, z = ~FTD, text = ~UniProt) %>%
  add_markers(color = ~Significance, 
              hovertemplate = paste('<b>UniProt</b>: %{text}',
                                    '<br><b>ALZ p value</b>: %{x}',
                                    '<br><b>DLB p value</b>: %{y}',
                                    '<br><b>FTD p value</b>: %{z}'),) %>%
  layout(title = 'P value of proteins for 3 diseases',
    scene = list(
      xaxis = list(title = "ALZ p value"),
      yaxis = list(title = "DLB p value"),
      zaxis = list(title = "FTD p value")
    )
  )

## TODO: Possibly add enum implementation from https://stackoverflow.com/questions/33838392/enum-like-arguments-in-r for significance_type
# Filters only the rows with the Significance matching significance_type and leaves only the UniProt column
extract_uniprots_according_to_significance_type <- function(protein_pvalues, significance_type) {
  protein_pvalues %>%
    dplyr::filter(Significance == significance_type) %>%
    dplyr::select(UniProt)
}

## TODO: delete commented out lines once the functionality of replacements is verified

# A list of UniProts significant for all 3 diseases
# biomarkers_for_all_dementias <- combined_3d_p_value_per_protein %>%
#   dplyr::filter(Significance == "Significant for all") %>%
#   dplyr::select(UniProt)
biomarkers_for_all_dementias <- extract_uniprots_according_to_significance_type(combined_3d_p_value_per_protein, "Significant for all")
# A list of UniProts significant for ALZ
# biomarkers_for_ALZ <- combined_3d_p_value_per_protein %>%
#   dplyr::filter(Significance == "Significant for ALZ") %>%
#   dplyr::select(UniProt)
biomarkers_for_ALZ <- extract_uniprots_according_to_significance_type(combined_3d_p_value_per_protein, "Significant for ALZ")
# A list of UniProts significant for DLB
# biomarkers_for_DLB <- combined_3d_p_value_per_protein %>%
#   dplyr::filter(Significance == "Significant for DLB") %>%
#   dplyr::select(UniProt)
biomarkers_for_DLB <- extract_uniprots_according_to_significance_type(combined_3d_p_value_per_protein, "Significant for DLB")

# A list of UniProts significant for FTD
# biomarkers_for_FTD <- combined_3d_p_value_per_protein %>%
#   dplyr::filter(Significance == "Significant for FTD") %>%
#   dplyr::select(UniProt)
biomarkers_for_FTD <- extract_uniprots_according_to_significance_type(combined_3d_p_value_per_protein, "Significant for FTD")

# Plot 4: 1-d curve of the 3 diseases
# we give the proteins an ordering, a panel of proteins according to the cleaned_VUMC_olink.tsv dataset
# vumc_ol needs to filter out the non important ones
# vumc_significants <- vumc_ol %>% 
#   dplyr::filter(UniProt %in% biomarkers_for_all_dementias$UniProt |
#                   UniProt %in% biomarkers_for_ALZ$UniProt |
#                   UniProt %in% biomarkers_for_DLB$UniProt |
#                   UniProt %in% biomarkers_for_FTD$UniProt)
# pval_vs_uniprot <- ggplot(vumc_significants, aes(x = UniProt, y = p.val, group = disease, color = disease)) +
#                             geom_line() + geom_point()
# pval_vs_uniprot
pval_vs_uniprot_significant <- ggplot(vumc_ol %>%
                                        dplyr::filter(UniProt %in% biomarkers_for_all_dementias$UniProt)
                                      , aes(x = UniProt, y = p.val, group = disease, color = disease)) +
  geom_line() + geom_point() +
  labs(title = "P values for proteins significant for all 3 diseases")
pval_vs_uniprot_significant
pval_vs_uniprot_ALZ_significant <- ggplot(vumc_ol %>%
                                        dplyr::filter(UniProt %in% biomarkers_for_ALZ$UniProt &
                                                        disease == "ALZ")
                                      , aes(x = UniProt, y = p.val, group = disease)) +
  geom_line() + geom_point() +
  labs(title = "P values for ALZ for proteins significant for ALZ")
pval_vs_uniprot_ALZ_significant
pval_vs_uniprot_DLB_significant <- ggplot(vumc_ol %>%
                                            dplyr::filter(UniProt %in% biomarkers_for_DLB$UniProt &
                                                            disease == "DLB")
                                          , aes(x = UniProt, y = p.val, group = disease)) +
  geom_line() + geom_point() +
  labs(title = "P values for DLB for proteins significant for DLB")
pval_vs_uniprot_DLB_significant
pval_vs_uniprot_FTD_significant <- ggplot(vumc_ol %>%
                                            dplyr::filter(UniProt %in% biomarkers_for_FTD$UniProt &
                                                            disease == "FTD")
                                          , aes(x = UniProt, y = p.val, group = disease)) +
  geom_line() + geom_point() +
  labs(title = "P values for FTD for proteins significant for FTD")
pval_vs_uniprot_FTD_significant

# Read the KTH
kth <- read.table("/Users/felicia.burtscher/Documents/UL/DATASETS/KTH/cleaned_KTH_entries.tsv", sep = "\t", header = T,
              quote = "")

# Preprocess datasets to have the unified format of UniProt and Pvalue, possibly after filtering
# some of them according to disease.
# datasets and should_filter_by_disease must be lists/vectors of the same length
# datasets HAS TO BE A LIST OF DATAFRAMES
# Any dataframe for which should_filter_by_disease is true needs to have a disease column
# If any of should_filter_by_disease is true, disease_to_filter_by should be a string that appears in the disease column
preprocess_datasets <- function(datasets, should_filter_by_disease, disease_to_filter_by = NULL) {
  datasets_length <- length(datasets)
  if (datasets_length != length(should_filter_by_disease)) {
    stop("datasets and should_filter_by_disease must be lists/vectors of the same length")
  }
  processed_datasets <- vector(mode = "list", length = datasets_length)
  ind <- 1
  for (dataset in datasets) {
    if(should_filter_by_disease[ind]) {
      dataset <- dataset %>% dplyr::filter(disease == disease_to_filter_by)
    }
    dataset <- dataset %>% dplyr::select(UniProt, p.val)
    processed_datasets[[ind]] <- dataset
    ind <- ind + 1
  }
  return(processed_datasets)
}

# Binds datasets together and produces a dataframe with UniProt p_values and a column for the source of that row
# TODO: fold the mutate into this function once we realize how to use case_when with expressions built
# on the fly with indices unknown in advance
bind_datasets_together <- function(processed_datasets, should_filter_by_disease, disease_to_filter_by = NULL) {
  # Unfortunately, failed to do the renaming of source here... might have to repeat it
  return(dplyr::bind_rows(processed_datasets, .id = "Source")
  )
}

# Redo the ALZ plot with KTH?
# alz_olink_and_kth <- dplyr::bind_rows(
#   vumc_ol %>%
#     dplyr::filter(disease == "ALZ") %>%
#     dplyr::select(UniProt, p.val),
#   kth %>% dplyr::select(UniProt, p.val), .id = "Source") %>%
#   dplyr::mutate(Source = case_when(Source == 1 ~ 'olink', Source == 2 ~ 'kth'))
alz_olink_and_kth <- bind_datasets_together(
  preprocess_datasets(
  datasets = list(vumc_ol,kth),
  should_filter_by_disease = c(TRUE,FALSE),
  disease_to_filter_by = "ALZ")) %>%
  dplyr::mutate(Source = case_when(Source == 1 ~ 'olink', Source == 2 ~ 'kth'))
  
  
pval_vs_uniprot_ALZ_kth <- ggplot(alz_olink_and_kth
                                          , aes(x = UniProt, y = p.val, color = Source)) +
  geom_line() + geom_point() +
  labs(title = "P values for ALZ including kth")
pval_vs_uniprot_ALZ_kth
pval_vs_uniprot_ALZ_kth_significant <- ggplot(alz_olink_and_kth %>%
                                                dplyr::filter(p.val <= 0.05)
                                  , aes(x = UniProt, y = p.val, color = Source)) +
  geom_vline(xintercept = kth$UniProt[which(kth$UniProt %in% vumc_ol$UniProt)]) + geom_point() +
  labs(title = "P values for ALZ for proteins significant for ALZ in olink and kth")
pval_vs_uniprot_ALZ_kth_significant

# This function takes a dataframe with uniprots and p values from more than one source
# filters and returns according to the list of limiting dataframes
filter_according_to_UniProt_masks <- function(data_frame, masking_frames_vector) {
  filtered_data_frame <- data_frame
  for (mask in masking_frames_vector) {
    filtered_data_frame <- filtered_data_frame %>%
      dplyr::filter(UniProt %in% mask$UniProt)
  }
  return(filtered_data_frame)
}

# alz_olink_and_kth_intersection_list <- alz_olink_and_kth %>%
#   dplyr::filter(UniProt %in% kth$UniProt & UniProt %in% vumc_ol$UniProt)
alz_olink_and_kth_intersection_list <- filter_according_to_UniProt_masks(alz_olink_and_kth, list(kth, vumc_ol))
pval_vs_uniprot_olink_kth_intersection <- ggplot(alz_olink_and_kth_intersection_list
                                              , aes(x = UniProt, y = p.val, color = Source)) + geom_point() +
  labs(title = "P values for ALZ for proteins for ALZ both in olink and kth")
pval_vs_uniprot_olink_kth_intersection

# Takes a datafrane with a UniProt column, and extracts the distinct UniProts
extract_unique_uniprots <- function(dataframe) {
  return(dataframe %>%
           dplyr::select(UniProt) %>% dplyr::distinct(UniProt))
}

alz_olink_and_kth_intersection_uniprots <- extract_unique_uniprots(alz_olink_and_kth_intersection_list)



# This function takes a dataframe with uniprots and p values from more than one source
# and filters according to the list of limiting dataframes. Also, takes a title
create_plot_of_significant_according_to_limiting_lists <- function(data_frame_to_be_plotted, masking_frames_vector, title) {
  return(ggplot(filter_according_to_UniProt_masks(data_frame_to_be_plotted, masking_frames_vector),
                  aes(x = UniProt, y = p.val, color = Source)) + geom_point() +
    labs(title = title))
}

# pval_vs_uniprot_olink_kth_intersection_ALZ_significant <- ggplot(alz_olink_and_kth_intersection_list %>%
#                                                                    dplyr::filter(UniProt %in% biomarkers_for_ALZ$UniProt)
#                                                  , aes(x = UniProt, y = p.val, color = Source)) + geom_point() +
#   labs(title = "P values for ALZ for significnat proteins for ALZ both in olink and kth")
pval_vs_uniprot_olink_kth_intersection_ALZ_significant <- create_plot_of_significant_according_to_limiting_lists(
  alz_olink_and_kth_intersection_list, list(biomarkers_for_ALZ), "P values for ALZ for significnat proteins for ALZ both in olink and kth"
)
pval_vs_uniprot_olink_kth_intersection_ALZ_significant
# only_olink_and_kth <- dplyr::bind_rows(
#   vumc_ol %>% 
#     dplyr::filter(UniProt %in% kth$UniProt &
#                     disease == "ALZ"),
#   kth %>% 
#     dplyr::filter(UniProt %in% vumc_ol$UniProt),  .id = "Source") %>%
#   dplyr::select(UniProt, p.val, Source) %>%
#   dplyr::mutate(Source = case_when(Source == 1 ~ 'olink', Source == 2 ~ 'kth'))
# ggplot(only_olink_and_kth %>%
#          dplyr::filter(p.val <= 0.05)
#        , aes(x = UniProt, y = p.val, color = Source)) + geom_point()
# dplyr::inner_join(kth, vumc_ol %>% dplyr::filter(disease == "ALZ"), by = "UniProt",
#                   suffix = c(".kth", ".olink")) %>%
#   dplyr::select(UniProt, p.val.kth, p.val.olink)



adni <- read.table("/Users/felicia.burtscher/Documents/UL/DATASETS/Brain 2020 AD/cleaned_ADNI_entries.tsv", sep = "\t", header = T,
                  quote = "")
emif <- read.table("/Users/felicia.burtscher/Documents/UL/DATASETS/Brain 2020 AD/cleaned_EMIF_entries.tsv", sep = "\t", header = T,
                   quote = "")
mspec <- read.table("/Users/felicia.burtscher/Documents/UL/DATASETS/VUMC/cleaned_VUMC_mspec.tsv", sep = "\t", header = T,
                   quote = "")

preprocessed_olink_and_adni <- preprocess_datasets(
  datasets = list(vumc_ol,adni),
  should_filter_by_disease = c(TRUE,FALSE),
  disease_to_filter_by = "ALZ")
# Now do the same thing we did with olink vs. kth with the rest
alz_olink_and_adni <- bind_datasets_together(preprocessed_olink_and_adni) %>%
  dplyr::mutate(Source = case_when(Source == 1 ~ 'olink', Source == 2 ~ 'adni'))

alz_olink_and_adni_intersection_list <- filter_according_to_UniProt_masks(alz_olink_and_adni,
                                                                          preprocessed_olink_and_adni)
alz_olink_and_adni_intersection_uniprots <- extract_unique_uniprots(alz_olink_and_adni_intersection_list)
alz_olink_and_adni_intersection_significant_in_both_datasets <- 
  dplyr::inner_join(preprocessed_olink_and_adni[[1]], preprocessed_olink_and_adni[[2]], by = "UniProt",
                   suffix = c("_olink", "_adni")) %>%
  dplyr::filter(p.val_olink <= 0.05 & p.val_adni <=0.05)
remove(preprocessed_olink_and_adni)
# pval_vs_uniprot_olink_adni_intersection_ALZ_significant <- ggplot(alz_olink_and_adni_intersection_list %>%
#                                                                    dplyr::filter(UniProt %in% biomarkers_for_ALZ$UniProt &
#                                                                                    UniProt %in% alz_olink_and_adni_intersection_significant_in_both_datasets$UniProt)
#                                                                  , aes(x = UniProt, y = p.val, color = Source)) + geom_point() +
#   labs(title = "P values for ALZ for significnat proteins for ALZ both in olink and adni")
pval_vs_uniprot_olink_adni_intersection_ALZ_significant <- create_plot_of_significant_according_to_limiting_lists(
  alz_olink_and_adni_intersection_list, list(biomarkers_for_ALZ, alz_olink_and_adni_intersection_significant_in_both_datasets),
  "P values for ALZ for significant proteins for ALZ both in olink and adni"
)
pval_vs_uniprot_olink_adni_intersection_ALZ_significant

preprocessed_olink_and_emif <- preprocess_datasets(
  datasets = list(vumc_ol,emif),
  should_filter_by_disease = c(TRUE,FALSE),
  disease_to_filter_by = "ALZ")

alz_olink_and_emif <- bind_datasets_together(preprocessed_olink_and_emif) %>%
  dplyr::mutate(Source = case_when(Source == 1 ~ 'olink', Source == 2 ~ 'emif'))

alz_olink_and_emif_intersection_list <- filter_according_to_UniProt_masks(alz_olink_and_emif,
                                                                          preprocessed_olink_and_emif)
alz_olink_and_emif_intersection_uniprots <- extract_unique_uniprots(alz_olink_and_emif_intersection_list)
alz_olink_and_emif_intersection_significant_in_both_datasets <- 
  dplyr::inner_join(preprocessed_olink_and_emif[[1]], preprocessed_olink_and_emif[[2]], by = "UniProt",
                    suffix = c("_olink", "_emif")) %>%
  dplyr::filter(p.val_olink <= 0.05 & p.val_emif <=0.05)
remove(preprocessed_olink_and_emif)

pval_vs_uniprot_olink_emif_intersection_ALZ_significant <- create_plot_of_significant_according_to_limiting_lists(
  alz_olink_and_emif_intersection_list, list(biomarkers_for_ALZ, alz_olink_and_emif_intersection_significant_in_both_datasets),
  "P values for ALZ for significant proteins for ALZ both in olink and emif"
)
pval_vs_uniprot_olink_emif_intersection_ALZ_significant

preprocessed_adni_and_emif <- preprocess_datasets(
  datasets = list(adni,emif),
  should_filter_by_disease = c(FALSE,FALSE))

alz_adni_and_emif <- bind_datasets_together(preprocessed_adni_and_emif) %>%
  dplyr::mutate(Source = case_when(Source == 1 ~ 'adni', Source == 2 ~ 'emif'))

alz_adni_and_emif_intersection_list <- filter_according_to_UniProt_masks(alz_adni_and_emif,
                                                                         preprocessed_adni_and_emif)
alz_adni_and_emif_intersection_uniprots <- extract_unique_uniprots(alz_adni_and_emif_intersection_list)
alz_adni_and_emif_intersection_significant_in_both_datasets <- 
  dplyr::inner_join(preprocessed_adni_and_emif[[1]], preprocessed_adni_and_emif[[2]], by = "UniProt",
                    suffix = c("_adni", "_emif")) %>%
  dplyr::filter(p.val_adni <= 0.05 & p.val_emif <=0.05)
remove(preprocessed_adni_and_emif)

pval_vs_uniprot_adni_emif_intersection_ALZ_significant <- create_plot_of_significant_according_to_limiting_lists(
  alz_adni_and_emif_intersection_list, list(biomarkers_for_ALZ, alz_adni_and_emif_intersection_significant_in_both_datasets),
  "P values for ALZ for significant proteins for ALZ both in adni and emif"
)
pval_vs_uniprot_adni_emif_intersection_ALZ_significant

# All 4 ALZ datasets
preprocessed_olink_adni_emif_kth <- preprocess_datasets(
  datasets = list(vumc_ol,adni,emif,kth),
  should_filter_by_disease = c(TRUE,FALSE,FALSE,FALSE),
  disease_to_filter_by = "ALZ")
alz_olink_adni_emif_kth <- bind_datasets_together(preprocessed_olink_adni_emif_kth) %>%
  dplyr::mutate(Source = case_when(Source == 1 ~ 'olink', Source == 2 ~ 'adni', Source == 3 ~ 'emif', Source == 4 ~ 'kth'))

alz_olink_adni_emif_kth_intersection_list <- filter_according_to_UniProt_masks(alz_olink_adni_emif_kth,
                                                                               preprocessed_olink_adni_emif_kth)
remove(preprocessed_olink_adni_emif_kth)

# Without KTH
preprocessed_olink_adni_emif <- preprocess_datasets(
  datasets = list(vumc_ol,adni,emif),
  should_filter_by_disease = c(TRUE,FALSE,FALSE),
  disease_to_filter_by = "ALZ")
alz_olink_adni_emif <- bind_datasets_together(preprocessed_olink_adni_emif) %>%
  dplyr::mutate(Source = case_when(Source == 1 ~ 'olink', Source == 2 ~ 'adni', Source == 3 ~ 'emif'))

alz_olink_adni_emif_intersection_list <- filter_according_to_UniProt_masks(alz_olink_adni_emif,
                                                                           preprocessed_olink_adni_emif)
alz_olink_adni_emif_intersection_uniprots <- extract_unique_uniprots(alz_olink_adni_emif_intersection_list)
remove(preprocessed_olink_adni_emif)

# Now Olink and mspec
preprocessed_olink_mspec_dlb <- preprocess_datasets(
  datasets = list(vumc_ol,mspec),
  should_filter_by_disease = c(TRUE,TRUE),
  disease_to_filter_by = "DLB")
dlb_olink_mspec <- bind_datasets_together(preprocessed_olink_mspec_dlb) %>%
  dplyr::mutate(Source = case_when(Source == 1 ~ 'olink', Source == 2 ~ 'mspec'))
dlb_olink_mspec_intersection_list <- filter_according_to_UniProt_masks(dlb_olink_mspec, preprocessed_olink_mspec_dlb)
dlb_olink_mspec_intersection_uniprots <- extract_unique_uniprots(dlb_olink_mspec_intersection_list)
dlb_olink_mspec_intersection_significant_in_both_datasets <- 
  dplyr::inner_join(preprocessed_olink_mspec_dlb[[1]], preprocessed_olink_mspec_dlb[[2]], by = "UniProt",
                    suffix = c("_olink", "_mspec")) %>%
  dplyr::filter(p.val_olink <= 0.05 & p.val_mspec <=0.05)
remove(preprocessed_olink_mspec_dlb)

preprocessed_olink_mspec_ftd <- preprocess_datasets(
  datasets = list(vumc_ol,mspec),
  should_filter_by_disease = c(TRUE,TRUE),
  disease_to_filter_by = "FTD")
ftd_olink_mspec <- bind_datasets_together(preprocessed_olink_mspec_ftd) %>%
  dplyr::mutate(Source = case_when(Source == 1 ~ 'olink', Source == 2 ~ 'mspec'))
ftd_olink_mspec_intersection_list <- filter_according_to_UniProt_masks(ftd_olink_mspec, preprocessed_olink_mspec_ftd)
ftd_olink_mspec_intersection_uniprots <- extract_unique_uniprots(ftd_olink_mspec_intersection_list)
ftd_olink_mspec_intersection_significant_in_both_datasets <- 
  dplyr::inner_join(preprocessed_olink_mspec_ftd[[1]], preprocessed_olink_mspec_ftd[[2]], by = "UniProt",
                    suffix = c("_olink", "_mspec")) %>%
  dplyr::filter(p.val_olink <= 0.05 & p.val_mspec <=0.05)
remove(preprocessed_olink_mspec_ftd)


# 15.11.2022
# visualising the original Olink data as box plots
olink_basic <- read.table("/Users/felicia.burtscher/Documents/UL/DATASETS/Olink/OLINK basic data files/protein_data_LODmax_10percent_outliers_rm_missing_imputed_16112020.txt", sep = "\t", header = T) 

create_box_plot_of_column <- function(dataframe, column) {
  column_symbol <- sym(column)
  return(
    ggplot(dataframe, aes(y = !!column_symbol)) +
      geom_boxplot() +
      scale_x_discrete(name = column) +
      labs(y = "values")
  )
}

create_box_plots_of_all_columns_starting_at_column_number <- function(dataframe, starting_column) {
  number_of_columns <- length(names(dataframe))
  column_names <- names(dataframe)[starting_column:number_of_columns]
  new_length <- number_of_columns - starting_column
  box_plots <- vector("list", length = new_length)
  for (i in 1:new_length) {
    box_plots[[i]] <- create_box_plot_of_column(dataframe, column_names[[i]])
  }
  return(box_plots)
}

# This ends up being ~4.7 GB
olink_basic_box_plots <- create_box_plots_of_all_columns_starting_at_column_number(olink_basic, 6)
# Use an index between 1 and 807 (for proteins)
## TODO: perhaps implementing a dictionary to call uniprot-ids (once converted to uniprot ids) or standard gene names
olink_basic_box_plots[[67]]
# Recommend to remove it when done using it
remove(olink_basic_box_plots)
