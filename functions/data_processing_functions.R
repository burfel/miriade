################################################################################
# Takes an interaction dataframe with a "dictionary" dataframe and replace
# the source and target columns using the translation from the dictionary.
# Optionally leaves behind only the translated source and target columns.
#
# @param interactions_df is the dataframe with the interaction data
# @param interaction_columns names of original source and target columns
# @param dictionary_df the dataframe with the tranlsation information
# @param dictionary_columns columns with the old type and the new type
#        to be translated to
# @param new_interaction_names names for the new source and target columns
# @param only_interaction_columns bool whether to only leave the interaction
#        columns in the new dataframe
#
# @return an interaction dataframe between the tranlsated types instead
#         of the original ones
################################################################################
translate_interactions_using_dictionary <- function(interactions_df,
                                                    interaction_columns,
                                                    dictionary_df,
                                                    dictionary_columns,
                                                    new_interaction_names,
                                                    only_interaction_columns = TRUE) {
  require(dplyr)
  new_df <- dplyr::inner_join(interactions_df, dictionary_df, 
                              by = setNames(dictionary_columns[[1]], interaction_columns[[1]])) %>%
    dplyr::rename(!!sym(new_interaction_names[[1]]) := dictionary_columns[[2]]) %>%
    dplyr::inner_join(dictionary_df, 
                      by = setNames(dictionary_columns[[1]], interaction_columns[[2]])) %>%
    dplyr::rename(!!sym(new_interaction_names[[2]]) := dictionary_columns[[2]])
  if (only_interaction_columns) {
    new_df <- dplyr::select(new_df, new_interaction_names[[1]], new_interaction_names[[2]]) 
  }
  return(dplyr::distinct(new_df))
}

################################################################################
# Takes a dataframe with objects of a certain type and the columns containing
# those objects to return a list of all objects of the type
#
# @param df is the dataframe with the desired objects
# @param object_columns names of all columns containing the desired objects
#
# @return a list with all distinct objects from the dataframe
################################################################################
extract_all_distinct_objects <- function(df, object_columns) {
  all_values <- union(df[[object_columns[[1]]]], df[[object_columns[[2]]]])
  length = length(object_columns)
  if (length > 2) {
    for (i in 3:length) {
      all_values <- union(all_values, df[[object_columns[[i]]]])
    }
  }
  return(na.omit(all_values))
}

################################################################################
#' Define data sets root
#'
#' This opens up a dialog to choose the root directory for data sets in one
#' of two cases:
#' 1) If the variable datasets_root_directory doesn't exist
#' 2) If it does exist, the user is prompted if they want to redefine it and
#'    they type in 'y' and press enter
#'
#' @return The path for the directory chosen as the root for the data sets
################################################################################
define_datasets_root <- function() {
  if(exists("datasets_root_directory")
     && !identical(readline(prompt = "datasets_root_directory already exists. If you want to redefine it, type 'y' and press enter: \n"), 'y')) {
    return(datasets_root_directory)
  }
  # Open a file dialog to store the root directory of all data sets
  return(choose.dir(caption = "Choose the root of all datasets"))
}

################################################################################
# Preprocess datasets to have the unified format of UniProt and Pvalue, possibly
# after filtering
# some of them according to disease.
# datasets and should_filter_by_disease must be lists/vectors of the same length
# datasets HAS TO BE A LIST OF DATAFRAMES
# Any dataframe for which should_filter_by_disease is true needs to have a
# disease column
# If any of should_filter_by_disease is true, disease_to_filter_by should be a
# string that appears in the disease column
################################################################################
preprocess_datasets <- function(datasets, should_filter_by_disease,
                                disease_to_filter_by = NULL) {
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

################################################################################
# Binds datasets together and produces a dataframe with UniProt p_values and a
# column for the source of that row
# TODO: fold the mutate into this function once we realize how to use case_when
# with expressions built
# on the fly with indices unknown in advance
################################################################################
bind_datasets_together <- function(processed_datasets) {
  # Unfortunately, failed to do the renaming of source here... might have to
  # repeat it
  return(dplyr::bind_rows(processed_datasets, .id = "Source")
  )
}

################################################################################
# This function takes a dataframe with uniprots and p values from more than one
# source
# filters and returns according to the list of limiting dataframes
################################################################################
filter_according_to_UniProt_masks <- function(data_frame, masking_frames_vector) {
  filtered_data_frame <- data_frame
  for (mask in masking_frames_vector) {
    filtered_data_frame <- filtered_data_frame %>%
      dplyr::filter(UniProt %in% mask$UniProt)
  }
  return(filtered_data_frame)
}

################################################################################
# Takes a datafrane with a UniProt column, and extracts the distinct UniProts
################################################################################
extract_unique_uniprots <- function(dataframe) {
  return(dataframe %>%
           dplyr::select(UniProt) %>% dplyr::distinct(UniProt))
}

################################################################################
#' Strip underscores and following chars from string
#'
#' @param str a string to strip
#'
#' @return the substring that comes before the underscore
################################################################################
take_only_pre_underscore_substring <- function(str) {
  return(strsplit(str, "_")[[1]][1])
}

################################################################################
#' Function to read password-protected xlsx file with a prompt
#'
#' @param filepath path to excel file
#' @param sheet index or name of the sheer
#'
#' @return A dataframe of the sheet
################################################################################

read_password_xlsx <- function(filepath, sheet) {
  wb <- XLConnect::loadWorkbook(filepath,
              password = rstudioapi::askForPassword("Enter password for file:"))

  data <- XLConnect::readWorksheet(wb, sheet = sheet)

  return(data)
}
