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
  all_values <- as.character(union(df[[object_columns[[1]]]], df[[object_columns[[2]]]]))
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
  # Depending on the OS, use appropriate method
  # For Windows and Mac we 
  # Open a file dialog to store the root directory of all data sets
  # For other unix the user needs to enter the path
  if (.Platform$OS.type == "windows") {
    folder <- choose.dir(caption = "Choose the root of all datasets")
  }
  else if (.Platform$OS.type == "unix") {
    folder <- rChoiceDialogs::rchoose.dir()
  }
  else {
    stop("Unsupported operating system")
  }
  
  return(folder)
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

read_password_xlsx <- function(filepath, sheet = 1) {
  wb <- XLConnect::loadWorkbook(filepath,
              password = rstudioapi::askForPassword("Enter password for file:"))

  data <- XLConnect::readWorksheet(wb, sheet = sheet)

  return(data)
}

################################################################################
#' Title
#'
#' @param table 
#' @param symbol_to_cut 
#' @param target_symbol 
#'
#' @return
#' @export
#'
#' @examples
################################################################################
convert_age_column_to_age_group_column <-
  function(table,
           symbol_to_cut,
           target_symbol,
           cutting_breaks = c(40, 66, 81, 116),
           cutting_labels = c("40-65", "66-80", "81-115")) {
  return(
    mutate(
      table,
      !!target_symbol := 
        cut(!!symbol_to_cut,
            breaks = cutting_breaks,
            labels = cutting_labels)))
}

################################################################################
#' Named list of collections
#' 
#' Creates a named list of collections
#'
#' @param names_vector vector of names to give the collections
#' @param ... variable number of collections
#'
#' @return a named list with the collections
#' @note Throws an exception if the number of names is not the same as
#'       the number of collections
################################################################################
named_list <- function(names_vector, ...) {
  collection_list <- list(...)

  # Check that the length of names is equal to the number of collections/dfs
  if (length(names_vector) != length(collection_list)) {
    stop("Number of names and collections must match")
  }

  names(collection_list) <- names_vector

  return(collection_list)
}

#' Read dataframe from a table file
#'
#' @param file_path Path to the file
#' @param header Whether the first row is a header
#'
#' @return Dataframe based on the file
read_df_from_file <- function(file_path, header = TRUE)
{
  if (endsWith(file_path, ".tsv")) {
    return(read.table(file = file_path, sep = "\t", header = header))
  } else if (endsWith(file_path, ".csv")) {
    return(tryCatch({
      read.table(file = file_path, sep = ",", header = header)
    }, error = function(e) {
      read.table(file = file_path, sep = ";", dec = ",", header = header)
    }))
  } else if (endsWith(file_path, ".xls")) {
    readxl::read_xls(file_path)
  } else if (endsWith(file_path, ".xlsx")) {
    return(tryCatch({
      readxl::read_xlsx(path = file_path)
    }, error = function(e) {
      read_password_xlsx(filepath = file_path)
    }))
  }
}

################################################################################
#' Meld together pieces of a dataset spread across several files
#'
#' @param dataset_paths 2 or 3 paths of dataset parts
#' @param by_columns A named vector where the names are the names of the columns
#'        on the "left" (i.e. lower indexed) file and the values are the names
#'        of the columns on the "right" (higher indexed) file
#' @param select_columns Either NULL or a vector of the columns to be selected.
#'        If explicitly set, will only return the specified cols
#' @param exclude_columns Either NULL or a list of vectors of col names with the
#'        same length as `dataset_paths`.
#'        This one does the exclusion before the join because same named columns
#'        will have their names suffixed after the join, possibly being missed
#'        by the exclusion as a result if it was done after the join.
#'        If explicitly set, will return all cols but those
#'
#' @return
#' @example If we have a data set with "cc1,d,e" columns, one with "c1,c2,c3"
#'          and a third one with "cc2,f,g", we could meld them like this:
#'          `melded <- meld_fractured_dataset(`
#'          ` dataset_paths = c("path1","path2","path3"),`
#'          ` by_columns(cc1 = "c1", c2 = "cc2"))`.
#'          To use exclusion in this case, here is an example:
#'          `exclude_columns = list(c("d"), c(), c("f","g"))`
#' @note If both `select_columns` and `exclude_columns` are set, `select_columns`
#'       will take precedence
################################################################################
meld_fractured_dataset <- function(
    dataset_paths,
    by_columns,
    select_columns = NULL,
    exclude_columns = NULL)
{
  # Dependencies
  require(dplyr)
  # Validation
  if(length(dataset_paths) < 2 | length(dataset_paths) > 3)
  {
    stop("You must only use this for 2 or 3 data files. No more, no less")
  }
  if (length(dataset_paths) - 1 != length(by_columns))
  {
    stop("The length of `by_columns` must be EXACTLY one less than that of `dataset_paths`")
  }
  if(!is.null(exclude_columns)) {
    # If both `select_columns` and `exclude_columns` are set NULLify
    # `exclude_columns` since `select_columns` takes precedence
    if(!is.null(select_columns)) {
      exclude_columns <- NULL
    } else if(length(dataset_paths) != length(exclude_columns)) {
      stop("The length of `exclude_columns` must be EXACTLY that of `dataset_paths`")
    }
  }
  # Melding
  df1 <- read_df_from_file(dataset_paths[[1]])
  df2 <- read_df_from_file(dataset_paths[[2]])
  if(!is.null(exclude_columns)) {
    df1 <- dplyr::select(df1, -exclude_columns[[1]])
    df2 <- dplyr::select(df2, -exclude_columns[[2]])
  }
  new_df <- dplyr::inner_join(df1, df2,
                              by = by_columns[1])
  if (length(dataset_paths) == 3) {
    df3 <- read_df_from_file(dataset_paths[[3]])
    if(!is.null(exclude_columns)) {
      df3 <- dplyr::select(df3, -exclude_columns[[3]])
    }
    new_df <- dplyr::inner_join(new_df, df3,
                                by = by_columns[2])
  }
  if(!is.null(select_columns)) {
    new_df <- dplyr::select(new_df, select_columns)
  }
  return(new_df)
}

################################################################################
#' Perform KW Wilcoxon process for grouping based p values
#'
#' @param melted_df The melted representation of the data
#' @param differentiating_feature_symbol a symbol of the column name representing
#'        the p value targets
#' @param grouping_feature_symbol a symbol representing the column name for the
#'        grouping to have the p values calculated for
#' @param measurement_symbol a symbol representing the column name for the values
#' @param significance_limit significance limit for the intermediate p values
#'
#' @return A df with rows for the values of the differentiating feature with
#'          values and adjusted pvalues for each one of the grouping comparisons
################################################################################
perform_kw_wilcoxon_according_to_grouping <- function(
    melted_df,
    differentiating_feature_symbol,
    grouping_feature_symbol,
    measurement_symbol,
    significance_limit
    )
{
  library(here)
  source(here("functions", "auxiliary_statistical_functions.R"))
  ### Run Kruskal-Wallis test for all biomarkers grouped by grouping_feature
  kw_grouped <-
    sapply(
      unique(melted_df[[as.character(differentiating_feature_symbol)]]),
      kruskal_wallis_test_for_sample_groups,
      melted_df,
      differentiating_feature_symbol,
      grouping_feature_symbol,
      measurement_symbol)

  significant_kw_grouped <-
    keep_only_significant_entries(
      melted_df,
      differentiating_feature_symbol = differentiating_feature_symbol,
      kw_per_differentiating_feature = kw_grouped,
      significance_limit = significance_limit)

  if (nrow(significant_kw_grouped) == 0) {
    message("No significant pvalues were found for the ",
            as.character(grouping_feature_symbol),
            " grouping")
    return(data.frame(Name = character(),
                      p.val = numeric(),
                      adj.p = numeric()))
  }

  kw_wilcoxon_grouped <-
    test_pairwise_wilcoxon(
      melted_df,
      differentiating_feature_symbol,
      significant_kw_grouped,
      grouping_feature_symbol,
      measurement_symbol)

  ### Adjust for measurements
  kw_wilcoxon_grouped_adj <- matrix(p.adjust(kw_wilcoxon_grouped, method = "BH"), nrow = nrow(kw_wilcoxon_grouped))
  rownames(kw_wilcoxon_grouped_adj) <- rownames(kw_wilcoxon_grouped)
  colnames(kw_wilcoxon_grouped_adj) <- colnames(kw_wilcoxon_grouped)
  
  # Use lapply() to create a list of data frames for each row in kw_wilcoxon_grouped
  dfs <- lapply(1:nrow(kw_wilcoxon_grouped), function(i) {
    data.frame(Name = colnames(kw_wilcoxon_grouped),
               p.val = t(kw_wilcoxon_grouped)[, i],
               adj.p = kw_wilcoxon_grouped_adj[i, ]) %>% 
      dplyr::rename(!!differentiating_feature_symbol := "Name")
  })
  
  # Name the list elements after the row names in kw_wilcoxon_grouped
  names(dfs) <- row.names(kw_wilcoxon_grouped)
  
  return(dfs)
}
