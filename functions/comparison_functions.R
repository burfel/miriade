#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#' Compare datasets and return list with the bound datasets, intersection list
#' and unique proteins
#'
#' @param datasets The datasets
#' @param dataset_names The names to give the datasets
#' @param should_filter_by_diseases logical vector of same length as `datasets`
#' @param additional_columns A single vetor, list of vectors of same length as
#'        `datasets` or NULL. If none desired it is NULL by default and can be ignored.
#'        For non NULL, it is used for column names each dataset should preserve
#'        in addition to Uniprot and p value.
#'        If it is a single vector, it is used across all datasets.
#'        If it is a list of vectors, they each match a dataset.
#' @param disease_to_filter_by If any of should_filter_by_disease is true,
#'        disease_to_filter_by should be a
#'        string that appears in the disease column
#'
#' @return A named list with the components of `raw` (the bound datasets),
#'         `intersection_list` and `unique_proteins`
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
compare_datasets <- function(datasets, dataset_names,
                             should_filter_by_disease,
                             additional_columns = NULL,
                             disease_to_filter_by = NULL)
{
  datasets_length <- length(datasets)
  if (datasets_length != length(dataset_names)) {
    stop("datasets and dataset_names must be the same length")
  }
  preprocessed_datasets <- preprocess_datasets(
    datasets = datasets,
    should_filter_by_disease = should_filter_by_disease,
    additional_columns = additional_columns,
    disease_to_filter_by = disease_to_filter_by)

  bound_datasets <- bind_datasets_together(preprocessed_datasets) %>%
    dplyr::mutate(Source = case_source(Source, dataset_names))

  intersection_list <- filter_according_to_UniProt_masks(bound_datasets,
                                                         preprocessed_datasets)

  # If intersection list is empty, no point in continuing
  if (nrow(intersection_list) == 0) {
    return(list(
      dataset_names = dataset_names,
      raw = bound_datasets,
      intersection_list = intersection_list))
  }
  unique_proteins <-
    extract_unique_uniprots(intersection_list,
                            additional_columns)

  significant_in_all <-
    join_datasets(preprocessed_datasets) %>%
    filter_joined_datasets(datasets_length) %>%
    rename_pvalue_columns(dataset_names)

  return(list(
    dataset_names = dataset_names,
    raw = bound_datasets,
    intersection_list = intersection_list,
    unique_proteins = unique_proteins,
    significant_in_all = significant_in_all))
}

# Auxiliary function for compare_datasets
case_source <- function(Source, dataset_names)
{
  if (length(dataset_names) == 2) {
    return(two_sources(Source, dataset_names))
  } else if(length(dataset_names) == 3) {
    return(three_sources(Source, dataset_names))
  } else {
    return(four_sources(Source, dataset_names))
  }
}

# Auxiliary function for compare_datasets
two_sources <- function(Source, dataset_names)
{
  return(case_when(
    Source == 1 ~ dataset_names[[1]],
    Source == 2 ~ dataset_names[[2]]))
}

# Auxiliary function for compare_datasets
three_sources <- function(Source, dataset_names)
{
  return(case_when(
    Source == 1 ~ dataset_names[[1]],
    Source == 2 ~ dataset_names[[2]],
    Source == 3 ~ dataset_names[[3]]))
}

# Auxiliary function for compare_datasets
four_sources <- function(Source, dataset_names)
{
  return(case_when(
    Source == 1 ~ dataset_names[[1]],
    Source == 2 ~ dataset_names[[2]],
    Source == 3 ~ dataset_names[[3]],
    Source == 4 ~ dataset_names[[4]]))
}

# Auxiliary function for compare_datasets
join_datasets <- function(datasets)
{
  datasets_length <- length(datasets)
  joined <- dplyr::inner_join(datasets[[1]], datasets[[2]],
                    by = "UniProt",
                    suffix = c("_1", "_2"))
  if(datasets_length > 2) {
    joined <- dplyr::inner_join(joined, datasets[[3]],
                                by = "UniProt") %>%
      dplyr::rename(!!sym("p.val_3") := "p.val")
  }
  if(datasets_length > 3) {
    joined <- dplyr::inner_join(joined, datasets[[4]],
                                by = "UniProt") %>%
      dplyr::rename(!!sym("p.val_4") := "p.val")
  }
  return(joined)
}

# Auxiliary function for compare_datasets
filter_joined_datasets <- function(joined_datasets, number_of_original_datasets)
{
  if (number_of_original_datasets == 2) {
    return(dplyr::filter(joined_datasets,
                         p.val_1 <= 0.05 & p.val_2 <=0.05))
  } else if (number_of_original_datasets == 3) {
    return(dplyr::filter(joined_datasets,
                         p.val_1 <= 0.05 & p.val_2 <=0.05 & p.val_3 <=0.05))
  } else if (number_of_original_datasets == 4) {
    return(dplyr::filter(joined_datasets,
                         p.val_1 <= 0.05 & p.val_2 <=0.05 & p.val_3 <=0.05 & p.val_4 <=0.05))
  }
}



# Auxiliary function for compare_datasets
rename_pvalue_columns <- function(joined_datasets, dataset_names)
{
  number_of_original_datasets <- length(dataset_names)
  renamed_df <- joined_datasets %>%
    dplyr::rename(!!suffixed_pval_symbol(dataset_names[[1]]) := "p.val_1") %>%
    dplyr::rename(!!suffixed_pval_symbol(dataset_names[[2]]) := "p.val_2")
  if (number_of_original_datasets == 3) {
    renamed_df <- renamed_df %>%
      dplyr::rename(!!suffixed_pval_symbol(dataset_names[[3]]) := "p.val_3")
  } else if (number_of_original_datasets == 4) {
    renamed_df <- renamed_df %>%
      dplyr::rename(!!suffixed_pval_symbol(dataset_names[[4]]) := "p.val_4")
  }
  return(renamed_df)
}

# Auxiliary function for rename_pvalue_columns
suffixed_pval_symbol <- function(suffix)
{
  return(sym(paste("p.val_", suffix, sep = "")))
}
