#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#' Kruskal-Wallis test for significant differentiating features
#' 
#' A helper function to calculate Kruskal-Wallis for all the grouping
#' features in the table, given the differentiating_feature_symbol
#' (e.g. HGNC_Symbol).
#' It means it checks the groups of samples per that grouping feature
#' (e.g. Diagnosis) and tests whether they are all from the same distribution.
#' If the resulting p value is below the chosen significance level, we reject
#' that hypothesis and say they are not from the same distribution.
#'
#' @param tested_differentiating_feature The specific value for the
#'        differentiating feature to be tested on.
#'        E.g. the HGNC_Symbol value "A1BG"
#' @param table The data table
#' @param differentiating_feature_symbol The column name symbol of the
#'        differentiating feature.
#'        E.g. sym("HGNC_Symbol")
#' @param grouping_feature_symbol The column name symbol of the feature by which
#'        the different groups are constructed.
#'        E.g. sym("Diagnosis")
#' @param measurement_symbol The column name symbol of the measurements of the
#'        differentiating feature.
#'        E.g. sym("median_ab_readout")
#'
#' @return a p-value to ascertain whether the samples of the differentiating
#'         feature from the different group come from the same distribution
#'         or whether the different groups have different distribution
#'
#' @examples
#' @example `sapply(unique(tab2$HGNC_Symbol), kruskal_wallis_test_for_sample_groups, tab2, sym("HGNC_Symbol"),`
#'                 `sym("Diagnosis"), sym("median_ab_readout"))`
#' @example `kruskal_wallis_test_for_sample_groups("A1BG", tab2, sym("HGNC_Symbol"),`
#'                  `sym("Diagnosis"), sym("median_ab_readout"))`
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
kruskal_wallis_test_for_sample_groups <- function(tested_differentiating_feature,
                                                  table,
                                                  differentiating_feature_symbol,
                                                  grouping_feature_symbol,
                                                  measurement_symbol) {
  ### This one-liner filters by differentiating_feature,
  ### splits measurement by grouping_feature
  kw_tibbles <- tibble::as_tibble(table) %>%
    dplyr::filter(!!differentiating_feature_symbol == tested_differentiating_feature) %>%
    dplyr::select(!!grouping_feature_symbol, !!measurement_symbol) %>%
    dplyr::group_split(!!grouping_feature_symbol, .keep = FALSE)
  ### We need to unlist the tibbles from the 'group_split' to run the kruskal test
  stats::kruskal.test(lapply(kw_tibbles, unlist))$p.value
}

keep_only_significant_entries <- function(melted_table,
                                          differentiating_feature_symbol,
                                          kw_per_differentiating_feature,
                                          significance_limit = 0.1) {
  ### Adjust for multiple hypothesis testing
  is_adjusted_pvalue_significant_predicate_vector <-
    is_adjusted_pvalue_significant_predicate(kw_per_differentiating_feature,
                                             significance_limit)
  unique_differentiating_values <-
    unique_table_column_values(melted_table, differentiating_feature_symbol)
  kw_biomarkers <-
    as.character(unique_differentiating_values[is_adjusted_pvalue_significant_predicate_vector])
  ### Keep only the significant entries
  return(dplyr::filter(melted_table, !!differentiating_feature_symbol %in% kw_biomarkers))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Wilcoxon pairwise comparisons for biomarkers passing the KW test
#
# A helper function to calculate pairwise Wilcoxon tests for the given list of
# group pairs. This function is specific for the variable/value names assigned
# in reshape2::melt, above
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
test_pairwise_wilcoxon <- function(melted_table,
                                   differentiating_feature_symbol,
                                   significant_kw,
                                   grouping_feature_symbol,
                                   measurement_symbol) {
  differentiating_name <- as.character(differentiating_feature_symbol)

  pairs <- unique_table_column_values(melted_table, grouping_feature_symbol) %>%
           combn(2, simplify = F) %>%
           purrr::set_names(purrr::map_chr(., ~ paste(., collapse = "_vs_")))

  unique_differentiator_values <- unique(significant_kw[[differentiating_name]])

  result <-
    sapply(
      unique_differentiator_values,
      collect_wilcoxon_pvalues_of_all_pairs_and_proteins,
      significant_kw,
      measurement_symbol,
      grouping_feature_symbol,
      differentiating_feature_symbol, pairs)
  if(any(class(result) == "numeric")) {
    result <- as.data.frame(t(result))
    rownames(result) <- names(pairs)
  }
  colnames(result) <-
    unique_differentiator_values

  return(result)
}

# AUXILIARY to test_pairwise_wilcoxon
is_adjusted_pvalue_significant_predicate <-
  function(kw_per_differentiating_feature, significance_limit) {
    return(p.adjust(kw_per_differentiating_feature, method = "BH") < significance_limit)
  }

# AUXILIARY to test_pairwise_wilcoxon
unique_table_column_values <- function(table, column_name_symbol) {
  return(unique(table[[as.character(column_name_symbol)]]))
}

# AUXILIARY to test_pairwise_wilcoxon
collect_wilcoxon_pvalues_of_all_pairs_and_proteins <-
  function(tested_differentiating_feature, df,
           measurement_symbol,
           grouping_feature_symbol, differentiating_feature_symbol, pairs) {
    differentiating_name <- as.character(differentiating_feature_symbol)
    filtered_df <- df[df[[differentiating_name]] == tested_differentiating_feature,]
    return(sapply(pairs, pvalue_from_wilcoxon_on_pair,
                  filtered_df,
                  measurement_symbol,
                  grouping_feature_symbol))
    
  }

# AUXILIARY to test_pairwise_wilcoxon
pvalue_from_wilcoxon_on_pair <- function(pair, df,
                                         measurement_symbol,
                                         grouping_feature_symbol) {
  return(perform_wilcoxon_on_pair(pair, df, measurement_symbol,
                           grouping_feature_symbol)$p.value)
}

# AUXILIARY to test_pairwise_wilcoxon
perform_wilcoxon_on_pair <- function(pair, df,
                                     measurement_symbol,
                                     grouping_feature_symbol) {
  measurement_name <- as.character(measurement_symbol)
  grouping_name <- as.character(grouping_feature_symbol)
  return(stats::wilcox.test(df[[measurement_name]][df[[grouping_feature_symbol]] == pair[1]],
                          df[[measurement_name]][df[[grouping_feature_symbol]] == pair[2]],
                          exact = F))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#' Prepare the cutting breaks and cutting labels for a conversion from discrete
#' values to groups based on quantiles
#'
#' @param df_column The column planned to be converted
#' @param desired_num_of_groups The number of groups desired in the final df
#'
#' @return A list with the `Breaks` and `Labels` to send to the `cut` method
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
prepare_breaks_and_labels_by_quantiles <- function(df_column,
                                                   desired_num_of_groups)
{
  quantiles <- quantile(df_column, probs = seq(0, 1, 1/desired_num_of_groups),
                        na.rm = TRUE)
  breaks <- c(quantiles[1:desired_num_of_groups] - 1,
              quantiles[desired_num_of_groups + 1])
  labels <- vector("character", desired_num_of_groups)
  for(i in 1:desired_num_of_groups) {
    end_value <- ifelse(i < desired_num_of_groups,
                        quantiles[[i+1]] - 1,
                        quantiles[[i+1]])
    labels[[i]] <- paste(quantiles[[i]], "-", end_value, sep = "")
  }
  return(list(Breaks = breaks, Labels = labels))
}
# TODO: TO BE CONTINUED
