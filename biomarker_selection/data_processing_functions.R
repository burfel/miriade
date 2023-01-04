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
