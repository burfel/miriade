#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#' Map uniprot to hgnc symbol and add it as column to the df
#'
#' @param df_with_uniprot_column
#' @param new_column_name
#'
#' @return The same `df_with_uniprot_column` with an added column for hgnc
#'         symbol given the name 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
map_uniprot_to_hgnc_and_cbind <- function(df_with_uniprot_column, new_column_name)
{
  require('biomaRt')
  if(!exists("ensembl_h")) {
    # Assign to global environment for caching purposes
    ensembl_h <<- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  }
  id_mapping <- getBM(attributes = c("uniprotswissprot", "hgnc_symbol"),
                      filters = "uniprotswissprot",
                      values = df_with_uniprot_column,
                      mart = ensembl_h)
  
  id_mapping <- dplyr::filter(id_mapping, hgnc_symbol != "")
  
  hgnc <- named_list(new_column_name, sapply(
    df_with_uniprot_column$UniProt,
    function(x) with(id_mapping, hgnc_symbol[uniprotswissprot == x][1])))
  
  return(cbind(hgnc, df_with_uniprot_column))
}
