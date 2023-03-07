# 15.11.2022
# visualising the original Olink data as box plots
library(here)
source(here("functions", "data_processing_functions.R"))

datasets_root_directory <- define_datasets_root()

olink_basic <- read.table(file.path(datasets_root_directory, "Olink/OLINK basic data files/protein_data_LODmax_10percent_outliers_rm_missing_imputed_16112020.txt"),
                          sep = "\t", header = T) 

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
