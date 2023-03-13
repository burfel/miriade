library(ggplot2)

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
  for (i in 1:new_length) { # TODO: Possibly replace with sapply
    box_plots[[i]] <- create_box_plot_of_column(dataframe, column_names[[i]])
  }
  return(box_plots)
}

create_box_plot <- function(dataframe, dataframe_name, x_column, y_column, with_points = TRUE) {
  x_column_symbol <- sym(x_column)
  y_column_symbol <- sym(y_column)
  return(
    ggplot(dataframe, aes(x = !!x_column_symbol, y = !!y_column_symbol)) +
      geom_boxplot() +
      scale_x_discrete() +
      labs(title = paste(dataframe_name, y_column)) +
      if(with_points) geom_point()
  )
}

create_box_plots_of_all_columns_starting_at_column_number <-
  function(dataframe, dataframe_name, x_column, y_starting_column, with_points = TRUE) {
    number_of_columns <- length(names(dataframe))
    y_column_names <- names(dataframe)[y_starting_column:number_of_columns]
    new_length <- number_of_columns - y_starting_column + 1
    box_plots <- vector("list", length = new_length)
    for (i in 1:new_length) { # TODO: Possibly replace with sapply
      box_plots[[i]] <- create_box_plot(dataframe, dataframe_name, x_column,
                                        y_column_names[[i]], with_points)
    }
    return(box_plots)
  }
