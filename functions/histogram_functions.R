################################################################################
# Function to create histogram plots
#
# @param data - dataframe to be plotted
# @param column_names - a list containing the column names from which the x axis,
#                fill and facet will take the data.
#                The fill is optional in case you want to split the data into
#                several groups.
#                The facet is optional in case you want to split into several
#                histograms each for a specific value from that column.
#                Should be of the format `list(x = ?, fill = ?, facet =?)`
# @param titles - a list containing the text for the titles of the x axis and
#                 fill legend. Same format as column_names (without facet).
# @param plot_title_suffix - the plot title will use the previous to populate a
#                 standard title format but will need a suffix to differentiate
#                 similar plots
# @param should_density - a boolean whether the density should be plotted.
#                         Default is FALSE. A count density requires binwidth
#                         to be defined..(?)
# @param type - count or frequency. Default is count.
# @param binwidth - the binwidth. The default is NULL, in which case this will
#                   default to whatever produces 30 bins.
#
# @return the histogram plot
################################################################################
make_histogram_plots <- function(data, column_names, titles, plot_title_suffix = NULL, should_density = FALSE, type = "count", binwidth = NULL) {
  are_there_facets <- "facet" %in% names(column_names)
  histogram <- generate_ggplot_object(data, column_names) %>%
    append_geom_histogram(data, column_names, type, binwidth) %>%
    {`if`(should_density, append_geom_density(., nrow(data), type, binwidth), .)} %>%
    {`if`(are_there_facets, append_geom_facet(., column_names[["facet"]]), .)} %>%
    append_histogram_labels(titles, plot_title_suffix, type, are_there_facets)
  return(histogram)
}

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
