library(ggplot2)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#' Plot correlation heatmap from (streched) correlation dataframe
#'
#' @param corr_df dataframe that is a result to a call to corrr's `correlate`
#'        followed by its `stretch`
#' @param dataset_name A string of the name of the dataset
#'
#' @return plot object for the correlation heatmap
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
plot_correlation_heatmap <- function(corr_df, dataset_name)
{
  corr_df %>%
    ggplot(aes(x, y, fill = r)) +
    geom_tile()  +
    labs(x = NULL, y = NULL,
         fill = "Pearson's\nCorrelation",
         title=paste("Correlations in", dataset_name)) +
    # map a red, white and blue color scale to correspond to -1:1 sequential gradient
    scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446",
                         limits=c(-1,1)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    # remove excess space on x and y axes
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0))
}