#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
## Project: MIRIADE
## Script purpose: Compare the non cleaned datasets
## Date: 13.03.2023
## Author: Felicia Burtscher
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
library(here)
library(dplyr)
library(limma)
source(here("functions", "mapping_functions.R"))
source(here("functions", "comparison_functions.R"))
source(here("biomarker_selection", "EDA", "prep_datasets_for_comparisons.R"))
# # Install
# if(!require(devtools)) install.packages("devtools")
# devtools::install_github("kassambara/ggpubr")
# install.packages("ggpubr")
# library("dplyr")
# library("ggpubr")
# ggqqplot(olink)
library(ggplot2)


###################################. NORMALITY TESTS #############################

###################################. OLINK #######################################
# # histogram <- ggplot(olink, aes(x=!!a, fill=!!b, color=!!b)) +
# #   geom_histogram(position="identity", alpha=0.5)
# histogram <- ggplot(olink +
#   geom_histogram(position="identity", alpha=0.5))

ggplot(olink, aes(sample = AA)) + geom_qq() +
geom_qq_line()

olink <- olink[order(olink$Diagnosis,decreasing=FALSE),]
olink_filtered <- olink[is.element(olink$Diagnosis, c('Control','AD','DLB', 'FTD')),]

olink_numbers <- olink_filtered
olink_numbers <- olink_numbers[, -c(1:5)]

# names(olink_numbers) <- NULL
# non_numeric_index <- which(!is.numeric(olink_numbers))

shapiro.test(as.numeric(olink_numbers[,1]))
shapiro_tests = lapply(seq(1,length(olink_numbers)),function(x){shapiro.test(as.numeric(olink_numbers[,x]))})
pvals = sapply(shapiro_tests, function(x){x$p.value})
number_of_nonsignificance <- length(pvals[pvals > 0.05])

###################################. LINEAR MODELS #############################
lm_res <- limma_gen(m, model.str = "~ AGE", coef.str = "AGE")
head(arrange(lm_res, adj.P.Val)) # top 6 rows sorted by adjusted p-value

# Conduct a simple linear regression between the exposure and the outcome
simple_model <- lm(emif$Diagosis ~ exposure, data = emif)

# Check for confounding by including a potential confounder in the model
adjusted_model <- lm(outcome ~ exposure + confounder, data = raw_data)

# Compare the coefficients for the exposure variable in both models
summary(simple_model)$coefficients
summary(adjusted_model)$coefficients


############# DOESN'T LOAD ###################
# Load the RforProteomics package
# 
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("RforProteomics")
library(RforProteomics)

# Create a histogram to visualize the data distribution
hist(emif)

# Create a Q-Q plot to compare the dataset to a normal distribution
qqnorm(emif)
qqline(emif, col = "red")

# Perform the Shapiro-Wilk test
shapiro.test(emif)

# Check the standardized residuals for normal distribution
check_normality(emif)
