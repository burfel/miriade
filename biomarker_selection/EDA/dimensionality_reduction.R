library(tidyverse)
library(missForest)
library(here)
source(here("functions", "data_processing_functions.R"))

datasets_root_directory <- define_datasets_root()

# Read data
kth <- read.table(file.path(datasets_root_directory, "KTH/KTH AD dataset for MIRIADE biomarker selection sample info and results.csv"),
                  sep = ";", header = T, quote = "")

# Remove non-numeric columns
kth_num <- kth %>%
  dplyr::select(-c(Class, Diagnosis, Gender, ApoE4))

# Impute missing values
kth_imp <- missForest(kth_num)

# Perform PCA
kth_pca <- prcomp(kth_imp$ximp)

# Visualize the results
biplot(kth_pca)


library(ggplot2)
# Scale the data
kth_scaled <- scale(kth_imp$ximp)

# Perform PCA
kth_pca <- prcomp(kth_scaled, center = TRUE, scale. = TRUE)

kth_twocomp_df <- cbind(as.data.frame(kth_pca[["x"]][,1:2]), Class = kth$Class)

# Visualize the results using ggbiplot
# plot PCA results
ggplot(kth_twocomp_df, aes(x = PC1, y = PC2, color = Class)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("blue", "red")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray35") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray35") +
  labs(x = paste0("PC1 (", round(kth_twocomp_df$PC1, 2) * 100, "%)"),
       y = paste0("PC2 (", round(kth_twocomp_df$PC2, 2) * 100, "%)"),
       title = "PCA plot using ggplot2",
       color = "Class") +
  theme_minimal()
