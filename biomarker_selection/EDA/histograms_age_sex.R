library(here)
source(here("functions", "histogram_functions.R"))

### read the VUMC Olink cohort file
vumc_olink_path <- "/Users/felicia.burtscher/Documents/UL/DATASETS/Olink/OLINK basic data files/clin_bl_basic_04022020.csv"
olink_cohort <- read.table(vumc_olink_path, sep = ",", header = T)

# select the three columns and filter out anything not the 3 diseases (AD, DLB or FTD)
age_sex_per_disease <- dplyr::select(olink_cohort, dx, age_gr, sex) %>%
  dplyr::filter(dx == "AD dementia" | dx == "DLB" | dx == "FTD") %>%
  dplyr::mutate(sex = case_when(sex == 1 ~ 'm', sex == 2 ~ 'f', TRUE ~ sex))

# Plot 1: distribution or histogram: age across diseases (x-axis: age, y-axis: #patients), per disease (diseases in different colour)
age_binwidth <- 5

age_freq_hist <- make_histogram_plots(age_sex_per_disease,
                                      column_names = list(x="age_gr", fill="dx"),
                                      titles = list(x="Age", fill="Disease"),
                                      plot_title_suffix = "for Olink",
                                      should_density = TRUE,
                                      type = "frequency",
                                      binwidth = age_binwidth)
age_freq_hist

age_count_hist <- make_histogram_plots(age_sex_per_disease,
                                       column_names = list(x="age_gr", fill="dx"),
                                       titles = list(x="Age", fill="Disease"),
                                       plot_title_suffix = "for Olink",
                                       should_density = TRUE,
                                       type = "count",
                                       binwidth = age_binwidth)
age_count_hist

# Plot 2: same as plot 1, but dissected according to sex, so 2 curves per disease
sex_hist <- make_histogram_plots(age_sex_per_disease,
                                 column_names = list(x="sex", fill="dx"),
                                 titles = list(x="Sex", fill="Disease"),
                                 plot_title_suffix = "for Olink")
sex_hist

# Plot 2.5: according to age and sex
age_and_sex_freq_hist <- make_histogram_plots(age_sex_per_disease,
                                              column_names = list(x="age_gr", fill="dx", facet="sex"),
                                              titles = list(x="Age", fill="Disease", facet="Sex"),
                                              plot_title_suffix = "for Olink",
                                              should_density = TRUE,
                                              type = "frequency",
                                              binwidth = age_binwidth)
age_and_sex_freq_hist

age_and_sex_count_hist <- make_histogram_plots(age_sex_per_disease,
                                               column_names = list(x="age_gr", fill="dx", facet="sex"),
                                               titles = list(x="Age", fill="Disease", facet="Sex"),
                                               plot_title_suffix = "for Olink",
                                               should_density = TRUE,
                                               type = "count",
                                               binwidth = age_binwidth)
age_and_sex_count_hist

# Samey for KTH
kth_path <- "/Users/felicia.burtscher/Documents/UL/DATASETS/KTH/KTH AD dataset for MIRIADE biomarker selection sample info and results.csv"
kth_age_sex_disease <- read.table(kth_path, sep = ";", fileEncoding = "UTF-8",
                                  header = T) %>%
  dplyr::select(Class, Diagnosis, Age, Gender, Tau, Ptau, Ab42, ApoE4, MMSE_score)

kth_age_count_hist <- make_histogram_plots(kth_age_sex_disease,
                                           column_names = list(x="Age", fill="Diagnosis"),
                                           titles = list(x="Age", fill="Diagnosis"),
                                           plot_title_suffix = "for KTH",
                                           should_density = TRUE,
                                           binwidth = age_binwidth)
kth_age_count_hist

kth_sex_hist <- make_histogram_plots(kth_age_sex_disease,
                                     column_names = list(x="Gender", fill="Diagnosis"),
                                     titles = list(x="Sex", fill="Diagnosis"),
                                     plot_title_suffix = "for KTH")
kth_sex_hist

kth_age_and_sex_count_hist <- make_histogram_plots(kth_age_sex_disease,
                                                   column_names = list(x="Age", fill="Diagnosis", facet="Gender"),
                                                   titles = list(x="Age", fill="Diagnosis", facet="Sex"),
                                                   plot_title_suffix = "for KTH",
                                                   should_density = TRUE,
                                                   type = "count",
                                                   binwidth = age_binwidth)
kth_age_and_sex_count_hist

#Samey for emif
emif_path <- "/Users/felicia.burtscher/Documents/UL/DATASETS/EMIF-AD MBD Study/Clinical_EMIFMBD_MIRIADE.csv"
emif_age_sex_disease <- read.table(emif_path, sep = ";", fileEncoding = "UTF-8",
                                   header = T, dec = ",") %>%
  dplyr::filter(Age != "" & Diagnosis !="") %>%
  dplyr::mutate(Gender = case_when(Gender == 0 ~ 'm', Gender == 1 ~ 'f', TRUE ~ as.character(Gender)))

emif_age_count_hist <- make_histogram_plots(emif_age_sex_disease ,
                                            column_names = list(x="Age", fill="Diagnosis"),
                                            titles = list(x="Age", fill="Diagnosis"),
                                            plot_title_suffix = "for Emif",
                                            should_density = TRUE,
                                            binwidth = age_binwidth)
emif_age_count_hist

emif_sex_hist <- make_histogram_plots(emif_age_sex_disease,
                                      column_names = list(x="Gender", fill="Diagnosis"),
                                      titles = list(x="Sex", fill="Diagnosis"),
                                      plot_title_suffix = "for Emif")
emif_sex_hist

emif_age_and_sex_count_hist <- make_histogram_plots(emif_age_sex_disease,
                                                    column_names = list(x="Age", fill="Diagnosis", facet="Gender"),
                                                    titles = list(x="Age", fill="Diagnosis", facet="Sex"),
                                                    plot_title_suffix = "for Emif",
                                                    should_density = TRUE,
                                                    type = "count",
                                                    binwidth = age_binwidth)
emif_age_and_sex_count_hist
