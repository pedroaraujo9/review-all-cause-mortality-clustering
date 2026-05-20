library(tidyverse)
library(patchwork)
library(clustmrev)

#### papers publish years ####
paper_year = c(
  2002, 
  2006, 
  2008, 
  2011, 
  2013, 
  2014, 
  2017, 
  2018, 
  2019, 2019, 
  2020, 
  2021, 2021, 2021, 2021, 
  2022, 2022, 2022, 
  2023, 2023, 2023,
  2024, 2024, 2024, 2024, 
  2025, 2025
)

length(paper_year)
median(paper_year)

#### read data ####
data_95 = readRDS("../data/both_1960_2010_0_95.rds")
data_110 = readRDS("../data/both_1960_2010_0_110.rds")

data_male_95 = readRDS("../data/male_1960_2010_0_95.rds")
data_male_110 = readRDS("../data/male_1960_2010_0_110.rds")

data_female_95 = readRDS("../data/female_1960_2010_0_95.rds")
data_female_110 = readRDS("../data/female_1960_2010_0_110.rds")

#### fits ####
review_95_fit = fit_review_methods(data_95)
review_110_fit = fit_review_methods(data_110)

review_fit = review_95_fit

review_male_95_fit = fit_review_methods(data_male_95)
review_female_95_fit = fit_review_methods(data_female_95)

#### compare fits ####
comp_95_110 = compare_fits(review_95_fit, review_110_fit)
comp_95_110$ARI_plot

comp_sexes = compare_fits(review_male_95_fit, review_female_95_fit)
comp_sexes$ARI_plot

#### internal quality metrics ####
review_female_95_fit$internal_metrics_plot
review_male_95_fit$internal_metrics_plot
review_fit_95$internal_metrics_plot

review_female_95_fit$internal_metrics_plot
review_male_95_fit$internal_metrics_plot

plot_clust_graph(review_fit_95, clust_number = c(2, 3, 4), seed = c(2, 11, 3))
plot_clust_graph(review_male_95_fit, clust_number = c(2, 3, 4), seed = c(1, 4, 4))
plot_clust_graph(review_female_95_fit, clust_number = c(2, 3, 4), seed = c(8, 2, 1))

#### analysis ####
devtools::load_all()
hc = analyse_hell_complete(review_fit = review_fit_95, clust_size = 3)
hc_male = analyse_hell_complete(review_fit = review_fit_95, clust_size = 3)
hc_female = analyse_hell_complete(review_fit = review_fit_95, clust_size = 3)

hc$d_cluster
hc$dendrogram

analyse_ilc_kmeans(review_fit = review_fit, clust_size = 3)
analyse_func_kmeans(review_fit = review_fit, clust_size = 2)
analyse_pca_fuzzy(review_fit = review_fit, clust_size = 2)

