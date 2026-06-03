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
data = readRDS("data/sex=both-minyear=1960-maxyear=2010-minage=0-maxage=90.rds")
data_male = readRDS("data/sex=male-minyear=1960-maxyear=2010-minage=0-maxage=90.rds")
data_female = readRDS("data/sex=female-minyear=1960-maxyear=2010-minage=0-maxage=90.rds")
data_90 = readRDS("data/sex=both-minyear=1990-maxyear=2010-minage=0-maxage=90.rds")
data_2019 = readRDS("data/sex=both-minyear=1990-maxyear=2019-minage=0-maxage=90.rds")

#### fits ####
review_fit = clustmrev::fit_review_methods(data)
review_male_fit = clustmrev::fit_review_methods(data_male)
review_female_fit = clustmrev::fit_review_methods(data_female)
review_90_fit = clustmrev::fit_review_methods(data_90)
review_2019_fit = clustmrev::fit_review_methods(data_2019)

#### compare fits ####
plot_clust_graph(review_fit, clust_number = c(2, 3, 4), seed = c(2, 11, 3))
plot_clust_graph(review_male_fit, clust_number = c(2, 3, 4), seed = c(15, 3, 1))
plot_clust_graph(review_female_fit, clust_number = c(2, 3, 4), seed = c(6, 2, 3))
plot_clust_graph(review_90_fit, clust_number = c(2, 3, 4), seed = c(2, 5, 2))
plot_clust_graph(review_2019_fit, clust_number = c(2, 3, 4), seed = c(4, 6, 5))

#### internal quality metrics ####
review_fit$internal_metrics_plot
ggsave("plots/internal-quality-metrics.pdf", width = 7, height = 8)
ggsave("~/Documents/GitHub/Review-mortality-clustering/ISR-sub-2/plots/results/internal-quality-metrics.pdf", width = 7, height = 8)

review_male_fit$internal_metrics_plot
review_female_fit$internal_metrics_plot

#### hellinger-complete ####
hc = analyse_hell_complete(review_fit = review_fit, clust_size = 3, new_label = c())

hc$d_cluster
ggsave("plots/H-complete-dx-curve-cluster.pdf", width = 5, height = 3)
ggsave("~/Documents/GitHub/Review-mortality-clustering/ISR-sub-2/plots/results/H-complete-dx-curve-cluster.pdf", width = 5, height = 3)

hc$dendrogram + 
  ggplot2::scale_y_continuous(limits = c(-0.11, 0.25), breaks = c(0.1, 0.25)) + 
  theme_bw()
ggsave("plots/H-complete-dendro.pdf", width = 5, height = 3)
ggsave("~/Documents/GitHub/Review-mortality-clustering/ISR-sub-2/plots/results/H-complete-dendro.pdf", width = 7, height = 4)

#### ILC-k-means ####
ilc = analyse_ilc_kmeans(
  review_fit = review_fit, clust_size = 2,
  ages = c(0, 15, 45, 80), periods = c(1960, 1995, 2010)
)
ilc$beta_plot
ggsave("plots/ILC-k-means-beta.pdf", width = 5, height = 3) 
ggsave("~/Documents/GitHub/Review-mortality-clustering/ISR-sub-2/plots/results/ILC-k-means-beta.pdf", width = 5, height = 3) 

ilc$mx_over_time_plot + theme(text = element_text(size = 15))
ggsave("plots/ILC-k-means-time.pdf", width = 10, height = 5)
ggsave("~/Documents/GitHub/Review-mortality-clustering/ISR-sub-2/plots/results/ILC-k-means-time.pdf", width = 10, height = 5)

ilc$mx_curve_plot + geom_point() + theme(text = element_text(size = 15))
ggsave("plots/ILC-k-means-curve.pdf", width = 10, height = 3)
ggsave("~/Documents/GitHub/Review-mortality-clustering/ISR-sub-2/plots/results/ILC-k-means-curve.pdf", width = 10, height = 3)

#### PCA-fuzzy #### 
country_order = data.frame(
  country = review_fit$data$country %>% unique(),
  cluster = review_fit$review_fit$pca_fuzzy$membership_level[[2]][, 2]
) %>%
  arrange(cluster) %>%
  .$country

pcaf = analyse_pca_fuzzy(review_fit = review_fit, clust_size = 2, country_order = country_order)
pcaf$prob_plot
ggsave("plots/PCA-fuzzy-membership-level.pdf", width = 7, height = 5)
ggsave("~/Documents/GitHub/Review-mortality-clustering/ISR-sub-2/plots/results/PCA-fuzzy-membership-level.pdf", width = 7, height = 5)

#### func-k-means #### 
fk = analyse_func_kmeans(review_fit = review_fit, clust_size = 2)
fk$ex_plot
ggsave("plots/funck-methods-single.pdf", width = 6, height = 3.5) 
ggsave("~/Documents/GitHub/Review-mortality-clustering/ISR-sub-2/plots/results/funck-methods-single.pdf", width = 6, height = 3.5) 

#### table #### 
data.frame(
  country = unique(review_fit$data$country),
  hc = hc$cluster, 
  ilc = ilc$cluster,
  pcaf = pcaf$cluster,
  fk = fk$cluster
) %>% 
  xtable::xtable() %>%
  print(include.rownames = F)

review_fit$review_fit$pca_fuzzy$membership_level[[2]]








