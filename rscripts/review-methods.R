library(tidyverse)
library(clustmrev)

lf_95 = readRDS("../data/both_1960_2010_0_95.rds")
lf_110 = readRDS("../data/both_1960_2010_0_110.rds")

devtools::load_all()

review_fit_95 = fit_review_methods(lf_95)
review_fit_110 = fit_review_methods(lf_110)

review_fit_95$hell_complete$dendrogram_plot
