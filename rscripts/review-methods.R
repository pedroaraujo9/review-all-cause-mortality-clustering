library(tidyverse)
library(NbClust)
library(cluster)
library(vegan)
library(dtw)
library(patchwork)
library(e1071)
library(vital)
library(fda)
library(ggraph)
library(tidygraph)
library(igraph)
library(mclust)
library(ggdendro)
source("rscripts/utils.R")
select = dplyr::select

#### data #### 
lt = readRDS("data/life_tables_5x1.rds")
lt %>% glimpse()

removed_countries = c(
  "Hong Kong", "Iceland", "Chile", "Croatia", 
  "Republic of Korea", "Luxembourg", "East Germany",
  "West Germany", "England and Wales (Total Population)",
  "England and Wales (Civilian Population)", "Scotland",
  "Northern Ireland", "New Zealand Maori", 
  "New Zealand Non-Maori"
)

period_range = c(1960, 2010)

lt = lt %>%
  dplyr::filter(year >= period_range[1], year <= period_range[2], 
                !(country %in% removed_countries)) 

n_years = lt$year %>% unique() %>% length()

# countries with all periods
countries_analyzed = lt %>%
  dplyr::select(country, year) %>% 
  distinct() %>%
  group_by(country) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  filter(n == n_years) %>%
  .$country

n_country = length(countries_analyzed)

lt = lt %>% filter(country %in% countries_analyzed)

countries_analyzed

age_map = data.frame(
  age_group_label = lt$age %>% unique() %>% sort(),
  age_group = 1:length(lt$age %>% unique() %>% sort())
)

period_map = data.frame(
  period_label = lt$year %>% unique() %>% sort(),
  period = 1:length(lt$year %>% unique() %>% sort())
)

ireland_log_mx_plot = lt %>%
  filter(country == "Ireland") %>%
  ggplot(aes(x=age, y=log(mx), color=year, group=year)) + 
  geom_line() +
  guides(color="none") + 
  viridis::scale_color_viridis() + 
  labs(x="Age group", 
       y=latex2exp::TeX("$\\log($${}_{5} m_{x})$ (log central mortaliry rate)"), 
       color="Year") + 
  theme(text = element_text(size = 12))

ireland_log_mx_plot
ggsave("plots/ireland_log_mx.pdf", width = 4.5, height = 3)

ireland_dx_plot = lt %>%
  filter(country == "Ireland") %>%
  ggplot(aes(x=age, y=dx, color=year, group=year)) + 
  geom_line() + 
  viridis::scale_color_viridis() + 
  labs(x="Age group", 
       y=latex2exp::TeX("${}_{5}d_{x}$ (death count distribution)"), 
       color="Year") + 
  theme(text = element_text(size = 12))

ireland_dx_plot
ggsave("plots/ireland_dx.pdf", width = 5, height = 3)


ireland_log_mx_plot + ireland_dx_plot
ggsave("plots/ireland_log_mx_dx.pdf", width = 9, height = 3)

mx_tidy = lt %>% dplyr::select(country, year, age, mx)
mx = mx_tidy %>% spread(age, mx)
mx_matrix = mx %>% dplyr::select(-country, -year) %>% as.matrix()

#### Helling-Ward distance ####
n = lt$country %>% unique() %>% length()
time_unique = lt$year %>% unique()
dH = matrix(0, nrow = n, ncol = n)

for(i in seq_along(time_unique)) {
  
  d_time = lt %>%
    filter(year == time_unique[i]) %>%
    mutate(dx_norm = dx/100000) %>% #divide by l0 = 100000
    select(country, age, dx_norm) %>%
    spread(age, dx_norm) %>%
    select(-country) %>%
    as.matrix() %>%
    sqrt() %>%
    dist() %>%
    as.matrix() %>%
    `/`(sqrt(2))
  
  dH = dH + d_time/length(time_unique)
  
}

colnames(dH) = rownames(dH) = lt$country %>% unique()
dH = as.dist(dH)

ward_fit = hclust(dH, method = "ward.D")
ggdendrogram(ward_fit)
ggsave("plots/dendro.pdf", width = 6, height = 3)

h_ward_fit = fit_compare(
  diss = dH, G = 2:10, seed = 1, method = "ward"
)

h_ward_fit$metrics_plot
ggsave("plots/h-ward-metrics.pdf", width = 8, height = 3)

h_ward_class_matrix = h_ward_fit$class_matrix
h_ward_class = h_ward_class_matrix[, 1]

lt %>%
  mutate(dx_norm = sqrt(dx/100000)) %>%
  select(country, year, age, dx_norm) %>%
  mutate(id = paste0(country, "-", year)) %>%
  mutate(class = h_ward_class[country]) %>%
  as_tibble() %>%
  group_by(class, age) %>%
  summarise(med = median(dx_norm),
            li = quantile(dx_norm, 0.025),
            ui = quantile(dx_norm, 0.975)) %>%
  ggplot(aes(x=age, y=med, color=factor(class))) + 
  geom_line(linewidth = 1) + 
  geom_ribbon(aes(x=age, ymin = li, ymax=ui, fill=factor(class)), 
              inherit.aes = F, alpha = 0.2) + 
  labs(x="Age group", y=latex2exp::TeX("Median $d_{x i}^{*}$"), 
       fill = "Cluster", color = "Cluster")

ggsave("plots/dx-curve-cluster.pdf", width = 5, height = 3)

#### ILC-k-means ####
ILC = lapply(countries_analyzed, function(country_name){
  fit = mx_tidy %>%
    filter(country == country_name) %>%
    dplyr::select(year, age, mx) %>%
    dplyr::arrange(age, year) %>%
    as_tibble() %>%
    as_vital(index = year, key = age, .age = "age") %>%
    model(lee_carter = LC(log(mx), scale = T)) 
  
  list(kt = fit %>% time_components() %>% .$kt, 
       bx = fit %>% age_components() %>% .$bx)
})

kappa_t = ILC %>% purrr::map(~{.x$kt}) %>% do.call(rbind, .) 
beta_x = ILC %>% purrr::map(~{.x$bx}) %>% do.call(rbind, .)

ILC_k_means_fit = fit_compare(
  data = beta_x, 
  diss = dist(beta_x), 
  G = 2:10, 
  n_start = 100, 
  n_iters = 1000, 
  method = "kmeans",
  seed = 1
)

ILC_k_means_fit$metrics_plot
ggsave("plots/ILC-kmeans-metrics.pdf", width = 8, height = 3)

ILC_kmeans_class_matrix = ILC_k_means_fit$class_matrix
ILC_kmeans_class = ifelse(ILC_kmeans_class_matrix[, 1] == 1, 2, 1)

beta_x %>%
  as.data.frame() %>%
  mutate(country = countries_analyzed, class = factor(ILC_kmeans_class)) %>%
  gather(age_group, bx, -country, -class) %>%
  as_tibble() %>%
  mutate(age_group = age_group %>% str_extract("\\d{1,10}") %>% as.numeric()) %>%
  left_join(age_map, by = "age_group") %>%
  ggplot(aes(x=age_group_label, y=bx, group=country, color=class)) + 
  geom_line() + 
  labs(x="Age group", y=latex2exp::TeX("$\\beta_{x i}$"), color="Cluster") + 
  geom_hline(yintercept = 0, linetype = "dashed") 

ggsave("plots/ILC-k-means-beta.pdf", width = 5, height = 3) 

kappa_t %>%
  as.data.frame() %>%
  mutate(country = countries_analyzed, class = factor(ILC_kmeans_class)) %>%
  gather(period, kappa, -country, -class) %>%
  as_tibble() %>%
  mutate(period = period %>% str_extract("\\d{1,10}") %>% as.numeric()) %>%
  left_join(period_map, by = c("period")) %>%
  ggplot(aes(x=period_label, y=kappa, group=country, color=class)) + 
  geom_line() + 
  labs(x="Period", y=expression(kappa[t]), color="Cluster") + 
  theme(text = element_text(size = 10))

ggsave("plots/ILC-k-means-kappa.pdf", width = 5, height = 3) 

mx_ilc =  mx_tidy %>%
  as_tibble() %>%
  left_join(
    data.frame(
      country = countries_analyzed, 
      class = factor(ILC_kmeans_class)
    ),
    by = "country"
  ) %>% 
  mutate(class = paste0("Cluster ", class)) 

mx_ilc %>%
  group_by(year, class, age) %>%
  summarise(mx = mean(log(mx))) %>%
  filter(age %in% c(0, 15, 45, 80)) %>%
  mutate(age = ifelse(age == 0, paste0("Age group [0, 1)"), 
                      paste0("Age group [", age, ", ", age + 5, ")"))) %>%
  mutate(class = class %>% str_remove("Cluster ")) %>%
  ggplot(aes(x=year, y=(mx), color=class)) + 
  geom_line(alpha = 0.8) + 
  facet_wrap(. ~ age, scales = "free") + 
  #scale_color_viridis() + 
  labs(x = "Period", y=latex2exp::TeX("Average $\\log({}_{5}m_{x})$"),
       color = "Cluster") + 
  theme(text = element_text(size = 17))

ggsave("plots/time-ILC-k-means.pdf", width = 10, height = 5)

#### PCA fuzzy method ####
qx_tidy = lt %>%
  dplyr::select(country, year, age, qx) %>%
  dplyr::filter(age < 110)

qx = qx_tidy %>%
  spread(age, qx)

logitqx = qx_tidy %>%
  as_tibble() %>%
  mutate(year_age = paste0("v_", year, "_", age)) %>%
  dplyr::select(country, year_age, qx) %>%
  mutate(qx = log(qx/(1-qx))) %>%
  spread(year_age, qx) 

qx_matrix = logitqx %>% dplyr::select(-country) %>% as.matrix() %>% scale()
dim(qx_matrix)

set.seed(1)
qx_eigen_dec = qx_matrix %>% cor() %>% eigen()
lambda_cumsum = cumsum(qx_eigen_dec$values)/sum((qx_eigen_dec$values))
n_dim = length(lambda_cumsum[lambda_cumsum < 0.9]) + 1

qxPCS = qx_matrix %*% qx_eigen_dec$vectors[, 1:n_dim]
dim(qxPCS)

PCA_fuzzy_fit = lapply(2:10, function(g){
  set.seed(1)
  cmeans(qxPCS, centers = g, m = 2)
})

PCA_fuzzy_class_matrix = lapply(PCA_fuzzy_fit, function(fit){
  fit$cluster
}) %>%
  do.call(cbind, .)

# metrics 
metrics = lapply(2:10, function(k){
  
  fcm_result = PCA_fuzzy_fit[[k-1]]
  
  sil = cluster::silhouette(fcm_result$cluster, dist(qxPCS), FUN = mean)
  avg_sil = mean(sil[, 3])
  
  partition_coefficient = sum(fcm_result$membership^2) / nrow(qxPCS)
  partition_entropy = -sum(fcm_result$membership * log(fcm_result$membership)) / nrow(qxPCS)
  
  min_intercluster_dist = min(dist(fcm_result$centers))^2
  
  xie_beni = sum(apply(fcm_result$membership^2 * rowSums((qxPCS - fcm_result$centers[fcm_result$cluster, ])^2), 1, sum)) /
    (nrow(qxPCS) * min_intercluster_dist)
  
  global_mean = colMeans(qxPCS)
  
  
  c(
    "sil" = avg_sil, 
    "PC" = partition_coefficient,
    "xie_beni" = xie_beni 
  )
  
}) %>%
  do.call(rbind, .) %>%
  as.data.frame() %>%
  mutate(K = 2:10)

metrics %>%
  gather(metric, val, -K) %>%
  mutate(metric = ifelse(metric == "xie_beni", "Xie-Beni", metric),
         metric = ifelse(metric == "sil", "Silhouette", metric)) %>%
  ggplot(aes(x=K, y=val)) + 
  geom_point() + 
  geom_line() + 
  facet_wrap(. ~ metric, scales = "free") + 
  scale_x_continuous(breaks = 2:10) + 
  labs(x="Number of clusters", y="Metric")

ggsave("plots/PCA-fuzzy-metrics.pdf", width = 8, height = 3)

PCA_fuzzy_class = PCA_fuzzy_class_matrix[, 1]

mb_class = PCA_fuzzy_fit[[1]]$membership
rownames(mb_class) = logitqx$country
lev = mb_class %>% apply(MARGIN = 1, FUN = max) %>% round(2)

lapply(1:2, function(g){
  lev[PCA_fuzzy_class == g]
  
}) 

PCA_fuzzy_class

data.frame(
  level = lev, 
  class = PCA_fuzzy_class,
  country = names(lev)
) %>%
  mutate(level = ifelse(class == 1, level, 1-level)) %>%
  arrange(desc(level)) %>%
  mutate(country = factor(country, levels = rev(country))) %>%
  data.frame() %>%
  ggplot(aes(y=country, x=level, fill=factor(class))) + 
  geom_bar(stat="identity") +
  labs(x="") + 
  geom_vline(xintercept = 0.5, linetype = "dashed") + 
  labs(x="Membership level for Cluster 1 (West)", fill="Cluster", y="Country")

ggsave("plots/membership-level.pdf", width = 7, height = 5)


qx_fuzzy =  qx_tidy %>%
  as_tibble() %>%
  left_join(
    data.frame(
      country = countries_analyzed, 
      class = factor(PCA_fuzzy_class)
    ),
    by = "country"
  ) %>% 
  mutate(class = paste0("Cluster ", class)) 

qx_fuzzy %>%
  group_by(year, class, age) %>%
  summarise(qx = mean(log(qx/(1-qx)))) %>%
  filter(age %in% c(0, 15, 45, 80)) %>%
  mutate(age = ifelse(age == 0, paste0("Age group [0, 1)"), paste0("Age group [", age, ", ", age + 5, ")"))) %>%
  mutate(class = class %>% str_remove("Cluster ")) %>%
  ggplot(aes(x=year, y=qx, color=class)) + 
  geom_line(alpha = 0.8) + 
  facet_wrap(. ~ age, scales = "free") + 
  #scale_color_viridis() + 
  labs(x = "Period", y=latex2exp::TeX("Average logit(${}_{5}q_{x}$)"),
       color = "Cluster") + 
  theme(text = element_text(size = 17))

ggsave("plots/time-qx-fuzzy.pdf", width = 10, height = 5)

#### Func-k-means ####
ex_tidy = lt %>% dplyr::select(country, year, age, ex)
ex = ex_tidy %>% spread(age, ex)

ex_0 = ex %>% 
  dplyr::select(country, year , `0`) %>%
  spread(year, `0`)

ex_0_matrix = ex_0 %>%
  select(-country) %>%
  as.matrix()

years = lt$year %>% unique()
basis = create.bspline.basis(
  rangeval = period_range, 
  nbasis = 25,  
  norder = 3
)

md_fd_obj = smooth.basis(
  argvals = years, 
  y = t(ex_0_matrix), 
  fdParobj = basis
)

plot(md_fd_obj)

set.seed(1)
func_kmeans_fit = fit_compare(
  data = t(md_fd_obj$fd$coefs), 
  diss = dist(t(md_fd_obj$fd$coefs)), 
  G = 2:10, 
  n_start = 100, 
  n_iters = 1000, 
  method = "kmeans", 
  seed = 1
)

func_kmeans_fit$metrics_plot
ggsave("plots/func-k-means-metrics.pdf", width = 8, height = 3)

func_kmeans_class_matrix = func_kmeans_fit$class_matrix
func_kmeans_class = func_kmeans_class_matrix[, 1]

ex_0 %>%
  mutate(class = factor(func_kmeans_class)) %>%
  gather(time, e0, -country, -class) %>%
  as_tibble() %>%
  mutate(time = time %>% str_extract("\\d{1,10}") %>% as.numeric()) %>%
  ggplot(aes(x=time, y=e0, group=country, color=class)) + 
  geom_line() + 
  labs(x="Period", y=expression(e[0]), color="Cluster")  + 
  theme(text = element_text(size = 13))

ggsave("plots/func-methods-single.pdf", width = 6, height = 3.5) 

#### Rand index ####
class_matrix_list = list(
  "Hellinger-Ward" = h_ward_class_matrix,
  "ILC-k-means" = ILC_kmeans_class_matrix,
  "PCA-fuzzy" = PCA_fuzzy_class_matrix, 
  "func-k-means" = func_kmeans_class_matrix
)


methods = names(class_matrix_list)
rand_df = data.frame(
  method_1 = NULL, method_2 = NULL, G = NULL, rand = NULL
)

for(method_1 in methods) {
  for(method_2 in methods) {
    for(g in 2:10) {
      rand = mclust::adjustedRandIndex(
        class_matrix_list[[method_1]][, g-1], 
        class_matrix_list[[method_2]][, g-1]
      )
      
      rand_df = bind_rows(
        rand_df, 
        data.frame(
          method_1 = method_1, method_2 = method_2, G = g, rand = rand
        )
      )
    }
  }
}

rand_df %>%
  filter(method_1 != method_2) %>%
  filter(G <= 5) %>%
  ggplot(aes(x=G, y=rand, color=method_2)) + 
  geom_point() + 
  geom_line() + 
  labs(x="Number of clusters", color="Method", y="Adjusted rand index") + 
  facet_grid(. ~ method_1) + 
  theme(legend.position = "top")

ggsave("plots/rand-index.pdf", width = 7, height = 4)

##### graph ####
selected_class = data.frame(
  country = countries_analyzed, 
  h_ward_class_matrix[, 2],
  ILC_kmeans_class_matrix[, 1], 
  PCA_fuzzy_class_matrix[, 1],
  func_kmeans_class_matrix[, 1]
)

selected_class %>% 
  xtable::xtable() %>%
  print(include.rownames = F)

class_2 = data.frame(
  country = countries_analyzed, 
  h_ward_class_matrix[, 1],
  ILC_kmeans_class_matrix[, 1], 
  PCA_fuzzy_class_matrix[, 1],
  func_kmeans_class_matrix[, 1]
)

class_3 = data.frame(
  country = countries_analyzed, 
  h_ward_class_matrix[, 2],
  ILC_kmeans_class_matrix[, 2], 
  PCA_fuzzy_class_matrix[, 2],
  func_kmeans_class_matrix[, 2]
)

class_4 = data.frame(
  country = countries_analyzed, 
  h_ward_class_matrix[, 3],
  ILC_kmeans_class_matrix[, 3], 
  PCA_fuzzy_class_matrix[, 3],
  func_kmeans_class_matrix[, 3]
)


class_2 %>% plot_graph(seed = 1)
class_3 %>% plot_graph(seed = 1)
class_4 %>% plot_graph(seed = 1)

selected_class %>% plot_graph(seed = 3)
ggsave("plots/review-country-graph.pdf", width = 5, height = 5)






