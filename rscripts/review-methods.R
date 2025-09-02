library(tidyverse)
library(NbClust)
library(cluster)
library(vegan)
library(dtw)
library(patchwork)
library(e1071)
library(vital)
library(fda)
library(funHDDC) # needs manual tar.gz installation 
library(ggraph)
library(tidygraph)
library(igraph)
library(mclust)
select = dplyr::select

nb_clust_kmeans = function(data, diss = NULL) {
  
  indices = c(
    "kl", "ch", "hartigan", "ccc", "scott", 
    "marriot", "trcovw", "tracew", "friedman", 
    "rubin", "cindex", "db", "silhouette", "duda", 
    "pseudot2", "beale", "ratkowsky", "ball", 
    "ptbiserial", "gap", "frey", "mcclain", "gamma", 
    "gplus", "tau", 
    "dunn", "hubert", "sdindex", "dindex", "sdbw"
  )
  
  if(is.null(diss)){
    
    cluster_fit = lapply(indices, function(index){
      
      pdf(NULL) 
      
      fit = safely(NbClust, otherwise = NA)(
        data = data, 
        method = "kmeans", 
        index = index, 
        distance = "euclidean"
      )
      
      dev.off()  
      
      par(mfrow=c(1,1))
      return(fit)
    })
    
  }else{
    cluster_fit = lapply(indices, function(index){
      
      pdf(NULL) 
      
      fit = safely(NbClust, otherwise = NA)(
        data = data, 
        diss = diss, 
        method = "kmeans", 
        index = index, 
        distance = NULL
      )
      
      dev.off()  
      
      return(fit)
    })
  }
  
  number_cluster = lapply(cluster_fit, function(fit){
    if(is.null(fit$error)){
      if(is.null(fit$result$Best.partition)) {
        return(NA)
      }else{
        return(length(unique(fit$result$Best.partition)))
      }
    }else{
      return(NA)
    }
  }) %>%
    do.call(c, .)
  
  
  table_number_cluster = number_cluster %>% table() 
  opt_n = which.max(table_number_cluster) %>% names() %>% as.numeric()
  
  class = cluster_fit[[which(number_cluster == opt_n)[1]]]$result$Best.partition
  return(class)
}


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

period_map


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

# mx data 
mx_tidy = lt %>%
  dplyr::select(country, year, age, mx)

mx = mx_tidy %>%
  spread(age, mx)

mx_matrix = mx %>% dplyr::select(-country, -year) %>% as.matrix()

#### PCA k-medoids ####
# reference PCA
X_base = mx_matrix[mx$country == countries_analyzed[1], ] %>% log() %>% scale()
eigen_decomp_reference = X_base %>% cor() %>% eigen()
v_base = cbind(eigen_decomp_reference$vectors[, 1])
PC_base = X_base %*% v_base

# PCA's
PCA1_list = lapply(countries_analyzed, function(country){
  X = mx_matrix[mx$country == country, ] %>% log() %>% scale()
  eigen_decomp = X %>% cor() %>% eigen()
  v = cbind(eigen_decomp$vectors[, 1])
  PC = (X %*% v)
  tr_PC = vegan::procrustes(PC_base, PC, scale = F)
  new_PC = PC %*% tr_PC$rotation
  new_v = v %*% tr_PC$rotation
  
  list(PC1 = new_PC, v = new_v)
}) 


PC1 = PCA1_list %>% 
  lapply(function(PCA){
    PCA$PC1 %>% t()
  }) %>%
  do.call(rbind, .)

DTW_dist = matrix(0, nrow = n_country, ncol = n_country)

for(i in 1:n_country) {
  for(j in 1:n_country) {
    dt = dtw::dtw(PC1[, i], PC1[, j])
    DTW_dist[i, j] = dt$normalizedDistance
  }
}

PCA_kmedoids_class = nb_clust_kmeans(PC1, DTW_dist)

p1 = PCA1_list %>% 
  lapply(function(PCA){
    PCA$v %>% t()
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame() %>%
  mutate(country = countries_analyzed, class = factor(PCA_kmedoids_class)) %>%
  gather(age_group, eigenvec, -country, -class) %>%
  mutate(age_group = age_group %>% str_extract("\\d{1,10}") %>% as.numeric()) %>%
  left_join(age_map, by = "age_group") %>%
  ggplot(aes(x=age_group_label, y=eigenvec, group=country, color=class)) + 
  geom_line() + 
  geom_point() + 
  guides(color="none") + 
  labs(x="Age group", y="Eigenvector value", color="Cluster") + 
  theme(text = element_text(size = 15))

p2 = PC1 %>%
  as.data.frame() %>%
  mutate(country = countries_analyzed, class = factor(PCA_kmedoids_class)) %>%
  gather(period, pc1, -country, -class) %>%
  mutate(period = period %>% str_extract("\\d{1,10}") %>% as.numeric()) %>%
  left_join(period_map, by = "period") %>%
  ggplot(aes(x=period_label, y=pc1, group=country, color=class)) + 
  geom_line() + 
  labs(x="Period", y="PC1", color="Cluster") + 
  theme(text = element_text(size = 15))

p1
ggsave("plots/pca-k-medoids-eigenvec.pdf", width = 5, height = 3)

p2
ggsave("plots/pca-k-medoids-scores.pdf", width = 5, height = 3)


p1 + p2 
ggsave("plots/pca-k-medoids.pdf", width = 10, height = 3)

mx_pca =  mx_tidy %>%
  as_tibble() %>%
  left_join(
    data.frame(
      country = countries_analyzed, 
      class = factor(PCA_kmedoids_class)
    ),
    by = "country"
  ) %>% 
  mutate(class = paste0("Cluster ", class)) 

mx_pca %>%
  group_by(year, class, age) %>%
  summarise(mx = mean(log(mx))) %>%
  mutate(id = paste0(year, "-", class)) %>%
  ggplot(aes(x=age, y=(mx), color=year, group=id)) + 
  geom_line(alpha = 0.8) + 
  facet_grid(. ~ class) + 
  scale_color_viridis() + 
  labs(x = "Age group", y=latex2exp::TeX("Average $\\log$ ${}_{5}m_{x}$")) + 
  theme(text = element_text(size = 17))

ggsave("plots/mx-pca-k-medoids.pdf", width = 10, height = 3)


PCA1_list %>% 
  lapply(function(PCA){
    PCA$v %>% t()
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame() %>%
  mutate(country = countries_analyzed, class = factor(PCA_kmedoids_class)) %>%
  gather(age_group, eigenvec, -country, -class) %>%
  mutate(age_group = age_group %>% str_extract("\\d{1,10}") %>% as.numeric()) %>%
  left_join(age_map, by = "age_group") %>%
  group_by(age_group_label) %>%
  summarise(cv = sd(eigenvec)/abs(mean(eigenvec))) %>%
  arrange(desc(cv)) %>%
  View()


mx_pca %>%
  group_by(year, class, age) %>%
  summarise(mx = mean(log(mx))) %>%
  filter(age %in% c(0, 15, 45, 80)) %>%
  mutate(age = ifelse(age == 0, paste0("Age group [0, 1)"), paste0("Age group [", age, ", ", age + 5, ")"))) %>%
  ggplot(aes(x=year, y=(mx), color=class)) + 
  geom_line(alpha = 0.8) + 
  facet_wrap(. ~ age, scales = "free") + 
  #scale_color_viridis() + 
  labs(x = "Period", y=latex2exp::TeX("Average $\\log$ ${}_{5}m_{x}$"),
       color = "Class") + 
  theme(text = element_text(size = 17))

ggsave("plots/time-mx-pca-k-medoids.pdf", width = 10, height = 5)


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

qx_eigen_dec = qx_matrix %>% cor() %>% eigen()
lambda_cumsum = cumsum(qx_eigen_dec$values)/sum((qx_eigen_dec$values))
n = length(lambda_cumsum[lambda_cumsum < 0.9]) + 1

qxPCS = qx_matrix %*% qx_eigen_dec$vectors[, 1:n]
dim(qxPCS)


# metrics 
metrics = lapply(2:10, function(k){
  fcm_result = cmeans(qxPCS, centers = k, m = 2)
  
  sil = cluster::silhouette(fcm_result$cluster, dist(qxPCS), FUN = mean)
  avg_sil = mean(sil[, 3])
  
  partition_coefficient = sum(fcm_result$membership^2) / nrow(qxPCS)
  partition_entropy = -sum(fcm_result$membership * log(fcm_result$membership)) / nrow(qxPCS)
  
  min_intercluster_dist = min(dist(fcm_result$centers))^2
  
  xie_beni = sum(apply(fcm_result$membership^2 * rowSums((qxPCS - fcm_result$centers[fcm_result$cluster, ])^2), 1, sum)) /
    (nrow(qxPCS) * min_intercluster_dist)
  
  global_mean = colMeans(qxPCS)
  
  compactness = sum(fcm_result$membership^2 * rowSums((qxPCS - fcm_result$centers[fcm_result$cluster, ])^2))
  separation = sum(fcm_result$membership^2 * rowSums((fcm_result$centers - global_mean)^2))
  fukuyama_sugeno = compactness - separation
  c(
    "sil" = avg_sil, 
    "PC" = partition_coefficient,
    #"PE" = -partition_entropy, 
    "xie_beni" = xie_beni 
    #"FS" = fukuyama_sugeno
  )
  
}) %>%
  do.call(rbind, .) %>%
  as.data.frame() %>%
  mutate(K = 2:10)

metrics %>%
  gather(metric, val, -K) %>%
  mutate(metric = ifelse(metric == "xie_beni", "Xie-Beni", metric),
         metric = ifelse(metric == "sil", "silhouette", metric)) %>%
  ggplot(aes(x=K, y=val)) + 
  geom_point() + 
  geom_line() + 
  facet_wrap(. ~ metric, scales = "free") + 
  scale_x_continuous(breaks = 2:10) + 
  labs(x="Number of clusters", y="Metric value") + 
  theme(text = element_text(size = 15))

ggsave("plots/fuzzy-metrics.pdf", height = 3, width = 10)

set.seed(1)
PCA_fuzzy_class = cmeans(qxPCS, centers = 4, m = 2)$cluster

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
  mutate(id = paste0(year, "-", class)) %>%
  ggplot(aes(x=age, y=qx, color=year, group=id)) + 
  geom_line(alpha = 0.8) + 
  facet_wrap(. ~ class, ncol = 2) + 
  scale_color_viridis() + 
  labs(x = "Age group", y=latex2exp::TeX("Average logit ${}_{5}q_{x}$")) + 
  theme(text = element_text(size = 15))

ggsave("plots/qx-fuzzy.pdf", width = 8, height = 5)

qx_fuzzy %>%
  group_by(year, class, age) %>%
  summarise(qx = mean(log(qx/(1-qx)))) %>%
  filter(age %in% c(0, 15, 45, 80)) %>%
  mutate(age = ifelse(age == 0, paste0("Age group [0, 1)"), paste0("Age group [", age, ", ", age + 5, ")"))) %>%
  ggplot(aes(x=year, y=qx, color=class)) + 
  geom_line(alpha = 0.8) + 
  facet_wrap(. ~ age, scales = "free") + 
  #scale_color_viridis() + 
  labs(x = "Period", y=latex2exp::TeX("Average logit ${}_{5}q_{x}$"),
       color = "Class") + 
  theme(text = element_text(size = 17))

ggsave("plots/time-qx-fuzzy.pdf", width = 10, height = 5)

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

##### kmeans #####
ILC_k_means_class = nb_clust_kmeans(beta_x)

p1 = kappa_t %>%
  as.data.frame() %>%
  mutate(country = countries_analyzed, class = factor(ILC_k_means_class)) %>%
  gather(period, kappa, -country, -class) %>%
  as_tibble() %>%
  mutate(period = period %>% str_extract("\\d{1,10}") %>% as.numeric()) %>%
  left_join(period_map, by = c("period")) %>%
  ggplot(aes(x=period_label, y=kappa, group=country, color=class)) + 
  geom_line() + 
  labs(x="Period", y=expression(kappa[t]), color="Cluster") + 
  theme(text = element_text(size = 15))


p2 = beta_x %>%
  as.data.frame() %>%
  mutate(country = countries_analyzed, class = factor(ILC_k_means_class)) %>%
  gather(age_group, bx, -country, -class) %>%
  as_tibble() %>%
  mutate(age_group = age_group %>% str_extract("\\d{1,10}") %>% as.numeric()) %>%
  left_join(age_map, by = "age_group") %>%
  ggplot(aes(x=age_group_label, y=bx, group=country, color=class)) + 
  geom_line() + 
  guides(color="none") + 
  labs(x="Age group", y=expression(beta[x]), color="Cluster") + 
  theme(text = element_text(size = 15))

p1 + p2

p2
ggsave("plots/ILC-k-means-beta.pdf", width = 5, height = 3) 

p1
ggsave("plots/ILC-k-means-kappa.pdf", width = 5, height = 3) 


p2 + p1
ggsave("plots/ILC-k-means.pdf", width = 10, height = 3) 

mx_ilc =  mx_tidy %>%
  as_tibble() %>%
  left_join(
    data.frame(
      country = countries_analyzed, 
      class = factor(ILC_k_means_class)
    ),
    by = "country"
  ) %>% 
  mutate(class = paste0("Cluster ", class)) 


mx_ilc %>%
  group_by(year, class, age) %>%
  summarise(mx = mean(log(mx))) %>%
  mutate(id = paste0(year, "-", class)) %>%
  ggplot(aes(x=age, y=(mx), color=year, group=id)) + 
  geom_line(alpha = 0.8) + 
  facet_grid(. ~ class) + 
  scale_color_viridis() + 
  labs(x = "Age group", y=latex2exp::TeX("Average $\\log$ ${}_{5}m_{x}$")) + 
  theme(text = element_text(size = 17))

ggsave("plots/mx-ILC-k-means.pdf", width = 8, height = 3)

beta_x %>%
  as.data.frame() %>%
  mutate(country = countries_analyzed, class = factor(ILC_k_means_class)) %>%
  gather(age_group, bx, -country, -class) %>%
  as_tibble() %>%
  mutate(age_group = age_group %>% str_extract("\\d{1,10}") %>% as.numeric()) %>%
  left_join(age_map, by = "age_group") %>%
  group_by(age_group_label) %>%
  summarise(cv = sd(bx)/abs(mean(bx))) %>%
  arrange(desc(cv)) %>%
  View()

mx_ilc %>%
  group_by(year, class, age) %>%
  summarise(mx = mean(log(mx))) %>%
  filter(age %in% c(0, 15, 45, 80)) %>%
  mutate(age = ifelse(age == 0, paste0("Age group [0, 1)"), paste0("Age group [", age, ", ", age + 5, ")"))) %>%
  ggplot(aes(x=year, y=(mx), color=class)) + 
  geom_line(alpha = 0.8) + 
  facet_wrap(. ~ age, scales = "free") + 
  #scale_color_viridis() + 
  labs(x = "Period", y=latex2exp::TeX("Average $\\log$ ${}_{5}m_{x}$"),
       color = "Class") + 
  theme(text = element_text(size = 17))

ggsave("plots/time-ILC-k-means.pdf", width = 10, height = 5)

##### functional model-based ####
basis = create.bspline.basis(rangeval = range(mx_tidy$age), 
                             nbasis = 23,  
                             norder = 3)

md_fd_obj = smooth.basis(argvals = mx_tidy$age %>% unique(), 
                         y = t(beta_x), 
                         fdParobj = basis)$fd

mdfda = funHDDC(md_fd_obj, K = c(2:10))

func_model_based_class = mdfda$class
ILC_k_means_class

mclust::adjustedRandIndex(
  ILC_k_means_class, 
  func_model_based_class
)

beta_x %>%
  as.data.frame() %>%
  mutate(country = countries_analyzed, class = factor(func_model_based_class)) %>%
  gather(age_group, bx, -country, -class) %>%
  as_tibble() %>%
  mutate(age_group = age_group %>% str_extract("\\d{1,10}") %>% as.numeric()) %>%
  left_join(age_map, by = "age_group") %>%
  ggplot(aes(x=age_group_label, y=bx, group=country, color=class)) + 
  geom_line() + 
  guides(color="none") + 
  labs(x="Age group", y=expression(beta[x]), color="Cluster") + 
  theme(text = element_text(size = 15))

#### Func-k-means ####
ex_tidy = lt %>% dplyr::select(country, year, age, ex)
ex = ex_tidy %>% spread(age, ex)

ex_0 = ex %>% 
  dplyr::select(country, year , `0`) %>%
  spread(year, `0`)

years = lt$year %>% unique()
basis = create.bspline.basis(rangeval = period_range, 
                             nbasis = diff(period_range) + 1,  
                             norder = 3)

fd_coefs = lapply(countries_analyzed, function(country_name){
  ex_vec = ex_0 %>%
    filter(country == country_name) %>%
    dplyr::select(-country) %>%
    unlist() 
  
  fd_obj = smooth.basis(argvals = years, 
                        y = ex_vec, 
                        fdParobj = basis)
  fd_obj$fd$coefs %>% as.numeric()
  
}) %>% do.call(rbind, .) %>%
  scale()

func_k_means_class = nb_clust_kmeans(fd_coefs)

p1 = ex_0 %>%
  mutate(class = factor(func_k_means_class)) %>%
  gather(time, e0, -country, -class) %>%
  as_tibble() %>%
  mutate(time = time %>% str_extract("\\d{1,10}") %>% as.numeric()) %>%
  ggplot(aes(x=time, y=e0, group=country, color=class)) + 
  geom_line() + 
  labs(x="Period", y=expression(e[0]), color="Clusters for Func-k-means")  + 
  theme(text = element_text(size = 13), legend.position = "top")

#### Func-model-based ####
ex_0_matrix = ex_0 %>% dplyr::select(-country) %>% as.matrix() 

md_fd_obj = smooth.basis(argvals = years, 
                         y = t(ex_0_matrix), 
                         fdParobj = basis)$fd

mdfda = funHDDC(md_fd_obj, K = c(2:10))

func_model_based_class = mdfda$class

p2 = ex_0 %>%
  mutate(class = factor(func_model_based_class)) %>%
  gather(time, e0, -country, -class) %>%
  as_tibble() %>%
  mutate(time = time %>% str_extract("\\d{1,10}") %>% as.numeric()) %>%
  mutate(class = factor(class, levels = c(2, 1), labels = c("1", "2"))) %>%
  ggplot(aes(x=time, y=e0, group=country, color=class)) + 
  geom_line() + 
  labs(x="Period", y=expression(e[0]), color="Clusters for Func-model-based")  + 
  theme(text = element_text(size = 13), legend.position = "top")

p1
ggsave("plots/func-k-means.pdf", width = 5, height = 3.5) 

p2
ggsave("plots/func-model-based.pdf", width = 5, height = 3.5) 


p1 + p2
ggsave("plots/func-methods.pdf", width = 9, height = 3.5) 


ex_0 %>%
  mutate(class = factor(func_model_based_class)) %>%
  gather(time, e0, -country, -class) %>%
  as_tibble() %>%
  filter(time == 2010) %>%
  filter(class == 1) %>%
  arrange(desc(e0))

#### class analysis ####
class_df = data.frame(
  country = countries_analyzed, 
  PCA_fuzzy_class, 
  ILC_k_means_class, 
  func_k_means_class, 
  func_model_based_class
)

class_df %>% 
  xtable::xtable() %>%
  print(include.rownames = F)

class_df %>%
  filter(func_k_means_class == 2 & func_model_based_class == 1)

class_df %>%
  filter(func_k_means_class == 2 & func_model_based_class == 2)

class_df %>%
  filter(func_k_means_class == 1 & func_model_based_class == 1)

class_df %>%
  filter(func_k_means_class == 1 & func_model_based_class == 2)


match_matrix = matrix(0, n_country, n_country)
rownames(match_matrix) = colnames(match_matrix) = class_df$country

for (i in 1:n_country) {
  for (j in 1:n_country) {
    match_matrix[i, j] = sum(class_df[i, -1] == class_df[j, -1])
  }
}

country_conn = match_matrix %>%
  as.data.frame() %>%
  mutate(to = rownames(.)) %>%
  gather(from, weight, -to) %>%
  filter(from != to) %>%
  filter(weight != 0)

country_conn$id = NA
for(i in 1:nrow(country_conn)) {
  country_conn$id[i] = sort(c(country_conn$from[i], country_conn$to[i])) %>% 
    str_flatten(collapse = "-")
}

country_conn = country_conn %>%
  distinct(id, .keep_all = T) %>%
  select(-id)

country_graph = graph_from_data_frame(country_conn, directed = FALSE)
country_tbl_graph = as_tbl_graph(country_graph)


set.seed(1)
ggraph(country_tbl_graph, layout = "fr") +  # Fruchterman-Reingold layout
  geom_edge_link(aes(alpha = weight), color = "grey60") +  
  geom_node_point() + 
  geom_node_text(aes(label = name), repel = TRUE, family = "serif") + 
  labs(alpha = "aa") + 
  theme_graph()

ggsave("plots/review-country-graph.pdf", width = 5, height = 5)

