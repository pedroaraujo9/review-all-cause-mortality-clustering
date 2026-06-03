library(tidyverse)
source("rscripts/utils.R")

format_data = function(data_path, 
                       period_range, 
                       age_range) {
  
  removed_small_countries = c(
    "Hong Kong", "Iceland", "Chile", "Croatia", 
    "South Korea", "Luxembourg", "East Germany",
    "West Germany", "England and Wales Total",
    "England and Wales Civilian", "Scotland",
    "Northern Ireland", "New Zealand Maori", 
    "New Zealand Non-Maori", "France Total"
  )
  
  lf = read_table(data_path, skip = 2)
  
  lf = lf %>% select(PopName, Year, Age, mx, qx, ex, dx)
  colnames(lf) = colnames(lf) %>% str_to_lower()
  
  lff = lf %>%
    mutate(age = age %>% str_remove("\\-\\d{1,10}|\\+") %>% as.integer()) %>%
    filter(age >= age_range[1], age <= age_range[2]) %>%
    mutate(country = hmd_country_names(popname)) %>%
    filter(!(country %in% removed_small_countries)) %>%
    filter(year >= period_range[1], year <= period_range[2])
  
  countries_full = lff %>%
    select(country, year) %>%
    distinct() %>%
    group_by(country) %>%
    summarise(n = n()) %>%
    filter(n == max(n)) %>%
    .$country
  
  lff = lff %>% 
    filter(country %in% countries_full) %>%
    select(country, year, age, mx, qx, ex, dx) %>%
    arrange(country, year, age)
  
  return(lff)
  
}


for(min_year in c(1960, 1990)) {
  
  for(max_year in c(2010, 2019)) {
    
    for(max_age in c(90, 110)) {
      
      for(sex in c("male", "female", "both")) {
        sex_min = substr(sex, 1, 1)
        data_path = paste0("data/lt_", sex, "/", sex_min, "ltper_5x1/", sex_min, "ltper_5x1.txt")
        rds_path = paste0(
          "data/sex=", sex, "-minyear=", min_year, "-maxyear=", max_year, 
          "-minage=0-", "maxage=", max_age, ".rds"
        )
        
        format_data(
          data_path = data_path, 
          period_range = c(min_year, max_year), 
          age_range = c(0, max_age)
        ) %>% saveRDS(rds_path)
      }
      
    }
  }
}



ct = format_data(
  data_path = paste0("data/lt_", "both", "/", "b", "ltper_5x1/", "b", "ltper_5x1.txt"), 
  period_range = c(1960, 2013), 
  age_range = c(0, max_age)
) %>%
  .$country

c("Russia", "Ukraine", "Belarus") %in% ct



