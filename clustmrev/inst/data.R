library(tidyverse)

lt = readRDS("../data/life_tables_5x1.rds")
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
  filter(n == n_years) %>%+
  .$country

n_country = length(countries_analyzed)

hmd_data = lt %>% filter(country %in% countries_analyzed)
usethis::use_data(hmd_data, overwrite = TRUE)

