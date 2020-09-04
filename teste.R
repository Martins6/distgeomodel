library(dplyr)

# From maps
data("ozone", package = 'maps')
# Reading data
oz <- ozone %>%
  as_tibble() %>%
  rename(lat = y,
         long = x,
         value = median)

distkrige::auto_fit_variogram_haversine(oz)
