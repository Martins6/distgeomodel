# distkrige

R Package In-Development for Spatial Regression (Kriging) with Non-Euclidean Distances. To this day, you can only use haversine distances with lat-long coordinates. 
The package can be used only by itself, however it is best used as a complement for *geoR* or *gstat* packages.

# Installation

As the package is still in development, the only way to install is through devtools with Github. Just run the following code:

```{r}
devtools::install_github('Martins6/distkrige')
```

# Usage

There are two main functions. They serve the purpose of saying which is the best model to use for kriging (or spatial regression) in a list of many possibilities.
Those functions returns a bunch of nice objects to be used. The difference between them is that one focus what the complete model and the other focus on just the variogram.
There are, however, many auxiliary functions that can be used. They will be better improved for future use.

```
library(dplyr)
library(maps)
library(distkrige)

# From maps
data("ozone", package = 'maps')
# Adjusting data
oz <- ozone %>%
  as_tibble() %>%
  rename(lat = y,
         long = x,
         value = median)

auto_fit_variogram_haversine(oz)

auto_krige_model_CV_haversine(oz)
```

# Contributions

Contributions are very much welcome! Please do open issues and pull requests.
Check out the [project](https://github.com/Martins6/distkrige/projects/1?add_cards_query=is%3Aopen) area for more specific objectives for the package.

# References

Two references were key: Model-based Geostatistics (Diggle and Ribeiro Jr., 2007) and Applied Spatial Data Analysis with R (Bivand and Pebesma, 2008).
