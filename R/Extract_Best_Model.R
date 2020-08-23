# Title: Criteria for automatic selection of autocorrelation function based ond the haversine variogram
# Author: Adriel Martins
# Date: 02/03/2020
################****************************
# Libraries
# Source the code for the functions, with the dependent libraries.
# source('my_codes/Kriging.R')
# source('my_codes/Variogram_Geodesic_Fitting_Function.R')
# source('my_codes/Cross_Validation_Geodesic.R')
####################### EXTRACTING THE BEST MODEL THROUGH MSE OR MAE ####################
extract_best_model <- function(autokrige, criteria = 'MSE', distance = 'haversine'){

  if(distance == 'haversine'){
    # Variogram with Geodesic Distance
    v <- variog.geodesic(dt.geo)
  }
  if(distance == 'euclidean'){
    aux.bin <- dt.geo$coords %>%
      dist(method = 'euclidean') %>%
      vectorizing_fun()

    max.dist <- max(aux.bin)

    aux.bin1 <- aux.bin[aux.bin < max.dist] %>%
      quantile(probs = seq(0,1, l = 12))
    # Variogram with Euclidean Distance
    v <- variog(geodata = dt.geo, uvec = aux.bin1)
  }

  column <- str_c(c('Krige.CV.', criteria))
  autokrige[autokrige$column == min(autokrige$column)]

  # Fitting the optimal variogram
  ols.fit <- variofit(v,
                      cov.model = as.character(best.model[1,1]),
                      weights = as.character(best.model[1,3]))

  # Changing from GeoR to GSTAT form
  # Adapting the name of the model so that it can fit in the gstat form
  model.gstat <- as.character(best.model[1,1]) %>%
    stringr::str_sub(1L, 3L) %>%
    str_to_title()

  var.fit <- vgm(psill = ols.fit$cov.pars[1],
                 model = model.gstat,
                 range = ols.fit$cov.pars[2],
                 nugget = ols.fit$nugget,
                 cutoff = ols.fit$max.dist)

  return(var.fit)
}
