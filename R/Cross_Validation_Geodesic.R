# Title: Exploring datasets with large distance betweeen them
# Author: Adriel Martins
# Date: 02/03/2020
################****************************
# THIS CODE IS DEPENDENT ON THE 'kriging.R'.

cv_geodesic <- function(dt, varg, distance = 'haversine'){
  # *********************************************************************
  # 'dt' is a Data Frame object;
  # based on the distance we will have different requirements for the dt:
  # Euclidean distance: x, y, value coordinates
  # Haversine distance: lat, long, value coordinates
  
  # 'varg' is the variogram function from the gstat library adjusted 
  # accoording to the distance
  # *********************************************************************
  
  aux <- tibble(cv.pred.krige = rep(NA_real_, nrow(dt)),
                cv.var.pred.krige = rep(NA_real_, nrow(dt)),
                cv.resid = rep(NA_real_, nrow(dt)),
                z.score = rep(NA_real_, nrow(dt))
                )
  
  if(distance == 'haversine'){
    dista <- 'haversine'
    coord.x <- 'lat'
    coord.y <- 'long'
  }
  
  if(distance == 'euclidean'){
    dista <- 'euclidean'
    coord.x <- 'x'
    coord.y <- 'y'
  }
  
  for(i in 1:nrow(dt)){
    # Train dataset
    train <- dt[-i,]
    # Test point coordinates
    test <- dt[i,]
    # Kriging the train dataset on the test point
    cv.fit.krige <- ordinary_kriging(train, varg, distance = dista,
                                  df.coord = test[c(coord.x, coord.y)])
    # Acessing and storing CV values in our auxiliary matrix.
    cv.pred.krige <- cv.fit.krige$pred.krige
    cv.var.pred.krige <- cv.fit.krige$var.pred.krige
    cv.sd.pred.krige <- cv.var.pred.krige %>% sqrt()
    cv.resid <- test$value - cv.pred.krige
    z.score <- cv.resid/cv.sd.pred.krige
    
    aux[i,] <- tibble_row(cv.pred.krige,
                 cv.var.pred.krige,
                 cv.resid,
                 z.score)
  }
  
  if(distance == 'haversine'){
    # Putting coordinates on each of our values
    aux <- aux %>% 
      mutate(lat = dt$lat,
             long = dt$long)
    
  }
  if(distance == 'haversine'){
    # Putting coordinates on each of our values
    aux <- aux %>% 
      mutate(x = dt$x,
             y = dt$y)
  }
  
  return(aux)
}

################# TEST
# # From maps
# data("ozone", package = 'maps')
# # Reading data
# oz <- ozone %>%
#   as_tibble() %>%
#   rename(lat = y,
#          long = x,
#          value = median)
# # Transforming to geodata 
# oz.geo <- oz %>%
#   as.geodata()
# # Variogram with Geodesic Distance
# v <- variog.geodesic(oz.geo)
# # Modelling Variogram
# ols.fit <- variofit(v, cov.model = 'spherical', weights="equal")
# 
# # CHANGING TO GSTAT VARIOGRAM FORM
# var.fit <- vgm(psill = ols.fit$cov.pars[1],
#                model = "Sph",
#                range = ols.fit$cov.pars[2],
#                nugget = ols.fit$nugget,
#                cutoff = ols.fit$max.dist)
# # Kriging
# aux <- ordinary_kriging(oz, var.fit,
#                  distance = 'haversine',
#                  df.coord = oz[c('lat', 'long')])
# # Cross-Validation
# cv <- cv_geodesic(oz, var.fit)
# # Distribution of our z-score: we hope it will be a N(0,1).
# cv$z.score %>% hist()
# # Finding outliers based on the quantile on the z.score
# fit.normal <- cv$z.score %>% MASS::fitdistr('normal')
# threshold <- qnorm(p = c(0.05), mean = 0, sd = sqrt(fit.normal$estimate[2]), lower.tail = T)
# 
# cv.outliers <- cv %>% filter(z.score <= threshold | z.score >= -threshold)
# # Plotting the outliers
# world <- ne_countries(scale = "medium", returnclass = "sf")
# ggplot() +
#   geom_sf(data = world) +
#   coord_sf(xlim = oz.gstat@bbox[2,], ylim = oz.gstat@bbox[1,], expand = TRUE) +
#   geom_point(data = oz, aes(y = lat, x = long, colour = value)) +
#   geom_point(data = cv.outliers, aes(y = lat, x = long), colour = 'red')
