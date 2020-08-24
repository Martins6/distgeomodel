# Title: Exploring datasets with large distance betweeen them
# Author: Adriel Martins
# Date: 02/03/2020
################****************************
# THIS CODE IS DEPENDENT ON THE 'kriging.R'.

cv_krige <- function(dt, varg, distance = 'haversine'){
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
    cv.sd.pred.krige <- cv.var.pred.krige %>% sqrt
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
