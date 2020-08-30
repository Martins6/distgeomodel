# Title: Exploring datasets with large distance betweeen them
# Author: Adriel Martins
# Date: 02/03/2020
################****************************

#' 'Leave-one-out' Cross-Validation method to measure model performance.
#'
#' @param df A dataframe that has 'value' and two coordinates columns 'lat', 'long' for haversine or 'x' and 'y' for euclidean.
#' @param distance Type of distance to be choosen.
#' @param varg A variogram object from the gstat library.
#' @return A tibble containing data from the Cross-Validation.
#' @export
cv_krige_haversine <- function(df, varg, distance = 'haversine'){
  # *********************************************************************
  # 'dt' is a Data Frame object;
  # based on the distance we will have different requirements for the dt:
  # Euclidean distance: x, y, value coordinates
  # Haversine distance: lat, long, value coordinates

  # 'varg' is the variogram function from the gstat library adjusted
  # accoording to the distance
  # *********************************************************************

  aux <- tibble::tibble(cv.pred.krige = rep(NA_real_, nrow(df)),
                cv.var.pred.krige = rep(NA_real_, nrow(df)),
                cv.resid = rep(NA_real_, nrow(df)),
                z.score = rep(NA_real_, nrow(df))
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

  for(i in 1:nrow(df)){
    # Train dataset
    train <- df[-i,]
    # Test point coordinates
    test <- df[i,]
    # Kriging the train dataset on the test point
    cv.fit.krige <- ordinary_kriging(train, varg, distance = dista,
                                  df.coord = test[c(coord.x, coord.y)])
    # Acessing and storing CV values in our auxiliary matrix.
    cv.pred.krige <- cv.fit.krige$pred.krige
    cv.var.pred.krige <- cv.fit.krige$var.pred.krige
    cv.sd.pred.krige <- cv.var.pred.krige %>% sqrt
    cv.resid <- test$value - cv.pred.krige
    z.score <- cv.resid/cv.sd.pred.krige

    aux[i,] <- tibble::tibble_row(cv.pred.krige,
                 cv.var.pred.krige,
                 cv.resid,
                 z.score)
  }

  if(distance == 'haversine'){
    # Putting coordinates on each of our values
    aux <- aux %>%
      dplyr::mutate(lat = df$lat,
             long = df$long)

  }
  if(distance == 'haversine'){
    # Putting coordinates on each of our values
    aux <- aux %>%
      dplyr::mutate(x = df$x,
             y = df$y)
  }

  return(aux)
}
