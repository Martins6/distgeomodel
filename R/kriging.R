# Title: Implementing Ordinary Kriging
# Author: Adriel Martins
# Date: 02/03/2020
# References: Model-based Geostatistics - Peter J. Diggle, Paulo Ribeiro Jr.


#' Fitting multiple variogram models using haversine distance (default) or euclidean.
#'
#' @param df A dataframe that has 'value' and two coordinates columns 'lat', 'long' for haversine or 'x' and 'y' for euclidean.
#' @param distance Type of distance to be choosen.
#' @param model_classes Classes of variogram models to be considered.
#' @param weight_types Weights considered in the fitting of the variogram models.
#' @param best_model_criteria The criteria to be taken for the choosing of the best model.
#' @return A list containing all the results of the combinations of the model and a gstat object of the best model.
#' @export
################# ************ Creating the function
ordinary_kriging <- function(df, varg, distance = 'euclidean', df.coord){
  # *********************************************************************
  # 'distance' is the type of distance to use in our calculations

  # 'varg' is the variogram function from the gstat library

  # 'df' is a Data Frame object;
  # based on the distance we will have different requirements for the dt:
  # Euclidean distance: x, y, value coordinates
  # Haversine distance: lat, long, value coordinates

  # 'df.coord' will be the same as dt, but without the value column.
  # *********************************************************************
  flag <- nrow(df.coord) == 1

  if(flag){
    df.coord <- rbind(df.coord, df.coord)
  }

  # Our parameters from the variogram function
  sigma2 <- varg[2,2]
  tau2 <- varg[1,2]

  #################### ******************** EUCLIDEAN DISTANCE ####################
  if(distance == 'euclidean'){

    # Transforming our data-set into a SpatialPointsDataFrame from the gstat package.
    df_sp <- df %>%
      sp::`coordinates<-`(c("x", "y"))

    # Euclidean distance matrix
    dist_matrix <- df_sp@coords %>%
      dist(method = 'euclidean', upper = T, diag = T) %>%
      as.matrix()
    ## Euclidean distance matrix with our dataframe coordinates
    aux <- df_sp@coords %>%
      rbind(., df.coord) %>%
      dist(method = 'euclidean', upper = T, diag = T) %>%
      as.matrix()
    ## We want just the distance with each new coordinate with old coordinates
    l_interest <- (nrow(df_sp@coords) + 1):(nrow(aux))
    c_interest <- 1:(nrow(df_sp@coords))
    pred_dist_matrix <- aux[l_interest, c_interest]
  }

  #################### ******************** HAVERSINE DISTANCE ####################
  if(distance == 'haversine'){

    # Transforming our data-set into a SpatialPointsDataFrame from the gstat package.
    df_sp <- df %>%
      sp::`coordinates<-`(c("lat", "long"))

    # Euclidean distance matrix
    dist_matrix <- df_sp@coords %>%
      geodist::geodist(measure = 'haversine')
    ## Euclidean distance matrix with our T point
    aux <- df_sp@coords %>%
      rbind(., df.coord) %>%
      geodist::geodist(measure = 'haversine')
    ## We want just the distance with each new coordinate with old coordinates
    l_interest <- (nrow(df_sp@coords) + 1):(nrow(aux))
    c_interest <- 1:(nrow(df_sp@coords))
    pred_dist_matrix <- aux[l_interest, c_interest]
  }

  # Calculating Covariance Matrix of the Observed Values for the Model.
  # Creating our auxiliary Matrix
  res <- rep(NA, nrow(dist_matrix))
  # Calculating the covariance values for each distance of the observed values
  for(i in 1:nrow(dist_matrix)){
    aux <- gstat::variogramLine(varg, maxdist = 3000,
                         covariance = T, dist_vector = dist_matrix[i,])[2]
    # anexing in our auxiliary matrix
    res <- rbind(res, t(aux))
  }
  # Our complete auxiliary matrix, taking out the NA
  res <- res[-1,]
  # Our complete 'pure correlation' "R" matrix for the observed values
  R <- cov2cor(res)
  # Our complete model variance V matrix
  V <- R + (tau2/sigma2)*diag(nrow(dist_matrix))
  # Our true model variance sigma2*V matrix
  s2V <- sigma2*V

  # Calculating Correlation Matrix of the Predicted Values with the Observed Values.
  # Creating our auxiliary Matrix with the covariances
  res <- rep(NA, nrow(pred_dist_matrix))
  # Calculating the covariance values for each distance of the observed values
  for(i in 1:nrow(pred_dist_matrix)){
    aux <- gstat::variogramLine(varg, maxdist = 3000,
                         covariance = T, dist_vector = pred_dist_matrix[i,])[2]
    # anexing in our auxiliary matrix
    res <- rbind(res, t(aux))
  }
  # Our complete auxiliary covariance matrix, taking out the NA
  res <- res[-1,]
  # Calculating the variance or gamma(0)
  gamma_0 <- gstat::variogramLine(varg, maxdist = 3000,
                           covariance = T, dist_vector = 0)[2] %>%
    dplyr::pull() %>%
    as.numeric()
  # Creating the correlation matrix
  r <- res/gamma_0

  # Calculating beforehand some matrices
  Vminus <- solve(V)
  one_vec_n <- rep(1, nrow(dist_matrix)) %>% as.matrix()
  one_vec_p <- rep(1, nrow(pred_dist_matrix)) %>% as.matrix()
  Y <- df$value %>% as.vector() %>% as.matrix()
  # Estimating mu by generalized least squares
  mu_h <- solve( t(one_vec_n) %*% Vminus %*% one_vec_n ) %*% t(one_vec_n) %*% Vminus %*% Y
  mu_h <- mu_h %>% as.numeric()
  # Estimating our mean point prediction
  T_h <- mu_h * one_vec_p + r %*% Vminus %*% (Y - mu_h * one_vec_n ) %>%
    t() %>% as.vector()
  # Estimating our variance point prediction
  var_T_h <- sigma2 * ( 1 - diag( r %*% Vminus %*% t(r) ) ) %>% as.vector()

  T_pred <- cbind(df.coord, T_h, var_T_h) %>%
    dplyr::as_tibble() %>%
    dplyr::rename(pred.krige = T_h,
           var.pred.krige = var_T_h)

  if(flag){
    T_pred <- T_pred[1,]
  }

  return(T_pred)
}

