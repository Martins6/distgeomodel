# Title: Criteria for automatic selection of autocorrelation function based ond the haversine variogram
# Author: Adriel Martins
# Date: 02/03/2020

#' Fitting multiple variogram models using haversine distance (default) or euclidean.
#'
#' @param df A dataframe that has 'value' and two coordinates columns 'lat', 'long' for haversine or 'x' and 'y' for euclidean.
#' @param distance Type of distance to be choosen.
#' @param max_dist Maximum distance to be considered in the fitting of variogram models.
#' @param model_classes Classes of variogram models to be considered.
#' @param weight_types Weights considered in the fitting of the variogram models.
#' @param best_model_criteria The criteria to be taken for the choosing of the best model.
#' @return A list containing all the results of the combinations of the model and a gstat object of the best model.
#' @export
auto_krige_model_CV_haversine <-
  function(df,
           distance = 'haversine',
           best_model_criteria = 'MSE',
           max_dist,
           model_classes = c("matern", "exponential", "gaussian", "spherical"),
           weights_types = c("npairs", "cressie", "equal")) {

  # If distance = 'haversine',
  # dt must have 'lat' and 'long' columns, with only other column, named 'value'.
  # If distance = 'euclidean",
  # dt must have 'x' and 'y' columns, with only other column, named 'value'.

  # Transforming to geodata
  df_geo <- df %>% geoR::as.geodata()

  # Modelling Variogram with the SSE minimization
  if(distance == 'haversine'){

    if (missing(max_dist)) {
      # Variogram with Geodesic Distance
      v <- variog_non_euclidean(df_geo)
    } else{
      v <- variog_non_euclidean(df_geo, max.dist = max_dist)
    }

  }
  if(distance == 'euclidean'){
    aux_bin <- df_geo$coords %>%
      dist(method = 'euclidean') %>%
      vectorizing_fun()

    if(missing(max_dist)){max_dist <- max(aux_bin)}

    aux_bin1 <- aux_bin[aux_bin < max_dist] %>%
      quantile(probs = seq(0,1, l = 12))
    # Variogram with Euclidean Distance
    v <- geoR::variog(geodata = df_geo, uvec = aux_bin1)
  }

  # Our loop to find the best model for the variogram or autocorrelation function based on the SSE
  res <- tibble::tibble(Vario.Model = NA, Vario.MSSErr = NA, Vario.Weights = NA,
                Krige.CV.MSE = NA, Krige.CV.MAE = NA, Kappa_Matern = NA)
  # Modelling Variogram
  for(j in weights_types){
    for(i in model_classes){
      if(i == 'matern'){
        for(kap in c(1, 1.5)){
          # fit the variogram
          ols_fit <- geoR::variofit(v, cov.model = i, weights = j,
                              fix.kappa = T, kappa = kap, messages = F)
          # Changing from GeoR to GSTAT form
          # Adapting the name of the model so that it can fit in the gstat form
          model_gstat <- as.character(i) %>%
            stringr::str_sub(1L, 3L) %>%
            stringr::str_to_title()

          if(ols_fit$cov.pars[2] == 0){
            warning(paste('Numeric Error in fitting of the following model:', i, kap, j))
            next
          }
          var_fit <- gstat::vgm(psill = ols_fit$cov.pars[1],
                         model = model_gstat,
                         range = ols_fit$cov.pars[2],
                         nugget = ols_fit$nugget,
                         cutoff = ols_fit$max.dist)

          # Calculating the diference
          fit_cov <- gstat::variogramLine(var_fit, max(v$u), dist_vector = v$u)$gamma
          diff_cov <- v$v - fit_cov
          MSSErr <- mean(diff_cov^2)

          if(distance == 'haversine'){
            cv <- cv_krige_haversine(df, var_fit)
          }
          if(distance == 'euclidean'){
            cv <- cv_krige_haversine(df, var_fit, distance = 'euclidean')
          }

          # MAE and MSE
          aux <- cv %>%
            dplyr::summarise(MAE = mean(abs(z.score)),
                      MSE = mean(z.score^2))
          # Add the parameters to our auxiliary tibble
          res <- res %>%
            tibble::add_row(tibble::tibble_row(Vario.Model = i, Vario.MSSErr = MSSErr), Vario.Weights = j,
                    Krige.CV.MSE = aux$MSE, Krige.CV.MAE = aux$MAE, Kappa_Matern = kap)

        }
        next
      }
      # fit the variogram
      ols_fit <- geoR::variofit(v, cov.model = i, weights = j, messages = F)
      # Changing from GeoR to GSTAT form
      # Adapting the name of the model so that it can fit in the gstat form
      model_gstat <- as.character(i) %>%
        stringr::str_sub(1L, 3L) %>%
        stringr::str_to_title()

      if(ols_fit$cov.pars[2] == 0){
        warning(paste('Numeric Error in fitting of the following model:', i, j))
        next
      }
      var_fit <- gstat::vgm(psill = ols_fit$cov.pars[1],
                     model = model_gstat,
                     range = ols_fit$cov.pars[2],
                     nugget = ols_fit$nugget,
                     cutoff = ols_fit$max.dist)

      # Calculating the diference
      fit_cov <- gstat::variogramLine(var_fit, max(v$u), dist_vector = v$u)$gamma
      diff_cov <- v$v - fit_cov
      MSSErr <- mean(diff_cov^2)

      if(distance == 'haversine'){
        cv <- cv_krige_haversine(df, var_fit)
      }
      if(distance == 'euclidean'){
        cv <- cv_krige_haversine(df, var_fit, distance = 'euclidean')
      }

      # MAE and MSE
      aux <- cv %>%
        dplyr::summarise(MAE = mean(abs(z.score)),
                  MSE = mean(z.score^2))
      # Add the parameters to our auxiliary tibble
      res <- res %>%
        tibble::add_row(tibble::tibble_row(Vario.Model = i, Vario.MSSErr = MSSErr), Vario.Weights = j,
                Krige.CV.MSE = aux$MSE, Krige.CV.MAE = aux$MAE, Kappa_Matern = NA)
    }
  }

  # Finding the model that minimizes the best_model_criteria
  column <- stringr::str_c(c('Krige.CV.', best_model_criteria), collapse = '')
  best_model <- res %>%
    dplyr::filter(!is.na(Vario.Model)) %>%
    dplyr::filter(.data[[column]] == min(.data[[column]]))

  # Fitting the optimal variogram
  ols_fit <- geoR::variofit(v,
                      cov.model = as.character(best_model[1,1]),
                      weights = as.character(best_model[1,3]),
                      messages = F)

  # Changing from GeoR to GSTAT form
  # Adapting the name of the model so that it can fit in the gstat form
  model_gstat <- as.character(best_model[1,1]) %>%
    stringr::str_sub(1L, 3L) %>%
    stringr::str_to_title()

  var_fit <- gstat::vgm(psill = ols_fit$cov.pars[1],
                 model = model_gstat,
                 range = ols_fit$cov.pars[2],
                 nugget = ols_fit$nugget,
                 cutoff = ols_fit$max.dist)

  res <- res %>%
    # Filtering the first row
    dplyr::filter(!is.na(Vario.Model)) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(Vario.Model = if_else(!is.na(Kappa_Matern),
                                 stringr::str_c(Vario.Model, ' (Phi = ', as.character(Kappa_Matern), ')'),
                                 Vario.Model)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-Kappa_Matern)

  return(list(Results = res, var_fit_best_model = var_fit))
}
