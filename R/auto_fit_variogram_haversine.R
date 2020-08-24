# Title: Criteria for automatic selection of autocorrelation function based ond the haversine variogram
# Author: Adriel Martins
# Date: 02/03/2020

#' Fitting multiple variogram models using haversine distance (default) or euclidean.
#'
#' @param df A dataframe that has 'value' and two coordinates columns 'lat', 'long' for haversine or 'x' and 'y' for euclidean.
#' @param distance Type of distance to be choosen.
#' @param max_dist Maximum distance to be considered in the fitting of variogram models.
#' @return A list containing first the plot, caracteristic and gstat object of the best model; also all the results of the other models.
#' @export
auto_fit_variogram_haversine <- function(df, distance = 'haversine', max_dist){

  # If distance = 'haversine',
  # dt must have 'lat' and 'long' columns, with only other column, named 'value'.
  # If distance = 'euclidean",
  # dt must have 'x' and 'y' columns, with only other column, named 'value'.

  # Transforming to geodata
  df_geo <- df %>%
    geoR::as.geodata()

  # Modelling Variogram with the SSE minimization
  # All the model classes
  model_classes <- c("matern", "exponential", "gaussian", "spherical")
  # All the Weight Types
  weights_types <- c("npairs", "cressie", "equal")

  if(distance == 'haversine'){

    if (missing(max_dist)) {
      # Variogram with Geodesic Distance
      v <- variog_geodesic(df_geo)
    } else{
      v <- variog_geodesic(df_geo, max.dist = max_dist)
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
  res <- tibble::tibble(Model = NA, MSSErr = NA, Weights = NA, Kappa_Matern = NA)
  # Modelling Variogram
  for(j in weights_types){
    for(i in model_classes){
      if(i == 'matern'){
        for(kap in c(1, 1.5)){
          # fit the variogram
          ols_fit <- geoR::variofit(v, cov.model = i, weights = j,
                              fix.kappa = T, kappa = kap, messages = F)
          # Transforming to gstat form
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
          # Add the parameters to our auxiliary tibble
          res <- res %>%
            dplyr::add_row(tibble::tibble_row(Model = i, MSSErr = MSSErr), Weights = j, Kappa_Matern = kap)
        }
        next
      }

      # fit the variogram
      ols_fit <- geoR::variofit(v, cov.model = i, weights = j, messages = F)
      # Transforming to gstat form
      model_gstat <- as.character(i) %>%
        stringr::str_sub(1L, 3L) %>%
        stringr::str_to_title()

      if(ols_fit$cov.pars[2] == 0){next}

      var_fit <- gstat::vgm(psill = ols_fit$cov.pars[1],
                     model = model_gstat,
                     range = ols_fit$cov.pars[2],
                     nugget = ols_fit$nugget,
                     cutoff = ols_fit$max.dist)

      # Calculating the diference
      fit_cov <- gstat::variogramLine(var_fit, max(v$u), dist_vector = v$u)$gamma
      diff_cov <- v$v - fit_cov
      MSSErr <- mean(diff_cov^2)
      # Add the parameters to our auxiliary tibble
      res <- res %>%
        dplyr::add_row(tibble::tibble_row(Model = i, MSSErr = MSSErr), Weights = j, Kappa_Matern = NA)
    }

    print(paste(i, j))
  }

  # Finding the model that minimizes the SSE
  best_model <- res %>%
    dplyr::filter(!is.na(Model)) %>%
    dplyr::filter(MSSErr == min(MSSErr))

  print(best_model)

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

  # Plotting the results
  a <- cbind(v$v, v$u, v$n) %>%
    dplyr::as_tibble(.name_repair = 'unique') %>%
    `colnames<-`(c('gamma', 'dist', 'n'))

  b <- gstat::variogramLine(var_fit, maxdist = max(v$u)) %>%
    dplyr::as_tibble()

  if(distance == 'haversine'){
    title_plot <- 'Modelo Ótimo de Autocorrelação com a distância de haversine'
  }else{
    title_plot <- 'Modelo Ótimo de Autocorrelação com a distância euclidiana'
  }

  plot_variog_fit <- ggplot2::ggplot() +
    # Variogram points
    ggplot2::geom_point(ggplot2::aes(x = dist, y = gamma), a) +
    # Labels for each of the variogram points indicating how much points there are on each distance class
    ggplot2::geom_text(ggplot2::aes(x = dist, y = gamma, label = n),
              data = a, nudge_x = 10000) +
    # Autocorrelation function
    ggplot2::geom_path(ggplot2::aes(x = dist, y = gamma), b) +
    # Best model annotation
    ggplot2::geom_text(ggplot2::aes(x = a$dist[length(a$dist)], y = a$gamma[1]),
              label = paste('MSQR:', round(best_model[1,2], 3))) +
    ggplot2::labs(x = 'Distância',
         y = 'Semivariância',
         title = title.plot,
         subtitle = paste('Família:', stringr::str_to_title(best_model[1,1]),
                          ' ',
                          'Peso:', stringr::str_to_title(best_model[1,3]))) +
    ggplot2::theme_bw()

  res <- res %>%
    # Filtering the first row
    dplyr::filter(!is.na(Model)) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(Model = if_else(!is.na(Kappa_Matern),
                           str_c(Model, ' (Phi = ', as.character(Kappa_Matern), ')'),
                           Model)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-Kappa_Matern)

  return(list(plot = plot_variog_fit, variogram = var_fit,
              best_model = best_model, Results = res))

}

