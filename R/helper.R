#' Vectorizing a Matrix.
#'
#' @param a A matrix
#' @return A vectorized matrix, a vector.
vectorizing_fun <- function(a){

  # 'a' is a matrix

  a[!lower.tri(a)] <- NA # Filling the upper part of the triangle with NA's
  a <- as.vector(a) # Making the matrix a vector
  a <- a[!is.na(a)] # Taking out the NA's
  return(a) # returning only the lower triangle part of the matrix
}
