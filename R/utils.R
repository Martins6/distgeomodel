#' Vectorizing a Matrix.
#'
#' @param a A matrix
#' @return A vectorized matrix, a vector.
vectorizing_fun <- function(x){

  # 'x' is a matrix

  x[!lower.tri(x)] <- NA # Filling the upper part of the triangle with NA's
  x <- as.vector(x) # Making the matrix a vector
  x <- x[!is.na(x)] # Taking out the NA's
  return(x) # returning only the lower triangle part of the matrix
}
