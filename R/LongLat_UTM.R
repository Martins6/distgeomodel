# Title: Long-Lat to UTM transformation
# Author: Adriel Martins
# Date: 10/09/19
# ************************************************************************************************* #
# libraries
library(sp)
#library(rgdal)
library(dplyr)
library(stringr)

# Function
LongLatToUTM <- function(long = -20.902735,lat = -6.050278, sp_dt = 1, zone = 25){
  
  # If we want to change a already SpatialDataFrame!
  if(typeof(sp_dt) != 'double'){
    
    proj4string(sp_dt) <- CRS("+proj=longlat +datum=WGS84")
    
    sp_dt <- spTransform(sp_dt, CRS(paste("+proj=utm +zone=",zone,"+datum=WGS84", sep = '')))
    return(sp_dt)
    
  }
  
  xy <- tibble(Eastings = long, Northings = lat)
  coordinates(xy) <- c("Eastings", "Northings")
  # Most used Datum is WGS84
  proj4string(xy) <- CRS("+proj=longlat +datum=WGS84")
  
  # Example longitude number, for we assume that the points are close enough by.s
  long <- long[1]
  
  # Longitude can say where is the Zone, and the sign of the Latitude,
  # says if it is North or South.
  zone <- toString(31 + long/6)
  ## we need only the first two numbers of the above numeric.
  zone <- as.integer(str_sub(zone, start = 1, end = 2))
  
  res <- spTransform(xy, CRS(paste("+proj=utm +zone=",zone," +ellps=WGS84", sep = '')))
  return(as.data.frame(res))
  
}

UTM_To_LongLat <- function(Eastings, Northings, zone){
  
  xy <- tibble(Long = Eastings, Lat = Northings)
  coordinates(xy) <- c("Long", "Lat")
  proj4string(xy) <- CRS(paste("+proj=utm +zone=",zone," +ellps=WGS84", sep = ''))
  
  res <- spTransform(xy, CRS("+proj=longlat +datum=WGS84"))
  return(as.data.frame(res))
  
}

# # Example
# home <- c(-34.902735, -8.050278)
# LongLatToUTM(home[1],home[2])
# UTM_To_LongLat(aux[1,1], aux[1,2], 25)


################################## Rijksdriehoek to LongLat ###########################

rd_to_wgs84 <- function(x, y) {
  if (!is.numeric(x)) stop("x needs to be numeric.")
  if (nargs() == 1  && length(x) == 2) {
    y <- x[2]
    x <- x[1]
  } 
  if (!is.numeric(y)) stop("y needs to be numeric.")
  if (length(x) != length(y)) stop("x and y need to have same length.")
  
  x0 <- 155000.00
  y0 <- 463000.00
  phi0 <- 52.15517440
  lambda0 <- 5.38720621
  
  k <- list( 
    list(p=0, q=1, k=3235.65389), 
    list(p=2, q=0, k= -32.58297), 
    list(p=0, q=2, k=  -0.24750), 
    list(p=2, q=1, k=  -0.84978), 
    list(p=0, q=3, k=  -0.06550), 
    list(p=2, q=2, k=  -0.01709), 
    list(p=1, q=0, k=  -0.00738), 
    list(p=4, q=0, k=   0.00530), 
    list(p=2, q=3, k=  -0.00039), 
    list(p=4, q=1, k=   0.00033), 
    list(p=1, q=1, k=  -0.00012))
  l <- list(
    list(p=1, q=0, l=5260.52916), 
    list(p=1, q=1, l= 105.94684), 
    list(p=1, q=2, l=   2.45656), 
    list(p=3, q=0, l=  -0.81885), 
    list(p=1, q=3, l=   0.05594), 
    list(p=3, q=1, l=  -0.05607), 
    list(p=0, q=1, l=   0.01199), 
    list(p=3, q=2, l=  -0.00256), 
    list(p=1, q=4, l=   0.00128))
  dx <- (x - x0)/1E5
  dy <- (y - y0)/1E5
  phi <- rep(phi0, length(x))
  lambda <- rep(lambda0, length(x))
  
  for (i in seq_along(k)) {
    p <- k[[i]][["p"]]
    q <- k[[i]][["q"]]
    ki <- k[[i]][["k"]]
    phi <- phi + (ki * (dx^p) * (dy^q))/3600
  }
  for (i in seq_along(l)) {
    p <- l[[i]][["p"]]
    q <- l[[i]][["q"]]
    li <- l[[i]][["l"]]
    lambda <- lambda + (li * (dx^p) * (dy^q))/3600
  }
  if (nargs() == 1) {
    c(lambda[1], phi[1])
  } else list(lambda=lambda, phi=phi)
}




