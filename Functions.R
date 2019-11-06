calc.latlong.dist<- function(xy1,xy2) {
  # uses spherical law of cosines to calculate distance between two lat/long
  # coordinates in decimal degrees
  R <- 6371 # Earths radius
  xy1 <- (pi * xy1)/180 # radians
  xy2 <- (pi * xy2)/180 
  D <- acos(sin(xy1[,1])*sin(xy2[,1]) + cos(xy1[,1])*cos(xy2[,1])*cos(xy2[,2]-xy1[,2]))  
  return(R*D)
}

inv.logit <- function(x) { exp(x) / (1+exp(x))}

inv.cloglog <- function(x) {1-exp(-exp(x))}
