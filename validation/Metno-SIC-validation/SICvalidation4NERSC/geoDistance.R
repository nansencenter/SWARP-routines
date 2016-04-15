
# Simple function to compute distance, ignores all effects of curvature
# Assumes that Earth's radius is 6371 km
# Value returned is in km
# NOTE! lon1,lat1 must be scalars, lon2,lat2 may be vectors/arrays

geoDistance <- function(lon1,lat1,lon2,lat2) {


# Check dimensions:
size11 <- length(lon1)
size12 <- length(lat1)
size21 <- length(lon2)
size22 <- length(lat2)
if (size11 != size12 | size11 != 1 | size21 != size22)
  rmsValue <- NA
else {

  earthRadius <- 6371.
  oneDeg      <- pi*earthRadius/180.

  dLon <- (lon2 - lon1)*cos(pi*(lat1 + lat2)/360.)
  dLat <-  lat2 - lat1

  distance <- sqrt(dLon*dLon + dLat*dLat)*oneDeg
}

return(distance)

}
