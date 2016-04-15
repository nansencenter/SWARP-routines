
RMSdistance <- function(SRCpath, array1,array2,lon,lat,size, catValues,values, category) {

# Get a local function
Rfunction <- paste(SRCpath, "geoDistance.R", sep="/")
source(Rfunction)


# Check dimensions:
size1 <- length(array1)
size2 <- length(array2)


# Set target index in 'values':
catIndx <- which(catValues == category)


if (size1 != size2 | length(catIndx) != 1)
  rmsValue <- NA
else {

  nx <- size[1]
  ny <- size[2]

  target <- values[catIndx]

  # Find sea ice edge nodes in 'array2' (model product):

  indx2  <- which(array2 == target)
  n2     <- 0

  if (length(indx2) == 0) {

    print("One of the ice classes is empty in the model.")

  } else {

    x      <- indx2%%nx
    y      <- (indx2-1)%/%nx + 1

    match2 <- 0
    n2     <- 0
    for (i in 1:length(indx2)) {
      edge <- FALSE
      if (x[i] > 1)
        if (!is.na(array2[x[i]-1,y[i]  ]))
          if (array2[x[i]-1,y[i]  ] < target) edge <- TRUE
      if (x[i] < nx)
        if (!is.na(array2[x[i]+1,y[i]  ]))
          if (array2[x[i]+1,y[i]  ] < target) edge <- TRUE
      if (y[i] > 1)
        if (!is.na(array2[x[i]  ,y[i]-1]))
          if (array2[x[i],  y[i]-1] < target) edge <- TRUE
      if (y[i] < ny)
        if (!is.na(array2[x[i]  ,y[i]+1]))
          if (array2[x[i],  y[i]+1] < target) edge <- TRUE
      if (edge) {
        n2         <- n2 + 1
        match2[n2] <- i
      }
    }
  }
  # ...sea ice edge nodes are now
  #    indx2[match2[1:n2]]


  # Find sea ice edge nodes in 'array1' (observations):

  indx1  <- which(array1 == target)
  n1     <- 0

  if (length(indx1) == 0) {

    print("One of the ice classes is empty in the observations.")

  } else {

    x      <- indx1%%nx
    y      <- (indx1-1)%/%nx + 1

    match1 <- 0
    for (i in 1:length(indx1)) {
      edge <- FALSE
      if (x[i] > 1)
        if (!is.na(array1[x[i]-1,y[i]  ]))
          if (array1[x[i]-1,y[i]  ] < target) edge <- TRUE
      if (x[i] < nx)
        if (!is.na(array1[x[i]+1,y[i]  ]))
          if (array1[x[i]+1,y[i]  ] < target) edge <- TRUE
      if (y[i] > 1)
        if (!is.na(array1[x[i]  ,y[i]-1]))
          if (array1[x[i],  y[i]-1] < target) edge <- TRUE
      if (y[i] < ny)
        if (!is.na(array1[x[i]  ,y[i]+1]))
          if (array1[x[i],  y[i]+1] < target) edge <- TRUE
      if (edge) {
        n1         <- n1 + 1
        match1[n1] <- i
      }
    }
  }
  # ...sea ice edge nodes are now
  #    indx1[match1[1:n1]]


  if (n1*n2 > 0) {

    # Loop model edge nodes,
    # find shortest distance to observed edge node:

    sum <- 0.
    for (i2 in 1:n2) {
      lon2 <- lon[indx2[match2[i2]]]
      lat2 <- lat[indx2[match2[i2]]]
      minDistance <- min(geoDistance(lon2,lat2,
                                     lon[indx1[match1]],
                                     lat[indx1[match1]]))
      sum  <- sum + minDistance*minDistance
    }

    rmsValue <- sqrt(sum/n2)

  } else {

    rmsValue <- NA

  }
}

return(rmsValue)

}
