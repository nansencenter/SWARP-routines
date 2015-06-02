
# On input, MIZgridX must be NA for out-of-bounds, 0 elsewhere
# On output, MIZgridX will be reset to
#  1 inside  the MIZ
#  2 outside the MIZ, but inside MIZ footprint, on the ice infested side
# -2 outside the MIZ, but inside MIZ footprint, on the open ocean   side
 


setMIZgridX <- function(limits,varSizes,anaVar,footprint,MIZgridX) {

# Create the output list:
setMIZgridX <- list(numeric="grid")

# (Re-)Initialize MIZgridX to 0 where not out-of-bounds (NA):
MIZgridX <- 0*MIZgridX

nLimits  <- length(limits)
anaLow   <- limits[1]
anaHigh  <- limits[nLimits]


# Set  MIZgridX  to 1 in MIZ
for (y in 1:varSizes[2]) {
  ii <- which(anaVar[,y] < anaHigh & anaVar[,y] > anaLow)
  if (length(ii) > 0) MIZgridX[ii,y] <- 1
  rm(ii)
}

# Set  MIZgridX  to 2 in MIZ footprint region
for (y in 1:varSizes[2])
  for (x in 1:varSizes[1])
    if (!is.na(MIZgridX[x,y])) if (MIZgridX[x,y] == 1) {

      xa  <- max(1,          x-MIZnx)
      xe  <- min(varSizes[1],x+MIZnx)
      ya  <- max(1,          y-MIZny)
      ye  <- min(varSizes[2],y+MIZny)

      for (y0 in ya:ye) {
        yr <- y0 - y + MIZny + 1
        for (x0 in xa:xe) {
          xr <- x0 - x + MIZnx + 1
          if (!is.na(MIZgridX[x0,y0]) & !is.na(anaVar[x0,y0]))
            if (footprint[xr,yr] == 1 & MIZgridX[x0,y0] == 0)
              if (anaVar[x0,y0] >= anaHigh) {
                MIZgridX[x0,y0] <-  2
              } else {
                MIZgridX[x0,y0] <- -2
              }
        }
      }

    }


setMIZgridX$grid <- MIZgridX


return(setMIZgridX)

}
