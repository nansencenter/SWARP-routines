# Import library for handling netCDF files
library(RNetCDF)


# Set radius for an extended MIZ grid in km, to be used in
# def. of class in confusion matrix
MIZgridR <- 200.


# Read input home directory, file names from the command line:
cargs      <- Sys.getenv(c('SRCpath', 'IOpath', 'WORKpath', 'bdate'))
SRCpath    <- cargs[1]
IOpath     <- cargs[2]
WORKpath   <- cargs[3]
bdate      <- cargs[4]
if (SRCpath=='' | IOpath=='' | WORKpath=='' | bdate=='') {
  print("Syntax e.g.:")
  print(" SRCpath='/home/arnem/myOcean/CalVal/R'; IOpath='/home/arnem/myOcean/CalVal/IO'; WORKpath='/disk1/tmp'; bdate='2010-02-10'")
  print(" export SRCpath IOpath WORKpath bdate")
  print(" R CMD BATCH SIvalidate.R")
  quit("no")
}


# Add time of day to bulletin date
bDate    <- paste(bdate,"12:00:00")

# Get some local functions
Rfunction <- paste(SRCpath, "mapVector.R",          sep="/")
source(Rfunction)
Rfunction <- paste(SRCpath, "ncVarSizes.R",         sep="/")
source(Rfunction)
Rfunction <- paste(SRCpath, "rootMeanSquare.R",     sep="/")
source(Rfunction)
Rfunction <- paste(SRCpath, "RMSdistance.R",        sep="/")
source(Rfunction)
Rfunction <- paste(SRCpath, "biasDistance.R",       sep="/")
source(Rfunction)
Rfunction <- paste(SRCpath, "categorize.R",         sep="/")
source(Rfunction)
Rfunction <- paste(SRCpath, "categorizeNA.R",       sep="/")
source(Rfunction)
Rfunction <- paste(SRCpath, "categoryOverlap.R",    sep="/")
source(Rfunction)
Rfunction <- paste(SRCpath, "MIZconfusionMatrix.R", sep="/")
source(Rfunction)
Rfunction <- paste(SRCpath, "setMIZgridX.R",        sep="/")
source(Rfunction)

# Set names of output ASCII files:
resultFile     <- paste(IOpath, "SIvalidationResult.dat",   sep="/")
resultFilePmod <- paste(IOpath, "SIvalidationResultPmod.dat", sep="/")
resultFilePana <- paste(IOpath, "SIvalidationResultPana.dat", sep="/")
confusionFile  <- paste(IOpath, "SIconfusionResultBE.dat",    sep="/")


# Set file names, pointers of input files:
anaFile <- paste(WORKpath,  "iceChartResults.nc", sep="/")
modFile <- paste(WORKpath,  "TOPAZresults.nc",    sep="/")
mapFile <- paste(IOpath,    "gridmap.nc",         sep="/")

# Open input files (netCDF):
anaID   <- open.nc(anaFile, write=FALSE)
modID   <- open.nc(modFile, write=FALSE)
mapID   <- open.nc(mapFile, write=FALSE)

# Set ice charts classes; NOTE! Hardcoded here!
iceLimits <- c(        NA,      0.1,     0.35,    0.625,   0.85,   NA)
iceValues <- c(            0.,      0.2,    0.5,     0.75,    0.95)
iceCatgrs <- length(iceValues)
catValues  <- c(   0,   1,    2,     3,    4)
NAvalue    <- c(-1)
cntrValues <- c(NAvalue-0.5, NAvalue+0.5, catValues+0.5)

XiceLimits   <- c(0.1, 0.35, 0.625, 0.85)
nXcategories <- length(XiceLimits) + 1

nNAfirst     <- 0


# Set 'edgeCategory' so that the sea ice edge is defined as the nodes which
# is of this category (cfr. def. of 'catValues' above) and has at least one
# lower category neighbor (neighbors: [x-1,y], [x,y-1], [x+1,y], [x,y+1])
edgeCategory <- 2


# Set netCDF variable names:

anaVarName  <- "ice_concentration"  # from anaID
modVarName  <- "fice"               # from modID
xVarName    <- "x"                  # from modID
yVarName    <- "y"                  # from modID
xMapVarName <- "xMap"               # from mapID
yMapVarName <- "yMap"               # from mapID


# Find dimensional sizes:
anaSizes <- ncVarSizes(anaID,anaVarName)
modSizes <- ncVarSizes(modID,modVarName)


# Model results are scaled & offset, read attributes:
modScale  <- att.get.nc(modID, modVarName, "scale_factor")
modOffset <- att.get.nc(modID, modVarName, "add_offset")


# Coordinates, in units of 100km:
xCoord <- var.get.nc(modID,xVarName)
yCoord <- var.get.nc(modID,yVarName)


# Set grid size in km, check that grid size is constant:

# ...x
dx  <- xCoord[2] - xCoord[1]
xOK <- TRUE
for (x in 3:modSizes[1])
  if (xOK & xCoord[x]-xCoord[x-1] != dx) {
    print(paste("Non-constant grid in x, first at",x))
    xOK <- FALSE
  }
# ...units of 100km to 1km:
dx <- 100.*dx

# ...y
dy  <- yCoord[2] - yCoord[1]
yOK <- TRUE
for (y in 3:modSizes[2])
  if (yOK & yCoord[y]-yCoord[y-1] != dy) {
    print(paste("Non-constant grid in y, first at",y))
    yOK <- FALSE
  }
# ...units of 100km to 1km:
dy <- 100.*dy

if (!xOK | !yOK) quit("no")


# Set footprint  MIZfootprint  of node with radius of MIZgridR
#  =1: inside  footprint,
#  =0: outside footprint
MIZnx        <- (MIZgridR/dx + 0.5)%/%1
MIZny        <- (MIZgridR/dy + 0.5)%/%1
MIZfootprint <- array( 1,dim=c(2*MIZnx+1,2*MIZny+1))
MIZxc        <- MIZnx + 1
MIZyc        <- MIZny + 1
for (y in 1:MIZny)
  for (x in 1:MIZnx) {
    dist <- sqrt(x*x*dx*dx + y*y*dy*dy)
    if (dist > MIZgridR) {
      MIZfootprint[MIZxc+x,MIZyc+y] <- 0
      MIZfootprint[MIZxc-x,MIZyc+y] <- 0
      MIZfootprint[MIZxc+x,MIZyc-y] <- 0
      MIZfootprint[MIZxc-x,MIZyc-y] <- 0
    }
  }


# Mapping:
xMapVar <- var.get.nc(mapID,xMapVarName)
yMapVar <- var.get.nc(mapID,yMapVarName)

# Nodes where validation should not be performed
#  (non-overlapping or masked/e.g. Baltic)
# is set to NA; set masking array accordingly:
badNode <- is.na(xMapVar) | is.na(yMapVar)


# Create arrays for counting & aggregated values:
ana2numCount <- array(dim=modSizes[1:2])
ana2numSum   <- array(dim=modSizes[1:2])


# Create array for storing the extended MIZ region
#  ...will hold = 1 if inside
#               = 0 if outside & wet
#               =NA if out-of-bounds or dry
MIZgridX <- array(0,dim=modSizes[1:2])


# Read dates

# ...w/ analysis
anaDateString <- att.get.nc(anaID,"time","units")
anaTime       <- var.get.nc(anaID, "time")
anaDates      <- utcal.nc(anaDateString, anaTime, type="s")

# ...w/ model results, best estimate
modDateString <- att.get.nc(modID,"time","units")
modTime       <- var.get.nc(modID, "time")
modDates      <- utcal.nc(modDateString, modTime, type="s")


# Read model grid positions:

modLon        <- var.get.nc(modID, "longitude")
modLat        <- var.get.nc(modID, "latitude" )


# Set array with all dates where a model result is available:
#  start with best estimate (type 'mod') if there are no leading
#  forecast (type 'for') only dates; supplement if the othre
#  source has trailing dates:
numDates  <- modDates
nNumDates <- length(numDates)


# Set a reference date/unit, and compute offsets:

refDateString <- "hours since 2000-01-01 00:00:00"

anaOffsets    <- utinvcal.nc(refDateString, anaDates)
modOffsets    <- utinvcal.nc(refDateString, modDates)
numOffsets    <- utinvcal.nc(refDateString, numDates)
bOffset       <- utinvcal.nc(refDateString, bDate)


# ...TOPAZ gives times at midnight, while results are
#    averaged over the following 24h. We must add an offset:
modOffsets    <- modOffsets + 12
numOffsets    <- numOffsets + 12


# Set a time window, in units from 'refDateString' above:
#  (fields within +-timeWindow will be treated as synoptic)
timeWindow <- 6

# Find appropriate dates, i.e. model dates which have an
# analysis sufficiently near in time (within 'timeWindow'):
anaEntry <- mapVector(numOffsets,anaOffsets,timeWindow)
modEntry <- mapVector(numOffsets,modOffsets,timeWindow)
bEntry   <- mapVector(numOffsets,bOffset,   timeWindow)


# Set bDateNo to the first (and only) element where it is non-zero
#  (this is the index of the bulletin date)
bDateNo <- 1
while (bEntry[bDateNo] == 0) bDateNo <- bDateNo+1


# Initialize arrays
#  *Frqcy        will count #grids with each ice class
#  rmsConValue*  will store model vs. ice chart concentration rms offsets
#  rmsMIZvalue*  similar to rmsConValue*, but for the extended MIZ region only
#  overlapFrqcy* will count #grids where ice classes overlap in model and ice charts)
#  biasEdgValue* will store ice edge bias (positive when extent is largest in model)
#  rmsEdgValue*  will store model vs. ice chart edge position rms offsets

anaFrqcy         <- array(dim=c(nNumDates,iceCatgrs))
modFrqcy         <- array(dim=c(nNumDates,iceCatgrs))
forFrqcy         <- array(NA,dim=c(nNumDates,iceCatgrs))

overlapFrqcyBE   <- array(dim=c(nNumDates,iceCatgrs))
overlapFrqcyPmod <- array(dim=c(nNumDates,iceCatgrs))
overlapFrqcyPana <- array(dim=c(nNumDates,iceCatgrs))

rmsConValueBE    <- array(dim=c(nNumDates))
rmsConValuePmod  <- array(dim=c(nNumDates))
rmsConValuePana  <- array(dim=c(nNumDates))

rmsMIZvalueBE    <- array(dim=c(nNumDates))
rmsMIZvaluePmod  <- array(dim=c(nNumDates))
rmsMIZvaluePana  <- array(dim=c(nNumDates))

nObsEdgValueBE   <- array(dim=c(nNumDates))
nObsEdgValuePmod <- array(dim=c(nNumDates))
nObsEdgValuePana <- array(dim=c(nNumDates))

nModEdgValueBE   <- array(dim=c(nNumDates))
nModEdgValuePmod <- array(dim=c(nNumDates))
nModEdgValuePana <- array(dim=c(nNumDates))

biasEdgValueBE   <- array(dim=c(nNumDates))
biasEdgValuePmod <- array(dim=c(nNumDates))
biasEdgValuePana <- array(dim=c(nNumDates))

rmsEdgValueBE    <- array(dim=c(nNumDates))
rmsEdgValuePmod  <- array(dim=c(nNumDates))
rmsEdgValuePana  <- array(dim=c(nNumDates))



# Loop model dates:
# ...NOTE: Assumes 'mod' ("best estimate") and 'for' ("forecast")
#          are provided for the same domain!
#          ...AND with identical land-sea masks!
numStart <- c(1,1,1)
numCount <- c(modSizes[1:2],1)
anaStart <- c(1,1,1)
anaCount <- c(anaSizes[1:2],1)

anaFirst  <- TRUE
MIZmasked <- FALSE


anaCatgry  <- list(numeric="values", numeric="distribution")
modCatgry  <- list(numeric="values", numeric="distribution")
plotCatgry <- list(numeric="values", numeric="distribution")


BEconfusion <- array( 0,dim=c(nXcategories+2,nXcategories,nNumDates))

for (numNo in 1:nNumDates) {

  anaNo <- anaEntry[numNo]
  modNo <- modEntry[numNo]

  if (anaNo != 0) {

    print(paste("Processing time step no. ",numNo))

    # (re-)initialize arrays:
    ana2numCount[,] <- 0
    ana2numSum  [,] <- 0


    # Read input:

    # ...analysis is concentration in %
    anaStart[3] <- anaNo
    anaVar <- 0.01*var.get.nc(anaID, anaVarName,
                              start=anaStart, count=anaCount)

    # ...model input is scaled & offset concentration (i.e. in [0,1])
    # "best estimate":
    if (modNo != 0 ) {
      numStart[3] <- modNo
      modVar <- modScale*var.get.nc(modID, modVarName,
                                    start=numStart, count=numCount) +
                modOffset
    }

    for (y in 1:anaSizes[2])
      for (x in 1:anaSizes[1])
        if (!badNode[x,y]) {
          ana2numCount[xMapVar[x,y],yMapVar[x,y]] <-
          ana2numCount[xMapVar[x,y],yMapVar[x,y]] + 1
          ana2numSum  [xMapVar[x,y],yMapVar[x,y]] <-
          ana2numSum  [xMapVar[x,y],yMapVar[x,y]] + anaVar[x,y]
        }


    ana2numVar <- ana2numSum/ana2numCount

    # Overlay masks
    # ...assume stationary masks!
    if (anaFirst | !MIZmasked) {
      maskInd  <- is.na(modVar) | is.na(ana2numVar)
      # ...mask dry, out-of-bounds grids in def. of extended MIZ
      is.na(MIZgridX) <- c(maskInd)
      MIZmasked       <- TRUE
      anaFirst        <- FALSE
    }
    if (modNo != 0) is.na(modVar)     <- c(maskInd)
                    is.na(ana2numVar) <- c(maskInd)

    anaCatgry  <- categorize(ana2numVar, iceLimits, iceValues)
    ana2numVar <- anaCatgry$values

    anaFrqcy[numNo,] <- anaCatgry$distribution

    # Set spec. of extended MIZ
    MIZgridX <- setMIZgridX(XiceLimits,modSizes,ana2numVar,
                            MIZfootprint,MIZgridX)$grid
    MIZindx  <- which(MIZgridX != 0)


    if (modNo != 0) {
      modCatgry <- categorize(modVar, iceLimits, iceValues)
      modVar    <- modCatgry$values
      biasRslt  <- biasDistance(SRCpath, ana2numVar,modVar,
                                modLon,modLat, modSizes[1:2],
                                catValues,iceValues, edgeCategory)
      modFrqcy[numNo,]       <- modCatgry$distribution
      rmsConValueBE[numNo]   <- rootMeanSquare(ana2numVar,modVar)
      rmsMIZvalueBE[numNo]   <- rootMeanSquare(ana2numVar[MIZindx],
                                               modVar[MIZindx])
      nObsEdgValueBE[numNo]  <- biasRslt$nObs
      nModEdgValueBE[numNo]  <- biasRslt$nMod
      biasEdgValueBE[numNo]  <- biasRslt$biasValue
      rmsEdgValueBE[numNo]   <- RMSdistance(SRCpath, ana2numVar,modVar,
                                            modLon,modLat,modSizes[1:2],
                                            catValues,iceValues, edgeCategory)
      BEconfusion[,,numNo]   <- MIZconfusionMatrix(ana2numVar, modVar,
                                                   MIZgridX, XiceLimits)
      olapCount              <- categoryOverlap(ana2numVar, modVar, iceLimits)
      overlapFrqcyBE[numNo,] <- olapCount
    }


    rm(MIZindx)


  } else {     # belongs to:  if (anaNo != 0)

    print("No appropriate analysis is available.")

    # Assuming stationary land/sea masks!

    if (anaFirst) {

      # No analysis handled yet, none available for this date;
      # grab the first analysis to set the land-sea mask below
      #  [maskInd <- is.na(modVar) | is.na(ana2numVar) ...]

      # Initialize arrays:
      ana2numCount[,] <- 0
      ana2numSum  [,] <- 0

      # ...analysis is concentration in %
      anaStart[3] <- 1
      anaVar      <- 0.01*var.get.nc(anaID, anaVarName,
                                     start=anaStart, count=anaCount)

      for (y in 1:anaSizes[2])
        for (x in 1:anaSizes[1])
          if (!badNode[x,y]) {
            ana2numCount[xMapVar[x,y],yMapVar[x,y]] <-
            ana2numCount[xMapVar[x,y],yMapVar[x,y]] + 1
            ana2numSum  [xMapVar[x,y],yMapVar[x,y]] <-
            ana2numSum  [xMapVar[x,y],yMapVar[x,y]] + anaVar[x,y]
          }

      ana2numVar <- ana2numSum/ana2numCount

    }

    # Read input:
    # ...model input is scaled&offset abs. concentration
    if (modNo != 0) {
      numStart[3] <- modNo
      modVar      <- modScale*var.get.nc(modID, modVarName,
                                         start=numStart, count=numCount) +
                     modOffset
    }

    # Overlay masks
    if (anaFirst) {
      maskInd    <- is.na(modVar) | is.na(ana2numVar)
      anaFirst   <- FALSE
    }
    is.na(modVar) <- c(maskInd)


    # When no analysis is available, store #grids with each ice class
    # as negative numbers along the "diagonal" of the confusion matrix:

    if (modNo != 0) {
      modCatgry <- categorize(modVar, iceLimits, iceValues)
      modVar    <- modCatgry$values
      modFrqcy[numNo,] <- modCatgry$distribution
      for (ctg in 1:iceCatgrs) {
        if (ctg == 1)        xCtg <- 1             else
        if (ctg < iceCatgrs) xCtg <- ctg + 1       else
                             xCtg <- iceCatgrs + 2
        BEconfusion[xCtg,ctg,numNo] <- -modCatgry$distribution[ctg]
      }
    }

  }     # belongs to:  if (anaNo != 0)

}


# Validation of best estimate:

outList  <- c(numOffsets,1:nNumDates,anaEntry,rmsMIZvalueBE,rmsConValueBE,
              nObsEdgValueBE,nModEdgValueBE,biasEdgValueBE,rmsEdgValueBE,
              anaFrqcy,modFrqcy,overlapFrqcyBE)
tmpArray <- array(outList,dim=c(nNumDates,9+3*iceCatgrs))
outArray <- aperm(tmpArray,c(2,1))
write(outArray, file=resultFile,ncolumns=9+3*iceCatgrs, sep=" ")

# Write confusion matrix/best estimate vs. observations to file;
# output on confusionFile:
# (observations:o; best estimate:b; ice class:#)
#  o=0 & b=0 | o=1 & b=0 | o=2 & b=0 | o=3 & b=0 | o=4 & b=0 | o=5 & b=0 | o=6 & b=0
#  o=0 & b=1 | o=1 & b=1 | o=2 & b=1 | o=3 & b=1 | o=4 & b=1 | o=5 & b=1 | o=6 & b=1
#  o=0 & b=2 | o=1 & b=2 | o=2 & b=2 | o=3 & b=2 | o=4 & b=2 | o=5 & b=2 | o=6 & b=2
#  o=0 & b=3 | o=1 & b=3 | o=2 & b=3 | o=3 & b=3 | o=4 & b=3 | o=5 & b=3 | o=6 & b=3
#  o=0 & b=4 | o=1 & b=4 | o=2 & b=4 | o=3 & b=4 | o=4 & b=4 | o=5 & b=4 | o=6 & b=4
# ...
BEsize   <- dim(BEconfusion)
tmpArray <- array(BEconfusion,dim=c(BEsize[1],BEsize[2]*BEsize[3]))
write(tmpArray, file=confusionFile,ncolumns=BEsize[1], sep=" ")


quit("no")
