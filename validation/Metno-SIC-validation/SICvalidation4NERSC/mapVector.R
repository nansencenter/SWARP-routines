
# The function 'mapVector' maps the positions of similar values
# (identical values if 'window' is 0) in 

# Example 1:
# ...input:
# window     =  2
# fromVector =  8 12 16 20
# toVector   =  6  9 12 15
# ...output:
# vectorMap  =  1  3  4  0
# ...explanation:
#    algorithm returns first match in [value-window, value+window]

# Example 2 (same as 1, but with org. 'toVector[3]' value discarded):
# ...input:
# window     =  2
# fromVector =  8 12 16 20
# toVector   =  6  9 15
# ...output:
# vectorMap  =  1  0  3  0

# Example 3 (same as 1, but with 'window' = 1):
# ...input:
# window     =  1
# fromVector =  8 12 16 20
# toVector   =  6  9 12 15
# ...output:
# vectorMap  =  2  3  4  0

# Example 4 (same as 1, but with 'window' = 0):
# ...input:
# window     =  0
# fromVector =  8 12 16 20
# toVector   =  6  9 12 15
# ...output:
# vectorMap  =  0  3  0  0


mapVector <- function(fromVector, toVector, window) {


  fromSize <- length(fromVector)
  toSize   <- length(toVector)

  vectorMap <- rep(0, times=fromSize)


  # Check that all entries are ascending

  fromOk <- TRUE
  n      <- 1
  while (n < fromSize & fromOk) {
    n      <- n+1
    fromOk <- fromVector[n] > fromVector[n-1]
  }
  if (!fromOk) print("Error in function 'mapVector: 'fromVector' not ascending")

  toOk <- TRUE
  n    <- 1
  while (n < toSize & toOk) {
    n    <- n+1
    toOk <- toVector[n] > toVector[n-1]
  }
  if (!toOk) print("Error in function 'mapVector': 'toVector' not ascending")

  # 'window' must be non-negative, check:
  if ( window < 0) print("Error in function 'mapVector': input 'window'<0; invalid argument")

  # ranges must overlap, check:
  rangeOk <- toVector[1]      <= fromVector[fromSize]+window &
             toVector[toSize] >= fromVector[1]       -window

  if (!fromOk | !toOk | !rangeOk | window < 0) return(vectorMap)


  # Good to go!


  # Find first, last relevant element in 'toVector':
  toFirst <- 1
  while (toVector[toFirst] < fromVector[1]       -window) toFirst <- toFirst+1
  toLast  <- toSize
  while (toVector[toLast]  > fromVector[fromSize]+window) toLast  <- toLast -1

  # Find first, last relevant element in 'fromVector':
  fromFirst <- 1
  while (fromVector[fromFirst] < toVector[1]     -window) fromFirst <- fromFirst+1
  fromLast  <- fromSize
  while (fromVector[fromLast]  > toVector[toSize]+window) fromLast  <- fromLast -1

  #Map!
  if (fromLast >= fromFirst)
    for (fromN in fromFirst:fromLast) {
      toN <- toFirst
      while (fromVector[fromN] >  toVector[toN]+window) toN <- toN+1
      if    (fromVector[fromN] >= toVector[toN]-window) vectorMap[fromN] <- toN
    }

  return(vectorMap)

}

