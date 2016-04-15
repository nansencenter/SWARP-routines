#R The function 'MIZconfusionMatrix' is a variant of the function
#R 'confusionMatrix', which adds two categories for the observational
#R product (input 'values1'); these will be the 2nd category and the
#R 2nd-to-last category, which correspond to the first and last
#R categories, respectively, but is inside a footprint of grids with
#R values outside both the first and last categories. The mapping is
#R given by input 'MIZgrid', which will be -2 and 2 for grids in the
#R 2nd category and the 2nd-to-last category, respectively.

#R The no. of occurence of all vector element pairs is returned as
#R 'contingencyTable'. If one or both elements of a pair is NA, the
#R table is not updated.

#R For input 'values2', categories are defined by a vector 'limits'
#R which specifies the bounds for the various categories. Category
#R no. i is defined as <limits[i], limits[i+1]]
#R Categories for 'values1' id similar, but with two additional
#R category values as described above.


# Example:
# ...input:
# values1        =  8 12 NA  5 17 -1 22  1
# values2        = 10 10  8 NA 19  1 19  0
# MIZgrid        =  1  1 NA  1  1  0  2  1
# limits         =  0  5 10 15 20
# ...output:
#                =   0  1  0  0  0  0
#                    0  0  0  0  0  0
#                    1  0  0  0  0  0
#                    0  0  1  0  0  0
#                    0  0  1  0  0  0
#                    0  0  0  0  1  0
#                    0  0  0  0  1  0
#                    0  0  0  0  0  0
# ...explanation:
#    algorithm finds the index of the interval in 'limits'
#     (interval no. 1 is <-inf     , limits[1]],
#               no. 2 is <limits[1], limits[2]] etc.)
#     for each element pair in 'values1' and 'values2', and
#     increments the corresponding element in 'contingencyTable'
#     e.g. (values1[i],values2[i]) = (22,19) & MIZgrid[i] = 2
#         leads to an update
#     of 'contingencyTable[7,5]'


MIZconfusionMatrix <- function(values1, values2, MIZgrid, limits) {


# Check that limits are OK

nLimits     <- length(limits)
nLimVals    <- sum(!is.na(limits))

# 'limits' cannot include NA values, check:
okLimVals <- nLimits == nLimVals
if (!okLimVals) {
  print("Error in function 'confusionMatrix': input 'limits' vector contains NA entries")
} else {

  # 'limits' must be ascending, check:
  n <- 1
  ascending <- TRUE
  while (n < nLimits & ascending) {
    n <- n + 1
    ascending <- limits[n] > limits[n-1]
  }
  if (!ascending) print("Error in function 'categoryOverlap': input 'limits' vector is not ascending properly")

}


# Initialize distribution overlap count vector:
nCategories      <- nLimits + 1
contingencyTable <- array(data=0, dim=c(nCategories+2,nCategories))

v1size <- length(values1)
v2size <- length(values2)
if (v1size != v2size) {
  print("Warning from function 'confusionMatrix': input vectors have different sizes, only entries where both vectors have values will be considered")
  nPairs <- min(v1size,v2size)
} else {
  nPairs <- v1size
}

if (nPairs > 0) {
  v1 <- values1[1:nPairs]
  g1 <- MIZgrid[1:nPairs]
  v2 <- values2[1:nPairs]
  nGoodPairs <- sum(!is.na(v1+v2))
  if (nGoodPairs==0) print("Warning from function 'confusionMatrix': all pairs in input vectors include an NA value")
} else {
  nGoodPairs <- 0
}

if (!okLimVals | nGoodPairs == 0) return(contingencyTable)


# Got here? There is work to do!


# Find categories
category1 <- array(0,dim=length(v1))
category2 <- array(0,dim=length(v2))
category1[which(v1<=limits[1] & g1!=-2)] <- 1
category1[which(v1<=limits[1] & g1==-2)] <- 2
category2[which(v2<=limits[1])]          <- 1
if (nLimits > 1)
  for (n in 2:nLimits) {
    category1[which(v1>limits[n-1] & v1<=limits[n])] <- n + 1
    category2[which(v2>limits[n-1] & v2<=limits[n])] <- n
  }
category1[which(v1>limits[nLimits] & g1==2)] <- nCategories + 1
category1[which(v1>limits[nLimits] & g1!=2)] <- nCategories + 2
category2[which(v2>limits[nLimits])]         <- nCategories


# Produce 'contingencyTable'
for (n1 in 1:(nCategories+2))
  for (n2 in 1:nCategories)
    contingencyTable[n1,n2] <- sum(category1==n1 & category2==n2)


return(contingencyTable)

}
