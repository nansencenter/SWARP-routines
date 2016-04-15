
#R The function 'categoryOverlap' compares two input vectors
#R element-by-element, and checks if each element pair
#R belongs to the same category. The function returns this
#R count for all categories. Categories are defined by a
#R vector 'limits' which specify the bounds for the various
#R categories. Category no. i is defined as
#R  <limits[i], limits[i+1]]


# Example 1:
# ...input:
# values1        =  8 12 NA  5 17 -1 22  0
# values2        = 10 10  8 NA 19  1 19  5
# limits         =  0  5 10 15 20
# ...output:
#                =   0  1  0  1
# ...explanation:
#    algorithm finds the index of the interval in 'limits'
#     (interval 1 is  <limits[1], limits[2]] etc.)
#     for each element in 'values1', and increments the
#     corresponding element in the output vector if the
#     same entry in the other input vector ('values2')
#     belongs to the same 'limits' category


# Example 2 (same as 1, but with 'limits[1]' = 'limits[5]' = NA):
# ...input:
# values1        =  8 12 NA  5 17 -1 22  0
# values2        = 10 10  8 NA 19  1 19  5
# limits         = NA  5 10 15 NA
# ...output:
#                =   2  1  0  2


# Example 3 (same as 1, but with 'limits[1]' = NA):
# ...input:
# values1        =  8 12 NA  5 17 -1 22  0
# values2        = 10 10  8 NA 19  1 19  5
# limits         = NA  5 10 15 20
# ...output:
#                =   2  1  0  1


# Example 4 (same as 1, but with 'limits[5]' = NA):
# ...input:
# values1        =  8 12 NA  5 17 -1 22  0
# values2        = 10 10  8 NA 19  1 19  5
# limits         =  0  5 10 15 NA
# ...output:
#                =   0  1  0  2


# Example 5: counting all matching values <=10
# ...input:
# values1        =  8 12 NA  5 17 -1 22  0
# values2        = 10 10  8 NA 19  1 19  5
# limits         =  NA 10
# ...output:
#                =   3


categoryOverlap <- function(values1, values2, limits) {


# Check that limits are OK

nLimits     <- length(limits)
nLimVals    <- sum(!is.na(limits))

# 'limits' must have at least two entries, check:
okLimits <- nLimits > 1
if (!okLimits) print("Error in function 'categoryOverlap': input 'limits' vector has less than two values")

# 'limits' must have at least one actual value, check:
okLimVals <- nLimVals > 0
if (!okLimVals) print("Error in function 'categoryOverlap': input 'limits' vector only has NA entries")

# 'limits' must be ascending, check:
first <- 1
last  <- nLimits
if (is.na(limits[first])) first <- first + 1
if (is.na(limits[last] )) last  <- last  - 1
n <- first
ascending <- !is.na(limits[n])
while (n < last & ascending) {
  n <- n + 1
  ascending <- limits[n] > limits[n-1] & !is.na(limits[n])
}
if (!ascending) print("Error in function 'categoryOverlap': input 'limits' vector is not ascending properly")


# Initialize distribution overlap count vector:
nCategories <- nLimits - 1
overlap     <- array(data=0, dim=c(nCategories))


nGoodValues1 <- sum(!is.na(values1))
nValues1     <- length(values1)
nGoodValues2 <- sum(!is.na(values2))
nValues2     <- length(values2)
nValues      <- min(nValues1,nValues2)

if (!okLimits | !okLimVals | nGoodValues1 == 0 | nGoodValues2 == 0)
  return(overlap)


# There is work to do!

if (!is.na(limits[1]) & !is.na(limits[nLimits])) {

  # CASE 1: no NA in 'limits'

  for (n in 1:nValues)
    if (!is.na(values1[n])            & !is.na(values2[n])      &
        values1[n] >  limits[1]       & values2[n] >  limits[1] &
        values1[n] <= limits[nLimits] & values2[n] <= limits[nLimits]) {
      i <- 1
      while (values1[n] > limits[i+1]) i <- i + 1
      if (values2[n] > limits[i] & values2[n] <= limits[i+1])
        overlap[i] <- overlap[i] + 1
    }

}
else if  (!is.na(limits[1]) & is.na(limits[nLimits])) {

  # CASE 2: NA in upper 'limits' only

  for (n in 1:nValues)
    if (!is.na(values1[n]) & values1[n] >  limits[1] &
        !is.na(values2[n]) & values2[n] >  limits[1]) {
      if (values1[n] > limits[nLimits-1] &
          values2[n] > limits[nLimits-1])
        overlap[nLimits-1] <- overlap[nLimits-1] + 1
      else if (values1[n] <= limits[nLimits-1]) {
        i <- 1
        while (values1[n] > limits[i+1]) i <- i + 1
        if (values2[n] > limits[i] & values2[n] <= limits[i+1])
          overlap[i] <- overlap[i] + 1
      }
    }

}
else if  (is.na(limits[1]) & !is.na(limits[nLimits])) {

  # CASE 3: NA in lower 'limits' only

  for (n in 1:nValues)
    if (!is.na(values1[n]) & values1[n] <= limits[nLimits] &
        !is.na(values2[n]) & values2[n] <= limits[nLimits]) {
      if (values1[n] <= limits[2] & values2[n] <= limits[2])
        overlap[1] <- overlap[1] + 1
      else if (values1[n] > limits[2]) {
        i <- 2
        while (values1[n] > limits[i+1]) i <- i + 1
        if (values2[n] > limits[i] & values2[n] <= limits[i+1])
          overlap[i] <- overlap[i] + 1
      }
    }

}
else {

  # CASE 4 (remaining): NA in both 'limits'

  for (n in 1:nValues)
    if (!is.na(values1[n]) & !is.na(values2[n])) {
      if (values1[n] <= limits[2] & values2[n] <= limits[2])
        overlap[1] <- overlap[1] + 1
      else if (values1[n] > limits[nLimits-1] &
               values2[n] > limits[nLimits-1])
        overlap[nLimits-1] <- overlap[nLimits-1] + 1
      else if (values1[n] >  limits[2] &
               values1[n] <= limits[nLimits-1]) {
        i <- 2
        while (values1[n] > limits[i+1]) i <- i + 1
        if (values2[n] > limits[i] & values2[n] <= limits[i+1])
          overlap[i] <- overlap[i] + 1
      }
    }
}


return(overlap)

}
