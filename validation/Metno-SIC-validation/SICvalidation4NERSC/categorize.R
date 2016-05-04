
#R The function 'categorize' replaces entries in an input vector
#R with values provided for intervals in input values
#R (may e.g. be used to replace all values in <n*10,(n+1)*10]
#R  by n*10+5 for n= 0,1,...,10),
#R and returns a list variable with two vectors:
#R  $values:       a vector w/ the new values
#R                  (out-of-range values returned unaltered)
#R  $distribution: frequency used for each interval (count)

# Example 1:
# ...input:
# values         =  8 12 NA  5 20 -1 22  0
# limits         =  0  5 10 15 20
# categoryValues =   10 20 30 40
# ...output:
#  $values       = 20 30 NA 10 40 -1 22 0
#  $distribution =  1  1  1  1
# ...explanation:
#    algorithm finds the index of the interval in 'limits'
#     (interval 1 is  <limits[1], limits[2]] etc.)
#     and substitutes the value with the corrresponding
#     value in categoryValue


# Example 2 (same as 1, but with 'limits[1]' = 'limits[5]' = NA):
# ...input:
# values         =  8 12 NA  5 20 -1 22  0
# limits         = NA  5 10 15 NA
# categoryValues =   10 20 30 40
# ...output:
#  $values       = 20 30 NA 10 40 10 40 10
#  $distribution =  3  1  1  2


# Example 3 (same as 1, but with 'limits[1]' = NA):
# ...input:
# values         =  8 12 NA  5 20 -1 22  0
# limits         = NA  5 10 15 20
# categoryValues =   10 20 30 40
# ...output:
#  $values       = 20 30 NA 10 40 10 22 10
#  $distribution =  3  1  1  1


# Example 4 (same as 1, but with 'limits[5]' = NA):
# ...input:
# values         =  8 12 NA  5 20 -1 22  0
# limits         =  0  5 10 15 NA
# categoryValues =   10 20 30 40
# ...output:
#  $values       = 20 30 NA 10 40 -1 40  0
#  $distribution =  1  1  1  2


# Example 5: replacing all values <=10 with 0
# ...input:
# values         =  8 10 NA  5 20 -1 22  0
# limits         =  NA 10
# categoryValues =    0
# ...output:
#  $values       =  0  0 NA  0 20  0 22  0
#  $distribution =  5


# Example 6: same as 5, but replacing with NA
# ...input:
# values         =  8 10 NA  5 20 -1 22  0
# limits         =  NA 10
# categoryValues =   NA
# ...output:
#  $values       = NA NA NA NA 20 NA 22 NA
#  $distribution =  5


categorize <- function(values, limits, categoryValues) {


# Create the output list:
categorizeOut <- list(numeric="values", numeric="distribution")

# Check that limits, categoryValues are OK

nLimits     <- length(limits)
nLimVals    <- sum(!is.na(limits))
nCategories <- length(categoryValues)

# Initialize distribution count vector:
distribution <- array(data=0, dim=c(nCategories))


# 'limits' must have at least two entries, check:
okLimits <- nLimits > 1
if (!okLimits) print("Error in function 'categorize': input 'limits' vector has less than two values")

# 'limits' must have at least one actual value, check:
okLimVals <- nLimVals > 0
if (!okLimVals) print("Error in function 'categorize': input 'limits' vector only has NA entries")

# There must be one more limiting value than #categories, check:
okCategories <- nLimits - nCategories == 1
if (!okCategories) print("Error in function 'categorize': non-matching no. of categories and limits (latter must be one larger than former")

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
if (!ascending) print("Error in function 'categorize': input 'limits' vector is not ascending properly")


newValues   <- values
nGoodValues <- sum(!is.na(newValues))
nValues     <- length(newValues)

if (!okLimits | !okLimVals | !okCategories | nGoodValues == 0)
  return(newValues)


# There is work to do!

if (!is.na(limits[1]) & !is.na(limits[nLimits])) {

  # CASE 1: no NA in 'limits'

  for (n in 1:nValues)
    if (!is.na(newValues[n]) & newValues[n] >  limits[1] &
                               newValues[n] <= limits[nLimits]) {
      i <- 1
      while (newValues[n] > limits[i+1]) i <- i + 1
      newValues[n]    <- categoryValues[i]
      distribution[i] <- distribution[i] + 1
    }

}
else if  (!is.na(limits[1]) & is.na(limits[nLimits])) {

  # CASE 2: NA in upper 'limits' only

  for (n in 1:nValues)
    if (!is.na(newValues[n]) & newValues[n] >  limits[1]) {
      if (newValues[n] > limits[nLimits-1]) {
        newValues[n]            <- categoryValues[nLimits-1]
        distribution[nLimits-1] <- distribution[nLimits-1] + 1
      }
      else {
        i <- 1
        while (newValues[n] > limits[i+1]) i <- i + 1
        newValues[n]    <- categoryValues[i]
        distribution[i] <- distribution[i] + 1
      }
    }

}
else if  (is.na(limits[1]) & !is.na(limits[nLimits])) {

  # CASE 3: NA in lower 'limits' only

  for (n in 1:nValues)
    if (!is.na(newValues[n]) & newValues[n] <= limits[nLimits]) {
      if (newValues[n] <= limits[2]) {
        newValues[n]    <- categoryValues[1]
        distribution[1] <- distribution[1] + 1
      }
      else {
        i <- 2
        while (newValues[n] > limits[i+1]) i <- i + 1
        newValues[n]    <- categoryValues[i]
        distribution[i] <- distribution[i] + 1
      }
    }

}
else {

  # CASE 4 (remaining): NA in both 'limits'

  for (n in 1:nValues)
    if (!is.na(newValues[n])) {
      if (newValues[n] <= limits[2]) {
        newValues[n]    <- categoryValues[1]
        distribution[1] <- distribution[1] + 1
      }
      else if (newValues[n] > limits[nLimits-1]) {
        newValues[n]            <- categoryValues[nLimits-1]
        distribution[nLimits-1] <- distribution[nLimits-1] + 1
      }
      else {
        i <- 2
        while (newValues[n] > limits[i+1]) i <- i + 1
        newValues[n]    <- categoryValues[i]
        distribution[i] <- distribution[i] + 1
      }
    }
}

categorizeOut$values       <- newValues
categorizeOut$distribution <- distribution

return(categorizeOut)

}
