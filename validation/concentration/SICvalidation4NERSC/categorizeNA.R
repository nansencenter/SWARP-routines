
#R The function 'categorizeNA' is identical to 'categorize',
#R except that NA values are replaced by a user provided, 
#R specified value.


categorizeNA <- function(values, limits, categoryValues, NAvalue) {


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
nValues     <- length(newValues)

if (!okLimits | !okLimVals | !okCategories)
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
    else if (is.na(newValues[n])) newValues[n] <- NAvalue
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
    else if (is.na(newValues[n])) newValues[n] <- NAvalue

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
    else if (is.na(newValues[n])) newValues[n] <- NAvalue

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
    else if (is.na(newValues[n])) newValues[n] <- NAvalue

}

categorizeOut$values       <- newValues
categorizeOut$distribution <- distribution

return(categorizeOut)

}
