
rootMeanSquare <- function(array1, array2) {

# Check dimensions:
size1 <- length(array1)
size2 <- length(array2)

if (size1 != size2 | size1 == 0)
  rmsValue <- NA
else {
  n     <- 0
  sqSum <- 0
  for (i in 1:size1)
    if (!is.na(array1[i]) & !is.na(array2[i])) {
      n     <- n     + 1
      sqSum <- sqSum + (array2[i] - array1[i])**2
    }
  if (n == 0)
    rmsValue <- NA
  else
    rmsValue <- sqrt(sqSum/n)
}

return(rmsValue)

}
