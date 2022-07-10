# custom functions for use in ddPCR data analysis
# B. Layton



# a handy function for pretty figures and legends
unlog <- function(x) {round(10^x, digits = 0)}

# function for calculating standard error
se <- function(x) sd(x, na.rm = T)/sqrt(length(x))



#### Custom aggregation functions ####
# function to count positives in QC samples
blankMaxFun <- function (x){ # allows for corrections for re-runs when calling results for blank samples
  goodPlate <- max(x$Plate) # assumes most recent plate is one to use
  blankMax <- max(x$Positives[which(x$Plate == goodPlate)], na.rm = T)
  return(blankMax)
}

# numeric aggregation function to feed into dcast (needed bc mean will be skewed low if any rep is 0):
# this is no longer needed for most data bc of 1/2 LOD assignment, but comes in handy for maps and OVDL data
partialMean <- function(x) {
  if (sum(x, na.rm = T) == 0) {
    return(0) # if all are 0, return 0
  } else if (any(x == 0, na.rm = TRUE)) {
    x <- x[which(x != 0)]
    return(mean(x, na.rm = T)) # if some are 0, remove 0's and return the mean
  } else
    return(mean(x, na.rm = T)) # if no 0's, return the mean
}

# function for aggregating qualitative results - for use with dcast on QualCall variable only:
partialQualAggr <- function(x) { 
  if (all(is.na(x))) {
    return("NA") # if all are NAs - meaning all reps failed QC due to low droplets or all targets 0 - return NA
  } else if (all(str_detect(x, "Negative"), na.rm = TRUE)) {
    return("Negative") # if all reps are negative, return negative
  } else if (any(str_detect(x, "Positive"), na.rm = TRUE)) {
    return("Positive") # if any rep is positive, return positive
  } else {
    return("Inconclusive") # otherwise inconclusive
  }
}

# function for aggregating QCreason - for use with dcast on QCreason variable only:
aggrQCreason <- function(x) { 
  x <- str_replace_na(as.character(x))
  if (all(str_detect(x, "NA"))) {
    return("NA") # if all reps are NA, return NA
  } else {
    reasons <- str_c(x, collapse = ", ") # otherwise concatenate into a single string
    return(reasons) 
  }
}

# function for pulling out the lower 95% confidence interval from t.test 
lowerConfInf <- function(x) {
  if (all(is.na(x))) return(NA) # needed to avoid breaking the t.test function
  else {
    nx <- length(x)
    mx <- mean(x, na.rm = T)
    vx <- var(x, na.rm = T)
    if (is.na(vx)) return(NA) # needed to avoid breaking the t.test function
    else {
      stderr <- sqrt(vx/nx)
      if (stderr < 10 * .Machine$double.eps * abs(mx)) return(NA) # needed because negatives have all identical values for CalcCopies due to 1/2 LOD assignment
      else { 
        t <- t.test(x)
        return(t$conf.int[1])
      } 
    }
  }
}

# function for pulling out the upper 95% confidence interval from t.test 
upperConfInf <- function(x) {
  if (all(is.na(x))) return(NA) # needed to avoid breaking the t.test function
  else {
    nx <- length(x)
    mx <- mean(x, na.rm = T)
    vx <- var(x, na.rm = T)
    if (is.na(vx)) return(NA) # needed to avoid breaking the t.test function
    else {
      stderr <- sqrt(vx/nx)
      if (stderr < 10 * .Machine$double.eps * abs(mx)) return(NA) # needed because negatives have all identical values for CalcCopies due to 1/2 LOD assignment
      else {
        t <- t.test(x)
        return(t$conf.int[2])
      }
    }
  }
}

