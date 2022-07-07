# custom functions for use in covid ddPCR data analysis
# author: B. Layton unless otherwise noted



# a handy function for pretty figures and legends
unlog = function(x) {round(10^x, digits = 0)}

# function for calculating standard error
se = function(x) {sd(x, na.rm = T)/sqrt(length(x))}

#### Custom aggregation functions ####
# function for aggregating qualitative results - for use with dcast on QualCall variable only:
partialQualAggr = function(x) { 
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
aggrQCreason = function(x) { 
  x <- str_replace_na(as.character(x))
  if (all(str_detect(x, "NA"))) {
    return("NA") # if all reps are NA, return NA
  } else {
    reasons <- str_c(x, collapse = ", ") # otherwise concatenate into a single string
    return(reasons) 
  }
}

# function for pulling out the lower 95% confidence interval from t.test 
lowerConfInf = function(x) {
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
upperConfInf = function(x) {
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
