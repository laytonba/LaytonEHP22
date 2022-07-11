# COVID in wastewater triplex ddPCR assay data - data prep and QC
# started 6/26/2020 
# last edit July 10, 2022
# author: Blythe Layton

# script to create data frames used for analyses and figures related to 
# covid ddPCR data presented in Layton et al. EHP 2022
# https://ehp.niehs.nih.gov/doi/10.1289/EHP10289

# This is not the analysis script. This script generates aggregated, QC'd data
# ready for analysis. See https://github.com/laytonba/LaytonEHP22 for analysis
# script.



#### set up: import packages, functions, and data ####

library(ggplot2); library(reshape2); library(grid); 
library(chron); library(stringr); library(plyr); 
library(schoolmath); library(openxlsx); #library(lubridate)

theme_set(theme_bw()) # I prefer the bw background for graphics
Sys.setenv(TZ = 'GMT') # becomes important for time series plots

# load required functions to execute this script
source("https://raw.githubusercontent.com/laytonba/LaytonEHP22/main/ddPCR_functions.R")


## import data
# ddPCR data
data <- read.csv("https://raw.githubusercontent.com/laytonba/LaytonEHP22/main/ddPCR_data_for_github.csv")
# metadata
volFilt <- read.csv("https://raw.githubusercontent.com/laytonba/LaytonEHP22/main/mlFiltered.csv")
rna <- read.csv("https://raw.githubusercontent.com/laytonba/LaytonEHP22/main/RNAplates_github.csv")


#### Create main dataframe for analysis and QC ####

## create Site and Date variables from sample ID
splitCol <- as.data.frame(str_split_fixed(data$Sample, pattern = "-", 6))
data$Site <- str_c(splitCol$V1, splitCol$V2, sep = "-")
data$Date <- as.chron(str_c(splitCol$V3, splitCol$V4, splitCol$V5, sep = "/"))

## create Location factor variable
data$Location <- NA
data$Location[grep("Cor-", data$Sample, ignore.case = T)] <- "Corvallis"
data$Location[grep("EUG", data$Sample)] <- "Eugene"
data$Location[grep("Her", data$Sample, ignore.case = T)] <- "Hermiston"
data$Location[grep("NP-", data$Sample)] <- "Newport"
data$Location[grep("RED", data$Sample, ignore.case = T)] <- "Redmond"

## create SampleType factor variable 
data$SampleType <- NA
data$SampleType[grep("Inf", data$Sample, ignore.case = TRUE)] <- "Influent"
data$SampleType[which(data$Location == "Newport" & is.na(data$SampleType))] <- "Pump Station"
data$SampleType[grep("EB", data$Sample)] <- "Control"
data$SampleType[grep("NTC", data$Sample)] <- "Control"
data$SampleType[grep("NoTC", data$Sample)] <- "Control"
data$SampleType[grep("FB", data$Sample)] <- "Control"
data$SampleType[grep("Blank", data$Sample)] <- "Control"
data$SampleType[grep("Field", data$Sample, ignore.case = TRUE)] <- "Control"
data$SampleType[grep("BLNK", data$Sample)] <- "Control" 
data$SampleType[grep("Ctl", data$Sample)] <- "Control"
data$SampleType[grep("Ctrl", data$Sample)] <- "Control"
data$SampleType[grep("Pos", data$Sample)] <- "Control"
data$SampleType <- as.factor(data$SampleType)

## create ControlType factor variable 
data$ControlType <- NA
data$ControlType[grep("EB", data$Sample)] <- "EB" # extraction blank
data$ControlType[grep("FB", data$Sample)] <- "FB" # field blank
data$ControlType[grep("Blank", data$Sample, ignore.case = T)] <- "FB" 
data$ControlType[grep("Field", data$Sample, ignore.case = T)] <- "FB" 
data$ControlType[grep("BLNK", data$Sample, ignore.case = T)] <- "FB" 
data$ControlType[grep("Neg", data$Sample)] <- "Negative"
data$ControlType[grep("Pos", data$Sample)] <- "Positive"
data$ControlType[grep("NTC", data$Sample, ignore.case = T)] <- "NTC" # no template control
data$ControlType[grep("NoTC", data$Sample, ignore.case = T)] <- "NTC"
data$ControlType <- as.factor(data$ControlType)


## create Master Sample variable
# (allows pooling of filter replicates)
data$MasterID <- str_c(data$Site, data$Date, sep = "-")
data$MasterID <- str_replace_all(data$MasterID, "/", "-")

## create County variable 
data$County <- NA
data$County[which(data$Location == "Newport")] <- "Lincoln"
data$County[which(data$Location == "Bend" | data$Location == "Redmond")] <- "Deschutes"
data$County[which(data$Location == "Corvallis")] <- "Benton"
data$County[which(data$Location == "Eugene")] <- "Lane"
data$County[which(data$Location == "Hermiston")] <- "Umatilla"
data$County <- as.factor(data$County)



#### Attach metadata ####
# add ml filtered to main data frame
M <- match(data$Sample, volFilt$ID, nomatch = NA)
data$mlFiltered <- volFilt$Volume.filtered.mL[M]

# look for mismatches - indicates a data entry error somewhere
mismatch <- data[which(is.na(data$mlFiltered) & data$SampleType != "Control"), ]
print("Volume filtered mismatches:")
print(count(as.factor(mismatch$Sample))) # mismatches to fix....


##  RNA plate data  ##

rna$PCRplate <- str_remove_all(rna$PCRplate, " ")
rna$PCRplate <- strsplit(rna$PCRplate, split = ",") # converts to a list

data$RNAplate <- NA
data$RNAwell <- NA
data$RNA_ID <- NA

for (i in 1:nrow(rna)) { # from Blythe 6-22-21
  for (j in 1:length(rna$PCRplate[[i]])) { # using a nested loop to cycle through all elements of the PCRplate list
    d <- which(data$Sample == rna$Sample[i] & data$Plate == rna$PCRplate[[i]][j])
    data$RNAwell[d] <- rna$RNAwell[i]
    data$RNAplate[d] <- rna$RNAplate[i]
    data$RNA_ID[d] <- rna$RNA_ID[i]
  }
}

# check for mismatches
matchfail <- data[which(data$SampleType != "Control" & is.na(data$RNAwell) == T),]
print("RNA plate mismatches:")
print(count(matchfail$Sample))


## create Lysate Volume variable; indicates volume of lysate extracted 
m <- match(data$RNA_ID, rna$RNA_ID)
data$LysateVol <- rna$uLextracted[m] 
data$LysateVol[which(is.na(data$LysateVol) & data$Plate < 6)] <- 400 # a lot of controls have NAs bc of imperfect matching, or bc they weren't actually extracted (like NTCs)
data$LysateVol[which(is.na(data$LysateVol) & data$Plate >= 6)] <- 200


## create elution volume variable 
data$ElutionVol <- rna$uLeluted[m] # using same matching variable as previous section
data$ElutionVol[which(is.na(data$ElutionVol))] <- 30 # a lot of controls have NAs bc they weren't extracted. Only really affects Field blanks


#### QC ####
## score NTCs & negatives: 
# find the maximum number of positive droplets in N targets of NTCs/negative 
# controls on each ddPCR plate
data$NTCmax <- NA # note that this is in droplets
for (i in 1:max(data$Plate)) {
  if (length(which(data$Plate == i & data$ControlType == "NTC")) > 0) { 
    # this 'if/else' allows for some plates to not have NTCs, where the EB is "lab blank" for both
    data$NTCmax[which(data$Plate == i)] <- max(data$Positives[which((data$ControlType == "NTC" | data$ControlType == "Negative") & data$Target != "RP" & data$Plate == i)])
  }
  else data$NTCmax[which(data$Plate == i)] <- NA
}


## score EBs: 
# find the maximum number of positive droplets in the extraction blank for each 
# RNA plate, and attach this number to all samples on respective RNA plate
data$EBmax <- NA
ebs <- data[which(data$ControlType == "EB"), ]
ebMax <- ddply(ebs, ~ RNAplate, blankMaxFun)
for (i in 1:nrow(ebMax)) { 
  ebPCRplate <- which(data$RNAplate == ebMax$RNAplate[i])
  if (length(ebPCRplate) > 0) {
    data$EBmax[ebPCRplate] <- ebMax$V1[i]
  }
}
# if NTC amplified, EB result is inconclusive.
data$EBmax[which(data$NTCmax > 2 & data$EBmax > 2)] <- NA 


## score FBs: 
# find the maximum number of positive droplets in each field blank, and attach 
# this number to the samples associated with each field blank
# note that this is max of N targets only
data$FBmax <- NA
fbs <- data[which(data$ControlType == "FB" & data$Target != "RP"), ]
fbMax <- ddply(fbs, ~ Date + County + MasterID, blankMaxFun) 
for (i in 1:nrow(fbMax)) { 
  # using a date range of +/- 1 day from FB date. not perfect batching bc 
  # newport had fbs at almost every site for 2nd round of PS sampling
  fbPCRplate <- which(data$County == fbMax$County[i] &
                        (data$Date - 1 == fbMax$Date[i] |
                           data$Date  == fbMax$Date[i] |
                           data$Date + 1  == fbMax$Date[i]))
  data$FBmax[fbPCRplate] <- fbMax$V1[i]
}
# if EB amplified, FB result is inconclusive
data$FBmax[which(data$EBmax > 2 & data$FBmax > 2)] <- NA 
# if NTC amplified, FB result is inconclusive
data$FBmax[which(data$NTCmax > 2 & data$FBmax > 2)] <- NA 


## create QC pass/fail factor
# using a 3-droplet threshold for a positive call
data$QC <- "PASS"
# our ADG was routinely giving 8-12K droplets, so we set 6K as the threshold
data$QC[which(data$Accepted.Droplets < 6000)] <- "FAIL"
# amplification in any target of blank controls is fail
data$QC[which((data$ControlType == "NTC" |  data$ControlType == "EB") &
                data$Positives > 2)] <- "FAIL" 
# amplification in N targets of NegCtl or FB is fail
data$QC[which((data$ControlType == "Negative" | data$ControlType == "FB") & 
                data$Target != "RP" & 
                data$Positives > 2)] <- "FAIL" 
# note higher positive threshold for RP in FBs
data$QC[which(data$ControlType == "FB" & data$Target == "RP" &
                data$Positives > 5)] <- "FAIL"
# positive sample data are QC fail if QC samples failed
data$QC[which((data$SampleType != "Control" &
                 data$Positives > 2 &
                 data$Target != "RP") &
                (data$NTCmax > 2 | data$EBmax > 2 | data$FBmax > 2))] <- "FAIL" 


## create Qualitative Call variable
data$QualCall <- "Negative"
data$QualCall[which(data$Positives > 2)] <- "Positive" 
data$QualCall[which(data$QC == "FAIL")] <- "Inconclusive"  
data$QualCall[which(data$Accepted.Droplets < 6000)] <- NA
data$QualCall <- as.factor(data$QualCall)


## create QC Reason variable 
# (indicates why a sample failed QC)
# the ordering matters here...
data$QCreason <- NA 
data$QCreason[which(data$QualCall == "Inconclusive" & data$FBmax > 2)] <- "FB amplified"
data$QCreason[which(data$QualCall == "Inconclusive" & data$EBmax > 2)] <- "EB amplified"
data$QCreason[which(data$QualCall == "Inconclusive" & data$NTCmax > 2)] <- "NTC amplified"
data$QCreason[which(data$QC == "FAIL" & data$SampleType == "Control")] <- "target amplification"


## find samples that have 0 positive droplets for all three targets
# samples with 0 RP are an extraction failure; however if N targets amplify even
# in absence of RP it's still a valid result (per Bio-Rad protocol)
w <- dcast(subset(data, SampleType != "Control" & QC != "FAIL"), Sample + Plate + Well ~ Target, value.var = c("Positives"))
bad <- w[which(w$N1 == 0 & w$N2 == 0 & w$RP == 0), ]
for (i in 1:nrow(data)) { 
  if (is.na(match(data$Sample[i], bad$Sample)) == F) {
    badRow <- which(data$Sample[i] == bad$Sample & data$Plate[i] == bad$Plate & data$Well[i] == bad$Well)
    if (length(badRow) > 0) {
      data$QC[i] <- "FAIL"
      data$QualCall[i] <- NA
      data$QCreason[i] <- "all targets 0"
    }
  }
}

data$QCreason[which(data$Accepted.Droplets < 6000)] <- "low droplets"
data$QCreason <- as.factor(data$QCreason)
data$QC <- as.factor(data$QC)


#### Quantification ####

## calculate copies per reaction
uLRxn <- 20 # total volume of the ddPCR reaction. 22 is prepared, but ADG only aspirates 20
data$CopiesPerRxn <- data$CopiesPeruL * uLRxn 


## create clean copy number variable: 
# QC fails removed and inconclusive/BLOQ hits removed
data$CleanCopies <- data$CopiesPeruL
data$CleanCopies[which(data$QC == "FAIL")] <- NA # censor QC failed data
# set copy numbers to 0 for samples below positive threshold; this is changed to 1/2 LOD below
data$CleanCopies[which(data$QualCall == "Negative")] <- 0 


## calculate copies per L of sample
ShieldVol <- 1000 # uL. this was constant for data used in this study
data$CopiesPerL <- data$CleanCopies * uLRxn / data$TemplateVol * data$ElutionVol / data$LysateVol * ShieldVol / data$mlFiltered * 1000


## sample-specific LOD
# since samples had variable volumes of ml filtered, lysate extracted, and 
# RNA template analyzed, the theoretical LOD varies by sample
# using a theoretical limit of 3 copies/reaction, consistent with our 3 droplet threshold
data$SampleLOD <- 3 / data$TemplateVol * data$ElutionVol / data$LysateVol * ShieldVol / data$mlFiltered * 1000 

## assign 1/2 LOD values to non-detects
data$CalcCopies <- data$CopiesPerL # create new variable to be used in calculating aggregates
data$CalcCopies[which(data$QualCall == "Negative")] <- 0.5 * data$SampleLOD[which(data$QualCall == "Negative")]

## log-transform results
data$logCopies <- log10(data$CalcCopies) 


#### Aggregate results ####
## aggregate clean copy numbers of N1/N2 data, with bad QC data removed 
# note that we're using 1/2 LOD values for nondetects and these are included in 
# the means (logCopies comes from CalcCopies)

# now using MasterID to pool data across A and B filter replicates

s <- subset(data, Target != "RP" & 
              SampleType != "Control" & 
              QC != "FAIL")

samples.aggr <- ddply(s, ~ County + Location + Site + Date + SampleType + MasterID,
                            summarize, 
                      LogCopiesPerL = round(mean(logCopies, na.rm = T), digits = 3), 
                      upperCI = round(upperConfInf(logCopies), digits = 3), 
                      lowerCI = round(lowerConfInf(logCopies), digits = 3), 
                      StdErr = signif(se(logCopies), digits = 3), 
                      LogLOD = round(log10(mean(SampleLOD, na.rm = T)), digits = 3), 
                      QualCall = partialQualAggr(QualCall), 
                      aQCreason = aggrQCreason(QCreason), 
                      PCRplate = tail(Plate, n = 1))

# clean up for tidy plots and sharing
samples.aggr$upperCI[which(samples.aggr$QualCall == "Negative")] <- NA
samples.aggr$lowerCI[which(samples.aggr$QualCall == "Negative")] <- NA
samples.aggr$StdErr[which(samples.aggr$QualCall == "Negative")] <- NA
samples.aggr$QualCall[which(samples.aggr$QualCall == "Negative")] <- "Negative, 1/2 LOD assigned" 
samples.aggr$QualCall <- as.factor(samples.aggr$QualCall)
samples.aggr$aQCreason <- as.factor(samples.aggr$aQCreason)

head(samples.aggr)
print(tail(samples.aggr))
samples.aggr <- samples.aggr[, -14]
samples.aggr

# aggregate across technical replicates for each filter. 
# same as above except with Sample instead of MasterID
filters.aggr <- ddply(s, ~ County + Location + Site + Date + SampleType + Sample,
                            summarize, 
                      LogCopiesPerL = round(mean(logCopies, na.rm = T), digits = 3), 
                      upperCI = round(upperConfInf(logCopies), digits = 3), 
                      lowerCI = round(lowerConfInf(logCopies), digits = 3), 
                      StdErr = signif(se(logCopies), digits = 3), 
                      LogLOD = round(log10(mean(SampleLOD, na.rm = T)), digits = 3),
                      QualCall = partialQualAggr(QualCall), 
                      aQCreason = aggrQCreason(QCreason), 
                      PCRplate = tail(Plate, n = 1))

# clean up for tidy plots and sharing
filters.aggr$upperCI[which(filters.aggr$QualCall == "Negative")] <- NA
filters.aggr$lowerCI[which(filters.aggr$QualCall == "Negative")] <- NA
filters.aggr$StdErr[which(filters.aggr$QualCall == "Negative")] <- NA
filters.aggr$QualCall[which(filters.aggr$QualCall == "Negative")] <- "Negative, 1/2 LOD assigned" 
filters.aggr$QualCall <- as.factor(filters.aggr$QualCall)
filters.aggr$aQCreason <- as.factor(filters.aggr$aQCreason)

head(filters.aggr)
filters.aggr <- filters.aggr[, -14]
filters.aggr

# aggregate by RNA_ID and plate to compare reruns
RNA_ID.aggr <- ddply(s, ~ County + Location + Site + Date + SampleType + Sample + MasterID + RNA_ID + RNAwell + RNAplate, 
                     summarize, 
                     LogCopiesPerL = mean(logCopies, na.rm = T),  
                     upperCI = upperConfInf(logCopies), 
                     lowerCI = lowerConfInf(logCopies), 
                     StdErr = se(logCopies),
                     LogLOD = round(log10(mean(SampleLOD, na.rm = T)), digits = 3),
                     aQualCall = partialQualAggr(QualCall), 
                     aQCreason = aggrQCreason(QCreason), 
                     PCRplate = max(as.numeric(Plate)))
head(RNA_ID.aggr) 

# aggregated by target
samples <- subset(data, SampleType != "Control")
samples.by.target <- dcast(samples, County + Location + Site + Date + SampleType + MasterID ~ Target, 
                           mean, value.var = c("logCopies"), na.rm = T)
head(samples.by.target)


# aggregate data for EBs and FBs 
ef <- subset(data, Target != "RP" & (ControlType == "EB" | ControlType == "FB"))
blanks.agg <- ddply(ef, ~ Date + Plate + ControlType + Sample + NTCmax + EBmax, summarize, maxPositives = max(Positives, na.rm = T), aQualCall = partialQualAggr(QualCall), aQCreason = aggrQCreason(QCreason))
blanks.agg

# #### clean up workspace ####
rm("bad", 'ebMax', 'ebs', 'ef', "fbMax", 'fbs', 'matchfail', 'mismatch',
   "s", "splitCol","w", "badData", 'd', 'ebPCRplate', 'fbPCRplate',
    'goodPlate', 'i', 'j','m',"M", "QualColumn", "ShiedlVol", "uLRxn")

#### Export results ####

# write.csv(data, str_c("Layton_EHP22_ddPCRdata_QCd", Sys.Date(),".csv"))
# write.csv(samples.aggr, str_c("Layton_EHP22_ddPCR_aggregated", Sys.Date(),".csv"))

