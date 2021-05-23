#############################################################################
# MSc Earth Observation Assignment 06
# Robert Wilbrand - Group 1
#############################################################################

#############################################################################
# 0) Setup
#############################################################################

# Load necessary libraries
library(sf)
library(rgdal)
library(raster)
library(tidyverse)

# Choose working directory interactively or change to path name (preferable
# for purposes of reproducibility)
setwd("C:/MSc GCG/2021_SoSe/EO")

#############################################################################
# 1) Producing reference data 
#############################################################################

forestMap <- raster("Session05/randomForestPrediction.tif")
validationSample <- sampleStratified(forestMap, 25, na.rm = T, sp = T)
writeOGR(validationSample, "Session06", "validationSample",
         driver = 'ESRI Shapefile')

#############################################################################
# 2) Area-adjusted accuracy assessment
#############################################################################

classCounts <- freq(forestMap, useNA = "no")
sumNonNAPixels <- sum(classCounts[,2])

(classProportion <- as_tibble(classCounts) %>% 
  mutate(proportion = count/sumNonNAPixels))

validationReassembled <- st_read("Session06/validationReassembled.shp")

(confMatrix <- validationReassembled %>%
    st_set_geometry(NULL) %>% 
    dplyr::select(rndmFrP, Validation) %>%
    table %>% 
    as.matrix)

# Which class has the highest / lowest user큦 accuracy?
# -> Deciduous has the lowest UA (48%), Non-forest is highest (88%)

# Which class has the highest / lowest producer큦 accuracy?
# -> Before area adjustment, PA is lowest for Mixed (~47.1%) and highest for Non-forest(~91.7%)
# -> After area adjustment, it is lowest for Coniferous (~58.2%) and still highest for non-forest (~93.7%)

# How does the overall accuracy differ after area-adjustment? Why?
# -> It goes up from 70% to almost 75%. Correct classifications no longer receive equal weight.

# How do the map-based area estimates differ from those obtained using the reference data?
# -> Mixed stays roughly the same, deciduous and coniferous go up, non-forest goes down
# (see and of document for a direct number comparison)

#############################################################################
# 3) Knowledge transfer
#############################################################################

## a. Generate the confusion matrix containing probabilities

(confMatrixProbs <- 
  pull(classProportion, proportion)*confMatrix/rowSums(confMatrix)) %>% 
  round(4)

(probRSums <- rowSums(confMatrixProbs)) %>% round(4) # rowwise sums
(probCSums <- colSums(confMatrixProbs)) %>% round(4) # columnwise sums

## b. Calculate overall accuracy and class-wise user큦 and producer큦 accuracy from the confusion matrix

# Before area adjustment

(OA_unadjusted <- 100*sum(diag(confMatrix))/sum(confMatrix))
(PA_unadjusted <- 100*diag(confMatrix)/colSums(confMatrix)) %>% round(2)
(UA_unadjusted <- 100*diag(confMatrix)/rowSums(confMatrix))

# After area adjustment

(OA_adjusted <- 100*sum(diag(confMatrixProbs))) %>% round(2)
(PA_adjusted <- 100*diag(confMatrixProbs)/probCSums) %>% round(2)
(UA_adjusted <- 100*diag(confMatrixProbs)/probRSums)

## c. Produce error-adjusted area estimates from the confusion matrix

(adjAreas <- probCSums*sumNonNAPixels*30^2/10000)
(adjAreasPct <- 100*adjAreas/sum(adjAreas)) %>% round(1)

bind_cols(mapPct = 100*pull(classProportion, proportion),
          refPct = adjAreasPct)

# -> All of the above values are identical to those in the Excel document

# EOF
