#############################################################################
# MSc Earth Observation Assignment 5
# Robert Wilbrand - Group 1
#############################################################################

#############################################################################
# 0) Setup
#############################################################################

# Load necessary libraries
library(randomForest)
library(rgdal)
library(terra)
library(raster)
library(tidyverse)

# Set working directory
setwd("C:/MSc GCG/2021_SoSe/EO")

# Set seed for reproducibility
set.seed(462)

# Start timing
tic <- round_date(Sys.time())
print(paste("Start of script execution:", tic))

#############################################################################
# 1) Data preparation
#############################################################################

# Import data
doySlices <- dir(recursive = T, pattern = "DOY")
trainData<- readOGR("Session05/gcg_eo_s05/training", "train_merge")

# Extract DOYs from filenames
DOYs <- doySlices %>% str_match("(?<=DOY)[:digit:]{3}") %>% as.character
singleNames <- paste("DOY", DOYs, sep = "_")

# Define band names
bandNames <- c("blue", "green", "red", "nIR", "swIR1", "swIR2")

# Input raster filenames, output analysis-ready tibble
pixelValuesFromTrainData <- function(r){
  # Load rasters
  rs <- stack(r)
  
  # Extract number and paths of source rasters
  rsSources <- sources(rast(rs))$source %>% unique
  n <- length(rsSources)
  
  # Extract DOY from sources
  rsDates <- str_match(rsSources, "(?<=DOY)[:digit:]{3}") %>%
    as.character 

  # select layer 1-6 for each source raster
  rs <- subset(rs, which(1:nlayers(rs)%%8 %in% 1:6))
  
  # Change layer names
  lyrNames <- paste0(rep(bandNames, n), "_", rep(rsDates, each = 6))
  names(rs) <- lyrNames
  
  # Extract values from rasters at training data locations
  pixelValues <- raster::extract(rs, trainData, sp = T) %>%
    as_tibble %>% 
    dplyr::select(c(classID, all_of(lyrNames))) %>%
    mutate(across(classID, as.factor)) %>% 
    na.omit
  
  return(pixelValues)
}

# Create dataframes for all single rasters
sliceSingles <- map(doySlices, pixelValuesFromTrainData)
names(sliceSingles) <- singleNames

# List all possible combinations of 2 rasters
twoOfNine <- expand.grid(1:9, 1:9) %>%
  filter(Var1 < Var2) %>% 
  arrange(Var1, Var2)

# Define names for 2-raster combinations
doubleNames <- map2(twoOfNine$Var1, twoOfNine$Var2,
                    ~paste("DOY", DOYs[.x], DOYs[.y], sep = "_"))

# Create dataframes for all 2-raster combinations
sliceDoubles <- map2(twoOfNine$Var1, twoOfNine$Var2,
                     ~pixelValuesFromTrainData(doySlices[c(.x, .y)]))
names(sliceDoubles) <- doubleNames

#############################################################################
# 2) Training the model
#############################################################################

# Generate random forest models
rfSingle <- map(sliceSingles, ~randomForest(classID ~ ., data = .x))
rfDouble <- map(sliceDoubles, ~randomForest(classID ~ ., data = .x))
names(rfSingle) <- singleNames
names(rfDouble) <- doubleNames

# Extract OOB error estimates for final tree
(OOB_One <- map_dbl(rfSingle, ~tail(.x$err.rate[,1], 1)))
(OOB_Two <- map_dbl(rfDouble, ~tail(.x$err.rate[,1], 1)))

# Extract best-performing single and double raster models
(bestSingle <- rfSingle[[which.min(OOB_One)]])
OOB_One[which.min(OOB_One)]
(bestDouble <- rfDouble[[which.min(OOB_Two)]])
OOB_Two[which.min(OOB_Two)]

# -> The best-performing model uses DOY 196 (July 15) and 258 (September 15)
# A few other combinations (e.g. DOY 74 and 196) come close

tail(bestDouble$err.rate, 1)

# -> The OOB error estimate is highest for class 1 (~32.2%)

# Plot OOB error for best single raster model
bestSingle$err.rate %>%
  as_tibble %>%
  rowid_to_column(var = "ID") %>%
  pivot_longer(-ID) %>%
  ggplot(aes(x = ID, y = value, group = name, color = name)) +
  geom_line() +
  theme_bw() +
  labs(x = "Tree", y = "OOB error estimate", color = "Class ID",
       title = names(rfSingle)[which.min(OOB_One)]) +
  theme(plot.title = element_text(hjust = 0.5))

# Plot OOB error for best double raster model
bestDouble$err.rate %>%
  as_tibble %>%
  rowid_to_column(var = "ID") %>%
  pivot_longer(-ID) %>%
  ggplot(aes(x = ID, y = value, group = name, color = name)) +
  geom_line() +
  theme_bw() +
  labs(x = "Tree", y = "OOB error estimate", color = "Class ID",
       title = names(rfDouble)[which.min(OOB_Two)]) +
  theme(plot.title = element_text(hjust = 0.5))

# -> The OOB error estimate stabilizes around 150 trees

#############################################################################
# 3) Train final model  
#############################################################################

# Train final model
bestModelName <- names(sliceDoubles)[which.min(OOB_Two)]
(finalModel <- randomForest(classID ~ ., data = sliceDoubles[[bestModelName]],
                           ntree = 150))

# Plot OOB error estimates for the final model
finalModel$err.rate %>%
  as_tibble %>%
  rowid_to_column(var = "ID") %>%
  pivot_longer(-ID) %>%
  ggplot(aes(x = ID, y = value, group = name, color = name)) +
  geom_line() +
  theme_bw() +
  labs(x = "Tree", y = "OOB error estimate", color = "Class ID",
       title = names(rfDouble)[which.min(OOB_Two)]) +
  theme(plot.title = element_text(hjust = 0.5))

# Investigate variable importance
varImpPlot(finalModel)

# -> The most important predictor is the near infrared band for July 15

# Plot partial dependence plots
par(mfrow = c(2, 2))
partialPlot(finalModel, as.data.frame(sliceDoubles[[bestModelName]]),
            nIR_196, "1", col="red",
            lwd=2, main="classID = 1 (deciduous)")
partialPlot(finalModel, as.data.frame(sliceDoubles[[bestModelName]]),
            nIR_196, "2", col="orange",
            lwd=2, main="classID = 2 (mixed)")
partialPlot(finalModel, as.data.frame(sliceDoubles[[bestModelName]]),
            nIR_196, "3", col="limegreen",
            lwd=2, main="classID = 3 (coniferous)")
partialPlot(finalModel, as.data.frame(sliceDoubles[[bestModelName]]),
            nIR_196, "4", col="turquoise",
            lwd=2, main="classID = 4 (non-forest)")

#############################################################################
# 4) Create prediction raster  
#############################################################################

# Load validation raster
validationRaster <- stack(doySlices[c(5,7)]) %>% 
  dropLayer(c(7,8,15,16))
names(validationRaster) <- paste0(rep(bandNames, 2),
                                  "_", rep(c("196", "258"), each = 6))

# Create prediction
predictionRaster <- predict(validationRaster, finalModel)

# Plot result
layout(1)
customPalette <- c("limegreen", "goldenrod", "darkgreen", "grey70")
plot(rast(predictionRaster), col = customPalette)

# Write to disk
writeRaster(predictionRaster, "Session05/randomForestPrediction.tif",
            format = "GTiff", datatype = 'INT1U')

# Run predict to store RF probabilities for class 1-4
rfp <- predict(validationRaster, finalModel, type = "prob", index=1:4)

# RGB plot of classification probabilities
# red = mixed, green = coniferous, blue = deciduous
plotRGB(rfp*255, 2, 3, 1)

# Scale probabilities to integer values 0-100 and write to disk
writeRaster(rfp*100, "Session05/randomForestProbs.tif",
            datatype = 'INT1U', overwrite = T)

# Return script execution time
toc <- round_date(Sys.time())
print(paste("Script execution time:", as.period(toc-tic)))


# EOF