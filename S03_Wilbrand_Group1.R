#############################################################################
# MSc Earth Observation Assignment 03
# Robert Wilbrand - Group 1
#############################################################################

#############################################################################
# 0) Setup
#############################################################################
# Start timing of script-----------------------------------------------------
starttime <- Sys.time()

# Load necessary libraries---------------------------------------------------
library(rgdal)
library(raster)
library(terra) # Supposedly an evolution of the raster package, trying it out
library(tidyverse)
library(magrittr) # for pipe-friendly aliases to common operations

# Set working directory------------------------------------------------------
setwd("C:/MSc GCG/2021_SoSe/EO/gcg_eo_s03")

#############################################################################
# 1) NDVI and EVI calculation
#############################################################################

# Read in data---------------------------------------------------------------
pathToMarch <- dir(recursive = T, pattern = "201403")
marchStack <- rast(pathToMarch) # -> terra's SpatRaster format
marchStack[marchStack == -32768] <- NA
names(marchStack) <- c("blue", "green", "red", "nir", "swir1", "swir2")

# Define correction factors--------------------------------------------------
landsatCorrectionFactors <- c(2.5, 6, 7.5, 10000)

# NDVI functions-------------------------------------------------------------
normDiff <- function(a, b) (a - b) / (a + b)
NDVI <- function(ras, nIR, red){
  subset(ras, c(nIR, red)) %>%
    lapp(normDiff) # same as overlay() in raster package
}

# EVI function---------------------------------------------------------------
EVI <- function(ras, nIR, red, blue, correctionFactors){
  cf <- correctionFactors
  eviFunc <- function(n, r, b){
    # helper function that takes raw numbers to iterate over raster
    cf[1]*((n - r) / (n + cf[2]*r - cf[3]*b + cf[4]))
  }
  subset(ras, c(nIR, red, blue)) %>% 
    lapp(eviFunc) # same as overlay() in raster package
}

# Calculate EVI and NDVI for March image-------------------------------------
eviMarch <- EVI(marchStack, 4, 3, 1, landsatCorrectionFactors) %>% 
  clamp(-1, 1, values = F) %>%
  multiply_by(10000)
ndviMarch <- NDVI(marchStack, 4, 3) %>% 
  clamp(-1, 1, values = F) %>%
  multiply_by(10000)

# Write rasters to disk------------------------------------------------------
if (!dir.exists("spectral_indices")) dir.create("spectral_indices")
if (!file.exists("spectral_indices/EVI_March2014.tif")){
  writeRaster(eviMarch, "spectral_indices/EVI_March2014.tif",
              filetype = "GTiff", datatype = 'INT2S')}
if (!file.exists("spectral_indices/NDVI_March2014.tif")){
  writeRaster(ndviMarch, "spectral_indices/NDVI_March2014.tif",
              filetype = "GTiff", datatype = 'INT2S')}

#############################################################################
# 2) Tasseled Cap transformation
#############################################################################

# Define factors for Tasseled Cap transformation-----------------------------
tcc <- 
  cbind('brightness' = c( 0.2043,  0.4158,  0.5524, 0.5741,  0.3124,  0.2303),
        'greenness'  = c(-0.1603, -0.2819, -0.4934, 0.7940, -0.0002, -0.1446),
        'wetness'    = c( 0.0315,  0.2021,  0.3102, 0.1594, -0.6806, -0.6109))
rownames(tcc) <- c('blue', 'green', 'red', 'nIR', 'swIR1', 'swIR2')

# Set up helper function to multiply image stack by column in tcc------------
cwiseTransform <- function(n){ (marchStack*tcc[,n]) %>% app(sum) }

# Apply above function to each column----------------------------------------
# rast() turns list output into multilayer SpatRaster object
tcStack <- map(1:3, cwiseTransform) %>% rast
names(tcStack) <- c('Brightness', 'Greenness', 'Wetness')

# Write to disk--------------------------------------------------------------
if (!file.exists("spectral_indices/TasseledCap_March2014.tif")){
  writeRaster(tcStack, "spectral_indices/TasseledCap_March2014.tif",
              filetype = "GTiff", datatype = 'INT2S')}

# RGB plot of brightness - greenness - wetness-------------------------------
tcStack %>% stretch(minq = 0.01, maxq = 0.99) %>% plotRGB

#############################################################################
# 3) Explore  the training data
#############################################################################

# Read in training points----------------------------------------------------
trainConifer <- readOGR(getwd(), "Coniferous_training")

# Change SpatRaster to Raster for easier processing--------------------------
marchStack2 <- as(marchStack, "Raster")
names(marchStack2) <- names(marchStack)

# Extract raster values at training point locations--------------------------
trainValues <- raster::extract(marchStack2, trainConifer, sp = T)

# Convert to tibble and classID to factor------------------------------------
trainTibble <- 
  trainValues %>%
  as_tibble %>%
  mutate(across(classID, as.factor))

# Pivot to long format-------------------------------------------------------
trainLong <-
  trainTibble %>% 
  pivot_longer(blue:swir2,
               names_to = "band",
               values_to = "reflectance")

# Create boxplots of spectral bands per class--------------------------------
ggplot(trainLong, aes(band, reflectance, color=classID)) +
  geom_boxplot() +
  theme_bw()

# Convert raster to tibble, remove NAs, sample 100k points-------------------
set.seed(5678)
marchSample <-
  values(marchStack) %>% 
  as_tibble %>% 
  na.omit %>% 
  slice_sample(n = 1e5)
  
# Specify which bands to use for the x and y axis of the plot----------------
xband <- "red"
yband <- "nir"

# Create plot of band value density and training data------------------------
ggplot() + 
  geom_hex(marchSample, mapping = aes(get(xband), get(yband)), bins = 100) + 
  geom_point(trainTibble, mapping = aes(get(xband), get(yband),
                                        color=classID, shape=classID), 
             size = 2, inherit.aes = FALSE, alpha=1) + 
  
  scale_fill_gradientn(colours=c("black","white"), na.value=NA) + 
  scale_x_continuous(xband, limits=c(-10, quantile(marchSample[xband],
                                                   0.98, na.rm=T))) +
  scale_y_continuous(yband, limits=c(-10, quantile(marchSample[yband],
                                                   0.98, na.rm=T))) +
  scale_color_manual(values=c("red", "blue", "green", "purple")) +
  theme_bw()


# Print out total script execution time--------------------------------------
(elapsed <- Sys.time() - starttime)

# EOF