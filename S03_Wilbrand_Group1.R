#############################################################################
# MSc Earth Observation Assignment 03
# Robert Wilbrand - Group 1
#############################################################################

#############################################################################
# 0) Setup
#############################################################################
starttime <- Sys.time()

# Load necessary libraries---------------------------------------------------
library(raster)
library(terra) # Supposedly an evolution of the raster package, trying it out
library(tidyverse)

# Set working directory------------------------------------------------------
setwd("C:/MSc GCG/2021_SoSe/EO/gcg_eo_s03")

#############################################################################
# 1) NDVI and EVI calculation
#############################################################################

# Read in data---------------------------------------------------------------
pathToMarch <- dir(recursive = T, pattern = "201403")
marchStack <- rast(pathToMarch) # -> terra's SpatRaster format
marchStack[marchStack < 0] <- NA

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
eviMarch <- round(10000*EVI(marchStack, 4, 3, 1, landsatCorrectionFactors))
ndviMarch <- round(10000*NDVI(marchStack, 4, 3))

# Write rasters to disk------------------------------------------------------
if (!dir.exists("spectral_indices")) dir.create("spectral_indices")
if (!file.exists("spectral_indices/EVI_March2014.tif")){
  writeRaster(eviMarch, "spectral_indices/EVI_March2014.tif",
              filetype = "GTiff", datatype = 'INT2S')}
if (!file.exists("spectral_indices/EVI_March2014.tif")){
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
# 3)    
#############################################################################

(elapsed <- Sys.time() - starttime)

# EOF


