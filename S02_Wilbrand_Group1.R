#############################################################################
# MSc Earth Observation Assignment 02
# Robert Wilbrand - Group 1
#############################################################################

#############################################################################
# 0) Setup
#############################################################################

# Load necessary libraries
library(raster)
library(tidyverse)

# Change according to your own working directory
setwd("C:/MSc GCG/2021_SoSe/EO/gcg_eo_s02_sub")

#------ Definitions of custom functions -------------------------------------

# Helper-function to convert degrees to radians
deg2rad <- function(deg){ (deg * pi) / (180) }

#############################################################################
# 1.1)    
#############################################################################

# Read tifs, drop QA layer
l1tp_subset <-
  dir(path = "./DN",
      recursive = T,
      pattern = ".tif$",
      ignore.case = T,
      full.names = T) %>% 
  stack %>% 
  dropLayer(7)

# Plot each of the layers
map(1:6,  ~plot(l1tp_subset[[.x]]))

# Read metadata file
mtl_metadata <-
  dir("./DN",
      recursive = T,
      pattern = "MTL.txt$",
      full.names = T) %>% 
  read.delim(sep = "=",
             stringsAsFactors = F)

# Extract conversion values from metadata
# I prefer using dplyr operations to do this
REFLECTANCE_MULT_BAND <-
  mtl_metadata %>%
  filter(str_detect(GROUP, "REFLECTANCE_MULT_BAND")) %>%
  slice(2:7) %>%
  pull(2) %>%
  as.numeric

REFLECTANCE_ADD_BAND <-
  mtl_metadata %>%
  filter(str_detect(GROUP, "REFLECTANCE_ADD_BAND")) %>%
  slice(2:7) %>%
  pull(2) %>%
  as.numeric

SUN_ELEVATION <-
  mtl_metadata %>%
  filter(str_detect(GROUP, "SUN_ELEVATION")) %>%
  pull(2) %>%
  as.numeric %>% 
  deg2rad

### Copy-pasted instructions from script using grep() -----------------------
# 
# REFLECTANCE_MULT_BAND <-
#   as.numeric(mtl_metadata[grep('REFLECTANCE_MULT_BAND',
#                                mtl_metadata$GROUP),][2:7,2])
# REFLECTANCE_ADD_BAND <-
#   as.numeric(mtl_metadata[grep('REFLECTANCE_ADD_BAND',
#                                mtl_metadata$GROUP),][2:7,2])
# SUN_ELEVATION <-
#   as.numeric(mtl_metadata[grep('SUN_ELEVATION',
#                                mtl_metadata$GROUP),][2])
# ---------------------------------------------------------------------------

# Convert DNs to sun-angle corrected TOA reflectance
l1tp_converted <-
  ((REFLECTANCE_MULT_BAND*l1tp_subset +
     REFLECTANCE_ADD_BAND)/sin(SUN_ELEVATION)) %>% 
  magrittr::multiply_by(10000) %>% 
  round

# Write result to disk
writeRaster(l1tp_converted, "L1TP_Reflectance2.tif", format = "GTiff",
            datatype = 'INT2S')

#############################################################################
# 1.2)  
#############################################################################

# a)
# Forest:
# The red edge and peak in the nIR are visible for both datasets, but only the
# BOA data shows a green peak in the visible spectrum. TOA data in the visible
# part of the spectrum drops from blue to red

# Open soil:
# Both data sources show an increase from blue to nIR and a decline afterwards.
# The visible part of the spectrum is very flat for the TOA data.

# Water: For TOA, has the highest value in the blue band and drops off with
# each subsequent band. Similar for BOA, but with the difference that it peaks
# in the green band.

# b)
# A major difference is that TOA data shows consistently higher values in the
# blue band relative to other bands when compared with the BOA data. This is
# to be expected as a consequence of Rayleigh scattering.

#############################################################################
# 2.1)    
#############################################################################

# Load QA layer
l1tp_qa <-
  dir(path = "./DN",
      recursive = T,
      pattern = "BQA.tif$",
      ignore.case = T,
      full.names = T) %>% 
  raster

# Extract and print 3 most frequent values
(
mostFrequentValues <-
    freq(l1tp_qa) %>%
    as_tibble %>%
    arrange(desc(count)) %>%
    slice(1:3)
)

# Helper function for 16-bit values
intTo16Bits <- function(n) as.numeric(intToBits(n)[1:16])

# Display the bit values with corresponding bit numbers
map_dfc(mostFrequentValues$value, ~rev(intTo16Bits(.x))) %>% t %>%
  as.data.frame %>%
  magrittr::set_colnames(as.hexmode(15:0))

# What do the most frequent values mean?
# Decode each integer value into bits and describe their meaning here: 

# Most frequent value:
# -> low confidence each for cloud, cloud shadow, snow/ice and cirrus

# Second most frequent value:
# -> cloud, high cloud confidence, other confidences low as above

# Third most frequent value:
# -> high cloud shadow confidence, other confidences low as above

#############################################################################
# 2.2)    
#############################################################################

# Bit testing for high confidence cloud / cloud shadow mask
highConfFilter <- function(n){
  if (is.null(n) | is.na(n)) return(NA)
  pixel <- intTo16Bits(n) %>% as.logical
  if (pixel[1]) return(T)
  if (pixel[6] & pixel[7]) return(T)
  if (pixel[8] & pixel[9]) return(T)
  return(F) # only gets returned when previous checks fail
}

# Bit testing for medium OR high confidence requires only 3 bits
mediumHighConfFilter <- function(n){
  if (is.null(n) | is.na(n)) return(NA)
  pixel <- intTo16Bits(n) %>% as.logical
  if (pixel[1] | pixel[7] | pixel[9]) return(T)
  return(F)
}

# Create masks from the defined functions
highConfMask        <- calc(l1tp_qa, highConfFilter)
mediumHighConfMask  <- calc(l1tp_qa, mediumHighConfFilter)

writeRaster(highConfMask, "HighConfidenceMask.tif", format = "GTiff",
            datatype = "LOG1S")
# -> LOG1S would be sufficient, but is not available in GDAL
writeRaster(mediumHighConfMask, "MediumConfidenceMask.tif", format = "GTiff",
            datatype = "INT1U")

plot(highConfMask)
plot(mediumHighConfMask)
# Plot the difference for better visual interpretation
plot(mediumHighConfMask - highConfMask)

# Based on visually identifying and inspecting the largest patches of pixels
# that are unique to the medium confidence mask in QGIS, the high confidence
# mask appears more reliable

#############################################################################
# 2.3)    
#############################################################################

# Read in BOA rasters
l1tp_boa <-
  dir(path = "./SR",
      recursive = T,
      pattern = "band[2-7].tif$",
      ignore.case = T,
      full.names = T) %>% 
  stack

# Apply high confidence cloud mask
l1tp_boa_masked_high <-
  mask(l1tp_boa, highConfMask, maskvalue = 1, progress = "text")

# Apply medium confidence cloud mask
l1tp_boa_masked_medium <-
  mask(l1tp_boa, mediumHighConfMask, maskvalue = 1, progress = "text")

# Write to disk
writeRaster(l1tp_boa_masked_high, "BOA_masked_high.tif",
            format = "GTiff", datatype = 'INT2S')

writeRaster(l1tp_boa_masked_medium, "BOA_masked_medium.tif",
            format = "GTiff", datatype = 'INT2S')