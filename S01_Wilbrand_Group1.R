#############################################################################
# MSc Earth Observation Assignment 01
# Robert Wilbrand
#############################################################################

#############################################################################
# 1)    
#############################################################################

# Load necessary libraries
library(raster)
library(tidyverse)

# Change according to your own working directory
setwd("C:/MSc GCG/2021_SoSe/EO/gcg_eo_s01")

#############################################################################
# 2)    
#############################################################################

# Read and stack rasters from both subfolders
tif_list <- dir(recursive = T, pattern = "B[2-7].[Tt][Ii][Ff]$") 
march_stack <- stack(tif_list[1:6])
july_stack <- stack(tif_list[7:12])

#############################################################################
# 3)    
#############################################################################

projection(march_stack)
projection(july_stack)
# Both rasters are in the UTM Zone 34 projection

extent(march_stack)
extent(july_stack)

common_extent <- raster::intersect(extent(march_stack),
                                   extent(july_stack))

# Set custom extent and crop
roi.extent <- c(327945, 380325, 5472105, 5521095) %>% extent

march_cropped <- crop(march_stack, roi.extent)
july_cropped <- crop(july_stack, roi.extent)

#############################################################################
# 4)    
#############################################################################

# Write rasters to disk
writeRaster(march_cropped, "March_cropped.tif", format = "GTiff")
writeRaster(july_cropped, "July_cropped.tif", format = "GTiff")

# GTiff is a good default for rasters with a geographical interpretation,
# because it's compatible with all default GIS programs and other software

#############################################################################
# 5)    
#############################################################################

# True color plots
plotRGB(march_cropped, 3, 2, 1, scale = 65535, stretch = 'hist')
plotRGB(july_cropped, 3, 2, 1, scale = 65535, stretch = 'hist')

# False color plots
plotRGB(march_cropped, 4, 3, 2, scale = 65535, stretch = 'hist')
plotRGB(july_cropped, 4, 3, 2, scale = 65535, stretch = 'hist')

#############################################################################
# 6)    
#############################################################################

# Define coordinate (units = m)
coordinate <- data.frame('x' = 355623, 'y' = 5486216) 

# Extract cell values from image_stack
image_stack <- stack(march_cropped, july_cropped)
cellvals <- raster::extract(image_stack, coordinate)

#############################################################################
# 7)    
#############################################################################

# Create dataframe for plotting
profile <- data.frame('wavelength'=rep(c(482, 561, 655, 865, 1609, 2201),2), 
                      'date'=as.factor(c(rep('10 March', 6),
                                         rep('10 July', 6))),
                      'values'=as.numeric(cellvals))
print(profile)

# Plot results
ggplot(profile, aes(x=wavelength, y=values, group=date)) + 
  geom_line(aes(color=date)) + 
  # Add axis & legend names
  scale_y_continuous(name='DN') + 
  scale_x_continuous(name='Wavelength (nm)') + 
  scale_colour_manual(name='Acquisition DOY',
                      values=c('darkgreen', 'brown')) +
  theme_bw() # Chose black/white layout

# The spectral profile could be indicative of a forest area with no (or few)
# leaves in March (maybe picking up on ground vegetation with no interference
# from leaves), but which is green in July. An agricultural field is also
# conceivable if there is early vegetation in March, since the profile is not
# typical for bare soil