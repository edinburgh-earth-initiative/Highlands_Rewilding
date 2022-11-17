#Analysing Airborne Lidar data using lidR package 
#Author: Lucy Wells
#Date: 20/10/2022

###########Load required packages and set file paths######### 
#install.packages("lidR")
#install.packages("gstat")
library(lidR)
library(ggplot2)
library(gstat)
library(beepr)
library(raster)
library(sf)
library(future)
library(tidyr)
library(geometry)
get_lidr_threads() #see here for info on threads https://rdrr.io/cran/lidR/man/set_lidr_threads.html
set_lidr_threads(0)

path_to_data <- file.path("C:/Users/s2126368/OneDrive - University of Edinburgh/Bunloit/lidar") #change to pathway where data are stored on computer

##########Read lidar files. These are examples from a few locations in Scotland to have a look through#########

las_NH2746 <- readLAS(file.path(path_to_data, "NH2746_2PPM_LAS_PHASE1.laz"), select = "xyzinrc") #just north of Bunloit, near Meall na Faire. Mostly heather 
#^ use "select" to only load the attributes of interest and save memory. i = intensity, n = number of returns, r = return number, c = classification
summary(las_NH2746)#print summary of content e.g extent, location
las_check(las_NH2746) #validate data before analyzing

las_NM8428 <- readLAS(file.path(path_to_data, "NM8428_2PPM_LAS_PHASE1.laz"), select = "xyzinrc") #Just south of Oban, near Gallanach, mixture of improved grassland, broadleaved woodland, heather grassland, coniferous woodland
summary(las_NM8428)
las_check(las_NM8428)

las_NN55SW <- readLAS(file.path(path_to_data, "NN55SW_2PPM_LAS_PHASE2.laz"), select = "xyzinrc") #Black Wood of Rannoch, south side of Loch Rannoch, mixture of ancient woodland, heather, calcareous grassland
summary(las_NN55SW)
las_check(las_NN55SW)

las_NN55NW <- readLAS(file.path(path_to_data, "NN55NW_2PPM_LAS_PHASE2.laz"), select = "xyzinrc") #As above, adjoining tile
summary(las_NN55NW)
las_check(las_NN55NW)
plot(las_NN55NW)

#######Plot files#######
plot(las_NM8428) #basic plot just to have a look

plot(las_NM8428, color = "Intensity", bg = "white", axis = TRUE, legend = TRUE, breaks = "quantile") #Colour by attribute and add legend etc to improve display

######Plot cross section to look at vertical height######
#First create a function to combine the two steps of creating a cross-section and plotting it
plot_crossection <- function(las_NM8428,
                             p1 = c(min(las_NM8428@data$X), mean(las_NM8428@data$Y)),
                             p2 = c(max(las_NM8428@data$X), mean(las_NM8428@data$Y)),
                             width = 4, colour_by = NULL)
{
  colour_by <- enquo(colour_by)
  data_clip <- clip_transect(las_NM8428, p1, p2, width)
  p <- ggplot(data_clip@data, aes(X,Z)) + geom_point(size = 0.5) + coord_equal() + theme_minimal()
  
  if (!is.null(colour_by))
    p <- p + aes(color = !!colour_by) + labs(color = "")
  
  return(p)
}

#then use the function
plot_crossection(las_NM8428, colour_by = factor(Classification)) #Classification information here: https://desktop.arcgis.com/en/arcmap/latest/manage-data/las-dataset/lidar-point-classification.htm 

#^Note on classification. Classification of points as 'ground'/'not ground' can be done using a number of different algorithms, and its difficult to choose which one to use. 
#This R package provides a number of different methods, but advises that its best to use any classifications that are given by the data provider if possible, which is what I've done in these examples.

######Making a DTM#########
#DTMs (digital terrain models) represent the 'ground', and are used as the starting point for lots of analyses and 'normalizing' the data
#This package contains several different algorithms for creating them, which have different trade offs.

dtm_idw_NN55SW <- rasterize_terrain(las_NN55SW, algorithm = knnidw(k = 10L, p = 2));beep() #try using the invert distance weighting as test (medium in terms of edge effects, but has big steps in the dtm but so may not be suitable for bunloit considering the terrain of plateaus down to the loch)
plot_dtm3d(dtm_idw_NN55SW, bg = "white")

#######Use the DTM to normalize the point cloud#########
#Using Subtraction of the dtm
nlas_NN55SW_sub <- las_NN55SW - dtm_idw_NN55SW #this method is less accurate but simpler and less computationally heavy
plot(nlas_NN55SW_sub, bg = "white")

#Normalizing the point cloud through interpolation
nlas_NN55SW_int <- normalize_height(las_NN55NW, knnidw()) #more accurate method, but more computationally heavy and also suffers from edge effects

#######Make the canopy height model#########
#CHM is made using the normalized dtm, whereas the dsm is made without the normalized layer, and so represents the height of the canopy above sea level

#point-to-raster method
chm_ptr_NN55SW <- rasterize_canopy(las_NN55SW, res = 1, algorithm = p2r())
col <- height.colors(25)
plot(chm_ptr_NN55SW, col = col)

#Triangulation method
chm_tri_NN55SW <- rasterize_canopy(las_NN55SW, res = 1, algorithm = dsmtin(max_edge = 8));beep() #max_edge removes things like lochs, becuase it won't make triangles longer than 8
plot(chm_tri_NN55SW, col = col)

#############Derived metrics###########
#Derived metrics are scalar summaries of point distributions, that can be at different levels. For example, the standard deviation of point heights within a single tree crown is a metric calculated at the tree level. 
#Metrics can serve as a proxy for forest attributes, so are usually calculated based on the point heights (Z coordinates). They can be at the point cloud level, pixel level, tree crown level, inventory plot level

#Here, I have calculated metrics at the pixel level. This is because it can more easily link to the area-based approach which better links lidar data with field measured references to look at thinks like forest basal area
#The lidR package calculates the most commonly used metrics using the .stdmetrics function. The list of all metrics in this function can be found here: https://github.com/r-lidar/lidR/wiki/stdmetrics
#This produces rasters, which can easily be converted to a data frame to be used in further statistical analysis

metrics <- pixel_metrics(las_NN55NW, .stdmetrics, 10);beep() #calculate for pixel at 10 metres (can easily change)
plot(metrics, col = height.colors(50))
df <- as.data.frame(metrics, xy = TRUE) #convert spatraster object to a data frame
head(df)
write.csv(df, file = file.path(path_to_data, "df.csv"), row.names = FALSE)

######Using LAScatalog###############
#Data collected across the estates may be split into a number of point cloud files. The LAScatalog processing engine allows you to manage a number of files, rather than looping through each individual file
#In this example I have used the data sets from around the Black Wood of Rannoch: NN55SW, NN55NW, NN45NE, NN46SE, NN55SE, NN55NE, NN56SW

ctg <- readLAScatalog(file.path(path_to_data, "catalog"))
ctg
plot(ctg)
las_check(ctg)


##Here's the code for loading in and clipping the catalog to the habitat shapefiles for Bunloit##
poly <- shapefile("C:/Users/s2126368/OneDrive - University of Edinburgh/Bunloit/Shp files for Lucy/NVC.shp") #input the shapefile
opt_output_files(ctg) <- "C:/Users/s2126368/OneDrive - University of Edinburgh/Bunloit/lidar/plot_a_{ID}" #location you'd like the clipped LAS files to be saved 
subset <- clip_roi(ctg, poly) #clip to habitat shapefiles
##

#Dummy shapefiles for across the blackwood of rannoch here
poly_rannoch <- shapefile("C:/Users/s2126368/OneDrive - University of Edinburgh/Bunloit/Shp files for Lucy/rannoch_shapefile.shp")
plot(poly_rannoch)
opt_output_files(ctg) <- "C:/Users/s2126368/OneDrive - University of Edinburgh/Bunloit/lidar/rannoch_catalog/plot_a_{id}"

rannoch_subset <- clip_roi(ctg, poly_rannoch);beep()
plot(rannoch_subset)

rannoch_ctg <- readLAScatalog(file.path(path_to_data, "rannoch_catalog"))#re-read in catalog to do processing
opt_select(rannoch_ctg) <- "xyzinrc" #select parameters of interest (as above for the single file)

####Set up for parallel processing of the catalog - Doesn't seem to be working correclty at the moment (Stops DSM from calculating correctly)
#opt_independent_files(rannoch_ctg) <- TRUE #configures the processing for the case of independent files (here, habitat patches). Turns the engine intoa a simple loop rather than chunks covering the whole area wall-to-wall
#plan(multisession) #set up for parallel processing using the future package


#Making DTM using LAScatalog
dtm_rannoch <- rasterize_terrain(rannoch_ctg, 2, tin(), pkg = "terra", overwrite = TRUE) #resolution = 2, algorithm = tin()

#Use this DTM to normalize
opt_output_files(rannoch_ctg) <- paste0(tempdir(), "/{*}_norm")
rannoch_norm <- normalize_height(rannoch_ctg, dtm_rannoch)

#make CHM
rannoch_chm <- rasterize_canopy(rannoch_norm, res = 1, algorithm = p2r());beep()
writeRaster(rannoch_chm, filename = file.path(path_to_data, "rannoch_chm.tif"), overwrite = TRUE)

############Area Based Approach##################
rannoch_metrics <- grid_metrics(rannoch_ctg, .stdmetrics);beep() #wiki page here for definition of each metric https://github.com/r-lidar/lidR/wiki/stdmetrics 
plot(rannoch_metrics, col = height.colors(50))
df_rannoch <- as.data.frame(rannoch_metrics, xy = TRUE)
df_rannoch <- df_rannoch %>% drop_na()
head(df_rannoch)
 
write.csv(df_rannoch, file = file.path(path_to_data, "df_rannoch.csv"), row.names = FALSE)

############Non-standard metrics#################
###These can be used to give more complex information on woodland complexity. In particular 'entroy', which is a normalized Shannon vertical complexity index. This
#is a commonly used measure of woodland complexity. 
#These are not set up to run over LAScatalog files, so will need to do for each individual LAS file for each habitat
#swap these to rannoch_norm??

#Read in each individual LAS file representing a different habitat patch. Check them and plot them to have a look at each file
area_1 <- readLAS(file.path(path_to_data, "rannoch_catalog", "plot_a_1.las"), select = "xyzinrc")
summary(area_1)
las_check(area_1)
plot(area_1)

area_2 <- readLAS(file.path(path_to_data, "rannoch_catalog", "plot_a_2.las"), select = "xyzinrc")
summary(area_2)
las_check(area_2)
plot(area_2)

area_3 <- readLAS(file.path(path_to_data, "rannoch_catalog", "plot_a_3.las"), select = "xyzinrc")
summary(area_3)
las_check(area_3)
plot(area_3)

area_4 <- readLAS(file.path(path_to_data, "rannoch_catalog", "plot_a_4.las"), select = "xyzinrc")
summary(area_4)
las_check(area_4)
plot(area_4)

area_5 <- readLAS(file.path(path_to_data, "rannoch_catalog", "plot_a_5.las"), select = "xyzinrc")
summary(area_5)
las_check(area_5)
plot(area_5)

area_6 <- readLAS(file.path(path_to_data, "rannoch_catalog", "plot_a_6.las"), select = "xyzinrc")
summary(area_6)
las_check(area_6)
plot(area_6)

area_7 <- readLAS(file.path(path_to_data, "rannoch_catalog", "plot_a_7.las"), select = "xyzinrc")
summary(area_7)
las_check(area_7)
plot(area_7)

area_8 <- readLAS(file.path(path_to_data, "rannoch_catalog", "plot_a_8.las"), select = "xyzinrc")
summary(area_8)
las_check(area_8)
plot(area_8)

area_9 <- readLAS(file.path(path_to_data, "rannoch_catalog", "plot_a_9.las"), select = "xyzinrc")
summary(area_9)
las_check(area_9)
plot(area_9)

area_10 <- readLAS(file.path(path_to_data, "rannoch_catalog", "plot_a_10.las"), select = "xyzinrc")
summary(area_10)
las_check(area_10)
plot(area_10)

area_11 <- readLAS(file.path(path_to_data, "rannoch_catalog", "plot_a_11.las"), select = "xyzinrc")
summary(area_11)
las_check(area_11)
plot(area_11)

area_12 <- readLAS(file.path(path_to_data, "rannoch_catalog", "plot_a_12.las"), select = "xyzinrc")
summary(area_12)
las_check(area_12)
plot(area_12)

area_13 <- readLAS(file.path(path_to_data, "rannoch_catalog", "plot_a_13.las"), select = "xyzinrc")
summary(area_13)
las_check(area_13)
plot(area_13)

#Calculate the rumple_index. Computes the roughness of the surface. Higher Values mean a rougher surface

rumple_1 <- rumple_index(area_1$X, area_1$Y, area_1$Z)
rumple_2 <- rumple_index(area_2$X, area_2$Y, area_2$Z)
rumple_3 <- rumple_index(area_3$X, area_3$Y, area_3$Z)
rumple_4 <- rumple_index(area_4$X, area_4$Y, area_4$Z)
rumple_5 <- rumple_index(area_5$X, area_5$Y, area_5$Z)
rumple_6 <- rumple_index(area_6$X, area_6$Y, area_6$Z)
rumple_7 <- rumple_index(area_7$X, area_7$Y, area_7$Z)
rumple_8 <- rumple_index(area_8$X, area_8$Y, area_8$Z)
rumple_9 <- rumple_index(area_9$X, area_9$Y, area_9$Z)
rumple_10 <- rumple_index(area_10$X, area_10$Y, area_10$Z)
rumple_11 <- rumple_index(area_11$X, area_11$Y, area_11$Z)
rumple_12 <- rumple_index(area_12$X, area_12$Y, area_12$Z)
rumple_13 <- rumple_index(area_13$X, area_13$Y, area_13$Z)


#Entropy. Must give a value between 0 and 1. A normal distribution would give a value of 0.5, values closer to 1 have a more uniform distribution, values 
#closer to 0 indicate a unique possibility. 
entropy_1 <- entropy(area_1$Z, by = 1, zmax = NULL) #gives a number as an output
entropy_2 <- entropy(area_2$Z, by = 1, zmax = NULL) #gives a number as an output
entropy_3 <- entropy(area_3$Z, by = 1, zmax = NULL) #gives a number as an output
entropy_4 <- entropy(area_4$Z, by = 1, zmax = NULL) #gives a number as an output
entropy_5 <- entropy(area_5$Z, by = 1, zmax = NULL) #gives a number as an output
entropy_6 <- entropy(area_6$Z, by = 1, zmax = NULL) #gives a number as an output
entropy_7 <- entropy(area_7$Z, by = 1, zmax = NULL) #gives a number as an output
entropy_8 <- entropy(area_8$Z, by = 1, zmax = NULL) #gives a number as an output
entropy_9 <- entropy(area_9$Z, by = 1, zmax = NULL) #gives a number as an output
entropy_10 <- entropy(area_10$Z, by = 1, zmax = NULL) #gives a number as an output
entropy_11 <- entropy(area_11$Z, by = 1, zmax = NULL) #gives a number as an output
entropy_12 <- entropy(area_12$Z, by = 1, zmax = NULL) #gives a number as an output
entropy_13 <- entropy(area_13$Z, by = 1, zmax = NULL) #gives a number as an output


#put these into a data frame
Area <- c(1,2,3,4,5,6,7,8,9,10,11,12,13)
Habitat <- c("Forest/short veg mixture", "Forest/short veg mixture", "Short Vegetation", "Forest + river/gully", "Medium Vegetation", "River/Riparian", "House and Garden", "Short Vegetation", "Forest", "Forest", "Short Vegetation, isolated trees and power lines", "Short Vegetation", "Forest/Short veg mixture")
Complexity_Index <- c(entropy_1, entropy_2, entropy_3, entropy_4, entropy_5, entropy_6, entropy_7, entropy_8, entropy_9, entropy_10, entropy_11, entropy_12, entropy_13)
Rumple_Index <- c(rumple_1, rumple_2, rumple_3, rumple_4, rumple_5, rumple_6, rumple_7, rumple_8, rumple_9, rumple_10, rumple_11, rumple_12, rumple_13)
entropy_df <- data.frame(Area, Habitat, Complexity_Index, Rumple_Index)
entropy_df


#read in individual chm rasters for every habitat patch. Split using qgis
L <- list()
L$area_1_chm <- raster(file.path(path_to_data, "area_1_chm.tif"))
L$area_2_chm <- raster(file.path(path_to_data, "area_2_chm.tif"))
L$area_3_chm <- raster(file.path(path_to_data, "area_3_chm.tif"))
L$area_4_chm <- raster(file.path(path_to_data, "area_4_chm.tif"))
L$area_5_chm <- raster(file.path(path_to_data, "area_5_chm.tif"))
L$area_6_chm <- raster(file.path(path_to_data, "area_6_chm.tif"))
L$area_7_chm <- raster(file.path(path_to_data, "area_7_chm.tif"))
L$area_8_chm <- raster(file.path(path_to_data, "area_8_chm.tif"))
L$area_9_chm <- raster(file.path(path_to_data, "area_9_chm.tif"))
L$area_10_chm <- raster(file.path(path_to_data, "area_10_chm.tif"))
L$area_11_chm <- raster(file.path(path_to_data, "area_11_chm.tif"))
L$area_12_chm <- raster(file.path(path_to_data, "area_12_chm.tif"))
L$area_13_chm <- raster(file.path(path_to_data, "area_13_chm.tif"))

#Can turn each of these rasters into a data frame to do some stats on. For example here the standard deviation. This could give some information on canopy complexity
#e.g. a higher standard deviation suggests more variation in the canopy height
df_1 <- as.data.frame(L$area_1_chm, xy = TRUE) 
head(df_1)
df_1 <- df_1 %>% drop_na()
head(df_1)
sd(df_1$area_1_chm)
