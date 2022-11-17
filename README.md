# Highlands_Rewilding
Script for Highlands Rewilding using the lidR pacakge to assess woodland complexity 


I have used ALS data collected by the Scottish Government, which was commissioned for use for
flood risk and management, downloaded from here as LAS files. As such, the data was collected
mostly along rivers and lochs across Scotland. I have selected tiles that are from within the
Highlands, and include a mixture of habitats, including both deciduous and coniferous woodland.

For most of the analysis I have used the example of the Black Wood of Rannoch, which has a
mixture of habitats. These include some ancient Caledonian pinewood, areas that have short veg-
etation, and areas around houses and gardens. I have tiled these together in order to demonstrate
a workflow with multiple different LAS files, and selected a shapefile of habitat patches across the
area to imitate the shapefiles you have provided for Bunloit.

The software I have used is the package lidR in R, which is specifically for reading and analysing
lidar data in the form of LAS and LAZ files for forestry applications.
The github page for the package can be found here https://github.com/r-lidar/lidR and the book that I worked through to
create the script is here https://r-lidar.github.io/lidRbook/index.html
The overall aim of the script is to extract metrics describing the woodland structure, depending
on the habitat/location in the estate. The aim was to find a relatively simple way of describing
variation in structure, and as such the values are aimed at being used for comparison between
3 different areas, and potentially over time if repeat lidar surveys are carried out. The assumption
is that higher ’structural complexity’ = higher biodiversity, but this is not an exact measurement
in itself.
The script runs through loading LAS files into the package, and basic plotting to have a look
at the data visually. It then explains how to make DTMs (Digital Terrain Models) which are
used to normalize the point cloud and then create a canopy height model. More importantly, the
script then uses the pixel metrics function to calculate derived metrics at the pixel. This creates
a number of rasters which provide structural information including the standard deviation, mean
and height distribution of the point cloud in each pixel. The full list of calculated metrics can be
found on this wiki. These provide some basic but useful summary information, in particular the
standard deviation is often used in interpreting complexity.

Following on from this, I have provided an example of the same workflow, using the LAScatalog
processing engine, which is a more effective way of analysing multiple tiled LAS files. This was
done under the assumption that data collected over the Bunloit (and Beldorney) estates is likely
to be split into multiple files. Using this method would have the advantage of reducing edge effects
compared to looping over multiple data files in R to undertake the same calculations.
Finally, the script calculates some of the structural metrics described above. They are calculated
for each dummy ’habitat patch’ at the Blackwood of Rannoch site. These are placed in a data
frame so that they can be compared, and exported for further statistical analysis if needed.




