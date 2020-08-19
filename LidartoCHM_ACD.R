# GENERATION CANOPY HEIGHT MODEL and ABOVEGROUND CARBON DENSITY USING LIDR
# set up the directory
setwd("...")
# install Rtools
install.packages("sf") # for simple feature (vector data)
install.packages("rgl")
install.packages(c("lidR","rlas","raster","sp",'mapview','devtools'))

# Load packages
library(devtools)
devtools::install_github("gisma/uavRst", ref = "master")

library(rlas)         # read and Write 'las' and 'laz' Binary File Formats 
library(sp)           # represent spatial data
library(raster)       # analysis raster data
library(mapview)      # view map
library(rgl)          # plot 3D graphics
library(lidR)         # analysis lidar data 

# Piontcloudviewer
source("https://raw.githubusercontent.com/Jean-Romain/PointCloudViewer/master/sdl.R")
devtools::install_github("Jean-Romain/PointCloudViewer")


# TAPAJOS 2016
# TAPAJOS 2016 A02
# create a catalog on an interactive map
ctg1 <- catalog("./A2_A3_2016/TAP_A02_2016_LiDAR/TAP_A02_2016_laz")
#las2016a2 <- readLAS(ctg1)
#plot(las2016a2, backend = "pcv")

# generate the digital surface model(dsm)
thr <- c(0,2,5,10,15)
edg <- c(0, 1.5)
dsm1 <- grid_canopy(ctg1, 1, pitfree(thr, edg))
plot(dsm1)

# generate the digital terrain model(dtm)
dtm1 = lidR::grid_terrain(ctg1, res=1,  algorithm = lidR::knnidw(k=50, p=3))

# generate the canopy height model(chm)
chm1 = dsm1 - dtm1

# correct pixels artefacts with values <than 0
chm1[chm1 < 0]=NaN
plot(chm1)

#generate raster file
writeRaster(chm1,"chm_2016_a2.tif")


# TAPAJOS 2016 A03
# create a catalog on an interactive map
ctg2 <- catalog("./A2_A3_2016/TAP_A03_2016_LiDAR/TAP_A03_2016_laz")

# generate the digital surface model(dsm)
thr <- c(0,2,5,10,15)
edg <- c(0, 1.5)
dsm2 <- grid_canopy(ctg2, 1, pitfree(thr, edg))
plot(dsm2)

# generate the digital terrain model(dtm)
# spatial interpolation. Interpolation is done using a k-nearest neighbour (KNN) approach with an inverse-distance weighting (IDW).
dtm2 = lidR::grid_terrain(ctg2, res=1,  algorithm = lidR::knnidw(k=50, p=3)) 

# generate the canopy height model(chm)
chm2 = dsm2 - dtm2

# correct pixels artefacts with values <than 0
chm2[chm2 < 0]=NaN
plot(chm2)

#generate raster file
writeRaster(chm2,"chm_2016_a3.tif")

# merge two sites raster data in one
chm_2016 = merge(chm1, chm2, tolerance=0.05, filename="chm_2016.tif", overlap=TRUE, ext=NULL)
plot(chm_2016)



# TAPAJOS 2018
# TAPAJOS 2018 A02
# create a catalog on an interactive map
ctg3 <- catalog("./A2_A3_2018/TAP_A02_2018_LiDAR/TAP_A02_2018_LAS")
#las2018a2 <- readLAS(ctg3)
#plot(las2018a2, backend = "pcv")

# generate the digital surface model(dsm)
thr <- c(0,2,5,10,15)   # defult by the Khosravipour et al. (2014) algorithm
edg <- c(0, 1.5)        # max edge length, outer will be removed
dsm3 <- grid_canopy(ctg3, 1, pitfree(thr, edg))  #using pit free method to calculate the Canopy(dsm)
plot(dsm3)

# generate the digital terrain model(dtm)
dtm3 = lidR::grid_terrain(ctg3, res=1,  algorithm = lidR::knnidw(k=50, p=3))  #using knn IDW method to calculate the terrian(dtm)

# generate the canopy height model(chm)
chm3 = dsm3 - dtm3

# correct pixels artefacts with values <than 0
chm3[chm3 < 0]=NaN
plot(chm3)

#generate raster file
writeRaster(chm3, "chm_2018_a2.tif")

# TAPAJOS 2018 A03
# create a catalog on an interactive map
ctg4 <- catalog("./A2_A3_2018/TAP_A03_2018_LiDAR/TAP_A03_2018_LAS")

# generate the digital surface model(dsm)
thr <- c(0,2,5,10,15)
edg <- c(0, 1.5)
dsm4 <- grid_canopy(ctg4, 1, pitfree(thr, edg)) #(las, res, algorithm) #pitfree--This allows for virtual 'emulation' of the fact that a lidar point is not a point as such, but more a disc. This tweak densifies the point cloud and the resulting canopy model is smoother and contains fewer 'pits' and empty pixels.
plot(dsm4)

# generate the digital terrain model(dtm)
dtm4 = lidR::grid_terrain(ctg4, res=1,  algorithm = lidR::knnidw(k=50, p=3))

# generate the canopy height model(chm)
chm4 = dsm4 - dtm4

# correct pixels artefacts with values <than 0
chm4[chm4 < 0]=NaN
plot(chm4)

#generate raster file
writeRaster(chm4, "chm_2018_a3.tif")

#merge two sites raster data in one
chm_2018 = merge(chm3, chm4, tolerance=0.05, filename="chm_2018.tif", overlap=TRUE, ext=NULL)
plot(chm_2018)

# after resample data  grid size into 0.25ha = 2500m
chm_2016_50 <- aggregate(chm_2016,fact=50,fun=mean)
chm_2018_50 <- aggregate(chm_2018,fact=50,fun=mean)
writeRaster(chm_2016_50, "chm_2016_50b.tif")
writeRaster(chm_2018_50, "chm_2018_50b.tif")

# carbon stocks calculation
cellStats(chm_2016_50,mean)
plot(chm_2016_50)
crs(chm_2016_50)
# unite = kgC/m2
acd=0.025*(chm_2016_50 ^(1.99))
plot(acd)
writeRaster(acd,"acd_2016_50b.tif")

cellStats(chm_2018_50,mean)
acd=0.025*(chm_2018_50^(1.99))
plot(acd)
writeRaster(acd,"acd_2018_50b.tif")


