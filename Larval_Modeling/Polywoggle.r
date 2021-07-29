#!/bin/R
#Usage: Rscript Polywoggle.R <particle.data> <Output.txt>

#Get arguements
args <- commandArgs(trailingOnly=TRUE)
RELEASE <- args[1]
OUT <- args[2]

#Libraries
suppressMessages(library(raster))
suppressMessages(library(rgdal))
suppressMessages(library(sf))
suppressMessages(library(spatialEco))

#Read shapefile
suppressMessages(polys <- readOGR(dsn=getwd(), layer = "GOM_Polygons_all"))
polys <- spTransform(polys, CRS("+init=epsg:3160"))

#Make release file
release <- read.table(RELEASE,head=F)
colnames(release) <- c("Start", "Particle", "Time", "Longitude", "Latitude", "Depth", "Distance", "Exit_code", "Release_date")
if(release$Longitude[1] > 0){release$Longitude <- release$Longitude -360}

# prepare coordinates, data, and proj4string
coords <- release[ , c('Longitude','Latitude')]
dat   <- release[ , c("Start", "Particle", "Time", "Depth", "Distance", "Exit_code", "Release_date")]
crs    <- CRS("+init=epsg:3160")

# make the SpatialPointsDataFrame object
pts <- SpatialPointsDataFrame(coords = coords, data= dat, proj4string = crs)

#Informing user of progress
print("Data loaded scuccessfully")
print("Comparing files")

#Look for intersection
new_shape <- point.in.poly(pts, polys)

#Combine with location data
tmp.dat <- cbind(as.data.frame(pts@coords), new_shape@data)

#Output the results
write.table(tmp.dat, OUT, col.names=T, row.names=F, quote=F, sep="\t")
