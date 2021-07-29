#!/bin/R
#Usage: Rscript MultiPolywoggle.R <path to files> <particle.data prefix>

#Get arguements
args <- commandArgs(trailingOnly=TRUE)
PATH <- args[1]
PATTERN <- args[2]

#Libraries
suppressMessages(library(raster))
suppressMessages(library(rgdal))
suppressMessages(library(sf))
suppressMessages(library(spatialEco))
suppressMessages(library(dplyr))

#Read shapefile
polys <- readOGR(dsn="/home/afields/Workspace/snapper/support", layer = "GOM_Polygons_all")
polys <- spTransform(polys, CRS("+init=epsg:3160"))

#Get files list
files <- list.files(path=PATH,pattern=PATTERN)
files <- files[grep("intersect",files,invert=T)]

#Get Catch sites
dat.g <- read.table("/home/afields/Workspace/snapper/support/uniq_new_lats.txt", head=T)

print(c(PATH, PATTERN, files))

for(i in files){
print(paste("processing", i))
RELEASE <- i
if(!(file.exists(paste(PATH, RELEASE, sep="/")))){print(paste(i, "was not recognized by R")); next}
day <- matrix(unlist(strsplit(i,"_")))[3,1]

#Make release file
if(file.info(paste(PATH, RELEASE, sep="/"))$size==0){next}
release <- read.table(paste(PATH, RELEASE, sep="/"),head=F)
colnames(release) <- c("Start", "Particle", "Time", "Longitude", "Latitude", "Depth", "Distance", "Exit_code", "Release_date")
if(release$Longitude[1] > 0){release$Longitude <- release$Longitude -360}
release$Date <- rep(day,nrow(release))

#Prepare coordinates, data, and proj4string
coords <- release[ , c('Longitude','Latitude')]
dat   <- release[ , c("Start", "Particle", "Time", "Depth", "Distance", "Exit_code", "Release_date", "Date")]
crs    <- CRS("+init=epsg:3160")

#Make the SpatialPointsDataFrame object
pts <- SpatialPointsDataFrame(coords = coords, data= dat, proj4string = crs)

#Look for intersection
new_shape <- point.in.poly(pts, polys)

#Combine with location data
tmp.dat <- cbind(as.data.frame(pts@coords), new_shape@data)

cat.dat <- tmp.dat[tmp.dat$N %in% dat.g$N, ]

#Output the results
write.table(tmp.dat, paste(PATH,"/",i,".intersect",sep=""), col.names=T, row.names=F, quote=F, sep="\t")

if(exists("cat.dat")){
if(length(list.files(path=PATH,pattern="Particles_at_Catch.tab"))==0){write.table(cat.dat, paste(PATH,"Particles_at_Catch.tab",sep="/"), col.names=T, row.names=F, quote=F, sep="\t")
}else if(length(list.files(path=PATH,pattern="Particles_at_Catch.tab"))==1) {write.table(cat.dat, paste(PATH,"Particles_at_Catch.tab",sep="/"), col.names=F, row.names=F, quote=F, sep="\t", append=T)}
}
rm(cat.dat)
}
