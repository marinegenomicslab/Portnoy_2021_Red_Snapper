#Loading Libraries
library('devtools')
library('pegas')
library('adegenet')
library('vcfR') 
library('rospca')
library('dartR')
library('zvau')
library('geosphere')
library('stringr')
library('ggmap')
library('ggcompoplot')
library('vegan')
library('spdep')
library('adespatial')
library('igraph')
library('poppr') 
library('smatr')
library('radiator')
library('related')
library('ggcompoplot')
library('dartR')
library('scales')
library('akima')
library('spacemakeR')
library('hierfstat')
library('ggpubr')
library('PopGenReport')
library('raster')
library('rgdal')
library('sf')
library('spatialEco')

#Loading data
load("gen.net.gz")
load("Larval_conn_prob.gz")
Gulf_map <- ggmap(get_stamenmap(bbox = c(left = -100, bottom = 18, right =-74, top = 38.3), color="bw", force=T, zoom = 8, maptype=c("terrain-background")))

#Getting the Gulf samples
keep.ind <- gen.net@strata$INDV[gen.net2@strata$Region!="Atlantic"]
gen.Gulf <- gen.net[keep.ind, ]
setPop(gen.Gulf) <- ~ SubPOP
X.Gulf <- scaleGen(gen.Gulf, NA.method="mean", scale=F)

#Loading color schemes

#Getting the distance matrix
dist.mat <- distm(gen.Gulf@other$lat, fun=distHaversine)/1000
colnames(dist.mat) <- rownames(dist.mat) <- indNames(gen.Gulf)

#Listing the pairs to put points which need imaginary points between them
tmp.pairs <- matrix(c("FL_10.I2.Lib8","MXCII_30.I4.Lib9","CH_039.I1.Lib5","MXCII_31.I2.Lib9","MXCII_01.I7.Lib9","MXCII_06.I10.Lib9","MXCII_05.I7.Lib9","MXCII_08.I2.Lib9",
"MXCII_09.I4.Lib9","MXCII_10.I2.Lib9","MXCII_14.I4.Lib9","MXCII_32.I10.Lib9","MXCII_16.I4.Lib9","MXCII_26.I10.Lib9","MXCII_20.I7.Lib9","MXCII_31.I2.Lib9","MXCII_29.I10.Lib9",
"MXCII_26.I10.Lib9"),byrow=T,ncol=2)

#Getting the imaginary point locations
npoints <- data.frame(matrix(ncol=2))
colnames(npoints) <- c("Lon", "Lat")
for(i in 1:nrow(tmp.pairs)){
tmp.dist <- dist.mat[tmp.pairs[i,1],tmp.pairs[i,2]]
new.pts <- ceiling(tmp.dist/144)

Lon.diff <- abs(diff(as.numeric(as.matrix(gen.Gulf@strata$Lon[gen.Gulf@strata$INDV %in% tmp.pairs[i,]]))))
Lat.diff <- abs(diff(as.numeric(as.matrix(gen.Gulf@strata$Lat[gen.Gulf@strata$INDV %in% tmp.pairs[i,]]))))

Lon.start <- min(as.numeric(as.matrix(gen.Gulf@strata$Lon[gen.Gulf@strata$INDV %in% tmp.pairs[i,]])))
Lat.start <- min(as.numeric(as.matrix(gen.Gulf@strata$Lat[gen.Gulf@strata$INDV %in% tmp.pairs[i,]])))

for(j in 1:(new.pts-1)){
new.val <- c(Lon.start + (Lon.diff/(new.pts)) * j, Lat.start + (Lat.diff/(new.pts)) * j)
npoints <- rbind(npoints,new.val)
}
}

npoints <- npoints[-1,]
rownames(npoints) <- paste("Pt",seq(1:nrow(npoints)),sep="_")

# Getting Grids for new points to connect network
#Read shapefile
polys <- readOGR(dsn=getwd(), layer = "GOM_Polygons_all")
polys <- spTransform(polys, CRS("+init=epsg:3160"))

#Prepare coordinates, data, and proj4string
coords <- npoints
for(i in c('Lon','Lat')){coords[,i] <- as.numeric(as.matrix(coords[,i]))}
crs <- CRS("+init=epsg:3160")

# make the SpatialPointsDataFrame object
pts <- SpatialPointsDataFrame(coords = coords, data= npoints, proj4string = crs)

#Look for intersection
new_shape <- point.in.poly(pts, polys)

#View Results
npoints <- new_shape@data
npoints

#Getting Grids for Catch locations
#Prepare coordinates, data, and proj4string
coords <- gen.Gulf@strata[ , c('Lon','Lat')]
for(i in c('Lon','Lat')){coords[,i] <- as.numeric(as.matrix(coords[,i]))}
crs <- CRS("+init=epsg:3160")

# make the SpatialPointsDataFrame object
pts <- SpatialPointsDataFrame(coords = coords, data= gen.Gulf@strata, proj4string = crs)

#Look for intersection
new_shape <- point.in.poly(pts, polys)

#Output the results
strata(gen.Gulf) <- new_shape@data[match(indNames(gen.Gulf),new_shape@data$INDV),]
for(i in c("Lon", "Lat", "N")){gen.Gulf@strata[,i] <- as.numeric(as.matrix(gen.Gulf@strata[,i]))}
head(gen.Gulf@strata)

#Making an object incuding the imaginary points
which(npoints$N %in% gen.Gulf@strata$N)
new.lats <- rbind(gen.Gulf@strata[,c("Lon","Lat","N")], npoints)
rownames(new.lats) <- c(indNames(gen.Gulf),rownames(npoints))
write.table(new.lats, file="new.lats_v2", col.names=T, row.names=T, quote=F)

#Importing new.lats
new.lats <- read.table("new.lats_v2", head=T, row.names=1)

#Finding Location of points near the Campeche, MX point that was originally used
st_pt <- new.lats[which(rownames(new.lats) == "CH_002.I1.Lib5"),1:2]
end_pt <- as.matrix(new.lats[which(rownames(new.lats) == "CH_002.I1.Lib5"),1:2])
for(i in 1:10){
tmp_pt <- destPoint(destPoint(st_pt, 0, i*10000), 90, -i*10000)
end_pt <- rbind(end_pt, tmp_pt)
#end_pt <- rbind(end_pt, destPoint(st_pt, 45, -i*10000))
rownames(end_pt) <- c(rownames(end_pt)[1:i], paste("Pt",i,sep="_"))
}
end_pt

coords <- as.data.frame(end_pt)
for(i in c('Lon','Lat')){coords[,i] <- as.numeric(as.matrix(coords[,i]))}
crs <- CRS("+init=epsg:3160")
pts <- SpatialPointsDataFrame(coords = coords, data= coords, proj4string = crs)
new_shape <- point.in.poly(pts, polys)
new_shape@data

new.lats[new.lats$N == new_shape@data[1,3],1] <- new_shape@data[4,1]
new.lats[new.lats$N == new_shape@data[1,3],2] <- new_shape@data[4,2]
new.lats[new.lats$N == new_shape@data[1,3],3] <- new_shape@data[4,3]
write.table(new.lats, file="new.lats_v2", col.names=T, row.names=T, quote=F)


