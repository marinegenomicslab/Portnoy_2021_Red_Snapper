#!/bin/R

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
library("fitdistrplus")
library('CaDENCE')

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

##### CA Analysis with RDA #####
## Running analysis on the grids using all the imaginary points ##
# dbMEM #
#Run a CA on the genind object
setPop(gen.Gulf) <- ~N
obj <- genind2genpop(gen.Gulf)
ca1 <- dudi.coa(tab(obj),scannf=FALSE,nf=1000)
out.data <- ca1$li
rownames(out.data) <- levels(pop(gen.Gulf))

#Getting the unique grid locations
uniq.new.lats <-data.frame(matrix(ncol=3))
colnames(uniq.new.lats) <- colnames(new.lats)
for(i in unique(new.lats$N)){
tmp.v <- apply(new.lats[new.lats$N == i,], 2, mean)
uniq.new.lats <- rbind(uniq.new.lats,c(tmp.v,i))
}
uniq.new.lats <- uniq.new.lats[-1,]
}

#Getting distance matrix
dist.mat <- distm(uniq.new.lats[,1:2], fun=distHaversine)/1000
colnames(dist.mat) <- rownames(dist.mat) <- uniq.new.lats$N

#Defining the threshold
threshh <- 144

#Making matrix for the truncated network
nb.mat <- 1 - (dist.mat/(4 *threshh))^2
nb.mat[dist.mat > threshh] <- 0
diag(nb.mat) <- 0

#Making the network
nb <- mat2listw(nb.mat)

#Plot of connection network
tiff(filename ="dbMEM_Connection_network_CA.tiff", res=200, width=2000, height=2000)
par(mfrow=c(2,1))
plot(nb, uniq.new.lats, bg="red4", pch=21)
mtext(paste("dbMEM Nearest neighbors (thresh =", round(threshh,4), ")"), 3, cex=2, font=2, line=0.25)
grid(col="grey60")
box()
#Plot of distance vs spatial weight
plot(dist.mat[nb.mat!=0], nb.mat[nb.mat!=0], pch=21, bg="red4", ylab="dbMEM spatial weight", xlab="Distance (km)")
dev.off()

#Diagnolization of listw object
MEM.autocor <- "positive"
Dist_MEM <- adespatial::scores.listw(nb, MEM.autocor = MEM.autocor, store.listw = T)
dim(Dist_MEM)

save(Dist_MEM, file="dbMEMs.gz", compress=T)
#load("dbMEMs.gz")

#RDA
#Looking at VIF and correlation
Dist_MEM <- Dist_MEM[which(! uniq.new.lats$N %in% new.lats$N[grep("Pt_",rownames(new.lats))]),]
rda.dbMEM <- rda(out.data ~ ., Dist_MEM)
vif.dbMEM <- vif.cca(rda.dbMEM)
max(vif.dbMEM, na.rm=T)

cor.m <- cor(Dist_MEM)
for(i in 1:nrow(cor.m)){cor.m[i,i] <- 0}
max(cor.m)

#Removing multicollinearity within dbMEM data
new.data <- Dist_MEM
PC.rm.list <- vector()

#Removes the highest Axis until VIF drops below 3
while(max(vif.dbMEM) > 3){
MAX <- max(vif.dbMEM)
vif.test <- vif.dbMEM[which(vif.dbMEM==MAX)]
for(i in 1:length(vif.test)){
tmp.dat <- new.data[ , which(!(names(new.data) %in% names(vif.test)[i]))]
rda.dbMEM <- rda(out.data ~ ., tmp.dat)
vif.dbMEM <- vif.cca(rda.dbMEM)
if(max(vif.dbMEM)<MAX){
PC.rm.list <- c(PC.rm.list,names(vif.test)[i])
print(paste(names(vif.test)[i],": Max VIF = ", round(max(vif.dbMEM),3)))
break}
}
if(max(vif.dbMEM)>=MAX){print("All terms checked and none are colinear with max VIF"); break}
new.data <- tmp.dat
}
Dist_MEM <- new.data

#RDA dbMEM
{```{R}```
X.Gulf <- scaleGen(gen.Gulf, NA.method="mean", scale=F, center=T)
m1<-rda(out.data ~ ., Dist_MEM)
m0<-rda(out.data ~ 1, Dist_MEM)

set.seed(1235)
m.ord_dbMEM <- ordiR2step(m0, scope = formula(m1), Pin=0.05, Pout=0.1, permutations = 999, parallel=25)
m.ord_dbMEM

m.ord_dbMEM$anova

RsquareAdj(m.ord_dbMEM)

m.ord_dbMEM$CCA$eig[1:7]/m.ord_dbMEM$tot.chi

anova(m.ord_dbMEM, parallel=20)

tmp_anova <- anova(m.ord_dbMEM, by="axis", parallel=30)
tmp_anova
tmp_anova$Variance[1:4]/sum(tmp_anova$Variance)

anova_MARG <- anova(m.ord_dbMEM, by="margin", parallel=30)
anova_MARG
cbind(rownames(anova_MARG),round(anova_MARG$Variance[1:length(anova_MARG$Variance)]/sum(anova_MARG$Variance),6))

#Plotting MEMs
setPop(gen.Gulf) <- ~SubPOP
POP <- vector()
for(i in uniq.new.lats$N[!(uniq.new.lats$N %in% new.lats$N[grep("Pt",rownames(new.lats))])]){
if(i == 25742){POP <- c(POP,"CH_MX");next}
POP <- c(POP,as.character(head(pop(gen.Gulf)[gen.Gulf@strata$N == i],n=1)))}

tmp_data <- data.frame(MEM=Dist_MEM, POP=POP, Coords=uniq.new.lats[which(! uniq.new.lats$N %in% new.lats$N[grep("Pt_",rownames(new.lats))]),1:2,], RDA=m.ord_dbMEM$CCA$wa)
tmp_data$POP <- factor(tmp_data$POP, levels = c("DT", "Gulf_FL_S", "Gulf_FL_N", "PC", "AL", "E_LA", "W_LA", "TX", "N_MX", "VC_MX", "CP_MX"))

dbMEM20.map <- Gulf_map + geom_point(aes(x = Coords.Lon, y = Coords.Lat, fill=MEM.MEM20), data=tmp_data, shape=21, color="black", stroke=0.1, size =2.5) + scale_fill_gradient(low="black", high="red") + labs(fill = "dbMEM 20")
dbMEM20.map
 
RDA1.map <- Gulf_map + geom_point(aes(x = Coords.Lon, y = Coords.Lat, fill=RDA1), data=tmp_data, shape=21, color="black", stroke=0.1, size =2.5) + scale_fill_gradient(low="black", high="red") + labs(fill = "dbMEM RDA 1")
RDA1.map

ggsave(dbMEM20.map, file="dbMEM_CA_MEM20.tif", device="tiff")
ggsave(RDA1.map, file="dbMEM_CA_RDA1.tif", device="tiff")

# Adult MEMs #
# Run a CA on the genind object
setPop(gen.Gulf) <- ~N
obj <- genind2genpop(gen.Gulf)
ca1 <- dudi.coa(tab(obj),scannf=FALSE,nf=1000)
out.data <- ca1$li
rownames(out.data) <- levels(pop(gen.Gulf))

#Getting the unique grid locations
uniq.new.lats <-data.frame(matrix(ncol=3))
colnames(uniq.new.lats) <- colnames(new.lats)
for(i in unique(new.lats$N)){
tmp.v <- apply(new.lats[new.lats$N == i,], 2, mean)
uniq.new.lats <- rbind(uniq.new.lats,c(tmp.v,i))
uniq.new.lats <- uniq.new.lats[-1,]

#Getting Adult model
Adult.list <- list()
dat <- read.table("Adults_movement.txt", head=T)
flz <- fitdist(dat$Dispersion_m.day, "blnorm", start=list(prob=0.65, meanlog=1.5, sdlog=2.1))
flz

#Getting Adult distance
Adult.mat <- matrix(nrow=nrow(uniq.new.lats), ncol=nrow(uniq.new.lats))
for(i in 1:nrow(uniq.new.lats)){
tmp.v <- vector()
for(j in 1:nrow(uniq.new.lats)){
ROW <- uniq.new.lats[i, 1:2]
COL <- uniq.new.lats[j, 1:2]
tmp.geo <- distm(ROW, COL, fun=distHaversine)
if(!(exists('tmp.v'))){tmp.v <- tmp.geo} else {tmp.v <- c(tmp.v,tmp.geo)}
}
tmp.v2 <- tmp.v/6935	#19 years
tmp.v3 <- dblnorm(tmp.v2, flz$estimate[1], flz$estimate[2], flz$estimate[3])
Adult.mat[i,] <- tmp.v3
}
colnames(Adult.mat) <- rownames(Adult.mat) <- uniq.new.lats$N

Adult.dist <- sqrt(1-Adult.mat)

#Defining the threshold
threshh <- 0.996504515

#Making matrix for the truncated network
nb.mat <- 1 - (Adult.dist/(4 *threshh))^2
nb.mat[Adult.dist > threshh] <- 0
diag(nb.mat) <- 0

#Making the network
nb <- mat2listw(nb.mat)

#Plot of connection network
tiff(filename ="Adult_MEM_Connection_network_CA.tiff", res=200, width=2000, height=2000)
par(mfrow=c(2,1))
plot(nb, uniq.new.lats[uniq.new.lats$N != 26342,], col="red4", pch=19)
mtext(paste("Nearest neighbors (thresh=", round(threshh,4), ")"), 3, cex=2, font=2, line=0.25)
grid(col="grey60")
box()
#Plotting distance vs weight
dist.mat <- distm(uniq.new.lats[,1:2], fun=distHaversine)/1000
plot(dist.mat[nb.mat!=0], nb.mat[nb.mat!=0], pch=21, bg="red4", ylab="dbMEM spatial weight", xlab="Distance (km)")
dev.off()

#Diagnolization of listw object
MEM.autocor="positive"
Adult_MEM <- adespatial::scores.listw(nb, MEM.autocor = MEM.autocor, store.listw = T)
dim(Adult_MEM)

#RDA Adult MEMs
#Looking at VIF and correlation
Adult_MEM <- Adult_MEM[which(! uniq.new.lats$N %in% new.lats$N[grep("Pt_",rownames(new.lats))]),]
rda.A_MEM <- rda(out.data ~ ., Adult_MEM)
vif.A_MEM <- vif.cca(rda.A_MEM)
max(vif.A_MEM, na.rm=T)

cor.m <- cor(Adult_MEM)
for(i in 1:nrow(cor.m)){cor.m[i,i] <- 0}
max(cor.m)

#Removing multicollinearity within Adult MEM data
new.data <- Adult_MEM
PC.rm.list <- vector()

#Removes the highest Axis until VIF drops
while(max(vif.A_MEM) > 3){
MAX <- max(vif.A_MEM)
vif.test <- vif.A_MEM[which(vif.A_MEM==MAX)]
for(i in 1:length(vif.test)){
tmp.dat <- new.data[ , which(!(names(new.data) %in% names(vif.test)[i]))]
rda.A_MEM <- rda(out.data ~ ., tmp.dat)
vif.A_MEM <- vif.cca(rda.A_MEM)
if(max(vif.A_MEM)<MAX){
PC.rm.list <- c(PC.rm.list,names(vif.test)[i])
print(paste(names(vif.test)[i],": Max VIF = ", round(max(vif.A_MEM),3)))
break}
}
if(max(vif.A_MEM)>=MAX){print("All terms checked and none are colinear with max VIF"); break}
new.data <- tmp.dat
}
Adult_MEM <- new.data

#RDA Adult_MEM
m1<-rda(out.data ~ ., Adult_MEM)
m0<-rda(out.data ~ 1, Adult_MEM)

set.seed(1235)
m.ord_Adult <- ordiR2step(m0, scope = formula(m1), Pin=0.05, Pout=0.1, permutations = 999, parallel=25)
m.ord_Adult

m.ord_Adult$anova

RsquareAdj(m.ord_Adult)

m.ord_Adult$CCA$eig[1:7]/m.ord_Adult$tot.chi

anova(m.ord_Adult, parallel=20)

tmp_anova <- anova(m.ord_Adult, by="axis", parallel=30)
tmp_anova

tmp_anova$Variance[1:4]/sum(tmp_anova$Variance)

anova_MARG <- anova(m.ord_Adult, by="margin", parallel=30)
anova_MARG

cbind(rownames(anova_MARG),round(anova_MARG$Variance[1:length(anova_MARG$Variance)]/sum(anova_MARG$Variance),6))

#Plotting MEMs (distance)
setPop(gen.Gulf) <- ~SubPOP
POP <- vector()
for(i in uniq.new.lats$N[!(uniq.new.lats$N %in% new.lats$N[grep("Pt",rownames(new.lats))])]){
if(i == 25742){POP <- c(POP,"CP_MX");next}
POP <- c(POP,as.character(head(pop(gen.Gulf)[gen.Gulf@strata$N == i],n=1)))}

tmp_data <- data.frame(MEM=Adult_MEM, POP=POP, Coords=uniq.new.lats[which(! uniq.new.lats$N %in% new.lats$N[grep("Pt",rownames(new.lats))]),1:2], RDA=m.ord_Adult$CCA$wa)
tmp_data$POP <- factor(tmp_data$POP, levels = c("DT", "Gulf_FL_S", "Gulf_FL_N", "PC", "AL", "E_LA", "W_LA", "TX", "N_MX", "VC_MX", "CP_MX"))

A_MEM11.map <- Gulf_map + geom_point(aes(x = Coords.Lon, y = Coords.Lat, fill=MEM.MEM11), data=tmp_data, shape=21, color="black", stroke=0.1, size =2.5) + scale_fill_gradient(low="black", high="red") + labs(fill = "Adult MEM 11")
A_MEM20.map <- Gulf_map + geom_point(aes(x = Coords.Lon, y = Coords.Lat, fill=MEM.MEM20), data=tmp_data, shape=21, color="black", stroke=0.1, size =2.5) + scale_fill_gradient(low="black", high="red") + labs(fill = "Adult MEM 20")

multiplot(A_MEM11.map, A_MEM20.map)
ggsave(multiplot(A_MEM11.map, A_MEM20.map), file="Adult_CA_MEMs_map.tif", device="tiff")

RDA1.map <- Gulf_map + geom_point(aes(x = Coords.Lon, y = Coords.Lat, fill=RDA.RDA1), data=tmp_data, shape=21, color="black", stroke=0.1, size =2.5) + scale_fill_gradient(low="black", high="red") + labs(fill = "Adult RDA 1")
RDA2.map <- Gulf_map + geom_point(aes(x = Coords.Lon, y = Coords.Lat, fill=RDA.RDA2), data=tmp_data, shape=21, color="black", stroke=0.1, size =2.5) + scale_fill_gradient(low="black", high="red") + labs(fill = "Adult RDA 2")

multiplot(RDA1.map, RDA2.map)
ggsave(multiplot(RDA1.map, RDA2.map), file="Adult_CA_RDA_map.tif", device="tiff")

# Larval MEMs #
# Run a CA on the genind object
setPop(gen.Gulf) <- ~N
obj <- genind2genpop(gen.Gulf)
ca1 <- dudi.coa(tab(obj),scannf=FALSE,nf=1000)
s.label(ca1$li, sub="CA 1-2",csub=2)
out.data <- ca1$li
rownames(out.data) <- levels(pop(gen.Gulf))

#Putting regions and colors to the CA
N_regions <- unique(gen.Gulf@strata[,c("N","SubPOP")])
N_regions$N[N_regions$N == 26342] <- 25742
#Correcting for adjusting CP_MX near shore site
c_N <- as.character(N_regions$SubPOP)
c_N[c_N=="DT"] <- c_RDA[1]
c_N[c_N=="Gulf_FL_S"] <- c_RDA[2]
c_N[c_N=="Gulf_FL_N"] <- c_RDA[3]
c_N[c_N=="PC"] <- c_RDA[4]
c_N[c_N=="AL"] <- c_RDA[5]
c_N[c_N=="E_LA"] <- c_RDA[6]
c_N[c_N=="W_LA"] <- c_RDA[7]
c_N[c_N=="TX"] <- c_RDA[8]
c_N[c_N=="N_MX"] <- c_RDA[9]
c_N[c_N=="VC_MX"] <- c_RDA[10]
c_N[c_N=="CP_MX"] <- c_RDA[11]

#Getting the unique grid locations
uniq.new.lats <-data.frame(matrix(ncol=3))
colnames(uniq.new.lats) <- colnames(new.lats)
for(i in unique(new.lats$N)){
tmp.v <- apply(new.lats[new.lats$N == i,], 2, mean)
uniq.new.lats <- rbind(uniq.new.lats,c(tmp.v,i))
}
uniq.new.lats <- uniq.new.lats[-1,]

#On DT; Getting the matrix with the new points
Model.dat <- read.table("All_link.txt", head=T, row.names=1)
colnames(Model.dat) <- rownames(Model.dat)

Larv.mat <- matrix(nrow=nrow(uniq.new.lats), ncol=nrow(uniq.new.lats))
for(i in 1:nrow(uniq.new.lats)){
tmp.v <- vector()
for(j in as.character(uniq.new.lats$N)){
ROW <- as.character(uniq.new.lats$N)[i]
tmp.v <- c(tmp.v,Model.dat[ROW,j])
}
Larv.mat[i,] <- tmp.v
}
colnames(Larv.mat) <- rownames(Larv.mat) <- uniq.new.lats$N

#Averaging the similarity
Larv.avg <- matrix(nrow=nrow(uniq.new.lats), ncol=nrow(uniq.new.lats))
colnames(Larv.avg) <- rownames(Larv.avg) <- uniq.new.lats$N
for(i in 1:nrow(uniq.new.lats)){
for(j in 1:nrow(uniq.new.lats)){
INDV1 <- rownames(Larv.avg)[i]
INDV2 <- colnames(Larv.avg)[j]
Larv.avg[i, j] <- Larv.avg[j, i] <- mean(Larv.mat[INDV1, INDV2], Larv.mat[INDV2, INDV1])
}
}

#Getting Larval distance
Larval.dist <- sqrt(1-Larv.avg)

#Defining the threshold
threshh <- give.thresh(as.dist(Larval.dist))
threshh

#Making the weighted network
nb.mat <- 1 - (Larval.dist/(4 *threshh))^2
nb.mat[Larval.dist > threshh] <- 0
diag(nb.mat) <- 0

#Making the network
nb <- mat2listw(nb.mat)

#Plot of connection network
tiff(filename ="Larval_MEM_Connection_network_CA.tiff", res=200, width=2000, height=2000)
par(mfrow=c(2,1))
plot(nb, uniq.new.lats, col="red4", pch=19)
mtext(paste("Nearest neighbors (thresh=", round(threshh,4), ")"), 3, cex=2, font=2, line=0.25)
grid(col="grey60")
box()
#Plot of distance vs spatial weight
dist.mat <- distm(uniq.new.lats[,1:2], fun=distHaversine)/1000
plot(dist.mat[nb.mat!=0], nb.mat[nb.mat!=0], pch=21, bg="red4", ylab="dbMEM spatial weight", xlab="Distance (km)")
dev.off()

#Diagnolization of listw object
MEM.autocor="positive"
Larval_MEM <- adespatial::scores.listw(nb, MEM.autocor = MEM.autocor, store.listw = T)
dim(Larval_MEM)

#RDA Larval MEMs
#Looking at VIF and correlation
Larval_MEM <- Larval_MEM[which(! uniq.new.lats$N %in% new.lats$N[grep("Pt_",rownames(new.lats))]),]
rda.L_MEM <- rda(out.data ~ ., Larval_MEM)
vif.L_MEM <- vif.cca(rda.L_MEM)
max(vif.L_MEM, na.rm=T)

cor.m <- cor(Larval_MEM)
for(i in 1:nrow(cor.m)){cor.m[i,i] <- 0}
max(cor.m)

# Removing multicollinearity within dbMEM data
new.data <- Larval_MEM
PC.rm.list <- vector()

#Removes the highest Axis until VIF drops
while(max(vif.L_MEM) > 3){
MAX <- max(vif.L_MEM)
vif.test <- vif.L_MEM[which(vif.L_MEM==MAX)]
for(i in 1:length(vif.test)){
tmp.dat <- new.data[ , which(!(names(new.data) %in% names(vif.test)[i]))]
rda.L_MEM <- rda(out.data ~ ., tmp.dat)
vif.L_MEM <- vif.cca(rda.L_MEM)
if(max(vif.L_MEM)<MAX){
PC.rm.list <- c(PC.rm.list,names(vif.test)[i])
print(paste(names(vif.test)[i],": Max VIF = ", round(max(vif.L_MEM),3)))
break}
}
if(max(vif.L_MEM)>=MAX){print("All terms checked and none are colinear with max VIF"); break}
new.data <- tmp.dat
}
Larval_MEM <- new.data

#RDA Larval MEMs
m1<-rda(out.data ~ ., Larval_MEM)
m0<-rda(out.data ~ 1, Larval_MEM)

set.seed(1235)
m.ord_Larval <- ordiR2step(m0, scope = formula(m1), Pin=0.05, Pout=0.1, permutations = 999, parallel=25)
m.ord_Larval

m.ord_Larval$anova

RsquareAdj(m.ord_Larval)

m.ord_Larval$CCA$eig[1:7]/m.ord_Larval$tot.chi

anova(m.ord_Larval, parallel=20)

tmp_anova <- anova(m.ord_Larval, by="axis", parallel=30)
tmp_anova
tmp_anova$Variance[1:4]/sum(tmp_anova$Variance)

anova_MARG <- anova(m.ord_Larval, by="margin", parallel=30)
anova_MARG

cbind(rownames(anova_MARG),round(anova_MARG$Variance[1:length(anova_MARG$Variance)]/sum(anova_MARG$Variance),6))

#Plotting Larval MEMs
setPop(gen.Gulf) <- ~SubPOP
POP <- vector()
for(i in uniq.new.lats$N[!(uniq.new.lats$N %in% new.lats$N[grep("Pt",rownames(new.lats))])]){
if(i == 25742){POP <- c(POP,"CP_MX");next}
POP <- c(POP,as.character(head(pop(gen.Gulf)[gen.Gulf@strata$N == i],n=1)))}

tmp_data <- data.frame(MEM=Larval_MEM, POP=POP, Coords=uniq.new.lats[which(! uniq.new.lats$N %in% new.lats$N[grep("Pt_",rownames(new.lats))]),1:2], RDA=m.ord_Larval$CCA$wa)
tmp_data$POP <- factor(tmp_data$POP, levels = c("DT", "Gulf_FL_S", "Gulf_FL_N", "PC", "AL", "E_LA", "W_LA", "TX", "N_MX", "VC_MX", "CP_MX"))

L_MEM10.map <- Gulf_map + geom_point(aes(x = Coords.Lon, y = Coords.Lat, fill=MEM.MEM10), data=tmp_data, shape=21, color="black", stroke=0.1, size =2.5) + scale_fill_gradient(low="black", high="red") + labs(fill = "Larval MEM 10")
L_MEM10.map

ggsave(L_MEM10.map, file="LMEM10_map.tif", device="tiff")

RDA1.map <- Gulf_map + geom_point(aes(x = Coords.Lon, y = Coords.Lat, fill=RDA1), data=tmp_data, shape=21, color="black", stroke=0.1, size =2.5) + scale_fill_gradient(low="black", high="red") + labs(fill = "Larval RDA 1")
RDA1.map

ggsave(RDA1.map, file="Larval_CA_RDA1.tif", device="tiff")

# Combo #
#Using just the significant MEMs
#Combining data
names(Adult_MEM) <- apply(data.frame(names(Adult_MEM)), 1, function(x) paste(x, "Adult", sep="_"))
names(Larval_MEM) <- apply(data.frame(names(Larval_MEM)), 1, function(x) paste(x, "Larval", sep="_"))
RDA_data <- cbind(Dist_MEM[,grepl("MEM20", colnames(Dist_MEM))], Adult_MEM[,grepl("MEM11", colnames(Adult_MEM)) | grepl("MEM20", colnames(Adult_MEM))], Larval_MEM[,grepl("MEM10", colnames(Larval_MEM))])
colnames(RDA_data) <- c("MEM20","MEM11_Adult","MEM20_Adult","MEM10_Larval")
dim(RDA_data)

#RDA
#Looking at VIF and correlation
rda.Comb <- rda(out.data ~ ., RDA_data)
vif.Comb <- vif.cca(rda.Comb)
max(vif.Comb)

cor.m <- cor(RDA_data)
for(i in 1:nrow(cor.m)){cor.m[i,i] <- 0}
max(cor.m)

cor.m

#RDA Combo_MEM
m1<-rda(out.data ~ ., RDA_data)
m0<-rda(out.data ~ 1, RDA_data)

set.seed(1235)
m.ord_Combine <- ordiR2step(m0, scope = formula(m1), Pin=0.05, Pout=0.1, permutations = 999, parallel=20, )
m.ord_Combine

m.ord_Combine$anova

RsquareAdj(m.ord_Combine)

m.ord_Combine$CCA$eig[1:7]/m.ord_Combine$tot.chi

anova(m.ord_Combine, parallel=20)

tmp_anova <- anova(m.ord_Combine, by="axis", parallel=30)
tmp_anova
tmp_anova$Variance[1:5]/sum(tmp_anova$Variance)

anova_MARG <- anova(m.ord_Combine, by="margin", parallel=30)
anova_MARG
cbind(rownames(anova_MARG),round(anova_MARG$Variance[1:length(anova_MARG$Variance)]/sum(anova_MARG$Variance),6))

#Constrained analyses of each set of MEMs
#dbMEMs
prda1 <- rda(out.data ~  MEM20 + Condition(MEM11_Adult + MEM20_Adult + MEM10_Larval), data = RDA_data)
prda1_anov <- anova(prda1, permutations = 999, parallel=20)
cbind(rownames(prda1_anov),round(prda1_anov$Variance[1:length(prda1_anov$Variance)]/sum(prda1_anov$Variance),6))
RsquareAdj(prda1)
anova(prda1, by = "margin", permutations = 999, parallel=20)

#Adult MEMs
prda2 <- rda(out.data ~  MEM11_Adult + MEM20_Adult + Condition(MEM10_Larval + MEM20), data = RDA_data)
prda2_anov <- anova(prda2, permutations = 999, parallel=20)
cbind(rownames(prda2_anov),round(prda2_anov$Variance[1:length(prda2_anov$Variance)]/sum(prda2_anov$Variance),6))
RsquareAdj(prda2)
anova(prda2, by = "margin", permutations = 999, parallel=20)

#Larval MEMs
prda3 <- rda(out.data ~ MEM10_Larval + Condition(MEM20 + MEM11_Adult + MEM20_Adult), data = RDA_data)
prda3_anov <- anova(prda3, permutations = 999, parallel=20)
cbind(rownames(prda3_anov),round(prda3_anov$Variance[1:length(prda3_anov$Variance)]/sum(prda3_anov$Variance),6))
RsquareAdj(prda3)
anova(prda3, by = "margin", permutations = 999, parallel=20)

#Plotting MEMs
POP <- vector()
for(i in uniq.new.lats$N[!(uniq.new.lats$N %in% new.lats$N[grep("Pt",rownames(new.lats))])]){
if(i == 25742){POP <- c(POP,"CP_MX");next}
POP <- c(POP,as.character(head(pop(gen.Gulf)[gen.Gulf@strata$N == i],n=1)))}

tmp_data <- data.frame(MEM=RDA_data, POP=POP, Coords=uniq.new.lats[which(! uniq.new.lats$N %in% new.lats$N[grep("Pt_",rownames(new.lats))]),1:2], RDA=m.ord_Combine$CCA$wa, MAL=mahalanobis(m.ord_Combine$CCA$wa[,1:2], center=F, cov(m.ord_Combine$CCA$wa[,1:2])))
tmp_data$POP <- factor(tmp_data$POP, levels = c("DT", "Gulf_FL_S", "Gulf_FL_N", "PC", "AL", "E_LA", "W_LA", "TX", "N_MX", "VC_MX", "CP_MX"))

MEM20.map <- Gulf_map + geom_point(aes(x = Coords.Lon, y = Coords.Lat, fill=MEM.MEM20), data=tmp_data, shape=21, color="black", stroke=0.1, size =2.5) + scale_fill_gradient(low="black", high="red") + labs(fill = "dbMEM 20")
AMEM11.map <- Gulf_map + geom_point(aes(x = Coords.Lon, y = Coords.Lat, fill=MEM.MEM11_Adult), data=tmp_data, shape=21, color="black", stroke=0.1, size =2.5) + scale_fill_gradient(low="black", high="red") + labs(fill = "Adult MEM 11")
AMEM20.map <- Gulf_map + geom_point(aes(x = Coords.Lon, y = Coords.Lat, fill=MEM.MEM20_Adult), data=tmp_data, shape=21, color="black", stroke=0.1, size =2.5) + scale_fill_gradient(low="black", high="red") + labs(fill = "Adult MEM 20")
LMEM10.map <- Gulf_map + geom_point(aes(x = Coords.Lon, y = Coords.Lat, fill=MEM.MEM10_Larval), data=tmp_data, shape=21, color="black", stroke=0.1, size =2.5) + scale_fill_gradient(low="black", high="red") + labs(fill = "Larval MEM 10")

multiplot(MEM20.map, AMEM11.map, LMEM10.map, AMEM20.map, cols=2)
ggsave(multiplot(MEM20.map, AMEM11.map, LMEM10.map, AMEM20.map, cols=2), file="Combine_CA_MEM_maps.tif", device="tiff")

RDA1.map <- Gulf_map + geom_point(aes(x = Coords.Lon, y = Coords.Lat, fill=RDA.RDA1), data=tmp_data, shape=21, color="black", stroke=0.1, size =2.5) + scale_fill_gradient(low="black", high="red") + labs(fill = "Full RDA 1")
RDA2.map <- Gulf_map + geom_point(aes(x = Coords.Lon, y = Coords.Lat, fill=RDA.RDA2), data=tmp_data, shape=21, color="black", stroke=0.1, size =2.5) + scale_fill_gradient(low="black", high="red") + labs(fill = "Full RDA 2")
RDA3.map <- Gulf_map + geom_point(aes(x = Coords.Lon, y = Coords.Lat, fill=RDA.RDA3), data=tmp_data, shape=21, color="black", stroke=0.1, size =2.5) + scale_fill_gradient(low="black", high="red") + labs(fill = "Full RDA 3")
RDA4.map <- Gulf_map + geom_point(aes(x = Coords.Lon, y = Coords.Lat, fill=RDA.RDA4), data=tmp_data, shape=21, color="black", stroke=0.1, size =2.5) + scale_fill_gradient(low="black", high="red") + labs(fill = "Full RDA 4")

multiplot(RDA1.map, RDA3.map, RDA2.map, RDA4.map, cols=2)
ggsave(multiplot(RDA1.map, RDA3.map, RDA2.map, RDA4.map, cols=2), file="Combine_CA_RDA_Axis_maps.tif", device="tiff")

rda1 <- ggplot(data=tmp_data, aes(x=RDA.RDA1, y=RDA.RDA2, color=POP.SubPOP, fill=POP.SubPOP)) + geom_point(shape=21) + scale_fill_manual(values=c_RDA) + 
scale_color_manual(values=c_RDA) + stat_ellipse(aes(group=POP.SubPOP), type = "t") + theme_classic() + labs(x="RDA1", y="RDA2", color="Populations", fill="Populations")

rda2 <- ggplot(data=tmp_data, aes(x=RDA.RDA3, y=RDA.RDA4, color=POP.SubPOP, fill=POP.SubPOP)) + geom_point(shape=21) + scale_fill_manual(values=c_RDA) + 
scale_color_manual(values=c_RDA) + stat_ellipse(aes(group=POP.SubPOP), type = "t") + theme_classic() + labs(x="RDA3", y="RDA4", color="Populations", fill="Populations")

multiplot(rda1, rda2)
ggsave(multiplot(rda1, rda2), file="RDA_bipots.tif", device="tiff")

#Making the sine wave for the MEMs
#Preparing igraph
dist.mat <- distm(uniq.new.lats[!(uniq.new.lats$N %in% npoints$N[1:4]),1:2], fun=distHaversine)/1000
colnames(dist.mat) <- rownames(dist.mat) <- uniq.new.lats$N[!(uniq.new.lats$N %in% npoints$N[1:4])]

threshh <- 144
nb.mat <- 1 - (dist.mat/(4 *threshh))^2
nb.mat[dist.mat > threshh] <- 0
diag(nb.mat) <- 0

tmp_graph <- graph.adjacency(nb.mat, mode="undirected", weighted=T)

#Getting distance along network
DIST_m <- NULL
for(i in uniq.new.lats$N[!(uniq.new.lats$N %in% npoints$N[1:4])]){
if(i == "14682"){next}
ROUTE <- attr(shortest_paths(tmp_graph, "14682", as.character(i))$vpath[[1]], "names")
DIST <- 0
for(j in 1:(length(ROUTE)-1)){DIST <- DIST + dist.mat[ROUTE[j], ROUTE[j+1]]}
DIST_m <- rbind(DIST_m, c(i,DIST))
}
DIST_m <- rbind(DIST_m, c("14682", 0))

#Setting up DIST object for sine wave modeling
tmp_df2 <- data.frame(Cells=rownames(Dist_MEM), MEM20=as.numeric(Dist_MEM$MEM20))
tmp_df2$DIST <- as.numeric(DIST_m[match(tmp_df2[,1], DIST_m[,1]),2])
tmp_df2 <- tmp_df2[order(tmp_df2$DIST),]

boxplot(tmp_df2$MEM20, plot=F)$stats[c(2,4)]

#Getting Groups for coloring
SUBPOPs <- unique(N_regions$SubPOP)
#           "Gulf_FL_N"    "Gulf_FL_S"         "DT"      "AL"     "W_LA"             " E_LA"      "CP_MX"      "TX"       "N_MX"  "VC_MX"      "PC"
c_MEM <- c("mediumblue", "lightslateblue", "cadetblue", "cyan", "darkolivegreen3", "darkgreen", "grey25", "seagreen4", "grey65", "grey40", "steelblue1")

#Making a dataframe
SUBPOPs <- data.frame(POP=as.character(as.matrix(SUBPOPs)), Col=c_MEM)

#Getting the Min and Max Distance values for each group
POP_RANGE <- NULL
for(i in as.character(SUBPOPs$POP)){
CELLS <- N_regions[N_regions$SubPOP == i, "N"]
POP_DIST <- tmp_df2[tmp_df2[,1] %in% CELLS, "DIST"]
POP_RANGE <- rbind(POP_RANGE, c(i,min(POP_DIST), max(POP_DIST)))
}
for(i in 2:3){POP_RANGE[,i] <- as.numeric(as.matrix(POP_RANGE[,i]))}
SUBPOPs$Min <- as.numeric(POP_RANGE[,2])
SUBPOPs$Max <- as.numeric(POP_RANGE[,3])
}

#Modeling sine wave for Dist MEM
AIC_results <- NULL
for(i in seq(1, 5000, 1)){
PER <- (2 * 3.1415)/i
AIC_results <- rbind(AIC_results, c(i, AIC(lm(MEM20 ~ sin(PER*DIST) + cos(PER*DIST), data=tmp_df2))))
}
par(mfrow=c(3,1))
hist(AIC_results[,2], breaks=100, col="red4")
plot(AIC_results, xlab="Period", ylab="AIC", pch=19)
AIC_results[AIC_results[,2] == min(AIC_results[,2]),]
AIC_results[AIC_results[,2] == min(AIC_results[10:1000,2]),]
PER <- (2 * 3.1415)/(AIC_results[which(AIC_results[,2] == min(AIC_results[,2])),1])
fit <- lm(MEM20 ~ sin(PER*DIST) + cos(PER*DIST), data=tmp_df2)
fit
plot(MEM20 ~ DIST, data=tmp_df2, pch=21, bg="red4", xlab="Distance (km)", ylab="Dist MEM20")
lines(seq(1,ceiling(max(tmp_df2$DIST))), predict.lm(fit, data.frame(DIST=seq(1,ceiling(max(tmp_df2$DIST))))), col="mediumblue")

PER <- (2 * 3.1415)/413
fit <- lm(MEM1 ~ sin(PER*DIST) + cos(PER*DIST), data=tmp_df2[which(tmp_df2$MEM20 > 0.5 | tmp_df2$MEM20 < 0),])
fit
lines(seq(1,ceiling(max(tmp_df2$DIST))), predict.lm(fit, data.frame(DIST=seq(1,ceiling(max(tmp_df2$DIST))))), col="darkgreen", lty=3)

tiff("Dist_MEM20.tiff", res=200, width=2500, height=2000)

par(mfrow=c(2,1))
#Best AIC
PER <- (2 * 3.1415)/(AIC_results[which(AIC_results[,2] == min(AIC_results[,2])),1])
fit <- lm(MEM20 ~ sin(PER*DIST) + cos(PER*DIST), data=tmp_df2)
fit
plot(MEM20 ~ DIST, data=tmp_df2, pch=21, bg="red4", xlab="Distance (km)", ylab="Dist MEM20")
lines(seq(1,ceiling(max(tmp_df2$DIST))), predict.lm(fit, data.frame(DIST=seq(1,ceiling(max(tmp_df2$DIST))))), col="mediumblue")

#Coloring POP
lim <- par("usr")
for(i in 1:nrow(SUBPOPs)){
rect(SUBPOPs$Min[i], lim[3], SUBPOPs$Max[i], lim[4], border = rgb(t(col2rgb(SUBPOPs$Col[i])/255), alpha=0.25), col = rgb(t(col2rgb(SUBPOPs$Col[i])/255), alpha=0.25))
}
axis(3, at=apply(SUBPOPs[,3:4],1,mean), labels=SUBPOPs$POP, tick=F, las=2)
mtext(paste("AIC =",round(min(AIC_results[,2]),2)), 1, adj=0.95, line=-2)
mtext(paste("Period =",AIC_results[AIC_results[,2] == min(AIC_results[,2]),1], "km"), 1, adj=0.95, line=-1)

#Second best AIC
PER <- (2 * 3.1415)/413
fit <- lm(MEM20 ~ sin(PER*DIST) + cos(PER*DIST), data=tmp_df2)
fit
plot(MEM20 ~ DIST, data=tmp_df2, pch=21, bg="red4", xlab="Distance (km)", ylab="Dist MEM20")
lines(seq(1,ceiling(max(tmp_df2$DIST))), predict.lm(fit, data.frame(DIST=seq(1,ceiling(max(tmp_df2$DIST))))), col="mediumblue")

#Coloring POP
lim <- par("usr")
for(i in 1:nrow(SUBPOPs)){
rect(SUBPOPs$Min[i], lim[3], SUBPOPs$Max[i], lim[4], border = rgb(t(col2rgb(SUBPOPs$Col[i])/255), alpha=0.25), col = rgb(t(col2rgb(SUBPOPs$Col[i])/255), alpha=0.25))
}
axis(3, at=apply(SUBPOPs[,3:4],1,mean), labels=SUBPOPs$POP, tick=F, las=2)
mtext(paste("AIC =",round(AIC_results[413,2],2)), 1, adj=0.95, line=-2)
mtext(paste("Period =",413, "km"), 1, adj=0.95, line=-1)

dev.off()

#Modeling sine wave for Adult MEMs
#Setting up object for sine wave modeling
tmp_df2 <- data.frame(Cells=rownames(Adult_MEM), MEM11=as.numeric(Adult_MEM$MEM11), MEM20=as.numeric(Adult_MEM$MEM20))
tmp_df2$DIST <- as.numeric(DIST_m[match(tmp_df2[,1], DIST_m[,1]),2])
tmp_df2 <- tmp_df2[order(tmp_df2$DIST),]

boxplot(tmp_df2$MEM11, plot=F)$stats[c(2,4)]
boxplot(tmp_df2$MEM20, plot=F)$stats[c(2,4)]
}

#Dist MEM 11
AIC_results <- NULL
for(i in seq(1, 5000, 1)){
PER <- (2 * 3.1415)/i
AIC_results <- rbind(AIC_results, c(i, AIC(lm(MEM11 ~ sin(PER*DIST) + cos(PER*DIST), data=tmp_df2))))
}
par(mfrow=c(3,1))
hist(AIC_results[,2], breaks=100, col="red4")
plot(AIC_results, xlab="Period", ylab="AIC", pch=19)
AIC_results[AIC_results[,2] == min(AIC_results[,2]),]
AIC_results[AIC_results[,2] == min(AIC_results[1:500,2]),]
AIC_results[AIC_results[,2] == min(AIC_results[550:1000,2]),]

#Plotting sine wave data for review
PER <- (2 * 3.1415)/(AIC_results[which(AIC_results[,2] == min(AIC_results[,2])),1])
fit <- lm(MEM11 ~ sin(PER*DIST) + cos(PER*DIST), data=tmp_df2)
fit
plot(MEM11 ~ DIST, data=tmp_df2, pch=21, bg="red4", xlab="Distance (km)", ylab="Adult MEM11")
lines(seq(1,ceiling(max(tmp_df2$DIST))), predict.lm(fit, data.frame(DIST=seq(1,ceiling(max(tmp_df2$DIST))))), col="mediumblue")

PER <- (2 * 3.1415)/134
fit <- lm(MEM11 ~ sin(PER*DIST) + cos(PER*DIST), data=tmp_df2)
fit
lines(seq(1,ceiling(max(tmp_df2$DIST))), predict.lm(fit, data.frame(DIST=seq(1,ceiling(max(tmp_df2$DIST))))), col="darkgreen", lty=3)

PER <- (2 * 3.1415)/986
fit <- lm(MEM11 ~ sin(PER*DIST) + cos(PER*DIST), data=tmp_df2)
fit
lines(seq(1,ceiling(max(tmp_df2$DIST))), predict.lm(fit, data.frame(DIST=seq(1,ceiling(max(tmp_df2$DIST))))), col="darkgreen", lty=3)

#Plotting sine wave data to save
tiff("Adult_MEM11.tiff", res=200, width=2500, height=2000)
par(mfrow=c(1,1))
#Best fit
PER <- (2 * 3.1415)/(AIC_results[which(AIC_results[,2] == min(AIC_results[,2])),1])
fit <- lm(MEM11 ~ sin(PER*DIST) + cos(PER*DIST), data=tmp_df2)
fit
plot(MEM11 ~ DIST, data=tmp_df2, pch=21, bg="red4", xlab="Distance (km)", ylab="Adult MEM11")
lines(seq(1,ceiling(max(tmp_df2$DIST))), predict.lm(fit, data.frame(DIST=seq(1,ceiling(max(tmp_df2$DIST))))), col="mediumblue")

#Plotting location colors
lim <- par("usr")
for(i in 1:nrow(SUBPOPs)){
rect(SUBPOPs$Min[i], lim[3], SUBPOPs$Max[i], lim[4], border = rgb(t(col2rgb(SUBPOPs$Col[i])/255), alpha=0.25), col = rgb(t(col2rgb(SUBPOPs$Col[i])/255), alpha=0.25))
}
axis(3, at=apply(SUBPOPs[,3:4],1,mean), labels=SUBPOPs$POP, tick=F, las=2)
mtext(paste("AIC =",round(min(AIC_results[,2]),2)), 1, adj=0.95, line=-2)
mtext(paste("Period =",AIC_results[AIC_results[,2] == min(AIC_results[,2]),1], "km"), 1, adj=0.95, line=-1)

dev.off()

#Dist MEM 20
AIC_results <- NULL
for(i in seq(1, 5000, 1)){
PER <- (2 * 3.1415)/i
AIC_results <- rbind(AIC_results, c(i, AIC(lm(MEM20 ~ sin(PER*DIST) + cos(PER*DIST), data=tmp_df2))))
}
par(mfrow=c(3,1))
hist(AIC_results[,2], breaks=100, col="red4")
plot(AIC_results, xlab="Period", ylab="AIC", pch=19)
AIC_results[AIC_results[,2] == min(AIC_results[,2]),]
AIC_results[AIC_results[,2] == min(AIC_results[100:1000,2]),]
AIC_results[AIC_results[,2] == min(AIC_results[400:1000,2]),]
PER <- (2 * 3.1415)/(AIC_results[which(AIC_results[,2] == min(AIC_results[,2])),1])
fit <- lm(MEM20 ~ sin(PER*DIST) + cos(PER*DIST), data=tmp_df2)
fit
plot(MEM20 ~ DIST, data=tmp_df2, pch=21, bg="red4", xlab="Distance (km)", ylab="Adult MEM20")
lines(seq(1,ceiling(max(tmp_df2$DIST))), predict.lm(fit, data.frame(DIST=seq(1,ceiling(max(tmp_df2$DIST))))), col="mediumblue")

PER <- (2 * 3.1415)/364
fit <- lm(MEM11 ~ sin(PER*DIST) + cos(PER*DIST), data=tmp_df2)
fit
lines(seq(1,ceiling(max(tmp_df2$DIST))), predict.lm(fit, data.frame(DIST=seq(1,ceiling(max(tmp_df2$DIST))))), col="darkgreen", lty=3)

PER <- (2 * 3.1415)/472
fit <- lm(MEM11 ~ sin(PER*DIST) + cos(PER*DIST), data=tmp_df2)
fit
lines(seq(1,ceiling(max(tmp_df2$DIST))), predict.lm(fit, data.frame(DIST=seq(1,ceiling(max(tmp_df2$DIST))))), col="darkgreen", lty=3)

tiff("Adult_MEM20_US.tiff", res=200, width=2500, height=2000)
par(mfrow=c(3,1))
#Plotting Best AIC
PER <- (2 * 3.1415)/(AIC_results[which(AIC_results[,2] == min(AIC_results[,2])),1])
fit <- lm(MEM20 ~ sin(PER*DIST) + cos(PER*DIST), data=tmp_df2)
fit
plot(MEM20 ~ DIST, data=tmp_df2, pch=21, bg="red4", xlab="Distance (km)", ylab="Adult MEM20")
lines(seq(1,ceiling(max(tmp_df2$DIST))), predict.lm(fit, data.frame(DIST=seq(1,ceiling(max(tmp_df2$DIST))))), col="mediumblue")

#Coloring by POP
lim <- par("usr")
for(i in 1:nrow(SUBPOPs)){
rect(SUBPOPs$Min[i], lim[3], SUBPOPs$Max[i], lim[4], border = rgb(t(col2rgb(SUBPOPs$Col[i])/255), alpha=0.25), col = rgb(t(col2rgb(SUBPOPs$Col[i])/255), alpha=0.25))
}
axis(3, at=apply(SUBPOPs[,3:4],1,mean), labels=SUBPOPs$POP, tick=F, las=2)
mtext(paste("AIC =",round(min(AIC_results[,2]),2)), 1, adj=0.95, line=-2)
mtext(paste("Period =",AIC_results[AIC_results[,2] == min(AIC_results[,2]),1], "km"), 1, adj=0.95, line=-1)

dev.off()

#Modeling sine wave for Larval MEMs
#Setting up object for sine wave modeling
tmp_df2 <- data.frame(Cells=rownames(Larval_MEM), MEM10=as.numeric(Larval_MEM$MEM10))
tmp_df2$DIST <- as.numeric(DIST_m[match(tmp_df2[,1], DIST_m[,1]),2])
tmp_df2 <- tmp_df2[order(tmp_df2$DIST),]

boxplot(tmp_df2$MEM10, plot=F)$stats[c(2,4)]

#Larval MEM 10
AIC_results <- NULL
for(i in seq(1, 5000, 1)){
PER <- (2 * 3.1415)/i
AIC_results <- rbind(AIC_results, c(i, AIC(lm(MEM10 ~ sin(PER*DIST) + cos(PER*DIST), data=tmp_df2))))
}
par(mfrow=c(3,1))
hist(AIC_results[,2], breaks=100, col="red4")
plot(AIC_results, xlab="Period", ylab="AIC", pch=19)
AIC_results[AIC_results[,2] == min(AIC_results[,2]),]
AIC_results[AIC_results[,2] == min(AIC_results[150:1000,2]),]

#Plotting sine wave data for review
PER <- (2 * 3.1415)/(AIC_results[which(AIC_results[,2] == min(AIC_results[,2])),1])
fit <- lm(MEM10 ~ sin(PER*DIST) + cos(PER*DIST), data=tmp_df2)
fit
plot(MEM10 ~ DIST, data=tmp_df2, pch=21, bg="red4", xlab="Distance (km)", ylab="Larval MEM10")
lines(seq(1,ceiling(max(tmp_df2$DIST))), predict.lm(fit, data.frame(DIST=seq(1,ceiling(max(tmp_df2$DIST))))), col="mediumblue")

PER <- (2 * 3.1415)/186
fit <- lm(MEM10 ~ sin(PER*DIST) + cos(PER*DIST), data=tmp_df2)
fit
lines(seq(1,ceiling(max(tmp_df2$DIST))), predict.lm(fit, data.frame(DIST=seq(1,ceiling(max(tmp_df2$DIST))))), col="darkgreen", lty=3)

#Plotting sine wave data to save
tiff("Larval_MEM10.tiff", res=200, width=2500, height=2000)
par(mfrow=c(2,1))
#Best fit
PER <- (2 * 3.1415)/(AIC_results[which(AIC_results[,2] == min(AIC_results[,2])),1])
fit <- lm(MEM10 ~ sin(PER*DIST) + cos(PER*DIST), data=tmp_df2)
fit
plot(MEM10 ~ DIST, data=tmp_df2, pch=21, bg="red4", xlab="Distance (km)", ylab="Larval MEM10")
lines(seq(1,ceiling(max(tmp_df2$DIST))), predict.lm(fit, data.frame(DIST=seq(1,ceiling(max(tmp_df2$DIST))))), col="mediumblue")

#Plotting location colors
lim <- par("usr")
for(i in 1:nrow(SUBPOPs)){
rect(SUBPOPs$Min[i], lim[3], SUBPOPs$Max[i], lim[4], border = rgb(t(col2rgb(SUBPOPs$Col[i])/255), alpha=0.25), col = rgb(t(col2rgb(SUBPOPs$Col[i])/255), alpha=0.25))
}
axis(3, at=apply(SUBPOPs[,3:4],1,mean), labels=SUBPOPs$POP, tick=F, las=2)
mtext(paste("AIC =",round(min(AIC_results[,2]),2)), 1, adj=0.95, line=-2)
mtext(paste("Period =",AIC_results[AIC_results[,2] == min(AIC_results[,2]),1], "km"), 1, adj=0.95, line=-1)

dev.off()
