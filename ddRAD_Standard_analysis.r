#!/bin/R

`# Red Snapper ddRAD analysis JP on new genome #`

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

#Add functions
#Plot multiple ggplot objects in the same window or save file
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
} #http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/

#Importing Data
strata <- read.csv(file = "Loci.TRSdp3g75MIg9p2.1Fdupmaf05.2020-12-15/indv.csv", header = TRUE)

vcf<-read.vcfR(file="Loci.TRSdp3g75MIg9p2.1Fdupmaf05.recode.vcf")
gen.vcf<-vcfR2genind(vcf)
strata(gen.vcf) <- strata[match(indNames(gen.vcf),strata$INDV),]
gen.vcf@strata$Lat <- as.numeric(as.matrix(gen.vcf@strata$Lat))
gen.vcf@strata$Lon <- as.numeric(as.matrix(gen.vcf@strata$Lon))
head(gen.vcf@strata)
rm(vcf)

gen <- read.genepop(file = "new.gen", ncode=3L, quiet = FALSE)
strata(gen) <- strata[match(indNames(gen),strata$INDV),]
gen@strata$Lat <- as.numeric(as.matrix(gen@strata$Lat))
gen@strata$Lon <- as.numeric(as.matrix(gen@strata$Lon))
head(gen@strata)

gen@other$lat<-gen@strata[,c("Lon","Lat")]
gen@other$lat[,1]<-as.numeric(as.character(gen@strata$Lon))
gen@other$lat[,2]<-as.numeric(as.character(gen@strata$Lat))	#As numeric
gen@other$xy<-gen@strata[,c("Lon","Lat")]	#As factor

#Adding related ID to strata
gen.indv <- substr(gen@strata$INDV, 1, 20)
gen.indv <- gsub("_","-",gen.indv)
gen.indv <- data.frame(cbind(as.character(gen@strata$INDV), gen.indv))
names(gen.indv) <- c("gen.ID", "relate.ID")

gen@strata$relate.ID <- gen.indv[match(indNames(gen),gen.indv$gen.ID),"relate.ID"]
head(strata(gen))

#Redefining Grouping
#Breaking some large groups
setPop(gen) <- ~POP
tmp_SubPOP <- as.character(pop(gen))
tmp_SubPOP[which(gen@strata$POP == "MX" & gen@strata$Lon > -93)] <- "CP_MX"
tmp_SubPOP[which(gen@strata$POP == "MX" & gen@strata$Lon < -93 & gen@strata$Lat < 21)] <- "VC_MX"
tmp_SubPOP[which(gen@strata$POP == "MX" & gen@strata$Lon < -93 & gen@strata$Lat > 21 & gen@strata$Lat < 25)] <- "N_MX"
tmp_SubPOP[which(gen@strata$POP == "MX" & gen@strata$Lat > 25)] <- "TX"
tmp_SubPOP[which(gen@strata$POP == "TX" & gen@strata$Lat > 27.5)] <- "W_LA"
tmp_SubPOP[which(gen@strata$POP == "LA" & gen@strata$Lon < -92)] <- "W_LA"
tmp_SubPOP[which(gen@strata$POP == "LA" & gen@strata$Lon > -92)] <- "E_LA"
tmp_SubPOP[which(gen@strata$POP == "LA" & gen@strata$Lon > -89 & gen@strata$Lon < -87.5)] <- "AL"
tmp_SubPOP[which(gen@strata$POP == "LA" & gen@strata$Lon > -87.5)] <- "PC"
tmp_SubPOP[which(gen@strata$POP == "PC" & gen@strata$Lon < -87.5)] <- "AL"
tmp_SubPOP[which(gen@strata$POP == "MG")] <- "Gulf_FL_N"
tmp_SubPOP[which(gen@strata$POP == "FL_Gulf" & gen@strata$Lat > 28)] <- "Gulf_FL_N"
tmp_SubPOP[which(gen@strata$POP == "FL_Gulf" & gen@strata$Lat < 28)] <- "Gulf_FL_S"
tmp_SubPOP[which(gen@strata$POP == "Keys")] <- "DT"
table(tmp_SubPOP)
gen@strata$SubPOP <- as.factor(tmp_SubPOP)

setPop(gen.vcf) <- ~POP
tmp_SubPOP <- as.character(pop(gen.vcf))
tmp_SubPOP[which(gen.vcf@strata$POP == "MX" & gen.vcf@strata$Lon > -93)] <- "CP_MX"
tmp_SubPOP[which(gen.vcf@strata$POP == "MX" & gen.vcf@strata$Lon < -93 & gen.vcf@strata$Lat < 21)] <- "VC_MX"
tmp_SubPOP[which(gen.vcf@strata$POP == "MX" & gen.vcf@strata$Lon < -93 & gen.vcf@strata$Lat > 21 & gen.vcf@strata$Lat < 25)] <- "N_MX"
tmp_SubPOP[which(gen.vcf@strata$POP == "MX" & gen.vcf@strata$Lat > 25)] <- "TX"
tmp_SubPOP[which(gen.vcf@strata$POP == "TX" & gen.vcf@strata$Lat > 27.5)] <- "W_LA"
tmp_SubPOP[which(gen.vcf@strata$POP == "LA" & gen.vcf@strata$Lon < -92)] <- "W_LA"
tmp_SubPOP[which(gen.vcf@strata$POP == "LA" & gen.vcf@strata$Lon > -92)] <- "E_LA"
tmp_SubPOP[which(gen.vcf@strata$POP == "LA" & gen.vcf@strata$Lon > -89 & gen.vcf@strata$Lon < -87.5)] <- "AL"
tmp_SubPOP[which(gen.vcf@strata$POP == "LA" & gen.vcf@strata$Lon > -89)] <- "PC"
tmp_SubPOP[which(gen.vcf@strata$POP == "PC" & gen.vcf@strata$Lon < -87.5)] <- "AL"
tmp_SubPOP[which(gen.vcf@strata$POP == "MG")] <- "Gulf_FL_N"
tmp_SubPOP[which(gen.vcf@strata$POP == "FL_Gulf" & gen.vcf@strata$Lat > 28)] <- "Gulf_FL_N"
tmp_SubPOP[which(gen.vcf@strata$POP == "FL_Gulf" & gen.vcf@strata$Lat < 28)] <- "Gulf_FL_S"
tmp_SubPOP[which(gen.vcf@strata$POP == "Keys")] <- "DT"
table(tmp_SubPOP)
gen.vcf@strata$SubPOP <- as.factor(tmp_SubPOP)

#Changing "Keys" to "DT"
tmp_data <- as.character(gen$strata$Region)
tmp_data[which(tmp_data=="Keys")] <- "DT"
gen$strata$Region <- as.factor(tmp_data)

#Changing "Keys" to "DT"
tmp_data <- as.character(gen$strata$Region2)
tmp_data[which(tmp_data=="Keys")] <- "DT"
gen$strata$Region2 <- as.factor(tmp_data)

#Adding a new category
setPop(gen) <- ~SubPOP
tmp_MEMPOP <- as.character(pop(gen))
tmp_MEMPOP[which(gen@strata$SubPOP == "CP_MX" & gen@strata$Lat > 21.3)] <- "AMEM11_High"
tmp_MEMPOP[which(gen@strata$SubPOP == "Gulf_FL_S" & gen@strata$Lat > 26)] <- "AMEM11_Low"
table(tmp_MEMPOP)
gen@strata$MEMPOP <- as.factor(tmp_MEMPOP)

write.table(gen@strata, file="gen_strata.csv", col.names=T, row.names=F, quote=F, sep=",")

#Setting color schemes

#Looking for relatedness between individuals
setPop(gen)<-~SubPOP
test.cov <- genomic_converter(gen, output = "related", parallel.core=10)
test.out <- coancestry(test.cov$related , dyadml =1)
save(test.out, file="Snapper.relatedness", compress=T)
#load("Snapper.relatedness")

#Removing one of the samples if there is a relationship between them
dups.rel2 <- test.out$relatedness[test.out$relatedness$dyadml>0.2,2:3]
rm.list2 <- list()

for(i in 1:nrow(dups.rel2)){
j <- as.character(dups.rel2[i,])
if(sum(rm.list2 %in% j)>0){print("skipping bad sample in comparison"); next}

tmp <- gen@strata[gen@strata$relate.ID %in% j, ]
tmp$Diff <- abs(as.numeric(as.character(tmp[ , "E.HOM."])) - as.numeric(as.character(tmp[ , "O.HOM."])))

Miss.tab <- table(as.numeric(as.character(tmp$N_MISS)))
Depth.tab <- table(as.numeric(as.character(tmp$MEAN_DEPTH)))

MIN.Miss <- min(as.numeric(names(Miss.tab)))
MAX.Depth <- max(as.numeric(names(Depth.tab)))

if(Miss.tab[names(Miss.tab)==MIN.Miss]==1 && Depth.tab[names(Depth.tab)==MAX.Depth]==1) {keep <- as.character(tmp[tmp$N_MISS==MIN.Miss,1])
} else if (length(which(tmp$Lib!="cat"))==1) {keep <- as.character(tmp$INDV[which(tmp$Lib!="cat")])
} else {keep <- as.character(tmp[tmp$Diff == min(tmp$Diff),1])}
print(keep)

rm.list2[[i]] <- as.character(tmp[grep(keep, tmp$INDV, invert=T),1])
rm(keep, tmp)
}

#Investigating samples which are related, but not dups
dups.rel2 <- test.out$relatedness[test.out$relatedness$dyadml>0.2,2:4]
table(dups.rel2$group)

dups.rel3 <- test.out$relatedness[which(test.out$relatedness$group %in% c("FLLA", "FLPC", "FLTX", "LAPC", "LATX", "PCTX") & test.out$relatedness$dyadml>0.2),c(2:3,11)]
dups.rel3

for(i in 1:nrow(dups.rel3)){print(gen@strata[gen@strata$relate.ID==dups.rel3[i,1] | gen@strata$relate.ID==dups.rel3[i,2],])}

rm.list2 <- c(rm.list2,"FL_19.I10.Lib8")

#Removing duplates
ind.count <- table(gen2@strata$Sample)
dup.names <- names(ind.count[which(ind.count>1)])
ind.count[which(ind.count>1)]

for(i in dup.names[1:4]){print(gen2@strata[gen2@strata$Sample==i,])}

test.out$relatedness[which(test.out$relatedness$ind2.id=="ATX-012.I4.Lib9" & test.out$relatedness$ind1.id=="ATX-012.I10.Lib6"), ]
test.out$relatedness[which(test.out$relatedness$ind1.id=="ATX-013.I10.Lib6" & test.out$relatedness$ind2.id=="ATX-013.I7.Lib9"), ]
test.out$relatedness[which(test.out$relatedness$ind1.id=="LA-001.I2.Lib1" & test.out$relatedness$ind2.id=="LA-001.I4.Lib2"), ]

rm.list2 <- c(rm.list2,"ATX-012.I4.Lib9","ATX-012.I10.Lib6", "ATX-013.I10.Lib6", "ATX-013.I7.Lib9", "LA-001.I2.Lib1", "LA-001.I4.Lib2")
write.table(rm.list2, "Related_samples_to_rm2.txt", col.names=F, row.names=F, quote=F)

#Removing duplicated individuals and Select only Adults
rm.list2 <- as.character(as.matrix(read.table("Related_samples_to_rm2.txt", head=F)))
Adults <- indNames(gen)[gen@strata$Class=="Adult"]

set.ind <- subset(indNames(gen), !indNames(gen) %in% rm.list2)
set.ind <- subset(set.ind, set.ind %in% Adults)
gen2 <- gen[set.ind, ]
loc2 <- gen2@other$lat

set.ind<-subset(indNames(gen.vcf), !indNames(gen.vcf) %in% rm.list2)
set.ind <- subset(set.ind, set.ind %in% Adults)
gen2.vcf <- gen.vcf[set.ind, ]

#Library Bias (Outflank)
setPop(gen2.vcf)<-~Lib
tmp<-gl.outflank(gen2.vcf, qthreshold = 0.1)	#Can use plot=F option if you want

outliers<-tmp$outflank$results[which(tmp$outflank$results[15]=="TRUE"),1]
out.list<-dput(matrix(unlist(strsplit(as.vector(outliers),"[.]")),ncol=2,byrow=T)[,1])
tmp<-unique(matrix(unlist(strsplit(as.vector(out.list),"[_]")),ncol=2,byrow=T))
tmp[,2] <- as.numeric(as.matrix(tmp[,2]))
out.list <- unique(tmp[,1])
gen.out.list<-(out.list)
write.table(gen.out.list, "Outflank_Library_Outliers.list", quote=F, col.names=F, row.names=F)

#Library Bias (Bayescan)
setPop(gen2)<-~Lib
writeGenPop(gen2, "SNP.TRS.Lib_all.gen", "Red Snapper by Library without dups")

#Converting to BS format
### {bash} ###
java -jar /usr/local/bin/PGDSpider_2.1.1.0/PGDSpider2-cli.jar -inputfile SNP.TRS.Lib_all.gen -inputformat GENEPOP -outputfile SNP.TRS.Lib_all.BS -outputformat BAYESCAN -spid /home/afields/bin/genepop_to_BS.spid

#Running Bayescan
### {bash} ###
mkdir Bayescan
bayescan SNP.TRS.Lib_all.BS -od ./Bayescan -all-trace -threads 20 -thin 100 -nbp 30 -pr_odds 100

#Analyzing results
### {bash} ###
cd Bayescan
head -n2 ../SNP.TRS.Lib_all.gen | tail -n 1 | sed 's/, /\n/g' > contigs
echo "Contigs" | cat - contigs > tmp; mv tmp contigs
paste contigs SNP.TRS.Lib_al_fst.txt | awk '{$2=""; print}' > fst.txt
awk 'NR==1{next;} $5>0.05{print $1}' fst.txt | sort > Bayescan_Outliers_Lib.list
cd ..

#Removing Library outliers
OF.out <- read.table("Outflank_Library_Outliers.list", head=F)
tmp.out <- read.csv("Bayescan/Bayescan_Outliers_Lib.list", head=F)
gen.out.list <- unique(c(as.character(as.matrix(OF.out)), as.character(as.matrix(tmp.out))))

set.loc<-locNames(gen2)[which(!(locNames(gen2) %in% gen.out.list))]
gen3<-gen2[, loc=set.loc]

tmp<-vector()
for(i in gen.out.list){tmp <- append(tmp,grep(i, locNames(gen2.vcf), fixed=F, value=T))}

set.loc<-subset(locNames(gen2.vcf), !locNames(gen2.vcf) %in% tmp)
gen3.vcf<-gen2.vcf[, loc=set.loc]

#Outlier Detection (Outflank)
setPop(gen3.vcf)<-~SubPOP
tmp<-gl.outflank(gen3.vcf, qthreshold = 0.1)
outliers<-tmp$outflank$results[which(tmp$outflank$results[15]=="TRUE"),1]
out.list<-dput(matrix(unlist(strsplit(as.vector(outliers),"[.]")),ncol=2,byrow=T)[,1])
tmp<-unique(matrix(unlist(strsplit(as.vector(out.list),"[_]")),ncol=2,byrow=T))
tmp[,2] <- as.numeric(as.matrix(tmp[,2]))
out.list <- unique(tmp[,1])
gen.out.list<-(out.list)
write.table(gen.out.list, "Outflank_SubPOP_Outliers.list", quote=F, col.names=F, row.names=F)

#Outlier Detection (Bayescan)
#Exporting data for Bayescan
setPop(gen3)<-~SubPOP
writeGenPop(gen3, "SNP.TRS.Final_all_SubPop.gen", "Red Snapper data without dups by SubPOP")

#Converting to BS format
### {bash} ###
java -jar /usr/local/bin/PGDSpider_2.1.1.0/PGDSpider2-cli.jar -inputfile SNP.TRS.Final_all_SubPop.gen -inputformat GENEPOP -outputfile SNP.TRS.Final_all_SubPop.BS -outputformat BAYESCAN -spid /home/afields/bin/genepop_to_BS.spid

#Running Bayescan
### {bash} ###
bayescan SNP.TRS.Final_all_SubPop.BS -od ./Bayescan -all-trace -threads 20 -thin 100 -nbp 30 -pr_odds 100

#Analyzing results
### {bash} ###
cd Bayescan
head -n2 ../SNP.TRS.Final_all_SubPop.gen | tail -n 1 | sed 's/, /\n/g' > contigs
echo "Contigs" | cat - contigs > tmp; mv tmp contigs
paste contigs SNP.TRS.Final_all_SubPo_fst.txt | awk '{$2=""; print}' > fst.txt
awk 'NR==1{next;} $5>0.05{print $1}' fst.txt | sort > Bayescan_Outliers_SubPOP.list
cd ..

#Finding loci with an extreme numbers of alleles
tmp_dat <- apply(gen3@tab, 2, function(x) sum(x,na.rm=T))
tmp_dat <- tmp_dat[tmp_dat > 0]
tmp_names <- matrix(unlist(strsplit(names(tmp_dat),"[.]")),ncol=2, byrow=T)[,1]
tmp_tab <- table(tmp_names)
Increased_SNPs <- names(boxplot(tmp_tab,plot=F, range=2)$out)

#Removing population outliers
OF.out <- read.table("Outflank_SubPOP_Outliers.list", head=F)
tmp.out <- read.csv("Bayescan/Bayescan_Outliers_SubPOP.list", head=F)
gen.out.list <- unique(c(as.character(as.matrix(OF.out)), as.character(as.matrix(tmp.out))))

set.loc<-locNames(gen3)[which(!(locNames(gen3) %in% gen.out.list) & !(locNames(gen3) %in% Increased_SNPs))]
gen.net <- gen3[, loc=set.loc]
set.loc<-locNames(gen3)[which(locNames(gen3) %in% gen.out.list | locNames(gen3) %in% Increased_SNPs)]
gen.out <- gen3[, loc=set.loc]

save(gen.net, file = "gen.net.gz", compress=T)
save(gen.out, file = "gen.out.gz", compress=T)

#Finding the "center" of each population
setPop(gen.net2) <- ~SubPOP
Avg_Coords <- data.frame(matrix(ncol=2, nrow=length(levels(pop(gen.net2)))))
colnames(Avg_Coords) <- c("Lat", "Lon")
rownames(Avg_Coords) <- levels(pop(gen.net2))
for(i in levels(pop(gen.net2))){
tmp_df <- gen.net2@strata[gen.net2@strata$SubPOP == i,c("Lat", "Lon")]
LAT <- mean(c(max(tmp_df$Lat), min(tmp_df$Lat)))
LON <- mean(c(max(tmp_df$Lon), min(tmp_df$Lon)))
Avg_Coords[i,]<- c(LAT, LON)
}
write.table(Avg_Coords, file="Avg_Pops.txt", sep="\t")

#PCA
{```{R}```
#Standard
X.net <- scaleGen(gen.net, NA.method="mean", scale=F)
pca1 <- dudi.pca(X.net,cent=FALSE,scale=FALSE,scannf=FALSE,nf=5000)
setPop(gen.net)<-~POP2
tiff("PCA-Neutral_data.tiff", res=300, width=3000, height=3000)
par(mfrow=c(1,1))
ade4::s.class(pca.net$li, pop(gen.net), col=c_POP5, cstar=0, axesell=F, clabel=0)
mtext(paste("x-axis variation: ",format(round(100*pca.net$eig[1]/sum(pca.net$eig),3),nsmall=3),"%",sep=""),1, line=0.8, adj=0.15, cex=1)
mtext(paste("y-axis variation: ",format(round(100*pca.net$eig[2]/sum(pca.net$eig),3),nsmall=3),"%",sep=""),1, line=2, adj=0.15, cex=1)
tmp_tab <- cbind(names(table(pop(gen.net))),table(pop(gen.net)))
tmp_leg <- paste(tmp_tab[,1], " (n=", tmp_tab[,2], ")", sep="")
legend(x=-5, y=-2.75,legend=tmp_leg, col=c_POP5, cex=1, pch=16, bty="n", ncol=2)
dev.off()

#Mapping samples
Gulf_map <- ggmap(get_stamenmap(bbox = c(left = -100, bottom = 18, right =-74, top = 38.3), color="bw", force=T, zoom = 8, maptype=c("terrain-background")))

p1.map <- Gulf_map + theme_classic() + geom_point(aes(x = Lon, y = Lat, fill=SubPOP, alpha =0.75), data=gen.net@strata, color="black", shape=21, stroke=0.2, size =2.5) + 
theme(legend.position=c(0.9,0.2)) + ggtitle("Adult Red Snapper Samples") + scale_fill_manual(values=c_POP4) + guides(alpha=FALSE) + labs(fill=NULL)
p1.map

ggsave("Filtered_Red_Snapper_samples_BW.tiff", p1.map, device="tiff")

#Rarefied Allelic Richness
setPop(gen.net2) <- ~SubPOP
#Calculated rarified allellic richness
Adult.ar <- allelic.richness(gen.net2)

#Make Allelic Richness table
Adult.ar_out <- Adult.ar$Ar
Adult.ar_out <- data.frame(rbind(Adult.ar_out, apply(Adult.ar_out[1:nrow(Adult.ar$Ar),], 2, sum), apply(Adult.ar_out[1:nrow(Adult.ar$Ar),], 2, mean), apply(Adult.ar_out[1:nrow(Adult.ar$Ar),], 2, sd)))
Adult.ar_out$Locus_Mean <- apply(Adult.ar_out[,1:ncol(Adult.ar$Ar)], 1, mean)
Adult.ar_out$Locus_SD <- apply(Adult.ar_out[,1:ncol(Adult.ar$Ar)], 1, sd)
rownames(Adult.ar_out)[(nrow(Adult.ar$Ar)+1):(nrow(Adult.ar$Ar)+3)] <- c("Pop_Sum", "Pop_Mean", "Pop_SD")

write.table(Adult.ar_out, file="Allelic_Richness.txt", col.names=T, row.names=T, quote=F)

## Exporting data for external analysis ##
setPop(gen.net)<-~SubPOP
writeGenPop(gen.net, "SNP.TRS.Final_neutral_SubPOP.gen", "Filtered Adult Red Snapper neutral loci data without dups by SubPOP")
