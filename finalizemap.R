opar<-par()
setwd("~/workdir")
library(qtl)
mapthis4<-read.cross("csv",file="seto_genmap2_8Jan16.csv",genotypes=c("AA","BB"),estimate.map=FALSE)
class(mapthis4)[1]<-"dh"

# Order by length and output
summaryMap(mapthis4)
names(mapthis4$geno)<-c("1","2","3","4","5","6","7","8","11","10","12","9","13","14","15","16","17","18","19","20","21")
mapthis4$geno<-mapthis4$geno[c(1,2,3,4,5,6,7,8,12,10,9,11,13,14,15,16,17,18,19,20,21)]
summary(mapthis4)
plot.map(mapthis4)

# Check recombination
mapthis4<- est.rf(mapthis4)
plot.rf(mapthis4, alternate.chrid=T)

# Check marker order
rip <- vector("list",nchr(mapthis4))
names(rip) <- names(mapthis4$geno)
for(i in names(mapthis4$geno)){
  rip[[i]] <- ripple(mapthis4, i, 7, verbose=FALSE)
}
dif.nxo <- sapply(rip, function(a) a[1,ncol(a)]-a[2,ncol(a)])
dif.nxo
any(dif.nxo>0)
# showing improvements for lgs 2,5,12,14

##LG2
plotMap(mapthis4, chr = 2, show.marker.names = TRUE)
pull.map(mapthis4, chr=2, as.table=T)
rip[[2]][2,]
lg2<-subset(mapthis4,chr=2)
nmissing(lg2, what="mar")
#drop markers with missing data
mapthis4<-drop.markers(mapthis4, c("S1_2650208", "S1_2650209", "S1_2650211"))
newmap2<-est.map(mapthis4, chr=2, error.prob=0.005)
plot.map(mapthis4, newmap2, chr=2)
mapthis4<-replace.map(mapthis4, newmap2)

##LG5
plotMap(mapthis4, chr = 5, show.marker.names = TRUE)
pull.map(mapthis4, chr=5, as.table=T)
rip[[5]][2,]
lg5<-subset(mapthis4,chr=5)
nmissing(lg5, what="mar")
#drop markers with missing data
mapthis4<-drop.markers(mapthis4, c("S20_180403", "S20_158194"))
newmap5<-est.map(mapthis4, chr=5, error.prob=0.005)
plot.map(mapthis4, newmap5, chr=5)
mapthis4<-replace.map(mapthis4, newmap5)

##LG12
plotMap(mapthis4, chr = 12, show.marker.names = TRUE)
pull.map(mapthis4, chr=12, as.table=T)
rip[[12]][2,]
lg12<-subset(mapthis4,chr=12)
nmissing(lg12, what="mar")
#drop markers with missing data
mapthis4<-drop.markers(mapthis4, "S12_634666")
newmap12<-est.map(mapthis4, chr=12, error.prob=0.005)
plot.map(mapthis4, newmap12, chr=12)
mapthis4<-replace.map(mapthis4, newmap12)

##LG14
plotMap(mapthis4, chr = 14, show.marker.names = TRUE)
pull.map(mapthis4, chr=14, as.table=T)
rip[[14]][2,]
lg14<-subset(mapthis4,chr=14)
nmissing(lg14, what="mar")
#drop markers with missing data
mapthis4<-drop.markers(mapthis4, c("S10_1650256", "S10_1650254"))
newmap14<-est.map(mapthis4, chr=14, error.prob=0.005)
plot.map(mapthis4, newmap14, chr=14)
mapthis4<-replace.map(mapthis4, newmap14)

# Check marker order 2
rip2 <- vector("list",nchr(mapthis4))
names(rip2) <- names(mapthis4$geno)
for(i in names(mapthis4$geno)){
  rip2[[i]] <- ripple(mapthis4, i, 7, verbose=FALSE)
}
dif.nxo2 <- sapply(rip2, function(a) a[1,ncol(a)]-a[2,ncol(a)])
dif.nxo2
any(dif.nxo2>0)

# Check marker order with likelihood
rip3 <- vector("list",nchr(mapthis4))
names(rip3) <- names(mapthis4$geno)
for(i in names(mapthis4$geno)){
  rip3[[i]] <- ripple(mapthis4, i, 3, method="likelihood", error.prob=0.001, verbose=F)
}
lod <- sapply(rip3, function(a) a[2,ncol(a)-1])
lod[lod>0]
lod[lod>1]
# LGs 1,2,3,4,5,6,7,8,9,10,11,12,13,17,19 have lod>0 only for  lg 1 LOD >1

##LG1
summary(rip3[[1]])
plotMap(mapthis4, chr = 1, show.marker.names = TRUE)
pull.map(mapthis4, chr=1, as.table=T)
rip3[[1]][2,]
lg1<-subset(mapthis4,chr=1)
nmissing(lg1, what="mar")
compareorder(mapthis4,chr=1,c(1:7,9,8,10:23),error.prob = 0.005)
mapthis4<-switch.order(mapthis4,chr=1,c(1:7,9,8,10:23),error.prob = 0.005)
newmap1<-est.map(mapthis4, chr=1, error.prob=0.005)
plot.map(mapthis4, newmap1, chr=1)
mapthis4<-replace.map(mapthis4, newmap1)

## Order by length and output
summaryMap(mapthis4)
plot.map(mapthis4)
names(mapthis4$geno)<-c("1","2","3","4","5","6","7","8","9","10","11","12","13","16","14","15","17","18","19","20","21")
mapthis4$geno<-mapthis4$geno[c(1,2,3,4,5,6,7,8,9,10,11,12,13,15,16,14,17,18,19,20,21)]
summaryMap(mapthis4)
plotMap(mapthis4)
write.cross(mapthis4, "csv", "seto_genmap2_26Jan16")
