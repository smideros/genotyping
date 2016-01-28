opar<-par()
setwd("~/workdir")
library(qtl)
mapthis2<-read.cross("csv",file="seto_genmap2_6Jan16.csv",genotypes=c("AA","BB"))
class(mapthis2)[1]<-"dh"
summary(mapthis2)

######################
# Map Characteristics 
summaryMap(mapthis2)
plotMap(mapthis2)
plotRF(mapthis2)
plotRF(mapthis2, chr=9:22)
plotRF(mapthis2, chr=1:8)

##LG8
plotMap(mapthis2, chr = 8, show.marker.names = TRUE)
plotRF(mapthis2, chr=8)
pull.map(mapthis2, chr=8)
compareorder(mapthis2,chr=8,c(9:17,1:8),error.prob = 0.01)
compareorder(mapthis2,chr=8,c(8,7,6,5,4,3,2,1,9:17),error.prob = 0.01)

##LG7
plotMap(mapthis2, chr = 7, show.marker.names = TRUE, show.centimorgans=TRUE)
plotRF(mapthis2, chr=7)
pull.map(mapthis2, chr=7)

##LG6
plotMap(mapthis2, chr = 6, show.marker.names = TRUE)
plotRF(mapthis2, chr=6)
pull.map(mapthis2, chr=6)

##LG5
plotMap(mapthis2, chr = 5, show.marker.names = TRUE)
plotRF(mapthis2, chr=5)
pull.map(mapthis2, chr=5)
compareorder(mapthis2,chr=5,c(1,2,6,5,4,3,7:26))

##LG4
plotMap(mapthis2, chr = 4, show.marker.names = TRUE)
plotRF(mapthis2, chr=4)
pull.map(mapthis2, chr=4)

##LG3
plotMap(mapthis2, chr = 3, show.marker.names = TRUE)
plotRF(mapthis2, chr=3)
pull.map(mapthis2, chr=3)

##LG2
plotMap(mapthis2, chr = 2, show.marker.names = TRUE)
plotRF(mapthis2, chr=2)
pull.map(mapthis2, chr=2)

###########################
##LG1

plotMap(mapthis2, chr = 1, show.marker.names = TRUE)
plotRF(mapthis2, chr=1)
pull.map(mapthis2, chr=1)

## scaffold2 shows up in 1,7,10,14
plotMap(mapthis2, chr = c(1,7,10,14), show.marker.names = TRUE)
# sc2 markers from 7 were deleted before.
# drop sc2 markers from lg1.  High rate of missing data.
mapthis2<-drop.markers(mapthis2,c("S2_935050","S2_2397749"))
lg1<-subset(mapthis2,chr=1)
nmissing(lg1, what="mar")
# move sc2 markers from lg14 into lg10
lg14<-(markernames(mapthis2, chr=14))
for(mar in lg14){
mapthis2<-movemarker(mapthis2, mar, 10)}
plotMap(mapthis2, chr = c(1,7,10,14), show.marker.names = TRUE)

## LG10
#mapthis2<-orderMarkers(mapthis2, chr=10)
pull.map(mapthis2, chr=10)
#map10<-unlist(pull.map(mapthis2, chr=10))
#write.csv(map10,file="map10.csv")
compareorder(mapthis2,chr=10,c(1:14,16,15,17:23),error.prob = 0.01)
compareorder(mapthis2,chr=10,c(23:17,15,16,14:1),error.prob = 0.01)
mapthis2<-switch.order(mapthis2,chr=10,c(23:17,15,16,14:1),error.prob = 0.005)
plotMap(mapthis2, chr = 10, show.marker.names = TRUE)
plotRF(mapthis2, chr=10)
pull.map(mapthis2, chr=10)

## scaffold10 shows up in 1, 13 and 18
plotMap(mapthis2, chr = c(1,13,18), show.marker.names = TRUE)
##scaffold1 shows up in 4,9,13
plotMap(mapthis2, chr = c(4,9,13), show.marker.names = TRUE)
## are scaffolds 1,6,24,10 and 12 a single chr?
plotMap(mapthis2, chr = c(1,4,9,13,18,20), show.marker.names = TRUE)
## move scaffold10 from lg1 to lg13 didn't work.
#sc10from1<-c("S10_857699", "S10_858308", "S10_858408", "S10_858557")
#for(mar in sc10from1){
#  mapthis2<-movemarker(mapthis2,mar,13)
#}
mapthis2<-drop.markers(mapthis2,c("S10_857699", "S10_858308", "S10_858408", "S10_858557"))
## move sc10 from lg18 to lg13 didn't work
#lg18<-(markernames(mapthis2, chr=18))
#for(mar in lg18){
#  mapthis2<-movemarker(mapthis2, mar, 13)
#}
plotMap(mapthis2, chr = c(1,4,9,13,18,20), show.marker.names = TRUE)

## clean LG13
#pull.map(mapthis2, chr=13)
#new13<-est.map(mapthis2, chr=13)
#plotMap(mapthis2, new13, chr = 13, show.marker.names = TRUE)
#plotRF(mapthis2, chr=13)
#pull.map(mapthis2, chr=13)
#compareorder(mapthis2,chr=13,c(1:7,17:20,8:16),error.prob = 0.01)
#compareorder(mapthis2,chr=13,c(1:9,17:20,10:16),error.prob = 0.01)
#mapthis<-orderMarkers(mapthis2,chr=13)
#lg13<-subset(mapthis2,chr=13)
#nmissing(lg13, what="mar")


## clean again LG1
#mapthis<-orderMarkers(mapthis2,chr=1)
pull.map(mapthis2, chr=1)
lg1markers<-markernames(mapthis2,chr=1)
mapthis2<-drop.markers(mapthis2,c("S8_200988", "S8_1996532"))
newlg14<-(lg1markers[24:33])
for(mar in newlg14){
  mapthis2<-movemarker(mapthis2, mar, 14)
}
plotMap(mapthis2, chr = 1, show.marker.names = TRUE)
plotRF(mapthis2, chr=1)
pull.map(mapthis2, chr=1)

## move 9 to 4 didn't work
#lg9<-markernames(mapthis2,chr=9)
#for(mar in lg9){
#  mapthis2<-movemarker(mapthis2, mar, 4)
#}
#plotMap(mapthis2, chr = c(4,9), show.marker.names = TRUE)
#pull.map(mapthis2, chr=4)
#compareorder(mapthis2,chr=4,c(29:35,1:28),error.prob = 0.01)
#mapthis2<-switch.order(mapthis2,chr=4,c(29:35,1:28),error.prob = 0.005)
#plotMap(mapthis2, chr = 4, show.marker.names = TRUE)
#plotRF(mapthis2, chr=4)
#pull.map(mapthis2, chr=4)

############
## Calculate new genetic map
newmap<-est.map(mapthis2)
plotMap(mapthis2, newmap)
mapthis2<-replace.map(mapthis2,newmap)

## Order by length and output
summaryMap(mapthis2)
names(mapthis2$geno)<-c("5","4","3","2","6","8","10","12","17","1","11","7","14","16","15","19","13","20","18","21","9")
mapthis2$geno<-mapthis2$geno[c(10,4,3,2,1,5,12,6,21,7,11,8,17,13,15,14,9,19,16,18,20)]
summaryMap(mapthis2)
plotMap(mapthis2)
write.cross(mapthis2, "csv", "seto_genmap2_7Jan16")