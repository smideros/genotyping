opar<-par()
setwd("~/workdir")
library(qtl)
mapthis<-read.cross("csv",file="mapthis.csv",genotypes=c("AA","BB"),estimate.map=FALSE)
class(mapthis)[1]<-"dh"
summary(mapthis)

##LG22
mapthis<-orderMarkers(mapthis,chr=22)
pull.map(mapthis, chr=22)
#rip22<-ripple(mapthis, chr=22, window=7)
#summary(rip22)
#rip22lik<-ripple(mapthis,chr=22, window=3,method="likelihood",error.prob = 0.005)
#summary(rip22lik)

##LG20
mapthis<-orderMarkers(mapthis,chr=20)
pull.map(mapthis, chr=20)
#rip20<-ripple(mapthis, chr=20, window=5)
#summary(rip20)
#rip20lik<-ripple(mapthis,chr=20, window=5,method="likelihood",error.prob = 0.005)
#summary(rip20lik)

##LG19
mapthis<-orderMarkers(mapthis,chr=19)
pull.map(mapthis, chr=19)
#rip19<-ripple(mapthis, chr=19, window=6)
#summary(rip19)
#rip19lik<-ripple(mapthis,chr=19, window=4,method="likelihood",error.prob = 0.005)
#summary(rip19lik)
#Note last two markers might be switched but likelihood is not higher
#compareorder(mapthis,chr=19,c(1:3,4,6,5),error.prob = 0.01)

##LG18
mapthis<-orderMarkers(mapthis,chr=18)
pull.map(mapthis, chr=18)
#rip18<-ripple(mapthis, chr=18, window=7)
#summary(rip18)
#rip18lik<-ripple(mapthis,chr=18, window=7,method="likelihood",error.prob = 0.005)
#summary(rip18lik)
#compareorder(mapthis,chr=18,c(2,3,4,5,1,6,7),error.prob = 0.01)
#compareorder(mapthis,chr=18,c(1,5,2,4,3,6,7),error.prob = 0.01)
compareorder(mapthis,chr=18,c(1,5,4,3,2,6,7),error.prob = 0.01)
mapthis<-switch.order(mapthis,chr=18,c(1,5,4,3,2,6,7),error.prob = 0.005)

##LG17
mapthis<-orderMarkers(mapthis,chr=17)
pull.map(mapthis, chr=17)
#rip17<-ripple(mapthis, chr=17, window=7)
#summary(rip17)
#rip17lik<-ripple(mapthis,chr=17, window=4,method="likelihood",error.prob = 0.005)
#summary(rip17lik)
compareorder(mapthis,chr=17,c(1,3,2,4:8),error.prob = 0.01)
mapthis<-switch.order(mapthis,chr=17,c(1,3,2,4:8),error.prob = 0.005)

##LG16
mapthis<-orderMarkers(mapthis,chr=16)
pull.map(mapthis, chr=16)
#rip16<-ripple(mapthis, chr=16, window=7)
#summary(rip16)
#rip16lik<-ripple(mapthis,chr=16, window=4,method="likelihood",error.prob = 0.005)
#summary(rip16lik)

##LG15
mapthis<-orderMarkers(mapthis,chr=15)
pull.map(mapthis, chr=15)
#rip15<-ripple(mapthis, chr=15, window=7)
#summary(rip15)
#rip15lik<-ripple(mapthis,chr=15, window=4,method="likelihood",error.prob = 0.005)
#summary(rip15lik)
compareorder(mapthis,chr=15,c(1,3,2,4:7,9,8),error.prob = 0.01)
mapthis<-switch.order(mapthis,chr=15,c(1,3,2,4:7,9,8),error.prob = 0.005)

##LG14
mapthis<-orderMarkers(mapthis,chr=14)
pull.map(mapthis, chr=14)
#rip14<-ripple(mapthis, chr=14, window=7)
#summary(rip14)
#rip14lik<-ripple(mapthis,chr=14, window=4,method="likelihood",error.prob = 0.005)
#summary(rip14lik)
compareorder(mapthis,chr=14,c(2,1,3:5,7,6,8,9),error.prob = 0.01)
mapthis<-switch.order(mapthis,chr=14,c(2,1,3:5,7,6,8,9),error.prob = 0.005)

##LG13
mapthis<-orderMarkers(mapthis,chr=13)
pull.map(mapthis, chr=13)
#rip13<-ripple(mapthis, chr=13, window=7)
#summary(rip13)
#rip13lik<-ripple(mapthis,chr=13, window=4,method="likelihood",error.prob = 0.005)
#summary(rip13lik)
compareorder(mapthis,chr=13,c(1,2,4,3,5:7,9,8),error.prob = 0.01)
mapthis<-switch.order(mapthis,chr=13,c(1,2,4,3,5:7,9,8),error.prob = 0.005)

##LG12
mapthis<-orderMarkers(mapthis,chr=12)
pull.map(mapthis, chr=12)
#rip12<-ripple(mapthis, chr=12, window=7)
#rip12lik<-ripple(mapthis,chr=12, window=4,method="likelihood",error.prob = 0.005)
#summary(rip12lik)
#mapthis<-switch.order(mapthis,chr=12,c(10,9,8,7,6,5,4,3,2,1),error.prob = 0.005)

##LG11
mapthis<-orderMarkers(mapthis,chr=11)
pull.map(mapthis, chr=11)
#rip11<-ripple(mapthis, chr=11, window=7)
#summary(rip11)
#rip11lik<-ripple(mapthis,chr=11, window=4,method="likelihood",error.prob = 0.005)
#summary(rip11lik)
compareorder(mapthis,chr=11,c(2,3,1,4:11),error.prob = 0.01)
mapthis<-switch.order(mapthis,chr=11,c(2,3,1,4:11),error.prob = 0.005)

##LG21
mapthis<-orderMarkers(mapthis,chr=21)
pull.map(mapthis, chr=21)
#rip21<-ripple(mapthis, chr=21, window=7)
#summary(rip21)
#rip21lik<-ripple(mapthis,chr=21, window=4,method="likelihood",error.prob = 0.005)
#summary(rip21lik)
compareorder(mapthis,chr=21,c(2,1,3:11),error.prob = 0.01)
#mapthis<-switch.order(mapthis,chr=21,c(2,1,3:11),error.prob = 0.005)
#LG21 is wrong
mapthis<-subset(mapthis, -(chr=21))

##LG10
mapthis<-orderMarkers(mapthis,chr=10)
pull.map(mapthis, chr=10)
#rip10<-ripple(mapthis, chr=10, window=7)
#summary(rip10)
#rip10lik<-ripple(mapthis,chr=10, window=4,method="likelihood",error.prob = 0.005)
#summary(rip10lik)

##LG9
mapthis<-orderMarkers(mapthis,chr=9)
pull.map(mapthis, chr=9)
#rip9<-ripple(mapthis, chr=9, window=7)
#summary(rip9)
#rip9lik<-ripple(mapthis,chr=9, window=4,method="likelihood",error.prob = 0.005)
#summary(rip9lik)
#lg9<-subset(mapthis,chr=9)
#nmissing(lg9, what="mar")
mapthis<-drop.markers(mapthis,c("S1_1086596","S1_1082696","S1_1082510","S1_1082831","S1_1083343","S1_1083391","S1_1082786","S1_1079949"))
#rip9<-ripple(mapthis, chr=9, window=7)
#summary(rip9)
#rip9lik<-ripple(mapthis,chr=9, window=4,method="likelihood",error.prob = 0.005)
#summary(rip9lik)
#compareorder(mapthis,chr=9,c(1,2,4,5,3,7,6),error.prob = 0.01)
mapthis<-switch.order(mapthis,chr=9,c(1,2,4,5,3,7,6),error.prob = 0.005)

##LG8
mapthis<-orderMarkers(mapthis,chr=8)
pull.map(mapthis, chr=8)
#map8<-unlist(pull.map(mapthis, chr=8))
#write.csv(map8,file="map8.csv")
#compareorder(mapthis,chr=8,c(1,2,4,5,3,6,7,9,8,10:15,17,16),error.prob = 0.01)
mapthis<-switch.order(mapthis,chr=8,c(1,2,4,5,3,6,7,9,8,10:15,17,16),error.prob = 0.005)
#rip8<-ripple(mapthis, chr=8, window=7)
#summary(rip8)
#rip8lik<-ripple(mapthis,chr=8, window=4,method="likelihood",error.prob = 0.005)
#summary(rip8lik)
#LG8 is inverted

##LG7
mapthis<-orderMarkers(mapthis,chr=7)
pull.map(mapthis, chr=7)
#map7<-unlist(pull.map(mapthis, chr=7))
#write.csv(map7,file="map7.csv")
#rip7lik<-ripple(mapthis,chr=7, window=4,method="likelihood",error.prob = 0.005)
#summary(rip7lik)
#rip7<-ripple(mapthis, chr=7, window=7)
#summary(rip7)
#compareorder(mapthis,chr=7,c(3,4,1,2,5:7,9,8,10,11,13,12,14:20),error.prob = 0.01)
#compareorder(mapthis,chr=7,c(3,4,2,1,5:7,9,8,10,11,13,12,14:20),error.prob = 0.01)
#compareorder(mapthis,chr=7,c(3,4,2,1,5:20),error.prob = 0.01)
#compareorder(mapthis,chr=7,c(3,4,1,2,5:20),error.prob = 0.01)
#compareorder(mapthis,chr=7,c(3,4,2,1,5:7,9,8,10:20),error.prob = 0.01)
#compareorder(mapthis,chr=7,c(3,4,2,1,5:7,9,8,10,11,13,12,14:20),error.prob = 0.01)
mapthis<-switch.order(mapthis,chr=7,c(3,4,2,1,5:7,9,8,10:20),error.prob = 0.005)

##LG6
mapthis<-drop.markers(mapthis,c("S2_3552054","S2_3491272","S2_3436629","S33_40146","S6_1972979","S6_1975615","S6_2139773","S6_2139975"))
mapthis<-orderMarkers(mapthis,chr=6)
pull.map(mapthis, chr=6)
#map6<-unlist(pull.map(mapthis, chr=6))
#write.csv(map6,file="map6.csv")
##rip6lik<-ripple(mapthis,chr=6, window=3,method="likelihood",error.prob = 0.005)
#summary(rip6lik)
#rip6<-ripple(mapthis, chr=6, window=7)
#summary(rip6)
#compareorder(mapthis,chr=6,c(1:9,11,10,14,12,13,15:18),error.prob = 0.01)
#compareorder(mapthis,chr=6,c(1:9,11,10,12:18),error.prob = 0.01)
#compareorder(mapthis,chr=6,c(1:9,11,10,13,14,12,15:18),error.prob = 0.01)
#compareorder(mapthis,chr=6,c(1:9,11,10,12,13,14,15:18),error.prob = 0.01)
#compareorder(mapthis,chr=6,c(1:11,13,14,12,15:18),error.prob = 0.01)
mapthis<-switch.order(mapthis,chr=6,c(1:11,13,14,12,15:18),error.prob = 0.005)

##LG5
mapthis<-orderMarkers(mapthis,chr=5)
pull.map(mapthis, chr=5)
#map5<-unlist(pull.map(mapthis, chr=5))
#write.csv(map5,file="map5.csv")
#rip5lik<-ripple(mapthis,chr=5, window=3,method="likelihood",error.prob = 0.005)
#summary(rip5lik)
#rip5<-ripple(mapthis, chr=5, window=7)
#summary(rip5)
#compareorder(mapthis,chr=5,c(1,2,4,5,3,6,8,7,9,19,20,24,25,26,23,22,21,18,17,16,15,14,13,12,10,11),error.prob = 0.01)
#compareorder(mapthis,chr=5,c(1:9,19,20,24,25,26,23,22,21,18,17,16,15,14,13,12,10,11),error.prob = 0.01)
#compareorder(mapthis,chr=5,c(1:9,26,25,24,23,22,21,20,19,18,17,16,15,14,13,12,10,11),error.prob = 0.01)
#compareorder(mapthis,chr=5,c(1,2,5,3,4,6,8,7,9,25,26,24,23,22,21,20,19,18,10,11,17,16,15,12,13,14),error.prob = 0.01)
mapthis<-switch.order(mapthis,chr=5,c(1,2,5,3,4,6,8,7,9,25,26,24,23,22,21,20,19,18,10,11,17,16,15,12,13,14),error.prob = 0.005)

##LG4
mapthis<-drop.markers(mapthis,c("S1_3743323"))
mapthis<-orderMarkers(mapthis,chr=4)
pull.map(mapthis, chr=4)
#map4<-unlist(pull.map(mapthis, chr=4))
#write.csv(map4,file="map4.csv")
#rip4lik<-ripple(mapthis,chr=4, window=3,method="likelihood",error.prob = 0.005)
#summary(rip4lik)
#rip4<-ripple(mapthis, chr=4, window=7)
#summary(rip4)
#compareorder(mapthis,chr=4,c(1:8,10,9,11:14,17,18,15,16,20,19,21:28),error.prob = 0.01)
#compareorder(mapthis,chr=4,c(1:8,10,9,11:28),error.prob = 0.01)
compareorder(mapthis,chr=4,c(1:8,10,9,11:18,20,19,21:28),error.prob = 0.01)
mapthis<-switch.order(mapthis,chr=4,c(1:8,10,9,11:18,20,19,21:28),error.prob = 0.005)
#mapthis<-switch.order(mapthis,chr=4,c(1:6,8,7,10,11,9,12:28),error.prob = 0.005)
#lg4<-subset(mapthis,chr=4)
#nmissing(lg4, what="mar")

##LG3
mapthis<-orderMarkers(mapthis,chr=3)
pull.map(mapthis, chr=3)
#map3<-unlist(pull.map(mapthis, chr=3))
#write.csv(map3,file="map3.csv")
#rip3lik<-ripple(mapthis,chr=3, window=3,method="likelihood",error.prob = 0.005)
#summary(rip3lik)
#rip3<-ripple(mapthis, chr=3, window=7)
#summary(rip3)
#compareorder(mapthis,chr=3,c(16:30,1:15),error.prob = 0.01)
#compareorder(mapthis,chr=3,c(16:19,21,20,22:30,1:15),error.prob = 0.01)
#compareorder(mapthis,chr=3,c(16:19,21,20,22,25,24,23,27,26,28:30,1:15),error.prob = 0.01)
#compareorder(mapthis,chr=3,c(16:19,21,20,22,25,24,23,27,26,28:30,1,2,4,3,7,5,6,8:15),error.prob = 0.01)
compareorder(mapthis,chr=3,c(16:19,21,20,22,25,24,23,27,26,28:30,1,2,4,3,7,5,6,8:10,12,11,13,15,14),error.prob = 0.01)
## genome order had higher probability and smaller length
mapthis<-switch.order(mapthis,chr=3,c(16:19,21,20,22,25,24,23,27,26,28:30,1,2,4,3,7,5,6,8:10,12,11,13,15,14),error.prob = 0.005)

##LG2
mapthis<-drop.markers(mapthis,c("S28_108792", "S20_158192", "S20_158193", "S20_652578"))
mapthis<-orderMarkers(mapthis,chr=2)
pull.map(mapthis, chr=2)
#map2<-unlist(pull.map(mapthis, chr=2))
#write.csv(map2,file="map2.csv")
#rip2lik<-ripple(mapthis,chr=2, window=3,method="likelihood",error.prob = 0.005)
#summary(rip2lik)
#rip2<-ripple(mapthis, chr=3, window=7)
#summary(rip2)
## formed by scaffolds 14 and then 20. delete marker S28_108792 placed in the middle of scaffold 14
#compareorder(mapthis,chr=2,c(1:15,39,38,37,36,35,33,32,31,30,16:29,34),error.prob = 0.01)
#compareorder(mapthis,chr=2,c(1:15,39,38,37,35,31,32,33,36,30,16:29,34),error.prob = 0.01)
#compareorder(mapthis,chr=2,c(3,4,5,2,1,6:11,13,14,12,15,39,38,37,35,31,32,33,36,30,16:29,34),error.prob = 0.01)
#compareorder(mapthis,chr=2,c(3,4,5,2,1,6:11,13,14,12,15,39,38,37,35,31,32,33,36,30,19,21,20,22,25,23,17,18,16,24,26:29,34),error.prob = 0.01)
## scaffold 14 arranges in genome order resulting in significantly smaller map but lower LOD
#mapthis<-switch.order(mapthis,chr=2,c(3,4,5,2,1,6:11,13,14,12,15,39,38,37,35,31,32,33,36,30,16:29,34),error.prob = 0.005)
#lg2<-subset(mapthis,chr=2)
#nmissing(lg2, what="mar")
## markers S20_158194, 92 and 93 are extremely close and have many missing data. delete the latter two
## marker S20_652578 has 129 missing genotypes and is grossly missplaced
#compareorder(mapthis,chr=2,c(1:15,36,35,34,33,32,31,30,29,28,17:22,16,23:27),error.prob = 0.01)
#compareorder(mapthis,chr=2,c(2,3:5,1,6:15,36,35,34,33,32,31,30,29,28,17:22,16,23:27),error.prob = 0.01)
compareorder(mapthis,chr=2,c(2,3:5,1,6:15,36,35,34,32,30,31,29,33,28,19,17,18,20:22,16,23:27),error.prob = 0.01)
## still get smaller map but lower LOD with genome order
mapthis<-switch.order(mapthis,chr=2,c(2,3:5,1,6:15,36,35,34,32,30,31,29,33,28,19,17,18,20:22,16,23:27),error.prob = 0.005)

##LG1
mapthis<-drop.markers(mapthis,c("S11_146039", "S5_1492149", "S208_1253", "S14_1286364", "S272_827"))
mapthis<-drop.markers(mapthis,c("S8_406421","S8_406629","S8_442250","S8_406756","S8_437498", "S8_1207988"))
mapthis<-drop.markers(mapthis,c("S21_648078", "S25_4333", "S25_5687", "S23_102252"))
mapthis<-drop.markers(mapthis,c("S8_1494654"))
mapthis<-drop.markers(mapthis,c("S26_170722", "S26_143929"))
mapthis<-drop.markers(mapthis,c("S47_4029"))
mapthis<-orderMarkers(mapthis,chr=1)
pull.map(mapthis, chr=1)
map1<-unlist(pull.map(mapthis, chr=1))
write.csv(map1,file="map1.csv")
#rip1lik<-ripple(mapthis,chr=1, window=3,method="likelihood",error.prob = 0.005)
#summary(rip1lik)
#rip1<-ripple(mapthis, chr=1, window=7)
#summary(rip1)
## forms lg with scaffolds 8 , 12 and others
## will delete some of the others: S11_146039, S5_1492149, S208_1253, S14_1286364, S272_827
## map expanded substantially 
#lg1<-subset(mapthis,chr=1)
#nmissing(lg1, what="mar")
# markers are very close to S8_406996: S8_406421   S8_406629   S8_442250 S8_406756   S8_437498
# markers are very close so delete S8_1207988
# markers have many missing data and are messsing up the map: S21_648078, S25_4333, S25_5687, S23_102252
# markers has significant missing data and might be causing wrong alignment: S8_1494654
## try order: S8 ascending, S26, S27, S12 descending, S2, S23, S10 descening
#compareorder(mapthis,chr=1,c(10,12,11,13,14,15,16,17,18,20,19,21:24,7,5,6,2,3,4,8,9,1,26,25,31,28,27,30,29,32:37,39,38,40,42,41),error.prob = 0.01)
#compareorder(mapthis,chr=1,c(10:24,7,6,5,4,3,2,1,8,9,25:37,39,38,40,42,41),error.prob = 0.01)
# delete S26 which splits S8
## try order: S8 ascending, S47, S23, S10 ascending, S2, S12 ascending
#compareorder(mapthis,chr=1,c(1:16,18:22,17,23:40),error.prob = 0.01)
#compareorder(mapthis,chr=1,c(1:16,18,20,19,21,22,17,23:40),error.prob = 0.01)
#compareorder(mapthis,chr=1,c(1:18,20,19,21,22,23:40),error.prob = 0.01)
# best is with S47 within S8, now delete S47
compareorder(mapthis,chr=1,c(9:23,8,6,7,3,5,4,2,1,38,39,37,35,36,34,33,32,31,28,29,26,27,30,24,25),error.prob = 0.01)
#compareorder(mapthis,chr=1,c(9,11,10,12,13:19,21,20,22,23,8,6,7,3,5,4,2,1,38,39,37,35,36,34,33,32,31,28,29,26,27,30,24,25),error.prob = 0.01)
mapthis<-switch.order(mapthis,chr=1,c(9:23,8,6,7,3,5,4,2,1,38,39,37,35,36,34,33,32,31,28,29,26,27,30,24,25),error.prob = 0.005)

######################
# Map Characteristics 
summaryMap(mapthis)
plotMap(mapthis)
write.cross(mapthis, "csv", "seto_genmap2")
#plotMap(mapthis, chr=4, show.marker.names=TRUE, show.centimorgans=TRUE)
#plotRF(mapthis, chr=4)
