opar<-par()
setwd("~/workdir")
library(qtl)
mapthis<-read.cross("csvr",,"populationfilt2.txt.csv",genotypes=c("A","B"),estimate.map=FALSE)
class(mapthis)[1]<-"dh"
#mapthis<-convert2riself(mapthis)
summary(mapthis)
#plotMissing(mapthis)
#mapthis<-subset(mapthis,ind=(ntyped(mapthis)>50))
#drops markers that have less than 20 individuals genotyped
nt.bymar<-ntyped(mapthis,"mar")
todrop<-names(nt.bymar[nt.bymar<20])
mapthis<-drop.markers(mapthis, todrop)
summary(mapthis)
par(mfrow=c(1,2), las=1)
plot(ntyped(mapthis), ylab="No. typed markers", main="No. genotypes by individual") 
plot(ntyped(mapthis, "mar"), ylab="No. typed individuals",
     main="No. genotypes by marker")
cg<-comparegeno(mapthis)
par(mfrow=c(1,1))
hist(cg[lower.tri(cg)], breaks=seq(0,1,len=101),xlab="No. matching genotypes")
rug(cg[lower.tri(cg)])
#Looks like one pair of markers share 0% of markers
#To find the pair above:
wh<-which(cg<0.1, arr=TRUE)
wh<-wh[wh[,1]<wh[,2],]
wh
g<-pull.geno(mapthis)
table(g[205,],g[206,])
#individuals 205 and 206 are the parents
##Now find and drop duplicated markers
print(dup<-findDupMarkers(mapthis,exact.only = FALSE))
mapthis<-drop.markers(mapthis,unlist(dup))

## Markers with segregation distortion
gt<-geno.table(mapthis)
gt2<-gt[c("AA","BB")]
#chisq.test(c(28,53))
#head(gt)
chiFunction<-function(Geno){
  observed<-data.matrix(Geno)
  return(chisq.test(observed)$p.value)
}
chi1to1<-as.data.frame(apply(gt2,1,chiFunction))
colnames(chi1to1)<-c("chi1t1")
gt<-merge(gt,chi1to1,by=0)
distorted<-gt[gt$chi1t1<0.05/totmar(mapthis),]
distorted[with(distorted, order(chi1t1)),]
dim(distorted)
mapthis<-drop.markers(mapthis, distorted) # dropping 153 markers with segregation distortion
# the previous command doesn't work. Should have been:
# mapthis<-drop.markers(mapthis, distorted$Row.names)
# didn't actually discarted the distorted markers

## Individual genotype frequencies
mapthis<-subset(mapthis, ind=c("-stny001","-st52B")) #discard parents
g<-pull.geno(mapthis)
gfreq<-apply(g,1,function(a) table(factor(a,levels=1:2)))
gfreq<-t(t(gfreq)/colSums(gfreq))
par(mfrow=c(1,2),las=1)
for(i in 1:2)
  plot(gfreq[i,],ylab="Genotype frequency", main=c("AA", "BB")[i], 
       ylim=c(0,1))
summary(mapthis)

##############################################
## Linkage Map
##############################################
# To estimate the recombination frequencies
mapthis<-est.rf(mapthis)
checkAlleles(mapthis, threshold=5)
rf<-pull.rf(mapthis)
lod<-pull.rf(mapthis, what="lod")
par(mfrow=c(1,1))
#plot(as.numeric(rf), as.numeric(lod), xlab="Recombination fraction", ylab="LOD score")
#now form the linkage groups
lg<-formLinkageGroups(mapthis, max.rf=0.45, min.lod=6)
table(lg[,2])
mapthis<-formLinkageGroups(mapthis, max.rf=0.45, min.lod=6,reorgMarkers = TRUE)
#plotRF(mapthis,alternate.chrid = TRUE)
plotRF(mapthis,chr=1:41, alternate.chrid = TRUE)
plotRF(mapthis,chr=20:27, alternate.chrid = TRUE)

# select unlinked and throw them into a 'un' group
unlinked <- markernames(mapthis, chr=-(1:41)) # markers on chromosomes > 41
for(mar in unlinked)
  mapthis <- movemarker(mapthis, mar, "un")

# select first marker of lg 21 and plot rf and lod against all markers up to 27
mn21<-markernames(mapthis, chr=21)
rf<-pull.rf(mapthis)
lod<-pull.rf(mapthis, what="lod")
par(mfrow=c(2,1))
plot(rf, chr=20:27, mn21[1], bandcol="gray70", ylim=c(0,1), alternate.chrid=TRUE)
abline(h=0.5, lty=2)
plot(lod, chr=20:27, mn21[1], bandcol="gray70", alternate.chrid=TRUE)

# appears that lg 21, 23, 24, 25 and 26 should be linked
geno.crosstab(mapthis, mn21[1], mn21[2])
mn22<-markernames(mapthis, chr=22)
mn23<-markernames(mapthis, chr=23)
geno.crosstab(mapthis, mn21[1], mn22[1])
geno.crosstab(mapthis, mn21[1], mn23[1])

#move lg 23,24,25 and 26 into 21
new21<-markernames(mapthis, chr = 23:26)
for(mar in new21)
  mapthis <- movemarker(mapthis, mar, 21)

#select first marker of lg 22 and plot rf and lod against all markers
mn22<-markernames(mapthis, chr=22)
rf<-pull.rf(mapthis)
lod<-pull.rf(mapthis, what="lod")
par(mfrow=c(2,1))
plot(rf, chr=1:41, mn22[1], bandcol="gray70", ylim=c(0,1), alternate.chrid=TRUE)
abline(h=0.5, lty=2)
plot(lod, chr=1:41, mn22[1], bandcol="gray70", alternate.chrid=TRUE)

#select first marker of lg 27 and plot rf and lod against all markers
mn27<-markernames(mapthis, chr=27)
rf<-pull.rf(mapthis)
lod<-pull.rf(mapthis, what="lod")
par(mfrow=c(2,1))
plot(rf, chr=1:41, mn27[1], bandcol="gray70", ylim=c(0,1), alternate.chrid=TRUE)
abline(h=0.5, lty=2)
plot(lod, chr=1:41, mn27[1], bandcol="gray70", alternate.chrid=TRUE)
par(opar)

##Not sure that more linkage groups can be joined so will work with groups 1 to 22
softlinked <- markernames(mapthis, chr=-(1:22)) # markers on chromosomes > 22
for(mar in softlinked)
  mapthis <- movemarker(mapthis, mar, "un")
summary(mapthis)
write.cross(mapthis, "csv", "mapthis_29jan16",c(1:22))
