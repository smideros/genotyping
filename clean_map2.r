opar<-par()
setwd("~/workdir")
library(qtl)
mapthis3<-read.cross("csv",file="seto_genmap2_7Jan16.csv",genotypes=c("AA","BB"),estimate.map=FALSE)
class(mapthis3)[1]<-"dh"
summary(mapthis3)
plot.map(mapthis3)

dropone<-droponemarker(mapthis3,error.prob=0.005)
par(mfrow=c(2,1))
plot(dropone,lod=1,ylim=c(-100,0))
plot(dropone,lod=2,ylab="Change in chromosome length")
summary(dropone, lod.column=2)

##LG7
par(opar)
plotMap(mapthis3, chr = 7, show.marker.names = TRUE)
mapthis3<-drop.markers(mapthis3, "S12_30970")

##LG9
plotMap(mapthis3, chr = 9, show.marker.names = TRUE)
lg9<-subset(mapthis3,chr=9)
nmissing(lg9, what="mar")
mapthis3<-drop.markers(mapthis3, "S12_798070")

##LG10
plotMap(mapthis3, chr = 10, show.marker.names = TRUE)
lg10<-subset(mapthis3,chr=10)
nmissing(lg10, what="mar")
mapthis3<-drop.markers(mapthis3, "S9_537806")

newmap<-est.map(mapthis3, error.prob=0.005)
mapthis3<-replace.map(mapthis3, newmap)
summaryMap(mapthis3)
plot.map(mapthis3)

## Count Crossovers
plot(countXO(mapthis3), ylab="Number of crossovers")
xomap3<-countXO(mapthis3)
#library("ggplot2")
#qplot(xomap3)
#ggplot(xomap3)
boxplot(xomap3)
# St136 has 33 crossovers. Median and mean are 15 with a min of 5
# delete St136
mapthis3<-subset(mapthis3, ind=(countXO(mapthis3)<30))
newmap2<-est.map(mapthis3, error.prob=0.005)
plotMap(mapthis3, newmap2)
mapthis3<-replace.map(mapthis3, newmap)

## Estimating genotyping error rate
loglik<-err<-c(0.001, 0.0025, 0.005, 0.0075, 0.01, 0.0125, 0.015, 0.0175, 0.02)
for(i in seq(along=err)){
  cat(i, "of", length(err), "\n")
  tempmap<-est.map(mapthis3, error.prob=err[i])
  loglik[i]<-sum(sapply(tempmap, attr, "loglik"))
}
lod<-(loglik-max(loglik))/log(10)
plot(err,lod,xlab="Genotyping error rate", xlim=c(0,0.02), ylab=expression(paste(log[10], " likelihood")))

## Genotyping errors
mapthis3<-calc.errorlod(mapthis3,error.prob = 0.0025)
print(toperr<-top.errorlod(mapthis3, cutoff = 9))
plotGeno(mapthis3, chr=2, ind=toperr$id[toperr$chr==2],cutoff = 9, include.xo=F)
plotGeno(mapthis3, chr=3, ind=toperr$id[toperr$chr==3],cutoff = 9, include.xo=F)
plotGeno(mapthis3, chr=4, ind=toperr$id[toperr$chr==4],cutoff = 9, include.xo=F)
plotGeno(mapthis3, chr=6, ind=toperr$id[toperr$chr==6],cutoff = 9, include.xo=F)
plotGeno(mapthis3, chr=8, ind=toperr$id[toperr$chr==8],cutoff = 9, include.xo=F)
plotGeno(mapthis3, chr=9, ind=toperr$id[toperr$chr==9],cutoff = 9, include.xo=F)
plotGeno(mapthis3, chr=10, ind=toperr$id[toperr$chr==10],cutoff = 9, include.xo=F)
plotGeno(mapthis3, chr=11, ind=toperr$id[toperr$chr==11],cutoff = 9, include.xo=F)
plotGeno(mapthis3, chr=12, ind=toperr$id[toperr$chr==12],cutoff = 9, include.xo=F)
plotGeno(mapthis3, chr=17, ind=toperr$id[toperr$chr==17],cutoff = 9, include.xo=F)
plotGeno(mapthis3, chr=19, ind=toperr$id[toperr$chr==19],cutoff = 9, include.xo=F)
#change 34 genotypies to NA for markers with genotype lod>9
mapthis4<-mapthis3
for(i in 1:nrow(toperr)){
  chr<-toperr$chr[i]
  id<-toperr$id[i]
  mar<-toperr$marker[i]
  mapthis4$geno[[chr]]$data[mapthis3$pheno$id==id,mar] <- NA
}
mapthis4<-calc.errorlod(mapthis4,error.prob = 0.0025)
print(toperr2<-top.errorlod(mapthis4, cutoff = 6))
#new map
newmap3<-est.map(mapthis4, error.prob=0.0025)
plotMap(mapthis4, newmap3)
mapthis4<-replace.map(mapthis4, newmap3)

## Revisit Segregation Distortion
gt<-geno.table(mapthis4, scanone.output = T)
#chisq.test(c(24,32))
thr<--(log10(0.05/totmar(mapthis4)))
par(mfrow=c(2,1))
plot(gt, ylab=expression(paste(-log[10]," P-value")))
abline(h=thr, lty=2, col="gray")
plot(gt, lod=3:4, ylab="Genotype frequency")
abline(h=c(0.25,0.5), lty=2, col="gray")
# list of distorted markers (n=9)
gt[gt$neglog10P>thr,]
#three markers in lg3 closely located are all distorted
plotMap(mapthis4, chr=3, show.marker.names = T)
par(mfrow=c(2,1))
plot(gt, ylab=expression(paste(-log[10]," P-value")), chr=3)
abline(h=thr, lty=2, col="gray")
plot(gt, lod=3:4, ylab="Genotype frequency", chr=3)
abline(h=c(0.25,0.5), lty=2, col="gray")
lg3<-subset(mapthis4,chr=3)
nmissing(lg3, what="mar")
pull.map(mapthis4, chr=3)
plot(gt, ylab=expression(paste(-log[10]," P-value")), chr=2)
abline(h=thr, lty=2, col="gray")
plot(gt, lod=3:4, ylab="Genotype frequency", chr=2)
abline(h=c(0.25,0.5), lty=2, col="gray")
par(opar)
#not conviced that segregation distoriton is not due to real lethal alleles

## Order by length and output
summaryMap(mapthis4)
names(mapthis4$geno)<-c("1","2","3","5","4","6","7","8","11","10","12","9","14","13","16","15","17","18","19","20","21")
mapthis4$geno<-mapthis4$geno[c(1,2,3,5,4,6,7,8,11,10,12,9,14,13,16,15,17,18,19,20,21)]
summaryMap(mapthis4)
plotMap(mapthis4)
write.cross(mapthis4, "csv", "seto_genmap2_8Jan16")