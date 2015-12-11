setwd("~/workdir")
library(qtl)
mapthis<-read.cross("csvr",,"populationfilt2.txt.csv",estimate.map=FALSE)
#summary(mapthis)
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
##Markers with segregation distortion
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
## Individual genotype frequencies
mapthis<-subset(mapthis, ind=c("-stny001","-st52B")) #discard parents
g<-pull.geno(mapthis)
gfreq<-apply(g,1,function(a) table(factor(a,levels=1:3)))
gfreq<-t(t(gfreq)/colSums(gfreq))
par(mfrow=c(1,3),las=1)
for(i in 1:3)
  plot(gfreq[i,],ylab="Genotype frequency", main=c("AA", "AB", "BB")[i], 
       ylim=c(0,1))
##############################################
# To estimate the recombination frequencies
mapthis<-est.rf(mapthis)
checkAlleles(mapthis, threshold=5)
rf<-pull.rf(mapthis)
lod<-pull.rf(mapthis, what="lod")
plot(as.numeric(rf), as.numeric(lod), xlab="Recombination fraction", ylab="LOD score")
#now form the linkage groups
lg<-formLinkageGroups(mapthis, max.rf=0.35, min.lod=3)
table(lg[,2])