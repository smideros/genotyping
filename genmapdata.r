opar<-par()
setwd("~/workdir")
library(qtl)
setogenmap2<-read.cross("csv",file="seto_genmap2_26Jan16.csv",
                        genotypes=c("AA","BB"),estimate.map=FALSE)
class(setogenmap2)[1]<-"dh"
summary(setogenmap2)
plot.map(setogenmap2)
summaryMap(setogenmap2)

## Crossover counts
plot(countXO(setogenmap2), ylab="Number of crossovers")

## Genotyping errors
setogenmap2<-calc.errorlod(setogenmap2,error.prob = 0.0025)
print(toperr<-top.errorlod(setogenmap2, cutoff = 3))
plotGeno(setogenmap2, chr=2, ind=toperr$id[toperr$chr==2],cutoff = 3, include.xo=F)
plotGeno(setogenmap2, chr=5, ind=toperr$id[toperr$chr==5],cutoff = 3, include.xo=F)
plotGeno(setogenmap2, chr=6, ind=toperr$id[toperr$chr==6],cutoff = 3, include.xo=F)
plotGeno(setogenmap2, chr=19, ind=toperr$id[toperr$chr==19],cutoff = 3, include.xo=F)

## Missing genotype information
z<-plot.info(setogenmap2, method=c("both"),col=c("blue", "red"), 
          ylab="Proportion Missing", main="")
hist(nmissing(setogenmap2, what="mar"), breaks=50, xlim=c(0,200), xlab="No. missing genotypes", ylim = c(0,20), main="")

## By LG
plot.map(setogenmap2, chr=c(1:4), show.marker.names = T)
