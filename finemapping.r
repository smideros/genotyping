# finemapping.r
# further analysis of avr gene regions
# Santiago Mideros
# 18 April 2016

opar<-par()
setwd("~/workdir")
library(qtl)
avr<-read.cross("csvr",,"finemapavr.csv",genotypes=c("A","B"),
                    estimate.map=FALSE)
class(avr)[1]<-"dh"
summary(avr)
plot.map(avr)
summaryMap(avr)

cr<-read.cross("csv",file="seto2016avr.csv",
               genotypes=c("AA","BB"),estimate.map=FALSE)
class(cr)[1]<-"dh"
head(cr$pheno)
head(avr$pheno)

plotMap(cr, chr = 1, show.marker.names = T)
plotMap(cr, chr = 4, show.marker.names = T)

plotMap(avr, show.marker.names = T)
plotMap(avr, chr = 8, show.marker.names = T)


avr$pheno <- merge(cr$pheno,avr$pheno, by=c("id"), all.y = T)

#####################
## Marker Regression 
par(mfrow=c(1,2))
avr2.mr <- scanone(avr, pheno.col=c(2,3), method = "mr")
summary(avr2.mr, threshold = 2.67, format="allpheno")
plot.pxg(avr, pheno.col="ht1", "S8_2099705", infer = FALSE)
plot.pxg(avr, pheno.col="ci", "S8_2099705", infer = FALSE)
plot.pxg(avr, pheno.col="ht1", "S2_3552054", infer = FALSE, xlab="a")
plot.pxg(avr, pheno.col="ci", "S8_2099705", infer = FALSE)
par(opar)
plot(avr2.mr, lodcolumn = c(1:2))
