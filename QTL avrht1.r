## QTL2.r
# Maps single QTL on seto2016avr.csv file
# Author: Santiago Mideros
# April 2016

opar<-par()
setwd("~/workdir")
library(qtl)
cr<-read.cross("csv",file="seto2016avr.csv",
                        genotypes=c("AA","BB"),estimate.map=FALSE)
class(cr)[1]<-"dh"
summary(cr)
plot.map(cr)
summaryMap(cr)
#plot(cr)
#plotInfo(cr) #missing information

####################
## Phenotypes
####################
par(mfrow=c(1,2))
for(i in 2:3)
  plot.pheno(cr, pheno.col = i)
pairs(jitter( as.matrix(cr$pheno[,2:3])), cex=0.6, las=1)
par(mfrow=c(1,2), las=1, cex=0.8)
means <- apply(cr$pheno[,2:3], 1, mean)
plot(means)
plot(sample(means), xlab= "Random index", ylab="means")

###################
## Map review
###################
pull.map(cr, chr = c(1,18))
plot.map(cr, chr = c(1,18), show.marker.names = T)
cr.rip<-ripple(cr, chr=1, window=7)
summary(cr.rip)
cr.ripl<-ripple(cr, chr=1, window=3, method="likelihood", error.prob = 0.0025)
summary(cr.ripl)
cr2<-orderMarkers(cr, chr=1)
pull.map(cr2, chr = 1)
# marker S2_3552054 does link  in LG1 at the end of the chromosome and not
# end of scaffold 2
compareorder(cr, chr = 1, c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,19,20,21,22,
                            23,24,18), error.prob = 0.0025)
cr.new <- switch.order(cr,chr = 1, c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,19,
                                     20,21,22,23,24,18), error.prob = 0.0025)
nm.cr<-est.map(cr.new, error.prob = 0.0025)
plot.map(nm.cr, cr)
summaryMap(nm.cr)
# with marker S2_3552054 at the beginning of LG1, the new map has a 1001.5 cM gap
cr18<-movemarker(cr, 'S2_3552054', 18, 50)
pull.map(cr18, chr = c(1,18))
cr18.rip<-ripple(cr18, chr=18, window=7)
summary(cr18.rip)
cr18.ripl<-ripple(cr18, chr=18, window=6, method="likelihood", error.prob = 0.0025)
summary(cr18.ripl)
compareorder(cr18, chr = 18, c(6,1,2,3,4,5), error.prob = 0.0001)
cr18.new <- switch.order(cr18,chr = 18, c(6,1,2,3,4,5), error.prob = 0.0025)
nm<-est.map(cr18.new, error.prob = 0.0025)
plot.map(cr18, nm)
summaryMap(cr18.new)
# Marker S2_3552054 in LG18 causes a gap of 230cM
pull.map(cr18.new, chr = 18)
plot.map(cr18.new, chr = c(1,18), show.marker.names = T)
# S2_3552054 links to either LG1 or LG18 but only after significant gaps

####################
## Interval Mapping
####################
p <- 100
cr <- calc.genoprob(cr, step=1, error.prob=0.0025)
cr.i <- sim.geno(cr, step=1, n.draws = 8192)
names(cr$pheno)
#all.perm <- scanone(cr.i, pheno.col = c(2:16), n.perm = p)
#summary(all.perm)
effect <- function(x, y) {
  effectplot(cr.i, pheno.col= x, mname1 = y)
  eff <- effectplot(cr.i, pheno.col= x, mname1 = y, draw= FALSE)
  eff[[1]][1] - eff[[1]][2]
} # trait, marker
vrtn <- function(x,y,z){
  qtl <- makeqtl(cr.i, c(x), c(y))
  m.qtl <- fitqtl(cr.i, qtl=qtl,pheno.col = z)
  summary(m.qtl)
} # chromosome, location, trait
par(mfrow=c(1,1))

## CI
ci <- pull.pheno(cr.i, "ci")
ci <- replace(ci, ci==0.5, NA)
cr$pheno <- cbind(cr$pheno, binary=ci)
ci.imp <- scanone(cr.i, pheno.col="ci", method = "imp")
ci.bin <- scanone(cr, pheno.col = "binary", model = "binary")
ht1.imp <- scanone(cr.i, pheno.col="ht1", method = "imp")
ht1.bin <- scanone(cr, pheno.col = "ht1", model = "binary")
summary(ci.imp, threshold = 2.83)
summary(ci.bin, threshold = 2.87)
summary(ht1.imp, threshold = 2.6)
summary(ht1.bin, threshold = 2.67)
plot(ci.imp, col = c("green"), ylab="lod", main = "avr genes", chr = c(1,4))
plot(ci.bin, col = c("blue"), add = T, lty=2, chr=c(1,4))
plot(ht1.imp, col = c("darkcyan"), add = T)
plot(ht1.bin, col = c("orange"), add = TRUE, lty=2)

## QTL 19
bayesint(ci.imp, 4)
lodint(ci.imp, 4)
effect("ci","4@145")
vrtn(4, 145, "ci")
## QTL 20
bayesint(ci.bin, 4)
lodint(ci.bin, 4)
effect("binary","4@142")
vrtn(4, 142, "binary")
## QTL 21
bayesint(ht1.imp, 1)
lodint(ht1.imp, 1)
effect("ht1","S2_3552054")
vrtn(1, 170, "ht1")
## QTL 22
bayesint(ht1.bin, 1)
lodint(ht1.bin, 1)
effect("ht1", "S2_3552054")
vrtn(1, 170, "ht1")

#####################
## Marker Regression 
#####################
par(mfrow=c(1,2))
plot.pxg(cr, pheno.col="ht1", "S2_3552054", infer = FALSE)
plot.pxg(cr, pheno.col="ci", "S2_3552054", infer = FALSE)
avr.mr <- scanone(cr, pheno.col=c(2,3), method = "mr")
summary(avr.mr, threshold = 2.67, format="allpheno")
plot.pxg(cr, pheno.col="ht1", "S8_2099705", infer = FALSE)
plot.pxg(cr, pheno.col="ci", "S8_2099705", infer = FALSE)
plot.pxg(cr, pheno.col="ht1", "S2_3552054", infer = FALSE, xlab="a")
plot.pxg(cr, pheno.col="ci", "S8_2099705", infer = FALSE)
