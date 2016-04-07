## QTL figures.r
# Maps single QTL on seto2016.csv file
# Author: Santiago Mideros
# March 2016

opar<-par()
setwd("~/workdir")
library(qtl)
cr<-read.cross("csv",file="seto2016.csv",
                        genotypes=c("AA","BB"),estimate.map=FALSE)
class(cr)[1]<-"dh"
summary(cr)
plot.map(cr)
summaryMap(cr)
#plot(cr)
#plotInfo(cr) #missing information

####################
## Interval Mapping
####################
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
par(opar)


## abundance and melanine
am.imp <- scanone(cr.i, pheno.col=c(2,6), method = "imp")

## diameter and sporulation
ds.imp <- scanone(cr.i, pheno.col=c(5,7), method = "imp")

## planta
pl.imp <- scanone(cr.i, pheno.col=c(11,10,14), method = "imp")

## Plot
plot(pl.imp, lodcolumn=c(1:3), col = c("blue", "gold", "darkcyan"),
     ylab="LOD", chr=c(4), xlab= NULL)
plot(am.imp, lodcolumn=c(1), col = c("orange"), add = TRUE,
     chr=c(4))
plot(ds.imp, lodcolumn=c(1), col = c("violetred"), add = TRUE,
     chr=c(4))
abline(h=2.83, lty=2)


plot(ds.imp, lodcolumn=c(1), col = c("violetred"), lty = 1,
     chr=c(4), ylab="LOD")
plot(am.imp, lodcolumn=c(1), col = c("orange"), lty = 1, chr=c(4), add = T)
abline(h=2.83, lty=2)
