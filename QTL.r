## QTL.r
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
## Phenotypes
####################
par(mfrow=c(3,2))
for(i in 2:7)
  plot.pheno(cr, pheno.col = i)
for(i in c(8:10, 12:14))
  plot.pheno(cr, pheno.col = i)
par(mfrow=c(1,2))
for(i in c(11, 15))
  plot.pheno(cr, pheno.col = i)
pairs(jitter( as.matrix(cr$pheno[,2:15])), cex=0.6, las=1)
par(mfrow=c(1,2), las=1, cex=0.8)
means <- apply(cr$pheno[,2:15], 1, mean)
plot(means)
plot(sample(means), xlab= "Random index", ylab="means")

####################
## Normality
####################
par(opar)
names(cr$pheno)

par(mfrow=c(3,2))
for(i in 2:7){
  qqnorm(cr$pheno[,i], main = names(cr$pheno[i]))
  qqline(cr$pheno[,i])}
for(i in c(8:10, 12:14)){
  qqnorm(cr$pheno[,i], main = names(cr$pheno[i]))
  qqline(cr$pheno[,i])}
par(mfrow=c(1,2))
for(i in c(11, 15)){
  qqnorm(cr$pheno[,i], main = names(cr$pheno[i]))
  qqline(cr$pheno[,i])}

shap <- lapply(cr$pheno[,2:15], shapiro.test)
shap <- sapply(shap, `[`, c("statistic","p.value"))
shap<-data.frame(t(shap))
shap
row.names(shap[ which(shap$p.value > 0.01),]) # only normal is di6

# Prepare transformations
pheno <- cr$pheno

# Add minimum value so that we don't have negatives on sp, dip, dld, ddl
mph <- sapply(pheno[,2:15], min)
mph <-mph[ which(mph <0)]
pheno$sp<-unlist(pheno["sp"]+abs(mph["sp"]))
pheno$dip<-unlist(pheno["dip"]+abs(mph["dip"]))
pheno$dld<-unlist(pheno["dld"]+abs(mph["dld"]))
pheno$ddl<-unlist(pheno["ddl"]+abs(mph["ddl"]))

# natural log transformation
nl<-log(pheno[,2:15]+1)
shapl <- lapply(nl, shapiro.test)
shapl <- sapply(shapl, `[`, c("statistic","p.value"))
shapl <- data.frame(t(shapl))
shapl
row.names(shapl[ which(shapl$p.value > 0.01),]) # didn't help

# square root transformation
sqr<-sqrt(pheno[,2:15])
shaps <- lapply(sqr, shapiro.test)
shaps <- sapply(shaps, `[`, c("statistic","p.value"))
shaps <- data.frame(t(shaps))
shaps
row.names(shaps[ which(shaps$p.value > 0.04),]) #improved sporulation to 0.0449
par(mfrow=c(2,2))
plot.pheno(cr, pheno.col = "sp")
hist(sqr$sp, main = "sqrt(sp)")
qqnorm(cr$pheno[,7], main = "sp")
qqline(cr$pheno[,7])
qqnorm(sqr$sp, main = "sqrt(sp)")
qqline(sqr$sp)

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
par(opar)


## abundance and melanine
am.imp <- scanone(cr.i, pheno.col=c(2:3,6), method = "imp")
plot(am.imp, lodcolumn=c(1,3), col = c("blue", "black"),
     ylab="LOD", chr=c(2,4,18), xlab="Linkage Group")
abline(h=2.83, col="blue",lty=2)
abline(h=2.67,lty=2)
summary(am.imp, threshold = 2.67, format="allpheno")
# QTL 1
bayesint(am.imp, 2, lodcolumn=3)
lodint(am.imp, 2, lodcolumn=3)
effect("mel","2@64")
vrtn(2,64,"mel")
# QTL 6
bayesint(am.imp, 18, lodcolumn=3)
lodint(am.imp, 18, lodcolumn=3)
effect("mel","S6_1000851")
vrtn(18,45.75962,"mel")
# QTL 2
bayesint(am.imp, 2, lodcolumn=1)
lodint(am.imp, 2, lodcolumn=1)
effect("ab3","2@69")
vrtn(2,69,"ab3")
# QTL 3
bayesint(am.imp, 4, lodcolumn=1)
lodint(am.imp, 4, lodcolumn=1)
effect("ab3","4@138")
vrtn(4,138,"ab3")
# QTL 4
plot(am.imp, lodcolumn=c(1:3), col = c("blue", "gold", "black"),
     ylab="LOD", main= "ab3 ab5 mel", chr = 18)
bayesint(am.imp, 18, lodcolumn=1)
lodint(am.imp, 18, lodcolumn = 1)
effect("ab3","S6_1000851")
vrtn(18, 45.75962, "ab3")
# QTL 5
bayesint(am.imp, 18, lodcolumn=2)
lodint(am.imp, 18, lodcolumn=2)
effect("ab5","S6_1000851")
vrtn(18, 45.75962, "ab5")

## non-parametric melanine
m.np <- scanone(cr.i, pheno.col = 6, model="np", ties.random = FALSE)
#m.np2 <- scanone(cr.i, pheno.col = 6, model="np", ties.random = TRUE)
#m.imp <- scanone(cr.i, pheno.col= 6, method = "imp")
#m.em <- scanone(cr, pheno.col= 6)
#plot(m.np, m.np2, m.imp, col = c("blue", "gold", "black"))
#plot(m.em, add= TRUE, lty=2)
#plot(m.np, m.np2, m.imp, chr = 2, col = c("blue", "gold", "black"))
#plot(m.em, add= TRUE, lty=2, chr= 2)
#plot(m.np, m.np2, m.imp, chr = 18, col = c("blue", "gold", "black"))
#plot(m.em, add= TRUE, lty=2, chr= 18)
#m.perm <- scanone(cr.i, pheno.col=c(6), n.perm=p)
summary(m.np, threshold = 2.67)
## QTL 7
bayesint(m.np, 2)
lodint(m.np,2)
effect("mel","S1_2655665")
vrtn(2, 63.79343, "mel")
## QTL 8
bayesint(m.np, 18)
lodint(m.np,18)
effect("mel","S6_1000851")
vrtn(18, 45.75962, "mel")

## diameter
d.imp <- scanone(cr.i, pheno.col=c(4,5), method = "imp")
plot(d.imp, lodcolumn=c(1:2), col = c("blue", "gold", "black"),
     ylab="lod", main= "di6 di13")
#d.perm <- scanone(cr.i, pheno.col=c(4,5), n.perm=p)
#summary(d.perm, alpha = 0.05)
summary(d.imp, threshold = 2.71, format="allpheno")
plot(d.imp, lodcolumn=c(1:2), col = c("blue", "gold", "black"),
     ylab="lod", main= "di6 di13", chr = 18)
# QTL 9
bayesint(d.imp, 1, lodcolumn=2)
lodint(d.imp, 1, lodcolumn=2)
effect("di13","1@133")
vrtn(1, 133, "di13")
# QTL 10
bayesint(d.imp, 4, lodcolumn=2)
lodint(d.imp, 4, lodcolumn=2)
effect("di13","S8_2095392")
vrtn(4, 143.1121, "di13")
# QTL 11
bayesint(d.imp, 18, lodcolumn=2)
lodint(d.imp, 18, lodcolumn=2)
effect("di13","S6_1000851")
vrtn(18, 45.759, "di13")
# QTL 12
bayesint(d.imp, 18, lodcolumn=1)
lodint(d.imp, 18, lodcolumn=1)
effect("di6","S6_1000851")
vrtn(18, 45.75962, "di6")

## Sporulation
s.imp <- scanone(cr.i, pheno.col=c(7,16), method = "imp")
plot(s.imp, lodcolumn=c(1:2), col = c("blue", "gold", "black"),
     ylab="lod", main= "sp sqrt(sp)")
abline(h=2.8, lty=2)
#s.perm <- scanone(cr.i, pheno.col=c(7,16), n.perm=p)
#summary(s.perm, alpha = 0.05)
summary(s.imp, threshold = 2.8, format="allpheno")
plot(s.imp, lodcolumn=c(1:2), col = c("blue", "gold", "black"),
     ylab="lod", main= "sp sqrt(sp)", chr = 18)
# QTL 13
bayesint(s.imp, 18, lodcolumn=1)
lodint(s.imp, 18, lodcolumn=1)
effect("sp","18@44")
vrtn(18, 44, "sp")
# QTL 14
bayesint(s.imp, 18, lodcolumn=2)
lodint(s.imp, 18, lodcolumn=2)
effect("ssp","18@43")
vrtn(18, 43, "ssp")

## IP
ip.imp <- scanone(cr.i, pheno.col=c(8:10), method = "imp")
plot(ip.imp, lodcolumn=c(1:3), col = c("blue", "gold", "black"),
     ylab="lod", main= "ipd, ips, dip")
abline(h=2.61, lty=2)
summary(ip.imp, threshold = 2.61, format="allpheno")
plot(ip.imp, lodcolumn=c(1,3), col = c("blue", "gold", "black"),
     ylab="lod", main= "ipd ips", chr = 4)
# QTL 15
bayesint(ip.imp, 4, lodcolumn=1)
lodint(ip.imp, 4, lodcolumn=1)
effect("ipd","4@145")
vrtn(4, 145, "ipd")
# QTL 16
bayesint(ip.imp, 4, lodcolumn=3)
lodint(ip.imp, 4, lodcolumn=3)
effect("dip","4@145")
vrtn(4, 145, "dip")

## DLA
# no QTL for dls
dla.imp <- scanone(cr.i, pheno.col=c(12,14:15), method = "imp")
plot(dla.imp, lodcolumn=c(1:3), col = c("blue", "gold", "black"), ylab="lod",
     main="dld ddl d814")
plot(dla.imp, lodcolumn=c(3), col = c("blue", "gold", "black"), ylab="lod",
     main="d814") # no QTL for dla814
summary(dla.imp, threshold = 2.73, format="allpheno")
# QTL 17
bayesint(dla.imp, 4, lodcolumn=1)
lodint(dla.imp, 4, lodcolumn=1)
effect("dld","4@146")
vrtn(4, 146, "dld")
# QTL 18
bayesint(dla.imp, 4, lodcolumn=2)
lodint(dla.imp, 4, lodcolumn=2)
effect("ddl","4@145")
vrtn(4, 145, "ddl")

## CI
ci.imp <- scanone(cr.i, pheno.col=c(11), method = "imp")
ci <- pull.pheno(cr.i, 11)
ci <- replace(ci, ci==0.5, NA)
cr$pheno <- cbind(cr$pheno, binary=ci)
ci.bin <- scanone(cr, pheno.col = "binary", model = "binary")
plot(ci.imp, col = c("green"), ylab="lod", main= "CI", chr = 4)
plot(ci.bin, col = c("blue"), add = TRUE, chr=4)
summary(ci.imp, threshold = 2.83)
summary(ci.bin, threshold = 2.87)
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
