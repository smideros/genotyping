## qtlfigures2.R
## input: cross file seto2016.csv and qtl file qtl.csv
## Santiago Mideros
## 31March 2016

opar<-par()
setwd("~/workdir")
library(qtl)
library(ggplot2)
cr<-read.cross("csv",file="seto2016.csv",
               genotypes=c("AA","BB"),estimate.map=FALSE)
class(cr)[1]<-"dh"
summary(cr)
plot.map(cr)
summaryMap(cr)
crr <- subset(cr, chr=c("1","2","4","18"))
toplot <- read.csv("qtl.csv")

#################################################
## below are minor modifications to plot.qtlci.r Version. 2.1 
## by JT Lovell 
## downloaded from github 29March16
#################################################

roundUp <- function(x, nice=c(1,2,4,5,6,8,10)) {
  if(length(x) != 1) stop("'x' must be of length 1")
  10^floor(log10(x)) * nice[[which(x <= 10^floor(log10(x)) * nice)[[1]]]]
}
library("scales")
number_ticks <- function(n) {function(limits) pretty(limits, n)}

colnames(toplot)<-c("phe","chr","pos.left","pos","pos.right","numnames", "clr")
phe1<-unique(toplot[,c('phe','numnames')])
phe1<-phe1[order(phe1$numnames),]
phe1<-phe1[,1]
num1<-1:length(unique(toplot$phe))

toplot$pos<-as.numeric(as.character(toplot$pos))
toplot$pos.left<-as.numeric(as.character(toplot$pos.left))
toplot$pos.right<-as.numeric(as.character(toplot$pos.right))

#make dataset for marker ticks
marker.info<-data.frame(pull.map(crr, as.table=T))[,1:2]
colnames(marker.info)<-c("chr","pos")
seg.end<-1-1/(max(toplot$numnames)*.25)
seg.start<-seg.end-4/(max(toplot$numnames))
marker.info$seg.end<-seg.end
marker.info$seg.start<-seg.start
marker.info$chr<-as.numeric(paste(marker.info$chr))

#make dataset for axes info
chr.ys<-length(num1)
chr.len<-rep(as.numeric(chrlen(crr)),length(num1))
chrs<-rep(as.numeric(chrnames(crr)),length(num1))
chr.beg<-rep(0,length(chrs))
chr.ys<-rep(num1,each=nchr(crr))
chr.info<-data.frame(chr.len,chrs,chr.beg, chr.ys)
colnames(chr.info)[2]<-"chr"
toplot<-toplot[complete.cases(toplot),]

figure<-ggplot(marker.info)+
  geom_rug(aes(x=pos),alpha=1)+ # rug of markers
  geom_segment(aes(x=chr.beg, xend = chr.len, y=0, yend =0),
               alpha=1, data=chr.info)+
  geom_segment(aes(x=chr.beg, xend = chr.len, y=chr.ys, yend =chr.ys),
               alpha=.2, data=chr.info)+ #faint lines across each phenotype
  geom_segment(aes(x=pos.left, xend=pos.right, yend=numnames,y=numnames), 
               data=toplot, size=2, alpha=.5)+ # confidence intervals
  geom_point(aes(x=pos,y=numnames, color=clr),size=2.5, 
             data=toplot)+ #points for each QTL estimate
  facet_grid(.~chr, scale="free_x",space="free_x")+ #split by chromosome
  
  theme(
    panel.background = element_rect(fill = "white")
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank()
    ,strip.background = element_blank()
    ,axis.text.x = element_blank()
    ,axis.ticks = element_blank()
    ,axis.text.y = element_text(colour="black")
  ) +
 theme(legend.position="none")+
  #ggtitle("none")+
  #change axes
  scale_y_continuous("Phenotypes",labels=phe1, breaks=num1)+
  scale_x_continuous("")

#################################################
## show and save figure

figure
tiff("figure.tiff", width = 11, height = 2.8, units='in', res = 300)
figure
dev.off()
