## IPandDLA.R ##
# 2 Mar 2016
# Santiago Mideros
# Consolidates IP and DLA data
# Output dip.csv and ddla.csv

setwd("~/workdir")
library(lme4)

ip = read.table("ip.txt", col.names=c("isolate", "subset", "block", 
                                          "plot", "plant","genotype", "phenotype")) #ip

dla = read.table("dla.txt", col.names=c("isolate", "subset", "block", 
                                          "genotype", "plot", "phenotype")) #dla

## Check to ensure data imported correctly
str(ip)
head(ip)
tail(ip)
# Making sure variables are correctly characterized
ip$phenotype = as.numeric(ip$phenotype)
ip$subset = as.factor(ip$subset)
ip$block = as.factor(ip$block)
ip$plot <- as.factor(ip$plot)
ip$plant<- as.factor(ip$plant)
summary(ip)
dla$phenotype = as.numeric(dla$phenotype)
dla$subset = as.factor(dla$subset)
dla$block = as.factor(dla$block)
dla$plot <- as.factor(dla$plot)
summary(dla)
head(dla)

## Means accross 5 plants
################################################################
ip.plot.mean <- aggregate(ip$phenotype ~ ip$isolate + ip$subset + 
                            ip$block + ip$genotype + ip$plot, list(ip$isolate), mean)
colnames(ip.plot.mean)<-c("isolate", "subset", "block", "genotype", "plot", "phenotype")

ip.plot.mean<-ip.plot.mean[order(ip.plot.mean$subset, ip.plot.mean$block, ip.plot.mean$plot, 
                                 ip.plot.mean$genotype),]
head(ip.plot.mean, 20)
dim(ip.plot.mean)

## IP Controls
#################################################################
controls<-subset(ip.plot.mean, isolate == "StNY001" | isolate == "St52B")
controls<-controls[order(controls$isolate, controls$subset, controls$block),]
controls
controls.mean <- aggregate(controls$phenotype ~ controls$isolate + controls$subset + 
                             controls$genotype, list(controls$isolate), mean)
controls.mean <- controls.mean[order(controls.mean$`controls$subset`, 
                                     controls.mean$`controls$isolate` ),]
controls.mean

## Cast
##################################################################
library(reshape2)
#str(ip.plot.mean)
ip.mean.c<-dcast(ip.plot.mean, isolate + subset + block
                        ~ genotype,  na.rm=T)
head(ip.mean.c, 20)
ip.mean.c

str(dla)
dla.c<-dcast(dla, isolate + subset + block
                 ~ genotype,  na.rm=T)
head(dla.c)
dla.c

## Compute Differences
##################################################################
ip.mean.c$phenotype <- ip.mean.c$DK - ip.mean.c$S
hist(ip.mean.c$phenotype, col="gold")
dla.c$phenotype <- dla.c$S - dla.c$DK
hist(dla.c$phenotype, col="gold")
## Save Output
##################################################################
write.csv(ip.mean.c, file="dip.csv")
write.csv(dla.c, file="ddla.csv")

