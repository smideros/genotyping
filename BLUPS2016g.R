## BLUPS2016g.R ##
# 17 Feb 2016
# Santiago Mideros
# Calculates BLUEs from CLCs data set

setwd("~/workdir")
library(lme4)

#raw.data = read.table("sp.csv",
#                      col.names=c("isolate", "subset", "block", "plot", "phenotype")) #mel
#col.names=c("isolate", "subset", "block", "plot", "plant","genotype", "phenotype")) #ip
#col.names=c("isolate", "subset", "block", "dpi", "plate","plot", "phenotype")) #diam
#col.names=c("isolate", "subset", "block", "plot", "rep","counts", "phenotype")) #sp

raw.data = read.csv("ddla.csv") #dip, ddla

# genotype to include:
raw.data <- subset(raw.data, genotype == "DK")
raw.data<- droplevels(raw.data)

# dpi to include:
raw.data <- subset(raw.data, dpi == "13")
raw.data<- droplevels(raw.data)


## Check to ensure data imported correctly
str(raw.data)
head(raw.data)
tail(raw.data)
# Making sure variables are correctly characterized
raw.data$phenotype = as.numeric(raw.data$phenotype)
raw.data$subset = as.factor(raw.data$subset)
raw.data$block = as.factor(raw.data$block)
raw.data$plot <- as.factor(raw.data$plot)
raw.data$plant<- as.factor(raw.data$plant)
raw.data$dpi<- as.factor(raw.data$dpi)
raw.data$plate<- as.factor(raw.data$plate)
raw.data$rep<- as.factor(raw.data$rep)
raw.data$counts<- as.factor(raw.data$counts)

raw.data<-raw.data[order(raw.data$subset, raw.data$block, raw.data$plot),]
raw.data<-raw.data[order(raw.data$subset, raw.data$block),]
head(raw.data,20)

summary(raw.data)
#################################################################
## Exploratory Data: Controls
#################################################################
controls<-subset(raw.data, isolate == "StNY001" | isolate == "St52B")
controls<-controls[order(controls$isolate, controls$subset, controls$block),]
controls
controls.mean <- aggregate(controls$phenotype~controls$isolate+controls$subset, 
                           list(controls$isolate), mean)
controls.mean
#controls.mean <- aggregate(controls$phenotype ~ controls$isolate+ controls$subset +
#                             controls$rep, list(controls$isolate), mean)

#################################################################
## Exploratory Data Distribution
#################################################################
hist(raw.data$phenotype, col="gold")
# plot by block
raw.data$block<-factor(raw.data$block, c("1","2","3","4", "5","6","7","8","9","10","11","12","13","14","15","16",
                       "17","18","19","20","21","22","23","24","33d","25","26","27","28","29",
                       "30","31","32"))
boxplot(raw.data$phenotype~raw.data$block, xlab="Block", ylab="raw.data", main="Trait by Block", col="lightblue")  
# Box plot by subset
boxplot(raw.data$phenotype~raw.data$subset, xlab="subset", ylab="raw.data", main="Trait by subset", col="lightblue")  

#################################################################
## Models
#################################################################
#raw.data.lm.ssbbp = lm(raw.data$phenotype~ raw.data$isolate + raw.data$subset + 
#                         raw.data$subset:raw.data$block + raw.data$subset:raw.data$block:raw.data$rep)
#anova(raw.data.lm.ssbbp)
#raw.data.lm.ssp = lm(raw.data$phenotype~ raw.data$isolate + raw.data$subset + raw.data$subset:raw.data$rep)
#anova(raw.data.lm.ssp)
raw.data.lm.b = lm(raw.data$phenotype~ raw.data$isolate + raw.data$block)
anova(raw.data.lm.b)
raw.data.lm.s = lm(raw.data$phenotype~ raw.data$isolate + raw.data$subset)
anova(raw.data.lm.s) 
AIC(raw.data.lm.b, raw.data.lm.s)
BIC(raw.data.lm.b, raw.data.lm.s)
#summary(raw.data.lm.s) ## get the estimates for the fixed effects

#################################################################
## BLUE
#################################################################
## Now make it into a mixed model
raw.data.nona<-raw.data[complete.cases(raw.data$phenotype),]
raw.datavarcomp = lmer(raw.data.nona$phenotype ~ raw.data.nona$isolate + 
#                         (1|raw.data.nona$subset))# abun, diam, mel, sp, dip 
                          (1|raw.data.nona$block))# ip, dla, ddla

summary(raw.datavarcomp)
blue = fixef(raw.datavarcomp)  ## get the fixed effects for isolates = BLUEs
par(mfrow=c(1,2))
hist(blue, main = "Raw BLUEs")  
## so the way that R does this is it assigns the first entry as the intercept 
# which is basically zero compared to everything else.  
# So to get "normal" values for the isolates you need to add back on the intercept.
## fix that here 
intercept = blue[1]
blue[1]=0
blue <- blue + intercept
hist(blue) ## vola

##################################################################
## BLUP
##################################################################
## now the model with all random effects
# Linear Model with random effects for variance components
raw.varcomp.rnd = lmer(raw.data.nona$phenotype~ (1|raw.data.nona$isolate) + 
                         (1|raw.data.nona$subset))# abun, diam, mel, sp, dip
                          (1|raw.data.nona$block))# ip, dla
# Extract variance components
summary(raw.varcomp.rnd)
# get BLUPS
blup = ranef(raw.varcomp.rnd)  
par(mfrow=c(1,1))
hist(blup$'raw.data.nona$isolate'[,1])
## these values are way too shrunk.  Can't use this model for BLUPs.  
## Need to stick with using fixed effects for Isolate.  
## Has to do with the way your experiment is set up I think.  
## For a full explanation talk with Peter Bradbury. 

##################################################################
## Comparison Plots
##################################################################
# Create a histogram with the BLUx 
par(mfrow=c(1,2))
hist(raw.data$phenotype, col="gold", main = "Phenotype", breaks=10)
hist(blue, col="blue", main = "BLUEs", breaks=10)

## Compare BLUP to line averages on a scatterplot
raw.data.mean <- tapply(raw.data$phenotype, raw.data$isolate, na.rm=T, mean)
isolates<-levels(raw.data$isolate)
blues<-data.frame(blue)
row.names(blues)<-isolates
par(mfrow=c(1,1))
plot(blues$blue, raw.data.mean, grid())
points((blues[c("StNY001"),]), (blues[c("StNY001"),]), pch="y", col="blue")
points((blues[c("St52B"),]), (blues[c("St52B"),]), pch="b", col="brown")

## Check for BLUP effects on the Block and Subset
raw.nocontrols<-subset(raw.data, isolate!='StNY001')
raw.nocontrols<-subset(raw.nocontrols, isolate!='St52B')
str(raw.nocontrols)
# average over plant, recast, and take plants out
library(reshape2)
raw.nocontrols.c<-dcast(raw.nocontrols, isolate + subset + block + plot
           ~ plate, mean, margins =T, na.rm=T)
raw.nocontrols.c<-dcast(raw.nocontrols, isolate + subset + block + genotype + plot
                        ~ plant, mean, margins =T, na.rm=T) # for ip
#raw.nocontrols.c<-dcast(raw.nocontrols, isolate + subset + block + plot
#                        ~ counts, mean, margins =T, na.rm=T) # for sp
raw.nocontrols.c<-subset(raw.nocontrols.c, plot!= '(all)')
raw.nocontrols.c<-raw.nocontrols.c[,-(5:8)]
raw.nocontrols.c<-raw.nocontrols.c[,-(6:10)] #for ip
head(raw.nocontrols.c)
#blues<-merge(raw.nocontrols.c, blues, by.x='isolate', by.y = "row.names", all=T)
blues<-merge(raw.nocontrols, blues, by.x='isolate', by.y = "row.names", all=T)

blues<-droplevels(blues)
blues$block<-factor(blues$block, c("1","2","3","4", "5","6","7","8","9","10","11","12",
                                   "13","14","15","16","17","18","19","20","21","22","23",
                                   "24","33d","25","26","27","28","29","30","31","32"))
par(mfrow=c(1,2))
boxplot(raw.data$phenotype~raw.data$block, xlab="Block", ylab="trait", main="Trait by Block", col="gold")  
boxplot(blues$blue~blues$block, xlab="Block", ylab="trait", 
        main="BLUE by Block", col="blue")  
boxplot(raw.data$phenotype~raw.data$subset, xlab="subset", ylab="trait", main="Trait by subset", col="gold") 
boxplot(blues$blue~blues$subset, xlab="subset", ylab="BLUEs", 
        main="BLUE by subset", col="blue")

###########################################################
## Save Output
tail(blues)
write.csv(blues, file="outblue.csv")
