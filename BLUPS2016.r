

setwd("~/workdir")
library(lme4)

dla814 = read.table("dla814.txt", 
  col.names=c("isolate", "subset", "block", "genotype", "plot", "phenotype"))

## Check to ensure data imported correctly
str(dla814)
head(dla814)
tail(dla814)
# Making sure variables are correctly characterized
dla814$phenotype = as.numeric(dla814$phenotype)
dla814$subset = as.factor(dla814$subset)
dla814$block = as.factor(dla814$block)
# Attach dataset
attach(dla814)

#################################################################
## Exploratory Data: Controls
#################################################################
controls<-subset(dla814, isolate == "StNY001" | isolate == "St52B")
controls$subset <- as.factor(controls$subset)
controls$block
controls<-controls[order(controls$subset, controls$block),]
controls
controls.mean <- aggregate(controls$phenotype~controls$isolate+controls$subset, 
                           list(controls$isolate), mean)
controls.mean


#################################################################
## Exploratory Data Distribution
#################################################################
hist(phenotype, col="gold")
# plot by block
block<-factor(block, c("1","2","3","4", "5","6","7","8","9","10","11","12","13","14","15","16",
                       "17","18","19","20","21","22","23","24","33d","25","26","27","28","29",
                       "30","31","32"))
boxplot(phenotype~block, xlab="Block", ylab="DLA814", main="DLA by Block", col="lightblue")  
# Box plot by subset
boxplot(phenotype~subset, xlab="subset", ylab="DLA814", main="DLA by subset", col="lightblue")  

#################################################################
## Models
#################################################################
dla814.lm.ssb = lm(phenotype~ isolate + subset + subset:block)
anova(dla814.lm.ssb)
dla814.lm.b = lm(phenotype~ isolate + block)
anova(dla814.lm.b)
dla814.lm.s = lm(phenotype~ isolate + subset)
anova(dla814.lm.s) 
AIC(dla814.lm.ssb, dla814.lm.b, dla814.lm.s)
BIC(dla814.lm.ssb, dla814.lm.b, dla814.lm.s)
#summary(dla814.lm.s) ## get the estimates for the fixed effects

#################################################################
## BLUE
#################################################################
## Now make it into a mixed model
dla814varcomp = lmer(phenotype ~ isolate + (1|block))
dla814blue = fixef(dla814varcomp)  ## get the fixed effects for isolates = BLUEs
par(mfrow=c(1,2))
hist(dla814blue, main = "Raw BLUEs")  
## so the way that R does this is it assigns the first entry as the intercept 
# which is basically zero compared to everything else.  
# So to get "normal" values for the isolates you need to add back on the intercept.
## fix that here 
intercept = dla814blue[1]
dla814blue[1]=0
dla814blue <- dla814blue + intercept
hist(dla814blue) ## vola

##################################################################
## BLUP
##################################################################
## now the model with all random effects
# Linear Model with random effects for variance components
dla814varcomp.rnd = lmer(phenotype~ (1|isolate) + (1|block))
# Extract variance components
summary(dla814varcomp.rnd)
# get BLUPS
dla814blup = ranef(dla814varcomp.rnd)  
par(mfrow=c(1,1))
hist(dla814blup$isolate[,1])
## these values are way too shrunk.  Can't use this model for BLUPs.  
## Need to stick with using fixed effects for Isolate.  
## Has to do with the way your experiment is set up I think.  
## For a full explanation talk with Peter Bradbury. 

##################################################################
## Comparison Plots
##################################################################
# Create a histogram with the BLUx 
par(mfrow=c(1,2))
hist(phenotype, col="gold", main = "Phenotype", breaks=9)
hist(dla814blue, col="blue", main = "BLUEs", breaks=9)

## Compare BLUP to line averages on a scatterplot
DLA814.mean <- tapply(phenotype, isolate, na.rm=T, mean)
isolates<-levels(isolate)
blues<-data.frame(dla814blue)
row.names(blues)<-isolates
par(mfrow=c(1,1))
plot(blues$dla814blue, DLA814.mean, grid())
points((blues[c("StNY001"),]), (blues[c("StNY001"),]), pch="y", col="blue")
points((blues[c("St52B"),]), (blues[c("St52B"),]), pch="b", col="brown")

## Check for BLUP effects on the Block and Subset
dla.nocontrols<-subset(dla814, isolate!='StNY001')
dla.nocontrols<-subset(dla.nocontrols, isolate!='St52B')
blues<-merge(dla.nocontrols, blues, by.x='isolate', by.y = "row.names", all=T)
blues$block<-factor(blues$block, c("1","2","3","4", "5","6","7","8","9","10","11","12",
                                   "13","14","15","16","17","18","19","20","21","22","23",
                                   "24","33d","25","26","27","28","29","30","31","32"))
par(mfrow=c(1,2))
boxplot(phenotype~block, xlab="Block", ylab="DLA814", main="DLA by Block", col="gold")  
boxplot(blues$dla814blue~blues$block, xlab="Block", ylab="DLA814", 
        main="BLUE by Block", col="blue")  
boxplot(phenotype~subset, xlab="subset", ylab="DLA814", main="DLA by subset", col="gold") 
boxplot(blues$dla814blue~blues$subset, xlab="subset", ylab="BLUEs", 
        main="BLUE by subset", col="blue")

###########################################################
## Save Output

write.csv(blues, file="DLA814blue.csv")
