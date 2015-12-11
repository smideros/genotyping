#Load libraries
#{r, echo=FALSE}
library(ggplot2)
library(lme4)
library(plyr)


#Import and check data import
setwd("~/workdir")
mydata <- read.csv("uoiclinic15.csv", as.is=T)
head(mydata)
tail(mydata)
str(mydata)

##Subset data
mydata.corn<-mydata[mydata$HOST.HABITAT.TAXONOMIC.NAME == c("Zea mays"),]
head(mydata.corn)

##Bar Charts
counts<-sort(table(mydata.corn$DIAGNOSIS.ID.COMMON.NAME),decreasing = T)
par(mar=c(7,4,2,2)+4)
barplot(counts, ylab="Number of samples",cex.names=0.75,las=3)

##Data check
pairs(mydata)
hist(tillers$ilength)


#Run model
out2 <- lm(ilength ~ planttype-1 + event, tillers)
summary(out2)
anova(out2)

hist(tillers$lrate)
lmrate <- lm(lrate ~ planttype-1 + event, tillers)
summary(lmrate)
#se <- sqrt(diag(vcov(lmrate)))
#capture.output(summary(out),file="tillersanova.csv")

#ANOVA and pairwise comparisons
anova(lmrate)
pairwise.t.test(tillers$lrate, tillers$planttype, p.adj = "none")




##Plots
stderr <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))
se <- tapply(X = tillers$lrate, INDEX = list(tillers$planttype), FUN = stderr)
mm <- ddply(tillers, "planttype", summarise, lrate = mean(lrate, na.rm = TRUE))
#limits <- aes(ymax = la + se, ymin=la - se)
p <- ggplot(data=mm, aes(x=planttype, y=lrate, fill=planttype)) + geom_bar(stat="identity", colour="black") + xlab("Plant Type") + ylab("Lesion Rate")
pp <- p + theme_bw() + geom_errorbar(aes(ymin=lrate-se, ymax=lrate+se), width=.1) + theme(legend.position="none") + scale_fill_brewer(palette="Dark2")
pp
pdf("tillers.pdf")
print(pp)
dev.off()



