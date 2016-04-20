## pairwiselinkageht1.R
## Santiago Mideros
## 13 April 2016

library("genetics")
setwd("~/workdir")

###################
input <- read.csv("ht1hap.csv")
S2_3552054 <- genotype(input$S2_3552054)
S2_3549698 <- genotype(input$S2_3549698)
LD(S2_3552054,S2_3549698)
addmargins(table(S2_3552054,S2_3549698))
chisq.test(c(5,54,48,0))
addmargins(table(input$Phenotype,input$S2_3549698))
addmargins(table(input$Phenotype,input$S2_3552054))

ht2 <- read.csv("avrHt2.csv")
tables<-list()
for (i in 4:15){
  name <- colnames(ht2)[i]
  tmp <- addmargins(table(ht2$phenotype,ht2[,i]))
  tables[[name]] <- tmp
}
tables


