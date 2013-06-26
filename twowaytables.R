setwd("~/workdir/Seto")
data=read.table('~/Dropbox/Seto/data & analysis/avrHt1/phen_gen.txt', header=T)
Disease.n=as.numeric(data[,2]) # makes disease into 1 or 2 for R and S
Disease=(data[,2])
S2_3549689=data[,3]
S2_3552054=data[,4]
ObsTab1=table(Disease,S2_3549689)
ObsTab1
ObsTab2=table(Disease,S2_3552054)
ObsTab2