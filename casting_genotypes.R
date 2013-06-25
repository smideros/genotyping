# Genotype casting


getwd()
setwd("~/workdir/Seto")
data=read.table('setogen.txt', header=T)
dim(data)
head(data)

library(reshape)

data.c=cast(data, sample_id ~ genotype_plate_id + call)
head(data.c)
dim(data.c)
write.table(data.c, file='setogen2.txt', col.names=T)