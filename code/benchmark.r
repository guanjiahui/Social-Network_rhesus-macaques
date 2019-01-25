##########
#benchmark Data Mechanics 
nc6 <- read.csv("~/Dropbox/Research/VET_SNH_Monkey/DATA/Benchmark/nc6_ecoli_bandmatching.csv")
nc8 <- read.csv("~/Dropbox/Research/VET_SNH_Monkey/DATA/Benchmark/nc8_ecoli_bandmatching.csv")
nc13<- read.csv("~/Dropbox/Research/VET_SNH_Monkey/DATA/Benchmark/nc13_ecoli_bandmatching.csv")

name6=nc6[,1]
Ben6=as.matrix(nc6)[,-1]

name8=nc8[,1]
Ben8=as.matrix(nc8)[,-1]

name13=nc13[,1]
Ben13=as.matrix(nc13)[,-1]