NC13_RU3_Grooming_Matrix <- read.csv("~/Dropbox/Research/SNH_health profile data for Fushing-selected/NC13_RU3_Grooming_Matrix.csv")
NC13GroomR3=as.matrix(NC13_RU3_Grooming_Matrix [,-1])
colnames(NC13GroomR3)=NC13_RU3_Grooming_Matrix[,1]
rownames(NC13GroomR3)=NC13_RU3_Grooming_Matrix[,1]
##############
library(Perc)
win3=conductance(NC13GroomR3,maxLength = 4)
win_prob3=win3$p.hat



temp=c(0.0045,0.07,0.15,0.3,0.8,1)
Ens13.groom3=Eigen.plot2(temp, selected.id=c(1,2,3,4,5,6),win_prob3)
DCG13.groom3=DCGtree.plot(num.clusters.selected=c(1,1,2,3,5,9),
                        "NC13GroomR3 tree",Ens13.groom3,temp)
plot(DCG13.groom3,hang=-1,main="NC13GroomR3 tree")


G3=cutree(DCG.groom3,k=6)

