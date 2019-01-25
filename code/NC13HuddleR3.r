NC13_RU3_Huddling_Matrix <- read.csv("~/Dropbox/Research/SNH_health profile data for Fushing-selected/NC13_RU3_Huddling_Matrix.csv")
NC13HuddleR3=as.matrix(NC13_RU3_Huddling_Matrix[,-1])
colnames(NC13HuddleR3)=NC13_RU3_Huddling_Matrix[,1]
rownames(NC13HuddleR3)=NC13_RU3_Huddling_Matrix[,1]
##############

HuddleR3=matrix(0,nrow(NC13HuddleR3),ncol(NC13HuddleR3))
for (i in 1:ncol(NC13HuddleR3)){
  for (j in 1:ncol(NC13HuddleR3)){
    HuddleR3[i,j]=NC13HuddleR3[i,j]+NC13HuddleR3[j,i]
  }
}    

NC13HuddleR3=HuddleR3

temp=c(0.2,0.3,3,100)
Ens.huddle3=Eigen.plot2(temp, selected.id=c(1,2,3,4),NC13HuddleR3)
DCG.huddle3=DCGtree.plot(num.clusters.selected=c(2,2,16,18),
                         "NC13HuddleR3 tree",Ens.huddle3,temp)
########
heatmap.2(NC13HuddleR3,Rowv=as.dendrogram(DCG.huddle3),Colv=as.dendrogram(DCG.huddle3),
          trace="none",col =colorRampPalette(c("white","green","green4","violet","purple"))(100))

heatmap.2(NC13HuddleR3,col =colorRampPalette(c("white","green","green4","violet","purple"))(100),
          trace="none")
#########
#double check the temperature selection 
temp=c(0.1,0.2,0.3,0.5,0.6,0.8,1,100)
Ens.h3=Eigen.plot2(temp, selected.id=c(1,2,3,4,5,6,7,8),Eheat3)

DCG.h3=DCGtree.plot(num.clusters.selected=c(2,5,10,10,19,20,22,24),
                       "NC13HuddleR1 tree",Ens.h3,temp)
########






#############
library(sparcl)
# colors the leaves of a dendrogram
y3 = cutree(DCG.huddle3, 28)
y32=cutree(DCG.huddle3,3)
ColorDendrogram(DCG.huddle3, y = y3+1,  main = "NC13HuddleR3 tree", xlab="",
                branchlength = 5)
#####################


#####################



Eheat3=NC13HuddleR3
small3=which(NC13HuddleR3<5,arr.ind=TRUE)
Eheat3[small3]=NC13HuddleR3[small3]/5
Eheat3[which(NC13HuddleR3>=5,arr.ind=TRUE)]=1
#######
temp=c(0.2,0.3,2,1000)
Ens.heat3=Eigen.plot2(temp, selected.id=c(1,2,3,4),Eheat3)
DCG.heat3=DCGtree.plot(num.clusters.selected=c(2,2,20,23),
                       "NC13HuddleR1 tree",Ens.heat3,temp)
########

y3 = cutree(DCG.heat3, 24)
y32=cutree(DCG.heat3,3)
ColorDendrogram(DCG.heat3, y = y3+1,  main = "NC13HuddleR3 tree", xlab="",
                branchlength = 5)
#####################
#visualize it against baseline
ColorDendrogram(DCG.huddle3, y = y1[-c(37,90)]+1,  main = "NC13HuddleR3 tree", xlab="",
                branchlength = 5)