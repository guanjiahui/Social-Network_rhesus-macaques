NC13_RU2_Huddling_Matrix <- read.csv("~/Dropbox/Research/SNH_health profile data for Fushing-selected/NC13_RU2_Huddling_Matrix.csv")
NC13HuddleR2=as.matrix(NC13_RU2_Huddling_Matrix[,-1])
colnames(NC13HuddleR2)=NC13_RU2_Huddling_Matrix[,1]
rownames(NC13HuddleR2)=NC13_RU2_Huddling_Matrix[,1]
##############

HuddleR2=matrix(0,nrow(NC13HuddleR2),ncol(NC13HuddleR2))
for (i in 1:ncol(NC13HuddleR2)){
  for (j in 1:ncol(NC13HuddleR2)){
    HuddleR2[i,j]=NC13HuddleR2[i,j]+NC13HuddleR2[j,i]
  }
}    

NC13HuddleR2=HuddleR2




##############
ibrary(igraph)
GHR2=graph.adjacency(NC13HuddleR2,mode="directed",weighted=TRUE)

plot(GHR2,edge.arrow.size=.1,edge.color="orange",
     vertex.color="orange", vertex.frame.color="#ffffff")

l <- layout.random(GHR2)
plot(GHR2,layout=l,edge.arrow.size=.1)
###################


#######################
temp=c(0.2,0.3,3,100)
Ens.huddle2=Eigen.plot2(temp, selected.id=c(1,2,3,4),NC13HuddleR2)
DCG.huddle2=DCGtree.plot(num.clusters.selected=c(2,2,17,19),
                         "NC13HuddleR2 tree",Ens.huddle2,temp)
########


heatmap.2(NC13HuddleR2,Rowv=as.dendrogram(DCG.huddle2),Colv=as.dendrogram(DCG.huddle2),
          trace="none",col =colorRampPalette(c("white","green","green4","violet","purple"))(100))

heatmap.2(NC13HuddleR2,col =colorRampPalette(c("white","green","green4","violet","purple"))(100),
          trace="none")
############
plot(DCG.huddle2,hang=-1)
library(sparcl)
# colors the leaves of a dendrogram
y2 = cutree(DCG.huddle2, 20)
y22=cutree(DCG.huddle2,3)
ColorDendrogram(DCG.huddle2, y = y2+1,  main = "NC13HuddleR2 tree", 
                branchlength = 5)
#####################



Eheat2=NC13HuddleR2
small2=which(NC13HuddleR2<5,arr.ind=TRUE)
Eheat2[small2]=NC13HuddleR2[small2]/5
Eheat2[which(NC13HuddleR2>=5,arr.ind=TRUE)]=1
#######


temp=c(0.1,0.2,0.3,0.5,0.8,1,2,2000)
Ens.heat2=Eigen.plot2(temp, selected.id=c(1,2,3,4,5,6,7,8),Eheat2)




DCG.heat2=DCGtree.plot(num.clusters.selected=c(2,2,17,20),
                       "NC13HuddleR1 tree",Ens.heat2,temp)

temp=c(0.2,0.3,2,1000)
Ens.h2=Eigen.plot2(temp, selected.id=c(1,2,3,4),Eheat2)
DCG.h2=DCGtree.plot(num.clusters.selected=c(2,2,17,20),
                       "NC13HuddleR1 tree",Ens.h2,temp)


########


plot(DCG.heat2,hang=-1)
library(sparcl)
# colors the leaves of a dendrogram
y2 = cutree(DCG.heat2, 18)
y22=cutree(DCG.heat2,3)
ColorDendrogram(DCG.heat2, y = y2+1,  main = "NC13HuddleR2 tree", xlab="",
                branchlength = 5)