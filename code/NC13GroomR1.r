NC13_RU1_Grooming_Matrix <- read.csv("~/Dropbox/Research/SNH_health profile data for Fushing-selected/NC13_RU1_Grooming_Matrix.csv")
NC13GroomR1=as.matrix(NC13_RU1_Grooming_Matrix [,-1])
colnames(NC13GroomR1)=NC13_RU1_Grooming_Matrix[,1]
rownames(NC13GroomR1)=NC13_RU1_Grooming_Matrix[,1]
##############
ibrary(igraph)
GGrR1=graph.adjacency(NC13GroomR1,mode="directed",weighted=TRUE)


plot(GGrR1,edge.arrow.size=.1,edge.color="orange",
     vertex.color="orange", vertex.frame.color="#ffffff")

l <- layout.random(GGrR1)
plot(GGrR1,layout=l,edge.arrow.size=.1)
###################


temp=c(1e-06,0.00001,1)
Ens.groom1=Eigen.plot(temp, selected.id=c( 1,2,3),NC13GroomR1)
DCG.groom1=DCGtree.plot(num.clusters.selected=c(1,2,2),
                         "NC13GroomR1 tree",Ens.groom1,temp)
####################
temp=c(0.2,0.5,10,100)
Ens.groom1=Eigen.plot2(temp, selected.id=c(1,2,3,4),NC13GroomR1)
DCG.groom1=DCGtree.plot(num.clusters.selected=c(2,2,5,7),
                        "NC13GroomR1 tree",Ens.groom1,temp)
plot(DCG.groom1,hang=-1)
#####
heatmap.2(NC13GroomR1,Rowv=as.dendrogram(DCG.groom1),Colv=as.dendrogram(DCG.groom1),
          trace="none",col =colorRampPalette(c("white","green","green4","violet","purple"))(100))
heatmap.2(NC13GroomR1,col =colorRampPalette(c("white","green","green4","violet","purple"))(100),
          trace="none")
##########
##########
plot(DCG.groom1,hang=-1)
#################




###################

s1=sum(NC13GroomR1)
s2=sum(NC13GroomR2)
s3=sum(NC13GroomR3)
NC13GroomR1.del=NC13GroomR1[-c(37,90),]
NC13GroomR1.del=NC13GroomR1.del[,-c(37,90)]
s1.del=sum(NC13GroomR1.del)

eachR1=colSums(NC13GroomR1)
eachR2=colSums(NC13GroomR2)
eachR3=colSums(NC13GroomR3)
eachR1.del=colSums(NC13GroomR1.del)
eachSumG=cbind(eachR1.del,eachR2,eachR3)

SeachR1=rowSums(NC13GroomR1)
SeachR2=rowSums(NC13GroomR2)
SeachR3=rowSums(NC13GroomR3)
SeachR1.del=rowSums(NC13GroomR1.del)
SeachSumG=cbind(SeachR1.del,SeachR2,SeachR3)


nameR1=colnames(NC13GroomR1)
nameR2=colnames(NC13GroomR2)
nameR3=colnames(NC13GroomR3)
nameR1.del=nameR1[-c(37,90)]

###############
plot(c(s1,s2,s3),type="b")

plot(c(s1.del,s2,s3),type="b",ylab="Groom",xaxt="n",main="NC13 Grooming Total")
axis(1, at=1:3, labels=bb2)

matplot(t(SeachSumG), t="l", lty=1, las=1, ylab="Groom",col=1:99,
        xlab="Time", xaxt="n",main="NC13 sending Grooming")
axis(1, at=1:3, labels=bb2)

matplot(t(eachSumG), t="l", lty=1, las=1, ylab="Groom",col=1:99,
        xlab="Time", xaxt="n",main="NC13 receiving Grooming")
axis(1, at=1:3, labels=bb2)
############

##################
plot(density(SeachR1.del),ylim=c(0,0.07),main="Grooming sending density",lwd=3)
lines(density(SeachR2),col=2,lwd=3)
lines(density(SeachR3),col=4,lwd=3)
abline(v=SeachR1[37],col="yellow",lwd=3,lty=2)
abline(v=SeachR1[90],col="green",lwd=3,lty=2)

legend("topright",c("baseline","perturbation","postperturbation"),
       lty=c(1,1,1),col=c(1,2,4),cex=1.3)
###################################



##############
par(mfrow=c(1,3))
#Groom=data.frame(cbind(eachSumG,SeachSumG))
#pairs(~eachR1.del+SeachR1.del+eachR2 + SeachR2+eachR3+ SeachR3,data=Groom, 
    #  main="Simple Scatterplot Matrix")
plot(eachR1,SeachR1,xlab="received",ylim=c(0,max(SeachR1)),ylab="sent",main="R1",xlim=c(0,max(eachR1)))
lines(1:100,1:100)
plot(eachR2,SeachR2,xlab="received",ylim=c(0,max(SeachR1)),ylab="sent",main="R2",xlim=c(0,max(eachR1)))
lines(1:100,1:100)
plot(eachR3,SeachR3,xlab="received",ylim=c(0,max(SeachR1)),ylab="sent",main="R3",xlim=c(0,max(eachR1)))
lines(1:100,1:100)

#########################
library(Perc)

win_prob=conductance(NC13GroomR1,maxLength = 4)
win_prob1=win_prob$p.hat

temp=c(0.006,0.07,0.15,0.3,0.8,1)
Ens.groom1=Eigen.plot2(temp, selected.id=c(1,2,3,4,5,6),win_prob1)

DCG13.groom1=DCGtree.plot(num.clusters.selected=c(1,1,2,3,4,6),
                        "NC13GroomR1 tree",Ens.groom1,temp)
plot(DCG13.groom1,hang=-1,main="NC13GroomR1 Tree")

#####
heatmap.2(NC13GroomR1,Rowv=as.dendrogram(DCG13.groom1),Colv=as.dendrogram(DCG13.groom1),
          trace="none",col =colorRampPalette(c("white","green","green4","violet","purple"))(100))
heatmap.2(NC13GroomR1,col =colorRampPalette(c("white","green","green4","violet","purple"))(100),
          trace="none")
##########
library(sparcl)
plot(DCG.groom1,hang = -1)
x1 = cutree(DCG.groom1, 14)
x12=cutree(DCG.groom1,2)
ColorDendrogram(DCG.groom1, y = x1+3,  main = "NC13 Groom R1 tree",xlab="", 
                branchlength = 10)

######################










################
#want to check the change of the number of each monkey's friends
#first normalize the data to make it to be binary 

temp1=which(NC13GroomR1>0,arr.ind=TRUE)
BinaryGroomR1=NC13GroomR1
BinaryGroomR1[temp1]=1

temp2=which(NC13GroomR2>0,arr.ind=TRUE)
BinaryGroomR2=NC13GroomR2
BinaryGroomR2[temp2]=1

temp3=which(NC13GroomR3>0,arr.ind=TRUE)
BinaryGroomR3=NC13GroomR3
BinaryGroomR3[temp3]=1

##########
#check the row and colum sum for each stage 
rowBinary=c(rowSums(BinaryGroomR1)[-c(37,90)], rowSums(BinaryGroomR2),rowSums(BinaryGroomR3))
colBinary=c(colSums(BinaryGroomR1)[-c(37,90)], colSums(BinaryGroomR2),colSums(BinaryGroomR3))

TreeRow=hclust(dist(rowBinary))
TreeCol=hclust(dist(colBinary))

plot(rowBinary, type="l",lwd=3,xlab="monkey ID",ylab="number of monkey each indivudal was groomed",main="Groom Received")
lines(100:198,rowBinary[100:198],col=2,lwd=3)
lines(199:297,rowBinary[199:297],col=4,lwd=3)
legend("topright",c("baseline","perturbation","post"),
       lty=c(1,1,1),col=c(1,2,4),cex=1.3)
###################################


plot(colBinary, type="l",main="Groom Sent",lwd=3,xlab="monkey ID",ylab="number of monkey each indivudal groomed")
lines(100:198,colBinary[100:198],col=2,lwd=3)
lines(199:297,colBinary[199:297],col=4,lwd=3)
legend("topright",c("baseline","perturbation","post"),
       lty=c(1,1,1),col=c(1,2,4),cex=1.3)
###################################


RowGroom=as.dendrogram(TreeRow)
ColGroom=as.dendrogram(TreeCol)

plot(TreeRow, main="Groom(friends) Received: 6 clusters",hang=-1)
plot(TreeCol, main="Groom(friends) sent: 6 clusters",hang=-1)

rowMember=cutree(TreeRow,k=6)
#get the average value of each cluster 
avgClusterRow=c()
for (i in 1:6)
{
  temp=which(rowMember==i)
  avgClusterRow[i]=mean(rowBinary[temp])
}

#sort the cluster 
rowO=order(avgClusterRow)
rowMember1=numeric(99)
for(i in 1:6)
{
  temp=which(rowMember==rowO[i])
  rowMember1[temp]=i
}

#################
MRchange=matrix(0,99,4)
for(i in 1:99){
  MRchange[i,1]=as.numeric(names(rowMember)[i])
  MRchange[i,2]=rowMember1[i]
  MRchange[i,3]=rowMember1[i+99]
  MRchange[i,4]=rowMember1[i+198]
}




matplot(t(MRchange[,2:4]), t="l", lty=1, las=1, ylab="cluster No",col=1:99,
        xlab="Time", xaxt="n",main="Groom Received")
axis(1, at=1:3, labels=bb2)

#################
colMember=cutree(TreeCol,k=6)
#get the average value of each cluster 
avgClusterCol=c()
for (i in 1:6)
{
  temp=which(colMember1==i)
  avgClusterCol[i]=mean(colBinary[temp])
}

colO=order(avgClusterCol)
colMember1=numeric(99)
for(i in 1:6)
{
  temp=which(colMember==colO[i])
  colMember1[temp]=i
}


CRchange=matrix(0,99,4)
for(i in 1:99){
  CRchange[i,1]=as.numeric(names(rowMember)[i])
  CRchange[i,2]=colMember1[i]
  CRchange[i,3]=colMember1[i+99]
  CRchange[i,4]=colMember1[i+198]
}

matplot(t(CRchange[,2:4]), t="l", lty=1, las=1, ylab="cluster No",col=1:99,
        xlab="Time", xaxt="n",main="Groom Sent")
axis(1, at=1:3, labels=bb2)

########################
#now visualize the results
library(plotrix)
#small example 
slices <- c(10, 12, 4, 16, 8) 
lbls <- c("US", "UK", "Australia", "Germany", "France")
pie3D(slices,labels=lbls,explode=0.1,
      main="Pie Chart of Countries ")


#now apply to our data 
tmp=c()
for(i in 1:6)
  tmp[i]=length(which(CRchange[,2]==i))
slices <-tmp
lbls <- c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5","Cluster 6")
pct <- round(slices/sum(slices)*100)
lbls <- paste(lbls, pct) # add percents to labels 
lbls <- paste(lbls,"%",sep="") # ad % to labels 
pie(slices,labels = lbls, col=rainbow(length(lbls)),
    main="Baseline Composition")

######################################
TrendPie=function(original,current)
{
  #first check how many different clusters the current time stage has 
  ClusterNo=unique(current)
  m=length(ClusterNo)
  slice=list()
  membership=list()
  size=c()
  
  for(i in 1:m)
  {
    location=which(current==i)
    member=original[location] #find the invidual that falls in the cluster in current stage
    
    size[i]=length(member)
    #check how many different clusters the original one comes from 
    membership[[i]]=unique(member)
    N=length(membership[[i]])
    temp=numeric(N)
    
    for(j in 1:N)
    {
      temp[j]=length(which(member==membership[[i]][j]))
    }
    
    slice[[i]]=temp
  }#end of for loop of i 
  
  return(list(clusterNo=ClusterNo,size=size, slices=slice,sliceCluster=membership))
}#TrendPie() function 
#####################################
SendPie1=TrendPie(CRchange[,2],CRchange[,3])
N=length(SendPie1$clusterNo)
#par(mfrow=c(N,1))

m <- matrix(c(1,2,3,4,5,6,7,7,7),nrow = 3,ncol = 3,byrow = TRUE)

layout(mat = m,heights = c(0.45,0.45,0.1))

for(i in c(5,4,3,2,1))
{
  par(mar = c(2,2,1,1))
  slices=SendPie1$slices[[i]]
  lbls=SendPie1$sliceCluster[[i]]
  #pie3D(slices,labels=lbls,explode=0.2,radius = SendPie1$size[i]/10,
     #   main=paste("Pie Chart of Perturbation of Cluster",i))
  pct <- round(slices/sum(slices)*100)
  lbls <-pct # add percents to labels 
  lbls <- paste(lbls,"%",sep="") # ad % to labels 
  pie(slices,labels = lbls, col=SendPie1$sliceCluster[[i]]+1,
        radius = SendPie1$size[i]/48,font=2,cex=1.5,
    main=paste("Perturbation Cluster",i))
}
par(mar = c(2,2,1,1))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")


legend(x = "top",inset = 0, c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5","Cluster 6"), 
       bty="n",col=2:7,lwd=5, cex=1.2, pt.cex=1.8,horiz = TRUE)






SendPie2=TrendPie(CRchange[,3],CRchange[,4])
N=length(SendPie2$clusterNo)

m <- matrix(c(1,2,3,4,5,6,7,7,7),nrow = 3,ncol = 3,byrow = TRUE)

layout(mat = m,heights = c(0.45,0.45,0.1))

for(i in c(6,5,4,3,2,1))
{
  par(mar = c(2,2,1,1))
  slices=SendPie2$slices[[i]]
  lbls=SendPie2$sliceCluster[[i]]
  #pie3D(slices,labels=lbls,explode=0.2,radius = SendPie1$size[i]/10,
  #   main=paste("Pie Chart of Perturbation of Cluster",i))
  pct <- round(slices/sum(slices)*100)
  lbls <-pct # add percents to labels 
  lbls <- paste(lbls,"%",sep="") # ad % to labels 
  pie(slices,labels = lbls, col=SendPie2$sliceCluster[[i]]+1,
      radius = SendPie1$size[i]/48,font=2,cex=1.5,
      main=paste("Post Cluster",i))
}
par(mar = c(2,2,1,1))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")


legend(x = "top",inset = 0, c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5","Cluster 6"), 
       bty="n",col=2:7,lwd=5, cex=1.2, pt.cex=1.8,horiz = TRUE)






SendPie3=TrendPie(MRchange[,2],MRchange[,3])
N=length(SendPie3$clusterNo)
#par(mfrow=c(N,1))

m <- matrix(c(1,2,3,4,5,6,7,7,7),nrow = 3,ncol = 3,byrow = TRUE)

layout(mat = m,heights = c(0.45,0.45,0.1))

for(i in c(5,4,3,2,1))
{
  par(mar = c(2,2,1,1))
  slices=SendPie3$slices[[i]]
  lbls=SendPie3$sliceCluster[[i]]
  #pie3D(slices,labels=lbls,explode=0.2,radius = SendPie1$size[i]/10,
  #   main=paste("Pie Chart of Perturbation of Cluster",i))
  pct <- round(slices/sum(slices)*100)
  lbls <-pct # add percents to labels 
  lbls <- paste(lbls,"%",sep="") # ad % to labels 
  pie(slices,labels = lbls, col=SendPie1$sliceCluster[[i]]+1,
      radius = SendPie3$size[i]/48,font=2,cex=1.5,
      main=paste("Perturbation Cluster",i))
}

plot(1, type = "n", axes=FALSE, xlab="", ylab="")

legend(x = "top",inset = 0, c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5","Cluster 6"), 
       bty="n",col=2:7,lwd=5, cex=1.2, pt.cex=1.8,horiz = TRUE)




#############

SendPie2=TrendPie(CRchange[,3],CRchange[,4])
N=length(SendPie2$clusterNo)

m <- matrix(c(1,2,3,4,5,6,7,7,7),nrow = 3,ncol = 3,byrow = TRUE)

layout(mat = m,heights = c(0.45,0.45,0.1))

for(i in c(5,4,3,2,1))
{
  par(mar = c(2,2,1,1))
  slices=SendPie2$slices[[i]]
  lbls=SendPie2$sliceCluster[[i]]
  #pie3D(slices,labels=lbls,explode=0.2,radius = SendPie1$size[i]/10,
  #   main=paste("Pie Chart of Perturbation of Cluster",i))
  pct <- round(slices/sum(slices)*100)
  lbls <-pct # add percents to labels 
  lbls <- paste(lbls,"%",sep="") # ad % to labels 
  pie(slices,labels = lbls, col=SendPie2$sliceCluster[[i]]+1,
      radius = SendPie1$size[i]/48,font=2,cex=1.5,
      main=paste("Post Cluster",i))
}
par(mar = c(2,2,1,1))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")

plot(1, type = "n", axes=FALSE, xlab="", ylab="")

legend(x = "top",inset = 0, c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5","Cluster 6"), 
       bty="n",col=2:7,lwd=5, cex=1.2, pt.cex=1.8,horiz = TRUE)




#############################################
#make it in a transition probability matrix 
A13=Transition(MRchange[,2],MRchange[,3])
B13=Transition(MRchange[,3],MRchange[,4])


heatmap.2(A13,Rowv=FALSE,Colv = FALSE, main="Pert vs Base",trace="none")
heatmap.2(B13,Rowv=FALSE,Colv = FALSE, main="Post vs Pert",trace="none")








#####################
#build the DCG tree 
###################



win1=conductance(NC13GroomR1,maxLength = 4)
win_prob1=win1$p.hat

temp=c(0.008,0.07,0.15,0.3,1)
Ens13.groom1=Eigen.plot(temp, selected.id=c(1,2,3,4,5),win_prob1)

DCG13.groom1=DCGtree.plot(num.clusters.selected=c(1,1,2,3,6),
                        "NC13GroomR1 tree",Ens13.groom1,temp)
plot(DCG13.groom1,hang=-1,main="NC13GroomR1 tree")


library(gplots)

##############
#compute the entropy 
###########
G1=cutree(DCG.groom1,k=4)
EntryG=list()
EntryG[[1]]=G1[-c(37,90)]
EntryG[[2]]=G2
EntryG[[3]]=G3

par(mfrow=c(2,3))
EnG=matrix(0,3,3)
for (i in 1:3){
  for (j in 1:3){
    if (i!=j){
      EnG[i,j]=Entropy_sequence(EntryG[[i]],EntryG[[j]],Name[[i]],Name[[j]])$Sum
      plot(Entropy_sequence(EntryG[[i]],EntryG[[j]],Name[[i]],Name[[j]])$Entropy, type="b",
           xlab="group", ylab="Entropy",main=paste("stage",i,"from","stage",j))
    }
  }
}

for (i in 1:3)
  for (j in 1:3)
    cat(EnG[i,j]+EnG[j,i],"i",i,"j",j,"\n")

####################
#save the gromming tree 
save(DCG.groom1,DCG.groom2,DCG.groom3,file = "GroomingTREE.RData")

##########

#compute the second order of the spaghetti plots 

###################
SecondOrder=function(x){
  Second =x[1]+x[3]-2*x[2]
  return(Second)
}
##############

CurveCheck=function(dataFrame)
{
  count.cancave=0
  count.convex=0
  count.Nocurve=0
  N=length(dataFrame[,1])
  
  for(i in 1:N)
  {
    temp=SecondOrder(dataFrame[i,])
    if(temp>0)
      count.cancave=count.cancave+1
    else if(temp<0)
      count.convex=count.convex+1
    else
      count.Nocurve=count.Nocurve+1
  }

  return(list(convex=count.convex,concave=count.cancave,noChange=count.Nocurve))
}


###############

SendCurve=CurveCheck(SeachSumG)


count.increase=0
count.decrease=0
IndIncrease=0
IndDecrease=0
#############
for(i in 1:99){
  
  if(SeachSumG[i,3]>=SeachSumG[i,1]){
    count.increase=count.increase+1
    IndIncrease=c(IndIncrease,i)
  }
  
  else if(SeachSumG[i,3]<SeachSumG[i,1]){
    count.decrease=count.decrease+1
    IndDecrease=c(IndDecrease,i)
  }

}

SendIncrease=CurveCheck(SeachSumG[IndIncrease[-1],])
SendDecrease=CurveCheck(SeachSumG[IndDecrease[-1],])


#####################
#Compare to the Rank Dominance Probability NC13Rank 
SNH02_13B_AllDomProbMatrices <- read.csv("~/Dropbox/Research/SNH_health profile data for Fushing-selected/DATA/NC13/SNH02_13B_AllDomProbMatrices.csv")
NC13Rank=SNH02_13B_AllDomProbMatrices
NC13Rank=NC13Rank[-1,]
rownames(NC13Rank)=NC13Rank[,2]
NC13Rank=NC13Rank[,-c(1,2)]
NC13Rank=as.matrix(NC13Rank)
colnames(NC13Rank)=rownames(NC13Rank)
#########
#################################
DistSpecial=function(Mat){
  m=nrow(Mat)
  n=ncol(Mat)
  V=matrix(0,m,n)
  
  for(i in 1:m)
  {
    for(j in 1:n)
    {
      temp=0
      for(k in 1:n)
      {
        if(k!=i || k!=j)
          temp=temp+(Mat[i,k]-Mat[j,k])^2
      }
      V[i,j]=sqrt(temp+(Mat[i,j]-Mat[j,i])^2)
      
    }
  }
  
  return(V)
}
###################################
DistDominance=DistSpecial(NC13Rank)
temp=c(0.07,0.1,0.15,0.4,1)

Ens.Do1=Eigen.plot2(temp, selected.id=c(1,2,3,4,5),DistDominance)
DCG.Do1=DCGtree.plot(num.clusters.selected=c(1,2,3,4,5),
                        "NC13 Dominance Prob tree",Ens.Do1,temp)

plot(DCG.Do1,hang=-1,main="NC13 Dominance Prob Tree")


#get the rank 
library(Perc)
Rank13order=simRankOrder(NC13Rank)

Order13=Rank13order$BestSimulatedRankOrder





###################################

#reconstruct DCG for Grooming baseline 
Grooming13=MatrixSym(NC13GroomR1)
GroomDist13=DistPairwise(Grooming13)


temp=c(1,4,7.5,15,40)
Ens13.groom1=Eigen.plot(temp, selected.id=c(1,2,3,4,5),GroomDist13)

DCG13.groom1=DCGtree.plot(num.clusters.selected=c(1,2,3,5,8),
                          "NC13GroomR1 tree",Ens13.groom1,temp)
plot(DCG13.groom1,hang=-1,main="Grooming NC13")
######################3
#cut the DCG tree to obtain the membership information 

yNC13G1=cutree(DCG13.groom1,k=8)
