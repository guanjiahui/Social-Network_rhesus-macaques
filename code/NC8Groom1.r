######################
matrixNoSym=function(Initiator, Recipient,interaction)
{
  ID=unique(Initiator)
  n=length(ID)
  N=length(Initiator)
  Mat=matrix(0,n,n)
 
  colnames(Mat)=ID
  rownames(Mat)=ID
  
  for(i in 1:N)
  {
    I=which(ID==Initiator[i])
    J=which(ID==Recipient[i])
    if(interaction[i]==1)
       Mat[I,J]=Mat[I,J]+1
  }
  
  return(Mat)
}
######################
NC8_BL_Grooming <- read.csv("C:/Users/jiahui/Dropbox/Research/SNH_health profile data for Fushing-selected/DATA/NC8/NC8_BL_Grooming.csv")
NC8Groom1=matrixNoSym(NC8_BL_Grooming$Initiator,NC8_BL_Grooming$Recipient,NC8_BL_Grooming$Interaction)

NC8_PERT_Grooming <- read.csv("C:/Users/jiahui/Dropbox/Research/SNH_health profile data for Fushing-selected/DATA/NC8/NC8_PERT_Grooming.csv")

NC8Groom2=matrixNoSym(NC8_PERT_Grooming$Initiator,NC8_PERT_Grooming$Recipient,NC8_PERT_Grooming$Interaction)

NC8_POST_Grooming <- read.csv("C:/Users/jiahui/Dropbox/Research/SNH_health profile data for Fushing-selected/DATA/NC8/NC8_POST_Grooming.csv")

NC8Groom3=matrixNoSym(NC8_POST_Grooming$Initiator,NC8_POST_Grooming$Recipient,NC8_POST_Grooming$Interaction)
############
#check recipient and intiator are the same (i,e, diagonal is nonzero)
for(i in 1:ncol(NC8Groom1))
{
  for(j in 1:ncol(NC8Groom1))
  {
    if(i==j)
      cat(NC8Groom1[i,j],i,"\n")
  }
}

NC8Groom1[49,49]=0
NC8Groom3[25,25]=0
#> dim(NC8Groom1)
#[1] 96 96
#> dim(NC8Groom2)
#[1] 85 85
#> dim(NC8Groom3)
#[1] 83 83
#########################

name8G=list()
name8G[[1]]=colnames(NC8Groom1)
name8G[[2]]=colnames(NC8Groom2)
name8G[[3]]=colnames(NC8Groom3)
###
bb=intersect(name8G[[1]],name8G[[2]])
common2=intersect(bb,name8G[[3]])

###########
Find2=function(ID,subset) #where subset is a subset if ID 
{
  loc=c()
  for(i in 1:length(subset))
  {
    temp=which(ID==subset[i])
    loc=c(temp,loc)
  }
  return(loc)
}

loc=list()
for(i in 1:3)
{
  loc[[i]]=Find2(name8G[[i]],common2)
}

################################
NC82=list()
NC82[[1]]=NC8Groom1[loc[[1]],loc[[1]]]
NC82[[2]]=NC8Groom2[loc[[2]],loc[[2]]]
NC82[[3]]=NC8Groom3[loc[[3]],loc[[3]]]

eachR1=colSums(NC82[[1]])
eachR2=colSums(NC82[[2]])
eachR3=colSums(NC82[[3]])
eachSum=cbind(eachR1,eachR2,eachR3)

SeachR1=rowSums(NC82[[1]])
SeachR2=rowSums(NC82[[2]])
SeachR3=rowSums(NC82[[3]])
SeachSumG=cbind(SeachR1,SeachR2,SeachR3)


S82=c()
for(i in 1:3)
  S82[i]=sum(NC82[[i]])


###############

plot(S82,type="b",ylab="Groom",xaxt="n",ylim=c(min(s3),max(s1)),col=2,lwd=3,pch=16,
     main="NC8 Grooming Total",xlab.font=2)


lines(c(s1.del,s2,s3),lty=2,col=3,lwd=3)
points(c(s1.del,s2,s3),lty=2,col=3,lwd=3,pch=15)
axis(1, at=1:3, labels=bb2)

legend("topright",c("NC8","NC13"),
       lty=c(1,2),col=c(2,3),cex=1.3)


















#######################

matplot(t(SeachSumG), t="l", lty=1, las=1, ylab="Groom",col=1:99,
        xlab="Time", xaxt="n",main="NC8sending Grooming")
axis(1, at=1:3, labels=bb2)

matplot(t(eachSumG), t="l", lty=1, las=1, ylab="Groom",col=1:99,
        xlab="Time", xaxt="n",main="NC8 receiving Grooming")
axis(1, at=1:3, labels=bb2)
############

##################
plot(density(SeachR1),ylim=c(0,0.05),main="Grooming sending density",lwd=3)
lines(density(SeachR2),col=2,lwd=3)
lines(density(SeachR3),col=4,lwd=3)
legend("topright",c("baseline","perturbation","postperturbation"),
       lty=c(1,1,1),col=c(1,2,4),cex=1.3)
###################################

##################
plot(density(eachR1),ylim=c(0,0.05),main="Grooming receiving density",lwd=3)
lines(density(eachR2),col=2,lwd=3)
lines(density(eachR3),col=4,lwd=3)
legend("topright",c("baseline","perturbation","postperturbation"),
       lty=c(1,1,1),col=c(1,2,4),cex=1.3)
###################################



















##############################
#check the friends' circle 
#first normalize the data to make it to be binary 

temp1=which(NC82[[1]]>0,arr.ind=TRUE)
BinaryGroomR1=NC82[[1]]
BinaryGroomR1[temp1]=1

temp2=which(NC82[[2]]>0,arr.ind=TRUE)
BinaryGroomR2=NC82[[2]]
BinaryGroomR2[temp2]=1

temp3=which(NC82[[3]]>0,arr.ind=TRUE)
BinaryGroomR3=NC82[[3]]
BinaryGroomR3[temp3]=1

##########
#check the row and colum sum for each stage 
rowBinary=c(rowSums(BinaryGroomR1), rowSums(BinaryGroomR2),rowSums(BinaryGroomR3))
colBinary=c(colSums(BinaryGroomR1), colSums(BinaryGroomR2),colSums(BinaryGroomR3))

TreeRow=hclust(dist(rowBinary))
TreeCol=hclust(dist(colBinary))

plot(rowBinary, type="l",lwd=3,xlab="monkey ID",
     ylab="number of monkey each indivudal was groomed",main="Groom Received")
lines(84:166,rowBinary[84:166],col=2,lwd=3)
lines(167:249,rowBinary[167:249],col=4,lwd=3)
legend("topright",c("baseline","perturbation","post"),
       lty=c(1,1,1),col=c(1,2,4),cex=1.3)
###################################


plot(colBinary, type="l",main="Groom Sent",lwd=3,xlab="monkey ID",ylab="number of monkey each indivudal groomed")
lines(84:166,colBinary[84:166],col=2,lwd=3)
lines(167:249,colBinary[167:249],col=4,lwd=3)
legend("topright",c("baseline","perturbation","post"),
       lty=c(1,1,1),col=c(1,2,4),cex=1.3)


####################


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
rowMember1=numeric(83)
for(i in 1:6)
{
  temp=which(rowMember==rowO[i])
  rowMember1[temp]=i
}

#################
MRchange=matrix(0,83,4)
for(i in 1:83){
  MRchange[i,1]=as.numeric(names(rowMember)[i])
  MRchange[i,2]=rowMember1[i]
  MRchange[i,3]=rowMember1[i+83]
  MRchange[i,4]=rowMember1[i+166]
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
colMember1=numeric(83)
for(i in 1:6)
{
  temp=which(colMember==colO[i])
  colMember1[temp]=i
}


CRchange=matrix(0,83,4)
for(i in 1:83){
  CRchange[i,1]=as.numeric(names(rowMember)[i])
  CRchange[i,2]=colMember1[i]
  CRchange[i,3]=colMember1[i+83]
  CRchange[i,4]=colMember1[i+166]
}

matplot(t(CRchange[,2:4]), t="l", lty=1, las=1, ylab="cluster No",col=1:99,
        xlab="Time", xaxt="n",main="Groom Sent")
axis(1, at=1:3, labels=bb2)

A8=Transition(MRchange[,2],MRchange[,3])
B8=Transition(MRchange[,3],MRchange[,4])



heatmap.2(A8,Rowv=FALSE,Colv = FALSE, main="Pert vs Base",trace="none")
heatmap.2(B8,Rowv=FALSE,Colv = FALSE, main="Post vs Pert",trace="none")

#######################











#########################
#create DCGtree
library(Perc)
win81=conductance(NC8Groom1,maxLength = 4)
win_prob81=win81$p.hat

temp=c(0.033,0.2,1)
Ens.groom81=Eigen.plot(temp, selected.id=c(1,2,3),win_prob81)
DCG.groom81=DCGtree.plot(num.clusters.selected=c(1,2,4),
                         "NC8GroomR1 tree",Ens.groom81,temp)

plot(DCG.groom81,hang=-1,main="NC8 GroomR1 tree")

############################################
win82=conductance(NC8Groom2,maxLength = 4)
win_prob82=win82$p.hat

temp=c(0.045,0.15,1)
Ens.groom82=Eigen.plot(temp, selected.id=c(1,2,3),win_prob82)
DCG.groom82=DCGtree.plot(num.clusters.selected=c(1,2,3),
                        "NC8GroomR2 tree",Ens.groom82,temp)

plot(DCG.groom82,hang=-1,main="NC8 GroomR2 tree")

################################
win83=conductance(NC8Groom3,maxLength = 4)
win_prob83=win83$p.hat

temp=c(0.065,0.15,1)
Ens.groom83=Eigen.plot(temp, selected.id=c(1,2,3),win_prob83)
DCG.groom83=DCGtree.plot(num.clusters.selected=c(1,2,3),
                         "NC8GroomR3 tree",Ens.groom83,temp)

plot(DCG.groom83,hang=-1,main="NC8 GroomR3 tree")


################

MatrixSym=function(Mat)
{
  n=nrow(Mat)
  Matrix=matrix(0,n,n)
  for(i in 1:n)
  {
    for(j in 1:n)
    {
      Matrix[i,j]=Mat[i,j]+Mat[j,i]
    }
  }
  
  return(Matrix)
}

DistPairwise=function(Mat)
{
  n=nrow(Mat)
  Dist=matrix(0,n,n)
  for(i in 1:n)
  {
    for(j in 1:n)
    {
      Dist[i,j]=sqrt(dist(rbind(Mat[i,],Mat[j,]))^2-(Mat[i,j])^2)
    }
  }
  
  return(Dist)
}
Grooming8=MatrixSym(NC8Groom1)

GroomDist8=DistPairwise(Grooming8)


########################3
temp=c(2,4,8,20,50)
Ens.groom81=Eigen.plot(temp, selected.id=c(1,2,3,4,5),GroomDist8)

DCG.groom81=DCGtree.plot(num.clusters.selected=c(1,2,3,4,5),
                         "NC8GroomR1 tree",Ens.groom81,temp)

plot(DCG.groom81,hang=-1,main="Grooming NC8")

################
#cut the three
yNC8G1=cutree(DCG.groom81,k=8)








#####################################
#combine the huddling and grooming 
#####################################

NC8Comb=NC8Groom1+NC8Huddel1

###########3
sortMonkey=function(Mat,name)
{
  n=nrow(Mat)
  Mat2=matrix(0,n,n)
  order(name)
}
