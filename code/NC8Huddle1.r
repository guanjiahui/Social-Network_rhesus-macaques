##########################################################
MatrixTransform=function(Initiator, Recipient)
{
  ID1=unique(Initiator)
  ID2=unique(Recipient)
  ID=union(ID1,ID2)
  
  n=length(ID)
  N=length(Initiator)
  
  Mat=matrix(0,n,n)
  Matrix=matrix(0,n,n)
  
  colnames(Matrix)=ID
  rownames(Matrix)=ID
  
  for(i in 1:N)
  {
    I=which(ID==Initiator[i])
    J=which(ID==Recipient[i])
    Mat[I,J]=Mat[I,J]+1
  }
  
  for(i in 1:n)
  {
    for(j in 1:n)
    {
      Matrix[i,j]=Mat[i,j]+Mat[j,i]
    }
  }
  
  return(Matrix)
  
}
#############

Mat2=as.conflictmat(NC8_POST_Huddling[-572,1:2])
MatR1=matrix(0,ncol(Mat2),ncol(Mat2))
for (i in 1:ncol(Mat2)){
  for (j in 1:ncol(Mat2)){
    MatR1[i,j]=Mat2[i,j]+Mat2[j,i]
  }
}  
###############
NC8Huddel1=MatrixTransform(NC8_BL_Huddling$Initiator,NC8_BL_Huddling$Recipient)
NC8Huddel2=MatrixTransform(NC8_PERT_Huddling$Initiator,NC8_PERT_Huddling$Recipient)
NC8Huddel3=MatrixTransform(NC8_POST_Huddling$Initiator,NC8_POST_Huddling$Recipient)
############
#check recipient and intiator are the same (i,e, diagonal is nonzero)
for(i in 1:ncol(NC8Huddel3))
{
  for(j in 1:ncol(NC8Huddel3))
  {
    if(i==j)
      cat(NC8Huddel3[i,j],i,"\n")
  }
}

NC8Huddel3[7,7]=0

###############
temp=c(0.3,1,6,100)
Ens8.huddle1=Eigen.plot2(temp, selected.id=c(1,2,3,4),NC8Huddel1)

DCG8.huddle1=DCGtree.plot(num.clusters.selected=c(1,1,4,4),
                          "NC8 HuddleR1 tree",Ens8.huddle1,temp)
########

temp=c(0.3,1,6,100)
Ens8.huddle2=Eigen.plot2(temp, selected.id=c(1,2,3,4),NC8Huddel2)

DCG8.huddle2=DCGtree.plot(num.clusters.selected=c(2,2,5,8),
                          "NC8 HuddleR1 tree",Ens8.huddle2,temp)
########


temp=c(0.3,0.8,6,120)
Ens8.huddle3=Eigen.plot2(temp, selected.id=c(1,2,3,4),NC8Huddel3)

DCG8.huddle3=DCGtree.plot(num.clusters.selected=c(2,2,3,7),
                          "NC8 HuddleR1 tree",Ens8.huddle3,temp)
########

plot(DCG8.huddle1,hang=-1,main= "NC8 HuddleR1 tree")
plot(DCG8.huddle2,hang=-1,main= "NC8 HuddleR2 tree")
plot(DCG8.huddle3,hang=-1,main= "NC8 HuddleR3 tree")

##############################

library(sparcl)
# colors the leaves of a dendrogram
y1 = cutree(DCG8.huddle1, 14)
ColorDendrogram(DCG8.huddle1, y = y1,  main = "NC8 Huddle R1 tree",xlab="", 
                branchlength = 0.2)

y2 = cutree(DCG8.huddle2, 12)
ColorDendrogram(DCG8.huddle2, y = y2,  main = "NC8 Huddle R2 tree",xlab="", 
                branchlength = 0.8)

y3 = cutree(DCG8.huddle3, 6)
ColorDendrogram(DCG8.huddle3, y = y3,  main = "NC8 Huddle R3 tree",xlab="", 
                branchlength = 0.8)

######################
#> dim(NC8Huddel1)
#[1] 96 96
#> dim(NC8Huddel2)
#[1] 84 84
#> dim(NC8Huddel3)
#[1] 83 83
#
name8=list()
name8[[1]]=colnames(NC8Huddel1)
name8[[2]]=colnames(NC8Huddel2)
name8[[3]]=colnames(NC8Huddel3)

aa=intersect(name8[[1]],name8[[2]])
common=intersect(aa,name8[[3]])

#####################################
Find=function(ID,subset) #where subset is a subset if ID 
{
  loc=c()
  for(i in 1:length(subset))
  {
    temp=which(ID==subset[i])
    loc=c(temp,loc)
  }
  return(loc)
}

####################################
loc=list()
for(i in 1:3)
{
  loc[[i]]=Find(name8[[i]],common)
}

################################
NC8=list()
NC8[[1]]=NC8Huddel1[loc[[1]],loc[[1]]]
NC8[[2]]=NC8Huddel2[loc[[2]],loc[[2]]]
NC8[[3]]=NC8Huddel3[loc[[3]],loc[[3]]]

eachR1=colSums(NC8[[1]])
eachR2=colSums(NC8[[2]])
eachR3=colSums(NC8[[3]])
eachSum=cbind(eachR1,eachR2,eachR3)

S8=c()
for(i in 1:3)
  S8[i]=sum(NC8[[i]])
###############

plot(S8,type="b",ylab="Huddle",ylim=c(min(S8[[2]]),max(s1)),col=2,lwd=3,pch=16,
     xaxt="n",main=" Huddling Total (in common)")

#####################
#compare with NC13
lines(c(s1.del,s2,s3),lty=2,col=3,lwd=3)
points(c(s1.del,s2,s3),lty=2,col=3,lwd=3,pch=15)
axis(1, at=1:3, labels=bb2)

legend("topright",c("NC8","NC13"),
       lty=c(1,2),col=c(2,3),cex=1.3)







matplot(t(eachSum), t="l", lty=1, las=1, ylab="Huddle",
        xlab="Time", xaxt="n",main="NC8 Huddling")
bb2=c("baseline","pertubation","postpertubation")
axis(1, at=1:3, labels=bb2)
############

#put them into 2 groups
matplot(t(eachSum), t="l", lty=1, las=1, ylab="Huddle",col= cutree(DCG.huddle1,2),
        xlab="Time", xaxt="n",main="NC13")
bb2=c("baseline","pertubation","postpertubation")
axis(1, at=1:3, labels=bb2)
############



plot(density(eachR1),ylim=c(0,0.045), xlim=c(-15,105),main="Huddling density",lwd=3)
lines(density(eachR2),col=2,lwd=3)
lines(density(eachR3),col=4,lwd=3)

legend("topright",c("baseline","perturbation","postperturbation"),
       lty=c(1,1,1),col=c(1,2,4),cex=1.3)
###################################





#################
#compute the mutual entropy
EntryE=list()
EntryE[[1]]=y1
EntryE[[2]]=y2
EntryE[[3]]=y3


par(mfrow=c(2,3))
EnE=matrix(0,3,3)
for (i in 1:3){
  for (j in 1:3){
    if (i!=j){
      EnE[i,j]=Entropy_sequence(EntryE[[i]],EntryE[[j]],name8[[i]],name8[[j]])$Sum
      plot(Entropy_sequence(EntryE[[i]],EntryE[[j]],name8[[i]],name8[[j]])$Entropy, type="b",
           xlab="group", ylab="Entropy",main=paste("stage",i,"from","stage",j))
    }
  }
}

for (i in 1:3)
  for (j in 1:3)
    cat(EnE[i,j]+EnE[j,i],"i",i,"j",j,"\n")

#################
