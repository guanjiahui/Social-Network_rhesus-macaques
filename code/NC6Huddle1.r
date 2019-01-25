NC6_BL_EColi_similarity <- read.csv("~/Dropbox/Research/SNH_health profile data for Fushing-selected/DATA/NC6/NC6_BL_EColi_similarity.csv", header=FALSE)
nameEcol61=NC6_BL_EColi_similarity[1,-1]
NC6Ecol1=NC6_BL_EColi_similarity[-1,]
NC6Ecol1=NC6Ecol1[,-1]
for(j in 1:nrow(NC6Ecol1))
{
  for(i in 1:j)
  {
    NC6Ecol1[i,j]=NC6Ecol1[j,i]
  }
}

NC6Ecol1=as.matrix(NC6Ecol1)

aaa=NC6Ecol1[which(NC6Ecol1>0)]
quantile(aaa,0.01)

#############
#standardize NC8Ecoli since its distribution are from (0,100). As a similarity matrix input, we need 
#to keep all the entries within 1. 
NC6Ecol1=NC6Ecol1/100
colnames(NC6Ecol1)=nameEcol61


temp=c(0.1,0.18,0.25,1)
Ens.6Ecoli1=Eigen.plot2(temp, selected.id=c(1,2,3,4),NC6Ecol1)
###########
DCG.6Ecoli1=DCGtree.plot(num.clusters.selected=c(1,2,5,12),
                         "NC6 Ecoli Baseline",Ens.6Ecoli1,temp)
plot(DCG.6Ecoli1,hang=-1,labels=colnames(NC6Ecol1),main="NC6 Ecoli Baseline")

heatmap.2(NC6Ecol1*100,Rowv=as.dendrogram(DCG.6Ecoli1),Colv = as.dendrogram(DCG.6Ecoli1),
          trace="none",main="NC6 Ecoli")
heatmap.2(NC6Ecol1)



##################################################
#check Membership Methods
SingleTemp=EstClust(GetSim(NC6Ecol1,0.1))

hist(SingleTemp,breaks=50)
heatmap.2(SingleTemp,Rowv = FALSE, Colv = FALSE,trace="none",
          main = "Singele Temperature Matrix")
heatmap.2(SingleTemp,trace="none",
          main = "Single Temperature Matrix")

ind1=which(SingleTemp<0.4,arr.ind = FALSE)
ind2=which(SingleTemp>=0.4,arr.ind = FALSE)

S1=SingleTemp
S1[ind1]=0
S1[ind2]=1

heatmap.2(S1,Rowv = FALSE, Colv = FALSE,trace="none",
          main = "Singele Temperature after")
heatmap.2(S1,trace="none",
          main = "Single Temperature after")

##################################################
N=nrow(SingleTemp)
Dinvsqrt <- diag(sapply(1:N, function(i) 1/sqrt(sum(SingleTemp[i,]))))
Lsym <- diag(N) - Dinvsqrt %*% SingleTemp %*% Dinvsqrt
Eigen <- eigen(Lsym)$values
Eigen <- sort(1 - Eigen/Eigen[1], decreasing=TRUE)
plot(Eigen[1:10],type="b")


##################################################
#Huddling Analysis
##################################################
NC6_BL_Huddling <- read.csv("~/Dropbox/Research/SNH_health profile data for Fushing-selected/DATA/NC6_BL_Huddling.csv")

Monkey6=unique(NC6_BL_Huddling$Initiator)
#create the huddle Tree 
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
####################################

NC6Huddel1=MatrixTransform(NC6_BL_Huddling$Initiator,NC6_BL_Huddling$Recipient)

temp=which(NC6Huddel1>10,arr.ind = TRUE)
NC6Huddel1=NC6Huddel1/10
NC6Huddel1[temp]=1

############
temp=c(0.18,0.25,2,50)
Ens6.huddle1=Eigen.plot2(temp, selected.id=c(1,2,3,4),NC6Huddel1)
DCG6.huddle1=DCGtree.plot(num.clusters.selected=c(1,2,7,13),
                         "NC6 HuddleR1 tree",Ens6.huddle1,temp)
########
plot(DCG6.huddle1,labels=colnames(NC6Huddel1),main="NC6 Huddle Baseline")

heatmap.2(NC6Huddel1,Rowv=as.dendrogram(DCG6.huddle1),Colv = as.dendrogram(DCG6.huddle1),
          trace="none",main="NC6 Ecoli")

#########################
heatmap.2(NC6Huddel1,Rowv=as.dendrogram(DCG6.huddle1),Colv = as.dendrogram(DCG6.huddle1),
          trace="none")



###############################
temp=c(0.18,0.25,2,50)
Ensemble6Huddle=EnsembleExtract(num.clusters.selected=c(1,2,7,13),
                             Ens6.huddle1,temp)


heatmap.2(Ensemble6Huddle,Rowv = FALSE, Colv = FALSE,trace="none",
          main = "Ensemble Matrix")

heatmap.2(Ensemble6Huddle,trace="none",
          main = "Ensemble Matrix")
hist(Ensemble6Huddle,breaks=100)
#############################
#Grooming Analysis 
######################
NC6_BL_Grooming <- read.csv("~/Dropbox/Research/SNH_health profile data for Fushing-selected/DATA/NC6/NC6_BL_Grooming.csv")
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
###############
NC6Groom1=matrixNoSym(NC6_BL_Grooming$Initiator,NC6_BL_Grooming$Recipient,NC6_BL_Grooming$Interaction)

for(i in 1:ncol(NC6Groom1))
{
  for(j in 1:ncol(NC6Groom1))
  {
    if(i==j)
      cat(NC6Groom1[i,j],i,"\n")
  }
}

#################
library(Perc)
win61=conductance(NC6Groom1,maxLength = 4)
win_prob61=win61$p.hat

temp=c(0.02,0.033,0.2,1)
Ens.groom61=Eigen.plot(temp, selected.id=c(1,2,3,4),win_prob61)
DCG.groom61=DCGtree.plot(num.clusters.selected=c(1,2,3,4),
                         "NC6GroomR1 tree",Ens.groom61,temp)

plot(DCG.groom61,labels=colnames(NC6Groom1),main="NC6 GroomR1 tree")

heatmap.2(NC6Groom1,Rowv=as.dendrogram(DCG.groom61),Colv = as.dendrogram(DCG.groom61), 
          col=rainbow(100),
          trace="none",main="NC6 Ecoli")

heatmap.2(win_prob61,Rowv=as.dendrogram(hclust(as.dist(win_prob61))),
          Colv = as.dendrogram(hclust(as.dist(win_prob61))),
          trace="none",main="NC6 Ecoli")
###############################









################################
#entropy Comparison 
##################################################
yNC6E1= cutree(DCG.6Ecoli1, 18)
yNC6H1=cutree(DCG6.huddle1,9)
yNC6G1=cutree(DCG.groom61,8)

nameE6=list()
nameE6[[1]]=colnames(NC6Ecol1)
nameE6[[2]]=colnames(NC6Huddel1)
nameE6[[3]]=colnames(NC6Groom1)
nameE6[[4]]=colnames(NC6_Agg)
#######################

NC6Comp1=TriEntropy( yNC6E1,yNC6H1,yNC6G1,nameE6)

################################


##########################
#Withn and Between Ecoli Effect 
##########################
Huddle6Effect=EcoliEffect(yNC6H1,nameE6[[2]],nameE6[[1]],NC6Ecol1*100)$E
Groom6Effect=EcoliEffect(yNC6G1,nameE6[[3]],nameE6[[1]],NC6Ecol1*100)$E
############
Tree6H1=as.dendrogram(DCG6.huddle1)
Tree6H1U=cut(Tree6H1,h=1)$upper
plot(Tree6H1U)

Tree6G1=as.dendrogram(DCG.groom61)
Tree6G1U=cut(Tree6G1,h=10)$upper
plot(Tree6G1U)


heatmap.2(Huddle6Effect,Rowv = Tree6H1U, Colv = Tree6H1U,trace="none",
          main="NC6 Baseline Huddling")

heatmap.2(Groom6Effect,Rowv =Tree6G1U, Colv =Tree6G1U,trace="none",
          main="NC6 Baseline Grooming")

####################



#######
#save in excel file 
###########
y6Ecoli1=yNC6E1[loc[[1]]]
y6Huddle1=yNC6H1[loc[[2]]]
y6Groom1=yNC6G1[loc[[3]]]
ID6=as.numeric(nameE[[1]][loc[[1]]])
y6Baseline=data.frame(cbind(ID6,y6Ecoli1,y6Huddle1,y6Groom1))

write.csv(y6Baseline,file="NC6BaselineMembership.csv")

a6=order(y6Baseline$y6Ecoli1)
y6BaseSorted=y6Baseline[a6,]


write.csv(y6BaseSorted,file="NC6Base_Sorted.csv")
#############
library(ROCR)





EcoliEffect =function(y1,name1,ID2,Mat)
{
  N=length(unique(y1))

  check=0
  
  for (i in 1:N)
  {
    loc=which(y1==i)
    ID=name1[loc]
    
    ID=intersect(ID,ID2)
    n=length(ID)
    within=0
    
    if (n==0)
      check=c(check,i)
  }
  
  cat(check,"\n")
  
  if(length(check)>1)
  {
    EcoliWB=matrix(0,N-length(check)+1,N-length(check)+1)
    With=numeric(N-length(check)+1)
    check=check[-1]
    newSeq=(1:N)[-check]
  }

  else
  {
    EcoliWB=matrix(0,N,N)
    With=numeric(N)
    newSeq=1:N
  }
    
  
  
  
  ###########################
  #first compute the diganol
  ###########################
  for(i in newSeq)
  {
    loc=which(y1==i)
    ID=name1[loc]
    
    ID=intersect(ID,ID2)
    n=length(ID)
    within=0
    
    if (n==0)
      check=i
    else{
      for(p in 1:n)
      {
        for(q in 1:n)
        {
          a=which(ID2==ID[p])
          b=which(ID2==ID[q])
          temp=Mat[a,b]
          within=temp+within
        }
      }
      With[i]=length(loc)
      EcoliWB[i,i]=within/(n^2)
    }


  }
  
  ##########################
  #compute the upper off-diaganol
  ##########################
  for(i in newSeq)
  {
    for(j in newSeq)
    {
      if(i!=j)
      {
        lo1=which(y1==i)
        lo2=which(y1==j)
        IDa=name1[lo1]
        IDb=name1[lo2]
        IDa=intersect(IDa,ID2)
        IDb=intersect(IDb,ID2)
        na=length(IDa)
        nb=length(IDb)
        between=0
        
        # cat(na,nb,"i",i,"j",j,"\n")
        
        for(p in 1:na)
        {
          for(q in 1:nb)
          {
            a=which(ID2==IDa[p])
            b=which(ID2==IDb[q])
            temp=Mat[a,b]
            between=temp+between
          }
        }
        
        EcoliWB[i,j]=between/(na*nb)
      }
      
    }
  }
  
  
  
  
  return(list(W=With, E=EcoliWB))
}













####################
#recompute the DCG for Grooming in Cage 6

Grooming6=MatrixSym(NC6Groom1)

GroomDist6=DistPairwise(Grooming6)

temp=c(1,2,4,8,20,50)
Ens.groom61=Eigen.plot(temp, selected.id=c(1,2,3,4,5,6),GroomDist6)
DCG.groom61=DCGtree.plot(num.clusters.selected=c(1,2,3,4,5,6),
                         "NC6GroomR1 tree",Ens.groom61,temp)

plot(DCG.groom61,labels=colnames(NC6Groom1),main="NC6 GroomR1 tree")

plot(DCG.groom61,hang=-1,main="Grooming NC6")

#############
#cut the tree 
yNC6G1=cutree(DCG.groom61,k=8)
