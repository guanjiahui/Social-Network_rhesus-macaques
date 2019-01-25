#############
#process the data 
NC13_RU1_Huddling_Matrix <- read.csv("~/Dropbox/Research/SNH_health profile data for Fushing-selected/NC13_RU1_Huddling_Matrix.csv", header=FALSE)
NC13HuddleR1=as.matrix(NC13_RU1_Huddling_Matrix[-1,-1])

colnames(NC13HuddleR1)=NC13_RU1_Huddling_Matrix[-1,1]
rownames(NC13HuddleR1)=NC13_RU1_Huddling_Matrix[-1,1]

HuddleR1=matrix(0,nrow(NC13HuddleR1),ncol(NC13HuddleR1))
for (i in 1:ncol(NC13HuddleR1)){
  for (j in 1:ncol(NC13HuddleR1)){
    HuddleR1[i,j]=NC13HuddleR1[i,j]+NC13HuddleR1[j,i]
  }
}    

NC13HuddleR1=HuddleR1


####################
library(network)
GhuddleR1=network(NC13HuddleR1,directed=FALSE,matrix.type="adjacency")
plot(GhuddleR1,vertex.col=2+y12,vertex.cex=1.5,interactive=TRUE)

library(igraph)
GHR1=graph.adjacency(NC13HuddleR1,mode="undirected")


plot(GHR1,vertex.color=y12, vertex.frame.color="#ffffff")

l <- layout.random(GHR1)
plot(GHR1,layout=l,edge.arrow.size=.1,vertex.color=y12)
l<-layout.circle(GHR1)
l<-layout.sphere(GHR1)
################


################
#now doing DCG

##################
#NC13HuddleR1_dist=dist(NC13HuddleR1,diag=TRUE,upper=TRUE)
#NC13HuddleR1_dist=as.matrix(NC13HuddleR1_dist)
############
temp=c(0.2,0.3,10,100)
Ens.huddle1=Eigen.plot2(temp, selected.id=c(1,2,3,4),NC13HuddleR1)
DCG.huddle1=DCGtree.plot(num.clusters.selected=c(2,2,17,19),
                         "NC13HuddleR1 tree",Ens.huddle1,temp)
########
heatmap.2(NC13HuddleR1,Rowv=as.dendrogram(DCG.huddle1),Colv=as.dendrogram(DCG.huddle1),
          trace="none",col =colorRampPalette(c("white","green","green4","violet","purple"))(100))

heatmap.2(NC13HuddleR1,col =colorRampPalette(c("white","green","green4","violet","purple"))(100),
          trace="none")
############
plot(DCG.huddle1,hang=-1)
library(sparcl)
# colors the leaves of a dendrogram
y1 = cutree(DCG.huddle1, 17)
y12=cutree(DCG.huddle1,2)
ColorDendrogram(DCG.huddle1, y = y1,  main = "NC13HuddleR1 tree",xlab="", 
                labels=nameR1,branchlength = 1)
#another way to color it 
#den.huddle1=as.dendrogram(DCG.huddle1)
#dend2=cut(den.huddle1,h=3)
#plot(den.huddle1,hang=-1,nodePar = list(col=1:2))
#####################


#########################
s1=sum(NC13HuddleR1)
s2=sum(NC13HuddleR2)
s3=sum(NC13HuddleR3)
NC13HuddleR1.del=NC13HuddleR1[-c(37,90),]
NC13HuddleR1.del=NC13HuddleR1.del[,-c(37,90)]
s1.del=sum(NC13HuddleR1.del)

eachR1=colSums(NC13HuddleR1)
eachR2=colSums(NC13HuddleR2)
eachR3=colSums(NC13HuddleR3)
eachR1.del=colSums(NC13HuddleR1.del)
eachSum=cbind(eachR1.del,eachR2,eachR3)

nameR1=colnames(NC13HuddleR1)
nameR2=colnames(NC13HuddleR2)
nameR3=colnames(NC13HuddleR3)
nameR1.del=nameR1[-c(37,90)]

###############
plot(c(s1,s2,s3),type="b")
plot(c(s1.del,s2,s3),type="b",ylab="Huddle",ylim=c(min(s3),max(s1)),
     xaxt="n",main="NC13 Huddling Total")
axis(1, at=1:3, labels=bb2)

matplot(t(eachSum), t="l", lty=1, las=1, ylab="Huddle",
        xlab="Time", xaxt="n",main="NC13 Huddling")
bb2=c("baseline","pertubation","postpertubation")
axis(1, at=1:3, labels=bb2)
############
#put them into 17 groups
matplot(t(eachSum), t="l", lty=1, las=1, ylab="Huddle",col=y1,
        xlab="Time", xaxt="n",main="NC13")
bb2=c("baseline","pertubation","postpertubation")
axis(1, at=1:3, labels=bb2)
############
#put them into 2 groups
matplot(t(eachSum), t="l", lty=1, las=1, ylab="Huddle",col= cutree(DCG.huddle1,2),
        xlab="Time", xaxt="n",main="NC13")
bb2=c("baseline","pertubation","postpertubation")
axis(1, at=1:3, labels=bb2)
############

Huddle=data.frame(cbind(eachSum,t(t(y2))))
pairs(~eachR1.del+SeachR1.del+eachR2 + SeachR2+eachR3+ SeachR3,data=Groom, 
      main="Simple Scatterplot Matrix")

plot(density(eachR1.del),ylim=c(0,0.04),main="Huddling density",lwd=3)
lines(density(eachR2),col=2,lwd=3)
lines(density(eachR3),col=4,lwd=3)
abline(v=eachR1[37],col="yellow",lwd=3,lty=2)
abline(v=eachR1[90],col="green",lwd=3,lty=2)

legend("topright",c("baseline","perturbation","postperturbation"),
       lty=c(1,1,1),col=c(1,2,4),cex=1.3)
###################################
#next to compute the entropy 



Entropy=function(m,N){
  p=m/N
  if (p==0){
    E=0
  }
  else
    E=-log(p)*p
  return(E)
}

##############
Entropy_sequence=function(Entry1,Entry2,name1,name2){
  #a function that assume the elements of each trees(network) has the same location
  #and thus has the same total number of elements. That is, each element is matched. 
  ###########
  
  a=length(unique(Entry1))
 # b=length(unique(Entry2))
  entropy=numeric(a)
  SumEntropy=0
  N_total=length(intersect(name1,name2))
  
  for (i in 1:a){
    location=which(Entry1==i)
    N=length(location)
    membership=Entry2[location]
    M=unique(membership)
    
    for (j in 1:length(M)){
      m=length(which(membership==M[j]))
      temp=Entropy(m,N)
      entropy[i]=temp+entropy[i]
    }#end for j
    #cat("member",membership, "entropy",entropy[i],"\n")
    SumEntropy=SumEntropy+(N/N_total)*entropy[i]
    #cat("N",N,"\n")
  }#end for i
  #cat(SumEntropy,sum(entropy),"\n")
  Entropy_bottom=Entropy_de(Entry2)$Bottom
  
  return(list(Entropy=entropy,Sum=SumEntropy/Entropy_bottom))
  
  
}#Entropy_sequence() function 
###############################
Entropy_de=function(Entry1){
  a=length(unique(Entry1))
  entropy=0
  for(i in 1:a){
    prob=length(which(Entry1==i))/length(Entry1)
    temp=-log(prob)*prob
    entropy=temp+entropy
  }
  
  return(list(Bottom=entropy))
}

name1=colnames(NC13HuddleR1)
name2=colnames(NC13HuddleR2)
name3=colnames(NC13HuddleR3)





#################################
Entry1=y1[-c(37,90)]
Entry2=y2
Entry3=y3
Entry12=y12[-c(37,90)]
Entry22=y22
Entry32=y32

SM12=Entropy_sequence(Entry1,Entry2)$Entropy
SM13=Entropy_sequence(Entry1,Entry3)$Entropy

SM21=Entropy_sequence(Entry2,Entry1)$Entropy
SM23=Entropy_sequence(Entry2,Entry3)$Entropy

SM31=Entropy_sequence(Entry3,Entry1)$Entropy
SM32=Entropy_sequence(Entry3,Entry2)$Entropy

########
par(mfrow=c(3,2))
plot(SM12,type="b",
     xlab="group", ylab="Entropy",main="Entropy of Baseline vs Perturbation")
plot(SM21,type="b",
     xlab="group", ylab="Entropy",main="Entropy of Perturbation vs Baseline")

plot(SM23,type="b",
     xlab="group", ylab="Entropy",main="Entropy of perturbation vs post-pert ")

plot(SM32,type="b",
     xlab="group", ylab="Entropy",main="Entropy of Post-pert vs perturbation")

plot(SM13,type="b",
     xlab="group", ylab="Entropy",main="Entropy of Baseline vs post-pert")

plot(SM31,type="b",
     xlab="group", ylab="Entropy",main="Entropy of Post-pert vs Baseline")
##############
Entrpy=list()
Entrpy[[1]]=y1[-c(37,90)]
Entrpy[[2]]=y2
Entrpy[[3]]=y3
Entrpy[[4]]=y12[-c(37,90)]
Entrpy[[5]]=y22
Entrpy[[6]]=y32

Name=list(name1,name2,name3)

En=matrix(0,3,3)
for (i in 1:3){
  for (j in 1:3){
    if (i!=j)
    En[i,j]=Entropy_sequence(Entrpy[[i]],Entrpy[[j]],Name[[i]],Name[[j]])$Sum
  }
}

for (i in 1:3)
  for (j in 1:3)
    cat(En[i,j]+En[j,i],"i",i,"j",j,"\n")
###############

#######
bSM12=Entropy_sequence(Entry12,Entry22)$Entropy
bSM13=Entropy_sequence(Entry12,Entry32)$Entropy

bSM21=Entropy_sequence(Entry22,Entry12)$Entropy
bSM23=Entropy_sequence(Entry22,Entry32)$Entropy

bSM31=Entropy_sequence(Entry32,Entry12)$Entropy
bSM32=Entropy_sequence(Entry32,Entry22)$Entropy

par(mfrow=c(3,2))
plot(bSM12,type="b",ylim=c(0,1.0),
     xlab="group", ylab="Entropy",main="Entropy of Baseline vs Perturbation")
plot(bSM21,type="b",ylim=c(0,1.0),
     xlab="group", ylab="Entropy",main="Entropy of Perturbation vs Baseline")

plot(bSM23,type="b",ylim=c(0,1.0),
     xlab="group", ylab="Entropy",main="Entropy of perturbation vs post-pert ")

plot(bSM32,type="b",ylim=c(0,1.0),
     xlab="group", ylab="Entropy",main="Entropy of Post-pert vs perturbation")

plot(bSM13,type="b",ylim=c(0,1.0),
     xlab="group", ylab="Entropy",main="Entropy of Baseline vs post-pert")

plot(bSM31,type="b",ylim=c(0,1.0),
     xlab="group", ylab="Entropy",main="Entropy of Post-pert vs Baseline")


En2=matrix(0,3,3)
for (i in 1:3){
  for (j in 1:3){
    if (i!=j)
      En2[i,j]=Entropy_sequence(Entrpy[[i+3]],Entrpy[[j+3]],Name[[i]],Name[[j]])$Sum
  }
}

for (i in 1:3)
  for (j in 1:3)
    cat(En2[i,j]+En2[j,i],"i",i,"j",j,"\n")
###################

















#####################
#normalize the adjacency matrix 

EN1=EN
EN1[which(EN>5)]=5

Eheat1=NC13HuddleR1
small1=which(NC13HuddleR1<5,arr.ind=TRUE)
Eheat1[small1]=NC13HuddleR1[small1]/5
Eheat1[which(NC13HuddleR1>=5,arr.ind=TRUE)]=1
#######
temp=c(0.2,0.3,2,100)
Ens.heat1=Eigen.plot2(temp, selected.id=c(1,2,3,4),Eheat1)
DCG.heat1=DCGtree.plot(num.clusters.selected=c(2,2,15,17),
                       "NC13HuddleR1 tree",Ens.heat1,temp)
########
heatmap.2(Eheat1,col =colorRampPalette(c("white","green","green4","violet","purple"))(100),
          trace="none")


heatmap.2(Eheat1,Rowv=as.dendrogram(DCG.heat1),Colv=as.dendrogram(DCG.heat1),
          trace="none",
          col =colorRampPalette(c("white","green","green4","violet","purple"))(100))
############################




#d=heatmap(Eheat1)
#HC1=Eheat1[d$rowInd,d$colInd]
#D1=Eheat1[DCG.heat1$order,DCG.heat1$order]
#GetBipEnergy(HC1)
#GetBipEnergy(D1)

plot(DCG.heat1,hang = -1)
y1 = cutree(DCG.huddle1, 15)
y12=cutree(DCG.huddle1,2)
ColorDendrogram(DCG.huddle1, y = y1,  main = "NC13HuddleR1 tree",xlab="", 
                branchlength = 5)


######
#double check Temperature selection
temp=c(0.1,0.2,0.3,0.5,0.8,1,2,2000)
Ens.h1=Eigen.plot2(temp, selected.id=c(1,2,3,4,5,6,7,8),Eheat1)



#################
#visualize the entropy 
insertE=function(ary,ind,value){
  afterInsert=ary
  for (i in 1:length(ind)){
    temp=afterInsert[ind[i]:length(afterInsert)]
    afterInsert[ind[i]]=value[i]
    afterInsert[(ind[i]+1):(length(afterInsert)+1)]=temp
  }
  return(afterInsert)
}

y2.append=insertE(y2,c(37,90),c(0,0))

y3.append=insertE(y3,c(37,90),c(0,0))

ColorDendrogram(DCG.huddle1, y = y3.append, 
                main = "NC13HuddleR1 tree (post)",xlab="", 
                branchlength = 3)

#####################

NameHuddleR1=colnames(NC13HuddleR1)
NameHuddleR2=colnames(NC13HuddleR2)
NameHuddleR3=colnames(NC13HuddleR3)

save(DCG.huddle1,DCG.huddle2,DCG.huddle3,
     NameHuddleR1,NameHuddleR2,NameHuddleR3,
     file = "HuddlingTREE.RData")




















##############################
Eigen.plot2=function(tempinv,selected.id,D){
  tempinv.selected <- tempinv[selected.id]
  ensM<- list()    # your ensemble matrices at each temperature.
  for ( i in 1:length(selected.id))
    ensM[[i]]=EstClust(GetSim2(D,tempinv.selected[i]), MaxIt=1000, m=5)
  
  #check eigenvalues
  par(mfrow=c(2,2))
  for (j in 1:length(selected.id)){
    Ens=ensM[[j]]
    N <- nrow(Ens)
    Dinvsqrt <- diag(sapply(1:N, function(i) 1/sqrt(sum(Ens[i,]))))
    Lsym <- diag(N) - Dinvsqrt %*% Ens %*% Dinvsqrt
    Eigen <- eigen(Lsym)$values
    Eigen <- sort(1 - Eigen/Eigen[1], decreasing=TRUE)
    #cat(Eigen[1:25],"\n")
   # cat("difference",diff(Eigen[1:20]),"\n")
   plot(Eigen[1:15],type="b",main=paste(j,"T=",tempinv.selected[j]))
   #barplot(Eigen[1:25],main=j)
   # plot(diff(Eigen[1:20]),type="b")
  }

  #for (j in 1:length(selected.id))
  #  heatmap(ensM[[j]],main=j)
  
  return(ensM)
}

