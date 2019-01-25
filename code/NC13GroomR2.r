NC13_RU2_Grooming_Matrix <- read.csv("~/Dropbox/Research/SNH_health profile data for Fushing-selected/NC13_RU2_Grooming_Matrix.csv")
NC13GroomR2=as.matrix(NC13_RU2_Grooming_Matrix [,-1])
colnames(NC13GroomR2)=NC13_RU2_Grooming_Matrix[,1]
rownames(NC13GroomR2)=NC13_RU2_Grooming_Matrix[,1]
##############

win2=conductance(NC13GroomR2,maxLength = 4)
win_prob2=win2$p.hat

temp=c(0.0045,0.07,0.15,0.3,0.8,1)
Ens13.groom2=Eigen.plot2(temp, selected.id=c(1,2,3,4,5,6),win_prob2)
DCG13.groom2=DCGtree.plot(num.clusters.selected=c(1,1,2,3,4,6),
                        "NC13GroomR2 tree",Ens13.groom2,temp)
plot(DCG13.groom2,hang=-1,main="NC13GroomR2 tree")



G2=cutree(DCG.groom2,k=5)



















############################
Eigen.plot2=function(tempinv,selected.id,D){
  tempinv.selected <- tempinv[selected.id]
  ensM<- list()    # your ensemble matrices at each temperature.
  for ( i in 1:length(selected.id))
    ensM[[i]]=EstClust(GetSim2(D,tempinv.selected[i]), MaxIt=1000, m=5)
  
  #check eigenvalues
  par(mfrow=c(2,3))
  for (j in 1:length(selected.id)){
    Ens=ensM[[j]]
    N <- nrow(Ens)
    Dinvsqrt <- diag(sapply(1:N, function(i) 1/sqrt(sum(Ens[i,]))))
    Lsym <- diag(N) - Dinvsqrt %*% Ens %*% Dinvsqrt
    Eigen <- eigen(Lsym)$values
    Eigen <- sort(1 - Eigen/Eigen[1], decreasing=TRUE)
    #cat(Eigen[1:25],"\n")
    cat("difference","\n",diff(Eigen[1:20]),"\n")
    plot(Eigen[1:20],type="b",main=j)
  }

  return(ensM)
}