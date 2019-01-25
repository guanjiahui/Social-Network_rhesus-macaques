NC8_BL_EColi_similarity <- read.csv("C:/Users/jiahui/Dropbox/Research/SNH_health profile data for Fushing-selected/DATA/NC8/NC8_BL_EColi_similarity.csv", header=FALSE)
nameEcol81=NC8_BL_EColi_similarity[1,-1]
NC8Ecol1=NC8_BL_EColi_similarity[-1,]
NC8Ecol1=NC8Ecol1[,-1]
for(j in 1:nrow(NC8Ecol1))
{
  for(i in 1:j)
  {
    NC8Ecol1[i,j]=NC8Ecol1[j,i]
  }
}

NC8Ecol1=as.matrix(NC8Ecol1)

aaa=NC8Ecol1[which(NC8Ecol1>0)]
quantile(aaa,0.01)

#############
#standardize NC8Ecoli since its distribution are from (0,100). As a similarity matrix input, we need 
#to keep all the entries within 1. 
NC8Ecol1=NC8Ecol1/100
colnames(NC8Ecol1)=nameEcol81


temp=c(0.13,0.25,1)
Ens.8Ecoli1=Eigen.plot2(temp, selected.id=c(1,2,3),NC8Ecol1)
###########
DCG.8Ecoli1=DCGtree.plot(num.clusters.selected=c(1,2,4),
                         "NC8 Ecoli Baseline",Ens.8Ecoli1,temp)
plot(DCG.8Ecoli1,hang=-1,main="NC8 Ecoli Baseline")

#heatmap(NC8Ecol1,Rowv=as.dendrogram(DCG.8Ecoli1),Colv = as.dendrogram(DCG.8Ecoli1),trace="none")
#heatmap(NC8Ecol1)

yNC8E1= cutree(DCG.8Ecoli1, 6)
ColorDendrogram(DCG.8Ecoli1, y = yNC8E1, main = "NC8 Ecoli Baseline", 
                branchlength = 80)


heatmap.2(NC8Ecol1, Rowv=as.dendrogram(DCG.8Ecoli1), 
          Colv = as.dendrogram(DCG.8Ecoli1), trace="none",
          main="NC8 Ecoli Baseline")

###################
#####
#baseline tree for Huddle and Groom
#DCG8.huddle1
#DCG.groom81
yNC8E1= cutree(DCG.8Ecoli1, 6)
yNC8H1=cutree(DCG8.huddle1,4)
yNC8G1=cutree(DCG.groom81,5)

nameE8=list()
nameE8[[1]]=colnames(NC8Ecol1)
nameE8[[2]]=colnames(NC8Huddel1)
nameE8[[3]]=colnames(NC8Groom1)
nameE8[[4]]=colnames(NC8_Agg)
#######################

NC8Comp1=TriEntropy(yNC8E1,yNC8H1,yNC8G1,nameE8)

###########
#save in excel file 

aa=intersect(nameE8[[1]],nameE8[[2]])
common=intersect(aa,nameE8[[3]])

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
  loc[[i]]=Find(nameE8[[i]],common)
}

y8Ecoli1=yNC8E1[loc[[1]]]
y8Huddle1=yNC8H1[loc[[2]]]
y8Groom1=yNC8G1[loc[[3]]]
ID8=as.numeric(nameE8[[1]][loc[[1]]])
y8Baseline=data.frame(cbind(ID8,y8Ecoli1,y8Huddle1,y8Groom1))

write.csv(y8Baseline,file="NC8BaselineMembership.csv")


a8=order(y8Baseline$y8Ecoli1)
y8BaseSorted=y8Baseline[a8,]


write.csv(y8BaseSorted,file="NC8Base_Sorted.csv")




######################
TriEntropy=function(y1,y2,y3,nameE)
{
  
  EntryE=list()
  EntryE[[1]]=y1
  EntryE[[2]]=y2
  EntryE[[3]]=y3
  
  
  par(mfrow=c(2,3))
  EnE=matrix(0,3,3)
  for (i in 1:3){
    for (j in 1:3){
      if (i!=j){
        EnE[i,j]=Entropy_general(EntryE[[i]],EntryE[[j]],NameE[[i]],NameE[[j]])$Sum
        plot(Entropy_general(EntryE[[i]],EntryE[[j]],NameE[[i]],NameE[[j]])$Entropy, type="b",
             xlab="group", ylab="Entropy",main=paste("stage",i,"from","stage",j))
      }
    }
  }
  
  for (i in 1:3)
    for (j in 1:3)
      cat(EnE[i,j]+EnE[j,i],"i",i,"j",j,"\n")
  
  return(EnE)
}
############################