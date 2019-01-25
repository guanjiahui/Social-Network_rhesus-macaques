library(readr)
NC6_BL_GroomHuddlesym <- read_csv("~/Dropbox/Research/VET_SNH_Monkey/DATA/NC6_BL_GroomHuddlesym.csv")
NC8_BL_GroomHuddlesym <- read_csv("~/Dropbox/Research/VET_SNH_Monkey/DATA/NC8_BL_GroomHuddlesym.csv")
NC13_BL_GroomHuddlesym <- read_csv("~/Dropbox/Research/VET_SNH_Monkey/DATA/NC13_BL_GroomHuddlesym.csv")

a=unique(NC6_BL_GroomHuddlesym$Initiator)
b=unique(NC6_BL_GroomHuddlesym$Recipient)
##########
Matrixize=function(dataframe)
{
  ID=unique(dataframe$Initiator)
  n=length(ID)
  N=length(dataframe$Initiator)
  Mat=matrix(0,n,n)
  for(i in 1:N)
  {
    ind1=which(ID==as.numeric(dataframe[i,1]))
    ind2=which(ID==as.numeric(dataframe[i,2]))
    Mat[ind1,ind2]=Mat[ind1,ind2]+1
  }
  colnames(Mat)=ID
  return (Mat)
}
##############
HudGroom6=Matrixize(NC6_BL_GroomHuddlesym)
HudGroom8=Matrixize(NC8_BL_GroomHuddlesym)
HudGroom13=Matrixize(NC13_BL_GroomHuddlesym)
#################

checkSYM=function(Mat)
{
  res=0
  for(i in 1:ncol(Mat))
  {
    for(j in 1:ncol(Mat))
    {
      if(Mat[i,j]!=Mat[j,i])
        res=1
    }
  }
  return(res)
}

###########
checkSYM(HudGroom6)
checkSYM(HudGroom8)
checkSYM(HudGroom13)
#all of them are symmetric, max=1
###############

#################
#now take the pairwise difference and create a new distance matrix for each i,j
#########
DistPairwise=function(Mat)
{
  n=nrow(Mat)
  Dist=matrix(0,n,n)
  for(i in 1:n)
  {
    for(j in 1:n)
    {
      for(k in 1:n)
      {
        temp=0
        if(k!=i &&k!=j)
          temp=(Mat[i,k]-Mat[j,k])^2+temp
      }
      Dist[i,j]=sqrt(temp+(Mat[i,j])^2)
    }
  }
  
  return(Dist)
}
######################
HG6=DistPairwise(HudGroom6)

HG8=DistPairwise(HudGroom8)
HG13=DistPairwise(HudGroom13)

########################3
#now build the DCG Tree 
temp=c(0.1,1,2,20)
Ens6HG=Eigen.plot(temp, selected.id=c(1,2,3,4),HG6)

DCG.HG6=DCGtree.plot(num.clusters.selected=c(1,2,3,4),
                         "NC6Huddle_Groom tree",Ens6HG,temp)

#############
temp=c(0.1,0.7,2,20)
Ens8HG=Eigen.plot(temp, selected.id=c(1,2,3,4),HG8)

DCG.HG8=DCGtree.plot(num.clusters.selected=c(1,2,3,4),
                     "NC8Huddle_Groom tree",Ens8HG,temp)
############
temp=c(0.1,0.7,5)
Ens13HG=Eigen.plot(temp, selected.id=c(1,2,3),HG13)

DCG.HG13=DCGtree.plot(num.clusters.selected=c(1,2,3),
                     "NC13Huddle_Groom tree",Ens13HG,temp)
################


library(gplots)

heatmap.2(HG6,Rowv = as.dendrogram(DCG.HG6),
          Colv = as.dendrogram(DCG.HG6),trace="none",main="Grooming+Huddling")

heatmap.2(HudGroom6,Rowv = as.dendrogram(DCG.HG6),
          Colv = as.dendrogram(DCG.HG6),trace="none")

heatmap.2(HudGroom6,trace="none")

heatmap.2(HudGroom8,Rowv = as.dendrogram(DCG.HG8),
          Colv = as.dendrogram(DCG.HG8),trace="none")

heatmap.2(HudGroom13,Rowv = as.dendrogram(DCG.HG13),
          Colv = as.dendrogram(DCG.HG13),trace="none")

##################
#cut the three
yHG6=cutree(DCG.HG6,k=4)
yHG8=cutree(DCG.HG8,k=6)
yHG13=cutree(DCG.HG13,k=3)

######################################################
#initialize the name in order to use the function
nameE6[[5]]=colnames(HudGroom6)
nameE8[[5]]=colnames(HudGroom8)
nameE13[[5]]=colnames(HudGroom13)
#############


#############################
#within and between Effect 
#############################
HG6Effect=EcoliEffect(yHG6,nameE6[[5]],nameE6[[1]],NC6Ecol1*100)$E
HG8Effect=EcoliEffect(yHG8,nameE8[[5]],nameE8[[1]],NC8Ecol1)$E
HG13Effect=EcoliEffect(yHG13,nameE13[[5]],nameE13[[1]],NC13B)$E


#scale plot (matrix with diferent sizes of columns)
library(ggdendro)
library(grid)
library(reshape)
scalePlot(yHG6, HG6Effect,"NC6 Baseline")
scalePlot(yHG8, HG8Effect,"NC8 Baseline")
scalePlot(yHG13, HG13Effect,"NC13 Baseline")

scalePlot=function(y,Mat,Title)
{
  limits=c(33,100)
  n=dim(Mat)[1]
  colnames(Mat)<-sapply(1:n,function(i) sum(colWidth(y)[1:i]))
  rownames(Mat)<-sapply(1:n,function(i) sum(colWidth(y)[1:i]))
  
  Mat<- melt(Mat)
  Mat <- as.data.frame( Mat)
  names( Mat) <- c("Var1", "Var2", "value")
  v1m <- unique(Mat$Var1)
  Mat$Var1.min <- rep(c(0, v1m[-length(v1m)]), length.out = length(v1m))
  v2m <-unique(Mat$Var2)
  Mat$Var2.min <- rep(c(0, v2m[-length(v2m)]), each = length(v2m))
  
  ggplot(data =  Mat, aes(fill = value)) + 
    geom_rect(aes(ymin = Var1.min, ymax = Var1, xmin = Var2.min, xmax = Var2),colour = "grey50")+ggtitle(Title)+
    scale_fill_continuous(limits=limits,low='white',high='steelblue')
  
}
################################
#Boxplot for within-between effect of Ecoli 
###################################
#Within 
W6HG=sapply(1:nrow(HG6Effect), function(i) HG6Effect[i,i])
W8HG=sapply(1:nrow(HG8Effect), function(i) HG8Effect[i,i])
W13HG=sapply(1:nrow(HG13Effect), function(i) HG13Effect[i,i])
#between
B6HG=BetExtract(HG6Effect)
B8HG=BetExtract(HG8Effect)
B13HG=BetExtract(HG13Effect)


n1=length(W6HG)
n2=length(W8HG)
n3=length(W13HG)

n4=length(B6HG)
n5=length(B8HG)
n6=length(B13HG)

BWHG=data.frame(Group=c(rep("Group I",n1+n4),rep("Group II",n2+n5),rep("Group III",n3+n6)), 
                    Cluster=c(rep("within",n1),rep("between",n4), rep("within",n2), rep("between",n5),rep("within",n3), rep("between",n6)),
                    Value=c(W6HG,B6HG,W8HG,B8HG,W13HG,B13HG))

library(ggplot2)
ggplot(BWHG, aes(x=Group, y=Value, fill=Cluster)) +
  facet_grid(.~Group, scales="free")+
  ggtitle("Ecoli Effect on Huddling+Grooming") +
  labs(y="Ecoli Similarity") +
  theme(legend.position="top")+
  theme(legend.text=element_text(size=15))+
  theme(plot.title = element_text(face="bold"))+
  geom_boxplot()


################################
#ks test 
################################
ks.test(as.numeric(W6HG),as.numeric(B6HG),alternative = "less")
ks.test(as.numeric(W8HG),as.numeric(B8HG),alternative = "less")
ks.test(as.numeric(W13HG),as.numeric(B13HG),alternative = "less")



##################################
#Ecoli + Environment 
EcoEnv8<- read_csv("~/Dropbox/Research/VET_SNH_Monkey/DATA/NC8_BL_Ecoli_indiv+env.csv")
EcoEnv13<- read_csv("~/Dropbox/Research/VET_SNH_Monkey/DATA/NC13_BL_Ecoli_indiv+env.csv")

EcoEnv8=EcoEnv8[,-1]
EcoEnv13=EcoEnv13[,-1]

Symmtriz_clean=function(Mat)
{
  for (i in 1:nrow(Mat))
  {
    for (j in i:ncol(Mat))
    {
      Mat[i,j]=Mat[j,i]
    }
  }
  return (Mat)
}
###############
EcoEnv8=Symmtriz_clean(EcoEnv8)
EcoEnv13=Symmtriz_clean(EcoEnv13)

EcoEnv8=as.matrix(EcoEnv8)
EcoEnv13=as.matrix(EcoEnv13)

#####################
#standardize NC8Ecoli since its distribution are from (0,100). As a similarity matrix input, we need 
#to keep all the entries within 1. 
EcoEnv8=EcoEnv8/100
EcoEnv13=EcoEnv13/100
#####################
#run DCG 
############
#NC8
temp=c(0.1,0.2,0.7)
Ens8EcoEn=Eigen.plot2(temp, selected.id=c(1,2,3),EcoEnv8)

DCG.EcoEn8=DCGtree.plot(num.clusters.selected=c(1,2,8),
                         "NC13Huddle_Groom tree",Ens8EcoEn,temp)
################
#NC13
temp=c(0.1,0.2,0.7)
Ens13EcoEn=Eigen.plot2(temp, selected.id=c(1,2,3),EcoEnv13)

DCG.EcoEn13=DCGtree.plot(num.clusters.selected=c(1,2,13),
                      "NC13Huddle_Groom tree",Ens13EcoEn,temp)
################



heatmap.2(EcoEnv8,Rowv = as.dendrogram(DCG.EcoEn8),main="NC8 Ecoli+environmental isolates",
          Colv = as.dendrogram(DCG.EcoEn8),trace="none")

heatmap.2(EcoEnv13,Rowv = as.dendrogram(DCG.EcoEn13),main="NC13 Ecoli+environmental isolates",
          Colv = as.dendrogram(DCG.EcoEn13),trace="none")

save(DCG.EcoEn8,DCG.EcoEn13,file="Ecoli_Env.RData")
