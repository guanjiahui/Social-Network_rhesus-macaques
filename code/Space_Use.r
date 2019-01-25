NC6_Spaceusefreq <- read_csv("~/Dropbox/Research/VET_SNH_Monkey/DATA/NC6_Spaceusefreq.csv")


NC6_Space=NC6_Spaceusefreq[,-1]
NC6_Space=as.matrix(NC6_Space)

for (i in 1:nrow(NC6_Space))
{
  if(NC6_Space[i,i]!=0)
    print("bad")
  print("good")
}

NC8_Spaceusefreq <- read_csv("~/Dropbox/Research/VET_SNH_Monkey/DATA/NC8_Spaceusefreq.csv")
NC13_Spaceusefreq <- read_csv("~/Dropbox/Research/VET_SNH_Monkey/DATA/NC13_Spaceusefreq.csv")

NC8_Space=NC8_Spaceusefreq[,-1]
NC8_Space=as.matrix(NC8_Space)
NC13_Space=NC13_Spaceusefreq[,-1]
NC13_Space=as.matrix(NC13_Space)

for (i in 1:nrow(NC8_Space))
{
  if(NC8_Space[i,i]!=0)
    print("bad")
}

for (i in 1:nrow(NC13_Space))
{
  if(NC13_Space[i,i]!=0)
    print("bad")
}

library(readr)
NC6_Spaceusefreq <- read_csv("~/Dropbox/Research/VET_SNH_Monkey/DATA/NC6_SpaceuseInd.csv")
NC8_Spaceusefreq <- read_csv("~/Dropbox/Research/VET_SNH_Monkey/DATA/NC8_SpaceuseInd.csv")
NC13_Spaceusefreq <- read_csv("~/Dropbox/Research/VET_SNH_Monkey/DATA/NC13_SpaceuseInd.csv")

NC6_Space=NC6_Spaceusefreq[,-1]
NC6_Space=as.matrix(NC6_Space)
NC8_Space=NC8_Spaceusefreq[,-1]
NC8_Space=as.matrix(NC8_Space)
NC13_Space=NC13_Spaceusefreq[,-1]
NC13_Space=as.matrix(NC13_Space)

######################################################3


################################
#NC6
################################
temp=c(2,10,50,100)
Ens.space6=Eigen.plot(temp, selected.id=c(1,2,3,4),NC6_Space)

DCG.space6=DCGtree.plot(num.clusters.selected=c(1,2,5,7),
                      "Aggression NC6",Ens.space6,temp)

plot(DCG.space6,hang=-1,main="Aggression NC6")

###################################
#NC8
##################################

temp=c(2,5,10,40,80)
Ens.space8=Eigen.plot(temp, selected.id=c(1,2,3,4,5),NC8_Space)

DCG.space8=DCGtree.plot(num.clusters.selected=c(1,2,3,4,5),
                      "Aggression NC8",Ens.space8,temp)

plot(DCG.space8,hang=-1,main="Aggression NC8")

############################3
#NC13
#########################3
temp=c(10,50,80,130,230)
Ens.space13=Eigen.plot(temp, selected.id=c(1,2,3,4,5),NC13_Space)

DCG.space13=DCGtree.plot(num.clusters.selected=c(1,2,3,5,7),
                       "Aggression NC13",Ens.space13,temp)

plot(DCG.space13,hang=-1,main="Aggression NC13")

##################33

#cut tree
#########
ySpace6=cutree(DCG.space6,k=6)
ySpace8=cutree(DCG.space8,k=6)
ySpace13=cutree(DCG.space13,k=5)
##################







################################
#NC6
################################
temp=c(2,10,50,100)
Ens.space62=Eigen.plot(temp, selected.id=c(1,2,3,4),NC6_Space)

DCG.space62=DCGtree.plot(num.clusters.selected=c(1,2,4,8),
                        "Aggression NC6",Ens.space62,temp)

plot(DCG.space62,hang=-1,main="Aggression NC6")

###################################
#NC8
##################################

temp=c(2,5,10,40,80)
Ens.space82=Eigen.plot(temp, selected.id=c(1,2,3,4,5),NC8_Space)

DCG.space82=DCGtree.plot(num.clusters.selected=c(1,2,3,5,6),
                        "Aggression NC8",Ens.space82,temp)

plot(DCG.space82,hang=-1,main="Aggression NC8")

############################3
#NC13
#########################3
temp=c(10,50,80,130,230)
Ens.space132=Eigen.plot(temp, selected.id=c(1,2,3,4,5),NC13_Space)

DCG.space132=DCGtree.plot(num.clusters.selected=c(1,2,3,5,7),
                         "Aggression NC13",Ens.space132,temp)

plot(DCG.space132,hang=-1,main="Aggression NC13")

#cut tree
#########
ySpace6=cutree(DCG.space62,k=3)
ySpace8=cutree(DCG.space82,k=5)
ySpace13=cutree(DCG.space132,k=4)
##################








################################
#assigning names 
nameE6[[6]]=colnames(NC6_Space)
nameE8[[6]]=colnames(NC8_Space)
nameE13[[6]]=colnames(NC13_Space)

##########################
#Withn and Between Ecoli Effect 
##########################
Space6Effect=EcoliEffect(ySpace6,nameE6[[6]],nameE6[[1]],NC6Ecol1*100)$E
Space8Effect=EcoliEffect(ySpace8,nameE8[[6]],nameE8[[1]],NC8Ecol1)$E
Space13Effect=EcoliEffect(ySpace13,nameE13[[6]],nameE13[[1]],NC13B)$E

##################
#visualize the irregular matrix heatmap 
################
scalePlot(ySpace6, Space6Effect,"NC6 Baseline Spaceression")
scalePlot(ySpace8, Space8Effect,"NC8 Baseline Spaceression")
scalePlot(ySpace13, Space13Effect,"NC13 Baseline Spaceression")


################################
#Boxplot for within-between effect of Ecoli 
###################################
#Within 
W6Space=sapply(1:nrow(Space6Effect), function(i) Space6Effect[i,i])
W8Space=sapply(1:nrow(Space8Effect), function(i) Space8Effect[i,i])
W13Space=sapply(1:nrow(Space13Effect), function(i) Space13Effect[i,i])

B6Space=BetExtract(Space6Effect)
B8Space=BetExtract(Space8Effect)
B13Space=BetExtract(Space13Effect)

############################
#Boxplot
#######################3
library(ggplot2)
n1=length(W6Space)
n2=length(W8Space)
n3=length(W13Space)

n4=length(B6Space)
n5=length(B8Space)
n6=length(B13Space)

BWSpace=data.frame(Group=c(rep("Group I",n1+n4),rep("Group II",n2+n5),rep("Group III",n3+n6)), 
                 Cluster=c(rep("within-clusters",n1),rep("between-clusters",n4), 
                           rep("within-clusters",n2), rep("between-clusters",n5),
                           rep("within-clusters",n3), rep("between-clusters",n6)),
                 Value=c(W6Space,B6Space,W8Space,B8Space,W13Space,B13Space))


ggplot(BWSpace, aes(x=Group, y=Value, fill=Cluster)) +
  facet_grid(.~Group, scales="free")+
  #ggtitle("Ecoli Effect on Space-Use") +
  labs(y="E. coli % Similarity") +
  theme(legend.position="top")+
  theme(legend.text=element_text(size=18))+
 # theme(plot.title = element_text(face="bold"))+
  theme(legend.title=element_blank())+
  theme(axis.text=element_text(size=14))+
  theme(axis.title.y=element_text(size=14))+
  theme(axis.title.x=element_text(size=14))+
  geom_boxplot()

############################



####################3
#K-S Test 
################
ks.test(as.numeric(W6Space),as.numeric(B6Space))
ks.test(as.numeric(W8Space),as.numeric(B8Space))
ks.test(as.numeric(W13Space),as.numeric(B13Space))
########





#########################################################################################################
#######################################################################################################
#SPace_Use Ind
library(readr)
NC6_SpaceuseInd <- read_csv("~/Dropbox/Research/VET_SNH_Monkey/DATA/Major Revision/NC6_SpaceuseInd.csv")
NC8_SpaceuseInd <- read_csv("~/Dropbox/Research/VET_SNH_Monkey/DATA/Major Revision/NC8_SpaceuseInd.csv")
NC13_SpaceuseInd <- read_csv("~/Dropbox/Research/VET_SNH_Monkey/DATA/Major Revision/NC13_SpaceuseInd.csv")


NC6_Space1=NC6_SpaceuseInd[,-1]
NC6_Space1=as.matrix(NC6_Space1)
NC8_Space1=NC8_SpaceuseInd[,-1]
NC8_Space1=as.matrix(NC8_Space1)
NC13_Space1=NC13_SpaceuseInd[,-1]
NC13_Space1=as.matrix(NC13_Space1)


for (i in 1:nrow(NC8_Space1))
{
  if(NC8_Space[i,i]!=0)
    print("bad")
}

for (i in 1:nrow(NC13_Space1))
{
  if(NC13_Space[i,i]!=0)
    print("bad")
}

colnames(NC6_Space1)=gsub("x","",colnames(NC6_Space1))
colnames(NC8_Space1)=gsub("x","",colnames(NC8_Space1))
colnames(NC13_Space1)=gsub("x","",colnames(NC13_Space1))
###################




################################
#NC6
################################
temp=c(2,10,50,100)
Ens.Sp6=Eigen.plot(temp, selected.id=c(1,2,3,4),NC6_Space1)

DCG.Sp6=DCGtree.plot(num.clusters.selected=c(1,2,4,8),
                         "Aggression NC6",Ens.Sp6,temp)

plot(DCG.Sp6,hang=-1,main="Aggression NC6")

###################################
#NC8
##################################

temp=c(2,5,10,40,80)
Ens.Sp8=Eigen.plot(temp, selected.id=c(1,2,3,4,5),NC8_Space1)

DCG.Sp8=DCGtree.plot(num.clusters.selected=c(1,2,3,5,8),
                         "Aggression NC8",Ens.Sp8,temp)

plot(DCG.Sp8,hang=-1,main="Aggression NC8")

############################3
#NC13
#########################3
temp=c(10,50,80,130,230)
Ens.Sp13=Eigen.plot(temp, selected.id=c(1,2,3,4,5),NC13_Space1)

DCG.Sp13=DCGtree.plot(num.clusters.selected=c(1,2,3,5,7),
                          "Aggression NC13",Ens.Sp13,temp)

plot(DCG.Sp13,hang=-1,main="Aggression NC13")

#cut tree
#########
ySp6=cutree(DCG.Sp6,k=5)
ySp8=cutree(DCG.Sp8,k=8)
ySp13=cutree(DCG.Sp13,k=4)
##################



################################
#assigning names 
nameE6[[7]]=colnames(NC6_Space1)
nameE8[[7]]=colnames(NC8_Space1)
nameE13[[7]]=colnames(NC13_Space1)

##########################
#Withn and Between Ecoli Effect 
##########################
Space6Effect1=EcoliEffect(ySp6,nameE6[[7]],nameE6[[1]],NC6Ecol1*100)$E
Space8Effect1=EcoliEffect(ySp8,nameE8[[7]],nameE8[[1]],NC8Ecol1)$E
Space13Effect1=EcoliEffect(ySp13,nameE13[[7]],nameE13[[1]],NC13B)$E

##################
#visualize the irregular matrix heatmap 
################
scalePlot(ySpace6, Space6Effect1,"NC6 Baseline Spaceression")
scalePlot(ySpace8, Space8Effect1,"NC8 Baseline Spaceression")
scalePlot(ySpace13, Space13Effect1,"NC13 Baseline Spaceression")


################################
#Boxplot for within-between effect of Ecoli 
###################################
#Within 
W6Sp=sapply(1:nrow(Space6Effect1), function(i) Space6Effect1[i,i])
W8Sp=sapply(1:nrow(Space8Effect1), function(i) Space8Effect1[i,i])
W13Sp=sapply(1:nrow(Space13Effect1), function(i) Space13Effect1[i,i])

B6Sp=BetExtract(Space6Effect1)
B8Sp=BetExtract(Space8Effect1)
B13Sp=BetExtract(Space13Effect1)

############################
#Boxplot
#######################3
library(ggplot2)
n1=length(W6Sp)
n2=length(W8Sp)
n3=length(W13Sp)

n4=length(B6Sp)
n5=length(B8Sp)
n6=length(B13Sp)

BWSp=data.frame(Group=c(rep("Group I",n1+n4),rep("Group II",n2+n5),rep("Group III",n3+n6)), 
                   Cluster=c(rep("within-clusters",n1),rep("between-clusters",n4), 
                             rep("within-clusters",n2), rep("between-clusters",n5),
                             rep("within-clusters",n3), rep("between-clusters",n6)),
                   Value=c(W6Sp,B6Sp,W8Sp,B8Sp,W13Sp,B13Sp))


ggplot(BWSp, aes(x=Group, y=Value, fill=Cluster)) +
  facet_grid(.~Group, scales="free")+
  #ggtitle("Ecoli Effect on Space-Use") +
  labs(y="E. coli % Similarity") +
  theme(legend.position="top")+
  theme(legend.text=element_text(size=18))+
  # theme(plot.title = element_text(face="bold"))+
  theme(legend.title=element_blank())+
  theme(axis.text=element_text(size=14))+
  theme(axis.title.y=element_text(size=14))+
  theme(axis.title.x=element_text(size=14))+
  geom_boxplot()

############################



####################3
#K-S Test 
################
wilcox.test(as.numeric(W6Sp),as.numeric(B6Sp))
wilcox.test(as.numeric(W8Sp),as.numeric(B8Sp))
wilcox.test(as.numeric(W13Space),as.numeric(B13Sp))
########

#z-Score 
z_wilcox_rank(as.numeric(W6Sp),as.numeric(B6Sp))
z_wilcox_rank(as.numeric(W8Sp),as.numeric(B8Sp))
z_wilcox_rank(as.numeric(W13Space),as.numeric(B13Sp))