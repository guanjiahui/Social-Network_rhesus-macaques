NC8_Agg <- read.csv("~/Dropbox/Research/VET_SNH_Monkey/DATA/NC8_Aggrsym_matrix.csv")
NC6_Agg <- read.csv("~/Dropbox/Research/VET_SNH_Monkey/DATA/NC6_Aggrsym_matrix.csv")
NC13_Agg <- read.csv("~/Dropbox/Research/VET_SNH_Monkey/DATA/NC13_Aggrsym_matrix.csv")

colnames(NC8_Agg)=NC8_Agg[,1]
NC8_Agg=NC8_Agg[,-1]

colnames(NC6_Agg)=NC6_Agg[,1]
NC6_Agg=NC6_Agg[,-1]

colnames(NC13_Agg)=NC13_Agg[,1]
NC13_Agg=NC13_Agg[,-1]
################3
DiagClean=function(Mat)
{
  n=nrow(Mat)
  for(i in 1:n)
  {
     Mat[i,i]=0
  }
  return(Mat)
}

###############################################
NC6_Agg=DiagClean(NC6_Agg)
NC8_Agg=DiagClean(NC8_Agg)
NC13_Agg=DiagClean(NC13_Agg)
################333
Agg6=DistPairwise(NC6_Agg)
Agg8=DistPairwise(NC8_Agg)
Agg13=DistPairwise(NC13_Agg)

####################3
#diagonal check 
DiagCheck=function(Mat)
{
  n=nrow(Mat)
  Di=0
  for(i in 1:n)
  {
    if(Mat[i,i]!=0)
      Di=1
  }
  return(Di)
}

################


################################
#NC6
################################
temp=c(2,5,10,50)
Ens.agg6=Eigen.plot(temp, selected.id=c(1,2,3,4),Agg6)

DCG.agg6=DCGtree.plot(num.clusters.selected=c(1,2,3,5),
                      "Aggression NC6",Ens.agg6,temp)

plot(DCG.agg6,hang=-1,main="Aggression NC6")

###################################
#NC8
##################################

temp=c(2,5,10,50)
Ens.agg8=Eigen.plot(temp, selected.id=c(1,2,3,4),Agg8)

DCG.agg8=DCGtree.plot(num.clusters.selected=c(1,2,3,5),
                      "Aggression NC8",Ens.agg8,temp)

plot(DCG.agg8,hang=-1,main="Aggression NC8")


############################3
#NC13
#########################3
temp=c(2,3,8,20,50)
Ens.agg13=Eigen.plot(temp, selected.id=c(1,2,3,4,5),Agg13)

DCG.agg13=DCGtree.plot(num.clusters.selected=c(1,2,3,4),
                      "Aggression NC13",Ens.agg13,temp)

plot(DCG.agg13,hang=-1,main="Aggression NC13")

##################33

#cut tree
#########
yAgg6=cutree(DCG.agg6,k=4)
yAgg8=cutree(DCG.agg8,k=4)
yAgg13=cutree(DCG.agg13,k=5)
##################

#########
#compute the within and between
NC6Comp1=TriEntropy( yNC6E1,yNC6H1,yNC6G1,nameE6)
NC8Comp1=TriEntropy(yNC8E1,yNC8H1,yNC8G1,nameE8)
NC13Comp1=TriEntropy( yNC13E1,yNC13H1,yNC13G1,nameE)

Entropy_general(yAgg6,yNC6E1,nameE6[[4]],nameE6[[1]])$Sum
Entropy_general(yAgg8,yNC8E1,nameE8[[4]],nameE8[[1]])$Sum
Entropy_general(yAgg13,yNC13E1,nameE[[4]],nameE[[1]])$Sum


################################


##########################
#Withn and Between Ecoli Effect 
##########################
AGG6Effect=EcoliEffect(yAgg6,nameE6[[4]],nameE6[[1]],NC6Ecol1*100)$E
AGG8Effect=EcoliEffect(yAgg8,nameE8[[4]],nameE8[[1]],NC8Ecol1)$E
AGG13Effect=EcoliEffect(yAgg13,nameE[[4]],nameE[[1]],NC13B)$E

##################
#visualize the irregular matrix heatmap 
################
scalePlot(yAgg6, AGG6Effect,"NC6 Baseline Aggression")
scalePlot(yAgg8, AGG8Effect,"NC8 Baseline Aggression")
scalePlot(yAgg13, AGG13Effect,"NC13 Baseline Aggression")

#######################3


################################
#Boxplot for within-between effect of Ecoli 
###################################
#Within 
W6AGG=sapply(1:nrow(AGG6Effect), function(i) AGG6Effect[i,i])
W8AGG=sapply(1:nrow(AGG8Effect), function(i) AGG8Effect[i,i])
W13AGG=sapply(1:nrow(AGG13Effect), function(i) AGG13Effect[i,i])

B6AGG=BetExtract(AGG6Effect)
B8AGG=BetExtract(AGG8Effect)
B13AGG=BetExtract(AGG13Effect)

############################
#Boxplot
#######################3
n1=length(W6AGG)
n2=length(W8AGG)
n3=length(W13AGG)

n4=length(B6AGG)
n5=length(B8AGG)
n6=length(B13AGG)

BWAGG=data.frame(Group=c(rep("Group I",n1+n4),rep("Group II",n2+n5),rep("Group III",n3+n6)), 
                    Cluster=c(rep("within",n1),rep("between",n4), rep("within",n2), rep("between",n5),rep("within",n3), rep("between",n6)),
                    Value=c(W6AGG,B6AGG,W8AGG,B8AGG,W13AGG,B13AGG))

##################
ggplot(BWAGG, aes(x=Group, y=Value, fill=Cluster)) +
  facet_grid(.~Group, scales="free")+
  ggtitle("Ecoli Effect on Aggression") +
  labs(y="Ecoli Similarity") +
  theme(legend.position="top")+
  theme(legend.text=element_text(size=15))+
  theme(plot.title = element_text(face="bold"))+
  geom_boxplot()
###################################

BWAGG=data.frame(Group=c(rep("Group I",n1+n4),rep("Group II",n2+n5),rep("Group III",n3+n6)), 
                   Cluster=c(rep("within-clusters",n1),rep("between-clusters",n4), 
                             rep("within-clusters",n2), rep("between-clusters",n5),
                             rep("within-clusters",n3), rep("between-clusters",n6)),
                 Value=c(W6AGG,B6AGG,W8AGG,B8AGG,W13AGG,B13AGG))


ggplot(BWAGG, aes(x=Group, y=Value, fill=Cluster)) +
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
ks.test(as.numeric(W6AGG),as.numeric(B6AGG))
ks.test(as.numeric(W8AGG),as.numeric(B8AGG))
ks.test(as.numeric(W13AGG),as.numeric(B13AGG),alternative = "less")
########

#################3
#update the entropy plot 
x1=c(0.6707,1.4747,1.0467,1.0004,0.9344)
x2=c(0.79686, 1.33614,1.1873,1.1367,1.0799)
x3=c(0.84265,0.9164 ,0.76205,0.617,0.499)
x4=c(0.93084,1.062176,0.9398,0,0)
S=cbind(x1,x2,x3,x4)

#only baseline comparison 
Sb=S[c(1,2,3),c(1,2,4)]

bb=c("Ecol vs Huddle","Ecoli vs Groom", "Ecoli vs Aggression")
plot(Sb[1,],type="b",ylab="Entropy",xaxt="n",main="Entropy Comparison",ylim=c(min(Sb),max(Sb)),lwd=3)
axis(1, at=1:3, labels=bb)
for(i in 2:3)
{
  lines(Sb[i,],col=i,lwd=3)
  points(Sb[i,],col=i,lwd=3,pch=13+i)
}


legend("topright",c("NC6","NC8", "NC13"), lty=rep(1,3),
       col=c(1,2,3),cex=1.3)

