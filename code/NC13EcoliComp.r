#> library(sparcl)
#> # colors the leaves of a dendrogram
#  > y1 = cutree(DCG.NCB, 17)
 # > ColorDendrogram(DCG.NCB, y = y1, main = "NC13BL-EColi tree", 
#                    +                 branchlength = 80)
#  > plot(DCG.NCPE,hang=-1)
#  > y2= cutree(DCG.NCPE, 13)
#  > ColorDendrogram(DCG.NCPE, y = y2, main = "NC13PERT_EColi tree", 
#                    +                 branchlength = 80)
 # > ############################
#    > y3 = cutree(DCG.NCPO, 9)
#    > ColorDendrogram(DCG.NCPO, y = y3,  main = "NC13_PORT_EColi tree", 
 #                     +                 branchlength = 80)
 
 
 ########

 #Ecoli Comparison 
 
 EntryE=list()
 EntryE[[1]]=y1
 EntryE[[2]]=y2
 EntryE[[3]]=y3
 
 NameE=list()
 NameE[[1]]=colnames(NC13B)
 NameE[[2]]=colnames(NC13PE)
 NameE[[3]]=colnames(NC13PO)
 
 par(mfrow=c(2,3))
 EnE=matrix(0,3,3)
 for (i in 1:3){
   for (j in 1:3){
     if (i!=j){
       EnE[i,j]=Entropy_sequence(EntryE[[i]],EntryE[[j]],NameE[[i]],NameE[[j]])$Sum
       plot(Entropy_sequence(EntryE[[i]],EntryE[[j]],NameE[[i]],NameE[[j]])$Entropy, type="b",
            xlab="group", ylab="Entropy",main=paste("stage",i,"from","stage",j))
     }
   }
 }
 
 for (i in 1:3)
   for (j in 1:3)
     cat(EnE[i,j]+EnE[j,i],"i",i,"j",j,"\n")
 
 NC13EcoliComp=TriEntropy( yNC13E1,yNC13E2,yNC13E3,nameE)
 #################
 ######
 #baseline Comparison 

 yNC13E1= cutree(DCG.NCB, 17)

 
 yNC13H1=cutree(DCG.heat1,15)
 yNC13G1=cutree(DCG13.groom1,12)
 
 nameE=list()
 nameE[[1]]=colnames(NC13B)
 nameE[[2]]=colnames(NC13HuddleR1)
 nameE[[3]]=colnames(NC13HuddleR1)
 nameE[[4]]=colnames(NC13_Agg)
 #######################
 
 NC13Comp1=TriEntropy( yNC13E1,yNC13H1,yNC13G1,nameE)
 
 
 #####
 
 
 aa=intersect(nameE[[1]],nameE[[2]])
 common=intersect(aa,nameE[[3]])
 
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
   loc[[i]]=Find(nameE[[i]],common)
 }
 
 y13Ecoli1=yNC13E1[loc[[1]]]
 y13Huddle1=yNC13H1[loc[[2]]]
 y13Groom1=yNC13G1[loc[[3]]]
 ID13=as.numeric(nameE[[1]][loc[[1]]])
 y13Baseline=data.frame(cbind(ID13,y13Ecoli1,y13Huddle1,y13Groom1))
 
 write.csv(y13Baseline,file="NC13BaselineMembership.csv")
 
 
 a13=order(y13Baseline$y13Ecoli1)
 y13BaseSorted=y13Baseline[a13,]
 
 
 write.csv(y13BaseSorted,file="NC13Base_Sorted.csv")
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ######################################
 #now do the second time stage 
 ######################################
 
 yNC13E2= cutree(DCG.NCPE, 13)
 yNC13H2=cutree(DCG.heat2,18)
 yNC13G2=cutree(DCG13.groom2,12)
 
 nameE2=list()
 nameE2[[1]]=colnames(NC13PE)
 nameE2[[2]]=colnames(NC13HuddleR2)
 nameE2[[3]]=colnames(NC13HuddleR2)
 #######################
 
 NC13Comp2=TriEntropy( yNC13E2,yNC13H2,yNC13G2,nameE2)
 
 #######
 aa=intersect(nameE2[[1]],nameE2[[2]])
 common=intersect(aa,nameE2[[3]])
 
 #####################################
 loc=list()
 for(i in 1:3)
 {
   loc[[i]]=Find(nameE2[[i]],common)
 }
 
 y13Ecoli2=yNC13E2[loc[[1]]]
 y13Huddle2=yNC13H2[loc[[2]]]
 y13Groom2=yNC13G2[loc[[3]]]
 
 
 
 ######################################
 #now do the third time stage 
 ######################################
 yNC13E3= cutree(DCG.NCPO, 12)
 yNC13H3=cutree(DCG.heat3,24)
 yNC13G3=cutree(DCG13.groom3,17)
 
 nameE3=list()
 nameE3[[1]]=colnames(NC13PO)
 nameE3[[2]]=colnames(NC13HuddleR3)
 nameE3[[3]]=colnames(NC13HuddleR3)
 #######################
 
 NC13Comp3=TriEntropy( yNC13E3,yNC13H3,yNC13G3,nameE3)
 
 
 
 
 #######################
 #within and between 
 
 #1 Huddling 
 
 ####################################################################
EcoliEffect=function(y1,name1,ID2,Mat)
 {
   N=length(unique(y1))
   EcoliWB=matrix(0,N,N)
   With=numeric(N)
   ###########################
   #first compute the diganol
   ###########################
   for(i in 1:N)
   {
     loc=which(y1==i)
     ID=name1[loc]
     
     ID=intersect(ID,ID2)
     n=length(ID)
     within=0
    
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
   
   ##########################
   #compute the upper off-diaganol
   ##########################
   for(i in 1:N)
   {
     for(j in i:N)
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
   
   
   ######################
   #copy the lower triangle
   ########################
   for(j in 1:N)
   {
     for(i in j:N)
     {
       EcoliWB[i,j]=EcoliWB[j,i]
     }
   }
   
   
   return(list(W=With, E=EcoliWB))
 }
 
 ####################################################################
 #EcoliEffect=function(y1,name1,ID2,Mat)
 ###############################
 Huddle13Effect=EcoliEffect(yNC13H1,nameE[[2]],nameE[[1]],NC13B)$E
 Groom13Effect=EcoliEffect(yNC13G1,nameE[[3]],nameE[[1]],NC13B)$E
 
 ###############################
 Huddle13Effect2=EcoliEffect(yNC13H2,nameE2[[2]],nameE2[[1]],NC13PE)
 Groom13Effect2=EcoliEffect(yNC13G2,nameE2[[3]],nameE2[[1]],NC13PE)
 
 Huddle13Effect3=EcoliEffect(yNC13H3,nameE3[[2]],nameE3[[1]],NC13PO)
 Groom13Effect3=EcoliEffect(yNC13G3,nameE3[[3]],nameE3[[1]],NC13PO)
 
 #################
 Huddle8Effect=EcoliEffect(yNC8H1,nameE8[[2]],nameE8[[1]],NC8Ecol1)$E
 Groom8Effect=EcoliEffect(yNC8G1,nameE8[[3]],nameE8[[1]],NC8Ecol1)$E
 
 
 ################
 #visualzie it 
 #################
 
colWidth=function(seq)
{
  N=length(unique(seq))
  W=c()
  
  for(i in 1:N)
  {
    W[i]=length(which(seq==i))
  }
  
  return(W)
}
 
colWidth(yNC13H1)
 
 
 
 
 ######################
 mybreaks <- seq(38, 100, length.out=10000)
 colfunc2 <- colorRampPalette(c("white", "black"))
 heatmap.2(Huddle13Effect,Rowv = Tree13H1U, Colv =Tree13H1U,breaks = mybreaks,col= colfunc2 ,
           trace="none",main="NC13 Baseline Huddling")
 
 heatmap.2(Groom13Effect,Rowv = Tree13G1U, Colv = Tree13G1U,trace="none",breaks = mybreaks,col= colfunc2 ,
           main="NC13 Baseline Grooming")
 ########
 
 heatmap.2(Huddle8Effect,Rowv = Tree8H1U, Colv = Tree8H1U,trace="none",breaks = mybreaks,col= colfunc2 ,
           main="NC8 Baseline Huddling")
 
 heatmap.2(Groom8Effect,Rowv =Tree8G1U, Colv =Tree8G1U,trace="none",breaks = mybreaks,col= colfunc2 ,
           main="NC8 Baseline Grooming")
 
 
 
 
 
 
 
 
 
 ###########################################
 library(ggplot2)
 library(reshape2) 
 
 X <- matrix(nrow=3, ncol=3)
 X[1,] <- c(0.3, 0.4, 0.45)
 X[2,] <- c(0.3, 0.7, 0.65)
 X[3,] <- c(0.3, 0.4, 0.45)
 
 
 colnames(X)<-c(1.5, 6, 6.1)
 rownames(X)<-c(1.5, 6, 6.1)
 X <- melt(X)
 X <- as.data.frame(X)
 names(X) <- c("Var1", "Var2", "value")
 v1m <- unique(X$Var1)
 X$Var1.min <- rep(c(0, v1m[-length(v1m)]), length.out = length(v1m))
 v2m <- unique(X$Var2)
 X$Var2.min <- rep(c(0, v2m[-length(v2m)]), each = length(v2m))
 
 ggplot(data = X, aes(fill = value)) + 
   geom_rect(aes(ymin = Var1.min, ymax = Var1, xmin = Var2.min, xmax = Var2))
 ############################
 
 
 library(ggdendro)
 library(grid)
 
 scalePlot=function(y,Mat,Title)
 {
  limits=c(38,100)
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
 #####################################
scalePlot(yNC13H1, Huddle13Effect,"NC13 Baseline Huddling")
scalePlot(yNC13G1,Groom13Effect, "NC13 Baseline Grooming") 

scalePlot(yNC8H1, Huddle8Effect,"NC8 Baseline Huddling")
scalePlot(yNC8G1,Groom8Effect, "NC8 Baseline Grooming") 

scalePlot(yNC6H1, Huddle6Effect,"NC6 Baseline Huddling")
scalePlot(yNC6G1,Groom6Effect, "NC6 Baseline Grooming") 


 
###############
 Tree13H1=as.dendrogram(DCG.heat1)
 Tree13H1U=cut(Tree13H1,h=2)$upper
 plot(Tree13H1U)
 
 Tree13G1=as.dendrogram(DCG13.groom1)
 Tree13G1U=cut(Tree13G1,h=2)$upper
 plot(Tree13G1U)
 
 ########
 Tree8H1=as.dendrogram(DCG8.huddle1)
 Tree8H1U=cut(Tree8H1,h=0.2)$upper
 plot(Tree8H1U)
 
 Tree8G1=as.dendrogram(DCG.groom81)
 Tree8G1U=cut(Tree8G1,h=3)$upper
 plot(Tree8G1U)
 
 
 
 
 
 
 
 ##########
 
par(mfrow=c(3,2))
hist((as.numeric(Huddle6Effect)),main="Huddle 6",prob=TRUE,xlim=c(30,100))
lines(density(as.numeric(Huddle6Effect)))
 
hist((as.numeric(Groom6Effect)),breaks=10,main="Groom 6",prob=TRUE,xlim=c(30,100))
lines(density(as.numeric(Groom6Effect)))

hist((as.numeric(Huddle8Effect)),breaks=10,main="Huddle 8",prob=TRUE,xlim=c(30,100))
lines(density(as.numeric(Huddle8Effect)))

hist((as.numeric(Groom8Effect)),main="Groom 8",prob=TRUE,xlim=c(30,100))
lines(density(as.numeric(Groom8Effect)))

hist((as.numeric(Huddle13Effect)),main="Huddle 13",prob=TRUE,xlim=c(30,100))
lines(density(as.numeric(Huddle13Effect)))

hist((as.numeric(Groom13Effect)),main="Groom 13",prob=TRUE,xlim=c(30,100))
lines(density(as.numeric(Groom13Effect)))
par(mfrow=c(1,1))
 
##########
#run kolmogorov-Smirnov Test 
ks.test(as.numeric(Huddle6Effect),as.numeric(Groom6Effect))
ks.test(as.numeric(Huddle8Effect),as.numeric(Groom8Effect))
ks.test(as.numeric(Huddle13Effect),as.numeric(Groom13Effect))
 
ks.test(as.numeric(Huddle6Effect),as.numeric(Huddle8Effect))
ks.test(as.numeric(Huddle6Effect),as.numeric(Huddle13Effect))
ks.test(as.numeric(Huddle8Effect),as.numeric(Huddle13Effect))


ks.test(as.numeric(Groom6Effect),as.numeric(Groom8Effect))
ks.test(as.numeric(Groom6Effect),as.numeric(Groom13Effect))
ks.test(as.numeric(Groom8Effect),as.numeric(Groom13Effect))

 
 
 
 
 
 
 
 #############
 library(ROCR)
 

x1=c(0.6707,1.4747,1.0467,1.0004,0.9344)
x2=c(0.79686, 1.33614,1.1873,1.1367,1.0799)
x3=c(0.84265,0.9164 ,0.76205,0.617,0.499)
S=cbind(x1,x2,x3)

bb=c("Ecol vs Huddle","Ecoli vs Groom", "Huddle vs Groom")
plot(S[1,],type="b",ylab="Entropy",xaxt="n",main="Entropy Comparison",ylim=c(min(S),max(S)),lwd=3)
axis(1, at=1:3, labels=bb)
for(i in 2:5)
{
  lines(S[i,],col=i,lwd=3)
  points(S[i,],col=i,lwd=3,pch=13+i)
}


legend("topright",c("NC6","NC8", "NC13 Base","NC13 Pert","NC13 Post"), lty=rep(1,5),
       col=c(1,2,3,4,5),cex=1.3)

####################################
#only baseline comparison 
bb=c("Ecol vs Huddle","Ecoli vs Groom", "Huddle vs Groom")
plot(S[1,],type="b",ylab="Entropy",xaxt="n",main="Entropy Comparison",ylim=c(min(S),max(S)),lwd=3)
axis(1, at=1:3, labels=bb)
for(i in 2:3)
{
  lines(S[i,],col=i,lwd=3)
  points(S[i,],col=i,lwd=3,pch=13+i)
}


legend("topright",c("NC6","NC8", "NC13"), lty=rep(1,3),
       col=c(1,2,3),cex=1.3)

























################################
#Boxplot for within-between effect of Ecoli 
###################################
#Within 
W6huddle=sapply(1:nrow(Huddle6Effect), function(i) Huddle6Effect[i,i])
W8huddle=sapply(1:nrow(Huddle8Effect), function(i) Huddle8Effect[i,i])
W13huddle=sapply(1:nrow(Huddle13Effect), function(i) Huddle13Effect[i,i])

W6groom=sapply(1:nrow(Groom6Effect), function(i) Groom6Effect[i,i])
W8groom=sapply(1:nrow(Groom8Effect), function(i) Groom8Effect[i,i])
W13groom=sapply(1:nrow(Groom13Effect), function(i) Groom13Effect[i,i])

#Between 
###################
BetExtract=function(Mat)
{
  n=nrow(Mat)
  B=0
  
  for(i in 1:n)
  {
    for (j in i:n)
    {
      if(i !=j)
        B=c(B,Mat[i,j])
    }
  }
  B=B[-1]
  return(B)
}
#####################






















B6huddle=BetExtract(Huddle6Effect)
B8huddle=BetExtract(Huddle8Effect)
B13huddle=BetExtract(Huddle13Effect)

B6groom=BetExtract(Groom6Effect)
B8groom=BetExtract(Groom8Effect)
B13groom=BetExtract(Groom13Effect)
##################
#visualize it 
par(mfrow=c(3,2))
boxplot(W6huddle,B6huddle,names = c("within","between"),main="NC6 Huddling of Ecoli Effect ",ylim=c(35,100))
boxplot(W6groom,B6groom,names = c("within","between"),main="NC6 Grooming of Ecoli Effect ",ylim=c(35,100))

boxplot(W8huddle,B8huddle,names = c("within","between"),main="NC8 Huddling of Ecoli Effect ",ylim=c(35,100))
boxplot(W8groom,B8groom,names = c("within","between"),main="NC8 Grooming of Ecoli Effect ",ylim=c(35,100))

boxplot(W13huddle,B13huddle,names = c("within","between"),main="NC13 Huddling of Ecoli Effect ",ylim=c(35,100))
boxplot(W13groom,B13groom,names = c("within","between"),main="NC13 Grooming of Ecoli Effect ",ylim=c(35,100))
par(mfrow=c(1,1))
###############

boxplot(W6huddle,B6huddle,W8huddle, B8huddle,W13huddle,B13huddle,main="Ecoli Effect on Huddling",ylim=c(35,100),
        col =c(6,4,6,4,6,4) )

legend("topright", inset=.02,
       c("Within","Between"), fill=c(6,4))


bb=c("NC6","NC8", "NC13")
axis(1, at=1:3, labels=bb)



n1=length(W6huddle)
n2=length(W8huddle)
n3=length(W13huddle)

n4=length(B6huddle)
n5=length(B8huddle)
n6=length(B13huddle)

BWhuddle=data.frame(Group=c(rep("Group I",n1+n4),rep("Group II",n2+n5),rep("Group III",n3+n6)), 
                    Cluster=c(rep("within",n1),rep("between",n4), rep("within",n2), rep("between",n5),rep("within",n3), rep("between",n6)),
                    Value=c(W6huddle,B6huddle,W8huddle,B8huddle,W13huddle,B13huddle))
  

ggplot(BWhuddle, aes(x=Group, y=Value, fill=Cluster)) +
  facet_grid(.~Group, scales="free")+
  ggtitle("Ecoli Effect on Huddling") +
  labs(y="Ecoli Similarity") +
  theme(legend.position="top")+
  theme(legend.text=element_text(size=15))+
  theme(plot.title = element_text(face="bold"))+
  geom_boxplot()

######################
BWhuddle=data.frame(Group=c(rep("Group I",n1+n4),rep("Group II",n2+n5),rep("Group III",n3+n6)), 
                   Cluster=c(rep("within-clusters",n1),rep("between-clusters",n4), 
                             rep("within-clusters",n2), rep("between-clusters",n5),
                             rep("within-clusters",n3), rep("between-clusters",n6)),
                   Value=c(W6huddle,B6huddle,W8huddle,B8huddle,W13huddle,B13huddle))


ggplot(BWhuddle, aes(x=Group, y=Value, fill=Cluster)) +
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






####################33
#now for the grooming 
n1=length(W6groom)
n2=length(W8groom)
n3=length(W13groom)

n4=length(B6groom)
n5=length(B8groom)
n6=length(B13groom)

BWgroom=data.frame(Group=c(rep("Group I",n1+n4),rep("Group II",n2+n5),rep("Group III",n3+n6)), 
                    Cluster=c(rep("within",n1),rep("between",n4), rep("within",n2), rep("between",n5),rep("within",n3), rep("between",n6)),
                    Value=c(W6groom,B6groom,W8groom,B8groom,W13groom,B13groom))


ggplot(BWgroom, aes(x=Group, y=Value, fill=Cluster)) +
  facet_grid(.~Group, scales="free")+
  ggtitle("Ecoli Effect on Grooming") +
  labs(y="Ecoli Similarity") +
  theme(legend.position="top")+
  theme(legend.text=element_text(size=15))+
  theme(plot.title = element_text(face="bold"))+
  geom_boxplot()

#####################

BWgroom=data.frame(Group=c(rep("Group I",n1+n4),rep("Group II",n2+n5),rep("Group III",n3+n6)), 
                 Cluster=c(rep("within-clusters",n1),rep("between-clusters",n4), 
                           rep("within-clusters",n2), rep("between-clusters",n5),
                           rep("within-clusters",n3), rep("between-clusters",n6)),
                 Value=c(W6groom,B6groom,W8groom,B8groom,W13groom,B13groom))


ggplot(BWgroom, aes(x=Group, y=Value, fill=Cluster)) +
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



























#####################
#dataframe output for krishna 
##########################
n1=length(yNC13H1)
n2=length(union(nameE8[[2]],nameE8[[3]]))

ind1=rep(nameE[[2]],n1)
ind2=rep(nameE[[2]],each=n1)
EcoliSim=numeric(length(ind1))
for(i in 1:length(ind1))
{
    temp1=which(nameE[[1]]==ind1[i])
    temp2=which(nameE[[1]]==ind2[i])
    if(length(temp1)==0)
      EcoliSim[i]="NA"
    else if(length(temp2)==0)
      EcoliSim[i]="NA"
    else
      EcoliSim[i]=NC13B[temp1,temp2]
}

HuddleBool=numeric(length(ind1))
GroomBool=numeric(length(ind1))
for(i in 1:length(ind1))
{
    temp1=which(nameE[[2]]==ind1[i])
    temp2=which(nameE[[2]]==ind2[i])
    y1=yNC13H1[temp1]
    y2=yNC13H1[temp2]
    if(y1==y2)
      HuddleBool[i]=y1
    
    g1=yNC13G1[temp1]
    g2=yNC13G1[temp2]
    if(g1==g2)
      GroomBool[i]=g1
}
##############
nameU8=union(nameE8[[2]],nameE8[[3]])
ind18=rep(nameU8,n2)
ind28=rep(nameU8,each=n2)
EcoliSim8=numeric(length(ind18))
for(i in 1:length(ind18))
{
  temp1=which(nameE8[[1]]==ind18[i])
  temp2=which(nameE8[[1]]==ind28[i])
  if(length(temp1)==0)
    EcoliSim8[i]="NA"
  else if(length(temp2)==0)
    EcoliSim8[i]="NA"
  else
    EcoliSim8[i]=NC8Ecol1[temp1,temp2]
}

HuddleBool8=numeric(length(ind18))
GroomBool8=numeric(length(ind18))
for(i in 1:length(ind18))
{
  temp1=which(nameE8[[2]]==ind18[i])
  temp2=which(nameE8[[2]]==ind28[i])
  if(length(temp1)==0)
    HuddleBool8[i]="NA"
  else if(length(temp2)==0)
    HuddleBool8[i]="NA"
  else{
    y1=yNC8H1[temp1]
    y2=yNC8H1[temp2]
    if(y1==y2)
      HuddleBool8[i]=y1
  }

  temp1=which(nameE8[[3]]==ind18[i])
  temp2=which(nameE8[[3]]==ind28[i])
  if(length(temp1)==0)
    GroomBool8[i]="NA"
  else if(length(temp2)==0)
    GroomBool8[i]="NA"
  else{
    g1=yNC8G1[temp1]
    g2=yNC8G1[temp2]
    if(g1==g2)
      GroomBool8[i]=g1
  }
}
################
nameU6=union(nameE6[[2]],nameE6[[3]])
n3=length(nameU6)
ind16=rep(nameU6,n3)
ind26=rep(nameU6,each=n3)
EcoliSim6=numeric(length(ind16))
HuddleBool6=numeric(length(ind16))
GroomBool6=numeric(length(ind16))

for(i in 1:length(ind16))
{
  temp1=which(nameE6[[1]]==ind16[i])
  temp2=which(nameE6[[1]]==ind26[i])
  if(length(temp1)==0)
    EcoliSim6[i]="NA"
  else if(length(temp2)==0)
    EcoliSim6[i]="NA"
  else
    EcoliSim6[i]=NC6Ecol1[temp1,temp2]*100
}

for(i in 1:length(ind16))
{
  temp1=which(nameE6[[2]]==ind16[i])
  temp2=which(nameE6[[2]]==ind26[i])
  if(length(temp1)==0)
    HuddleBool6[i]="NA"
  else if(length(temp2)==0)
    HuddleBool6[i]="NA"
  else{
    y1=yNC6H1[temp1]
    y2=yNC6H1[temp2]
    if(y1==y2)
      HuddleBool6[i]=y1
  }
  
  temp1=which(nameE6[[3]]==ind16[i])
  temp2=which(nameE6[[3]]==ind26[i])
  if(length(temp1)==0)
    GroomBool6[i]="NA"
  else if(length(temp2)==0)
    GroomBool6[i]="NA"
  else{
    g1=yNC6G1[temp1]
    g2=yNC6G1[temp2]
    if(g1==g2)
      GroomBool6[i]=g1
  }
}
###############
CageID=c(rep("Cage 13",n1^2), rep("Cage 8",n2^2),rep("Cage6",n3^2))
Individual1=c(ind1,ind18,ind16)
Individual2=c(ind2,ind28,ind26)
EcoliSimilarity=c(EcoliSim,EcoliSim8,EcoliSim6)
Huddling=c(HuddleBool,HuddleBool8,HuddleBool6)
Grooming=c(GroomBool,GroomBool8,GroomBool6)

StringData=data.frame(cbind(CageID,Individual1,Individual2,EcoliSimilarity,Huddling,Grooming))
write.csv(StringData,file="String_All_three_cages.csv")




#######################
#provide the Ecoli with-between Matrix entry for Krishna 
#################

write.csv(Groom13Effect,file="EcoliGroomCage13.csv")
write.csv(Huddle13Effect,file="EcoliHuddleCage13.csv")
write.csv(Groom8Effect,file="EcoliGroomCage8.csv")
write.csv(Huddle8Effect,file="EcoliHuddleCage8.csv")
write.csv(Groom6Effect,file="EcoliGroomCage6.csv")
write.csv(Huddle6Effect,file="EcoliHuddleCage6.csv")
###########################
name13=colnames(NC13B)
name8=colnames(NC8Ecol1)
name6=colnames(NC6Ecol1)

save(DCG.NCB,DCG.8Ecoli1, DCG.6Ecoli1, name13,name8,name6,
     file = "EcoliTree.RData")


#########################



##################
#ROC curve
#################
library(ROCR)
data(ROCR.simple)
pred <- prediction( ROCR.simple$predictions, ROCR.simple$labels)
perf <- performance(pred,"tpr","fpr")
plot(perf)

W6huddle
B6huddle

ks.test(as.numeric(W6huddle),as.numeric(B6huddle))
ks.test(as.numeric(W8huddle),as.numeric(B8huddle))
ks.test(as.numeric(W13huddle),as.numeric(B13huddle))

ks.test(as.numeric(W6groom),as.numeric(B6groom))
ks.test(as.numeric(W8groom),as.numeric(B8groom))
ks.test(as.numeric(W13groom),as.numeric(B13groom))
########









#######################3
EntropySampling=function(y1,y2,name1,name2, iter)
{
  N2=length(unique(y2))
  N1=length(unique(y1))
  n1=length(y1)
  n2=length(y2)
  EntropyDist=c()
  
  for(i in 1:iter)
  {
    indSample=sample(n2)
    y2Sample=y2[indSample]
    EntropyDist[i]=Entropy_general(y1,y2Sample,name1,name2)$Sum
  }
  EntropyObs=Entropy_general(y1,y2,name1,name2)$Sum
  
  return(list(Dist=EntropyDist,Obs=EntropyObs))
  
}
#############################################
m <- matrix(c(1,1,1,2,3,4,5,6,7,8,9,10),nrow = 4,ncol = 3,byrow = TRUE)
layout(mat = m,heights = c(0.1,0.3,0.3,0.3))

par(mar = c(2,2,1,1))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend(x = "top",inset = 0, c("Sampled ", "Observed"), 
       bty="n",col=1:2,lwd=5, cex=1.2, pt.cex=1.8,horiz = TRUE)


A=EntropySampling(yNC6E1,yNC6G1,nameE6[[1]],nameE6[[3]],50)
plot(density(A$Dist),xlab="Entropy", lwd=3,
     main="Group1: Ecoli vs Grooming ")
abline(v=A$Obs,col=2,lwd=3)

A=EntropySampling(yNC8E1,yNC8G1,nameE8[[1]],nameE8[[3]],50)
plot(density(A$Dist),xlab="Entropy", lwd=3,
     main="Group2: Ecoli vs Grooming ")
abline(v=A$Obs,col=2,lwd=3)

A=EntropySampling(yNC13E1,yNC13G1,nameE[[1]],nameE[[3]],50)
plot(density(A$Dist),xlab="Entropy", lwd=3,
     main="Group3: Ecoli vs Grooming ")
abline(v=A$Obs,col=2,lwd=3)

#############################################
A=EntropySampling(yNC6E1,yNC6H1,nameE6[[1]],nameE6[[2]],50)
plot(density(A$Dist),xlab="Entropy", lwd=3,
     main="Group1: Ecoli vs Huddling")
abline(v=A$Obs,col=2,lwd=3)

A=EntropySampling(yNC8E1,yNC8H1,nameE8[[1]],nameE8[[2]],50)
plot(density(A$Dist),xlab="Entropy", lwd=3,
     main="Group2: Ecoli vs Huddling")
abline(v=A$Obs,col=2,lwd=3)

A=EntropySampling(yNC13E1,yNC13H1,nameE[[1]],nameE[[2]],50)
plot(density(A$Dist),xlab="Entropy", lwd=3,
     main="Group3: Ecoli vs Huddling")
abline(v=A$Obs,col=2,lwd=3)
####################
A=EntropySampling(yNC6G1,yNC6H1,nameE6[[3]],nameE6[[2]],50)
plot(density(A$Dist),xlab="Entropy", lwd=3,
     main="Group1: Grooming vs Huddling")
abline(v=A$Obs,col=2,lwd=3)

A=EntropySampling(yNC8G1,yNC8H1,nameE8[[3]],nameE8[[2]],50)
plot(density(A$Dist),xlab="Entropy", lwd=3,
     main="Group2: Grooming vs Huddling")
abline(v=A$Obs,col=2,lwd=3)

A=EntropySampling(yNC13G1,yNC13H1,nameE[[3]],nameE[[2]],50)
plot(density(A$Dist),xlab="Entropy", lwd=3,
     main="Group3: Grooming vs Huddling")
abline(v=A$Obs,col=2,lwd=3)


save(yNC13H1, yNC8H1, yNC6H1, yNC13G1, yNC6G1,yNC8G1,yAgg6,yAgg8,yAgg13, nameE13,nameE6,nameE8,
     DCG.agg6,DCG.agg8,DCG.agg13,DCG13.groom1, DCG.groom81,DCG.groom61,file = "Tree_and_Membership.RData")












