#Ecoli Heatmap 
#################
library(gplots)

##############
#NC6
heatmap.2(NC6Ecol1, Rowv=as.dendrogram(DCG.6Ecoli1), 
          Colv = as.dendrogram(DCG.6Ecoli1), trace="none",
          density.info = 'none',key= FALSE,
          labRow =FALSE, labCol = FALSE)

#########
#NC8

heatmap.2(NC8Ecol1, Rowv=as.dendrogram(DCG.8Ecoli1), 
          Colv = as.dendrogram(DCG.8Ecoli1), trace="none",
          density.info = 'none',key= FALSE,
          labRow =FALSE, labCol = FALSE)

#########
#NC13
heatmap.2(NC13B, Rowv=as.dendrogram(DCG.NCB), 
          Colv = as.dendrogram(DCG.NCB), trace="none",
          density.info = 'none',key= FALSE,
          labRow =FALSE, labCol = FALSE)

###################################################
#NC6 
heatmap.2(NC6Huddel1*20,Rowv=as.dendrogram(DCG6.huddle1),Colv = as.dendrogram(DCG6.huddle1),
          trace="none", density.info = 'none',
          labRow =FALSE, labCol = FALSE)

heatmap.2(GroomDist6,Rowv=as.dendrogram(DCG.groom61),Colv = as.dendrogram(DCG.groom61),
          trace="none", density.info = 'none',
          labRow =FALSE, labCol = FALSE)

heatmap.2(Agg8,Rowv=as.dendrogram(DCG.agg8),Colv = as.dendrogram(DCG.agg8),
          trace="none", density.info = 'none',
          labRow =FALSE, labCol = FALSE)

plot(DCG6.huddle1,hang=-1,labels= FALSE,main= "Huddling",xlab = "")
abline(h=2,col="red",lwd=3)

plot(DCG.groom61,hang=-1,labels= FALSE,main= "Grooming" ,xlab = "" )
abline(h=0.2,col="red",lwd=3)

plot(DCG.agg8,hang=-1,labels= FALSE,main= "Aggression",xlab = ""  )
abline(h=0.15,col="red",lwd=3)

#########################################################
#construct a permutation test for the significane of DCG
########################################################

######
Avg_withCluster=function(y_member, avg)
{
  N=length(y_member)
  n_cluster=max(unique(y_member))
  Obs=0
  
  for(i in 1:n_cluster)
  {
    id=which(y_member==i)
    n_member=length(id)
    
    obs_mean=mean(avg[id])
    Obs=Obs+obs_mean
  }
  return(Obs)
}
###########################
Perm_DCG=function(y_member, Mat,sample_iter)
{
  N=length(y_member)
  col_avg=colMeans(Mat)
  Sample_mean=c()
  
  if(N!=dim(Mat)[1])
    return("Error, the dimension doesn't match!")
  
  else
  {
     Obs_mean=Avg_withCluster(y_member, col_avg)
     
     ##sampling
     for (i in 1:sample_iter)
     {
       seq=sample(y_member)
       Sample_mean[i]=Avg_withCluster(y_member, col_avg)
     }
    
     #calculate the p-value
     p_value=1-ecdf(Sample_mean)(Obs_mean)
  }
  
  return(p_value)
}
############################
#yNC8E1
#yNC6E1
#yNC13E1

p_NC6=Perm_DCG(yNC6E1,NC6Ecol1,1000)
p_NC8=Perm_DCG(yNC8E1,NC8Ecol1,1000)
p_NC13=Perm_DCG(yNC13E1,NC13B,1000)


#Permutation test for aggression data 
p_agg6=Perm_DCG(yAgg6,Agg6, 1000)
p_agg8=Perm_DCG(yAgg8,Agg8, 1000)
p_agg13=Perm_DCG(yAgg13,Agg13, 1000)


#Permutation test for huddling data 

p_huddle6=Perm_DCG(y6Huddle1)



#Permutation test for gromming data 

p_groom6=Perm_DCG(yNC6G1,GroomDist6,1000)
p_groom8=Perm_DCG(yNC8G1,GroomDist8,1000)
p_groom13=Perm_DCG(yNC13G1,GroomDist13,1000)
  
#Permutation test for space-use data 
p_space6=Perm_DCG(ySpace6,NC6_Space,1000)
p_space8=Perm_DCG(ySpace8,NC8_Space,1000)
p_space13=Perm_DCG(ySpace13,NC13_Space,1000)
