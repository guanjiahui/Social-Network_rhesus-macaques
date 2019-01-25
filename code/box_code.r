#####################
#boxplot code
#######################3
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


################
#wilcoxon rank sum test 
ks.test(as.numeric(W6AGG),as.numeric(B6AGG))
ks.test(as.numeric(W8AGG),as.numeric(B8AGG))
ks.test(as.numeric(W13AGG),as.numeric(B13AGG),alternative = "less")


############

#Aggression 
wilcox.test(as.numeric(W6AGG),as.numeric(B6AGG))
wilcox.test(as.numeric(W8AGG),as.numeric(B8AGG))
wilcox.test(as.numeric(W13AGG),as.numeric(B13AGG))


#################
#Huddling 

wilcox.test(as.numeric(W6huddle),as.numeric(B6huddle))
wilcox.test(as.numeric(W8huddle),as.numeric(B8huddle))
wilcox.test(as.numeric(W13huddle),as.numeric(B13huddle))

#Grooming
wilcox.test(as.numeric(W6groom),as.numeric(B6groom))
wilcox.test(as.numeric(W8groom),as.numeric(B8groom))
wilcox.test(as.numeric(W13groom),as.numeric(B13groom))
########


#########
#compute z-score 
###############
#write a function 
z_wilcox_rank=function(A,B){
  df=data.frame(value=c(A,B), factor=as.factor(c(rep(1,length(A)),rep(2,length(B)))))
  
  wilcoxsign_test(value~factor,data=df, distribution="exact")
}

#Aggression 
z_wilcox_rank(as.numeric(W6AGG),as.numeric(B6AGG))
z_wilcox_rank(as.numeric(W8AGG),as.numeric(B8AGG))
z_wilcox_rank(as.numeric(W13AGG),as.numeric(B13AGG))


#################
#Huddling 

z_wilcox_rank(as.numeric(W6huddle),as.numeric(B6huddle))
z_wilcox_rank(as.numeric(W8huddle),as.numeric(B8huddle))
z_wilcox_rank(as.numeric(W13huddle),as.numeric(B13huddle))

#Grooming
z_wilcox_rank(as.numeric(W6groom),as.numeric(B6groom))
z_wilcox_rank(as.numeric(W8groom),as.numeric(B8groom))
z_wilcox_rank(as.numeric(W13groom),as.numeric(B13groom))
########


##########################
#Just Cage NC8, do a Wilcoxon rank-sum TEst 
#for comparing mean %E coli similairty of monkey-monkey  islolates vs monkey-env isolates 
library(readr)
NC8_IndivEnvEcoliString <- read_csv("~/Dropbox/Research/VET_SNH_Monkey/DATA/Major Revision/NC8_IndivEnvEcoliString.csv")

MonkeyEnv=NC8_IndivEnvEcoliString$Ecoli

#make a Ecoli matrix into a 1d array 

MonkeyMonkey=0
for(i in 1:nrow(NC8Ecol1))
{
  for(j in 1:ncol(NC8Ecol1))
  {
    if (i>j)
      MonkeyMonkey=c(MonkeyMonkey,NC8Ecol1[i,j])
  }
}

MonkeyMonkey=MonkeyMonkey[-1]

#Make a comparison using Wilcoxon

wilcox.test(MonkeyMonkey,MonkeyEnv)

#within community monkey-monkey vs monkey-env
wilcox.test(as.numeric(W8huddle),MonkeyEnv)
wilcox.test(as.numeric(W8groom),MonkeyEnv)
wilcox.test(as.numeric(W8AGG),MonkeyEnv)

wilcox.test(MonkeyMonkey,MonkeyEnv,alternative = "greater")
wilcox.test(as.numeric(W8huddle),MonkeyEnv,alternative = "greater")
wilcox.test(as.numeric(W8groom),MonkeyEnv,alternative = "greater")
wilcox.test(as.numeric(W8AGG),MonkeyEnv,alternative = "greater")


##############
wilcox.test(MonkeyMonkey,MonkeyEnv,alternative = "greater")
wilcox.test(as.numeric(W8huddle),MonkeyEnv,alternative = "greater")
wilcox.test(as.numeric(W8groom),MonkeyEnv,alternative = "greater")
wilcox.test(as.numeric(W8AGG),MonkeyEnv,alternative = "greater")






##############
#z-score 
#using asymptotic 
z_wilcox_rank=function(A,B){
  df=data.frame(value=c(A,B), factor=as.factor(c(rep(1,length(A)),rep(2,length(B)))))
  
  wilcoxsign_test(value~factor,data=df)
}
###############
z_wilcox_rank(as.numeric(W8huddle),MonkeyEnv)
z_wilcox_rank(as.numeric(W8groom),MonkeyEnv)
z_wilcox_rank(as.numeric(W8AGG),MonkeyEnv)

z_wilcox_rank(MonkeyMonkey,MonkeyEnv)
z_wilcox_rank(as.numeric(W8huddle),MonkeyEnv)
z_wilcox_rank(as.numeric(W8groom),MonkeyEnv)
z_wilcox_rank(as.numeric(W8AGG),MonkeyEnv)


##############
z_wilcox_rank(MonkeyMonkey,MonkeyEnv)
z_wilcox_rank(as.numeric(W8huddle),MonkeyEnv)
z_wilcox_rank(as.numeric(W8groom),MonkeyEnv)
z_wilcox_rank(as.numeric(W8AGG),MonkeyEnv)






################
#"permutation"
n1=length(W8huddle)
n2=length(W8groom)
n3=length(W8AGG)
A=sample(MonkeyMonkey,n1)
boxplot(as.numeric(W8huddle),A)
wilcox.test(as.numeric(W8huddle),A,alternative = "greater")
wilcox.test(as.numeric(W8groom),sample(MonkeyMonkey,n2),alternative = "greater")
wilcox.test(as.numeric(W8AGG),sample(MonkeyMonkey,n3),alternative = "greater")

n1=length(W6huddle)
n2=length(W6groom)
n3=length(W6AGG)
A=sample(MonkeyMonkey,n1)
wilcox.test(as.numeric(W6huddle),A,alternative = "greater")
wilcox.test(as.numeric(W6groom),sample(MonkeyMonkey,n2),alternative = "greater")
wilcox.test(as.numeric(W6AGG),sample(MonkeyMonkey,n3),alternative = "greater",exact=FALSE)


n1=length(W13huddle)
n2=length(W13groom)
n3=length(W13AGG)
A=sample(MonkeyMonkey,n1)
wilcox.test(as.numeric(W13huddle),A,alternative = "greater")
wilcox.test(as.numeric(W13groom),sample(MonkeyMonkey,n2),alternative = "greater")
wilcox.test(as.numeric(W13AGG),sample(MonkeyMonkey,n3),alternative = "greater")


n1=length(W6Sp)
n2=length(W8Sp)
n3=length(W13Sp)

A=sample(MonkeyMonkey,n1)
wilcox.test(as.numeric(W6Sp),A,alternative = "greater")
wilcox.test(as.numeric(W8Sp),sample(MonkeyMonkey,n2),alternative = "greater")
wilcox.test(as.numeric(W13Sp),sample(MonkeyMonkey,n3),alternative = "greater")

