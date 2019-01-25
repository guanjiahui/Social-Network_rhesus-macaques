load("NC.RData")


temp=c(0.008,0.07,0.15,0.3,1)
Ens13.groom1=Eigen.plot(temp, selected.id=c(1,2,3,4,5),win_prob1)

DCG13.groom1=DCGtree.plot(num.clusters.selected=c(1,2,3,8,18),
                          "NC13GroomR1 tree",Ens13.groom1,temp)
plot(DCG13.groom1,hang=-1,main="NC13GroomR1 tree")
save(DCG13.groom1,file="NC13Groom1.RData")

