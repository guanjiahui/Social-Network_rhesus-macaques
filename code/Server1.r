load("NC.RData")


temp=c(0.13,0.25,1)
Ens.8Ecoli1=Eigen.plot2(temp, selected.id=c(1,2,3),NC8Ecol1)
###########
DCG.8Ecoli1=DCGtree.plot(num.clusters.selected=c(1,2,4),
                         "NC8 Ecoli Baseline",Ens.8Ecoli1,temp)
plot(DCG.8Ecoli1,hang=-1,main="NC8 Ecoli Baseline")

save(DCG.8Ecoli1,file="NC8Ecoli.RData")



