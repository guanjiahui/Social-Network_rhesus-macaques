EcoliEffect=function(y1,name1,ID2,Mat)
{
  N=length(unique(y1))
  
  check=0
  
  for (i in 1:N)
  {
    loc=which(y1==i)
    ID=name1[loc]
    
    ID=intersect(ID,ID2)
    n=length(ID)
    within=0
    
    if (n==0)
      check=c(check,i)
  }
  
  cat(check,"\n")
  
  if(length(check)>1)
  {
    EcoliWB=matrix(0,N-length(check)+1,N-length(check)+1)
    With=numeric(N-length(check)+1)
    check=check[-1]
    newSeq=(1:N)[-check]
  }
  
  else
  {
    EcoliWB=matrix(0,N,N)
    With=numeric(N)
    newSeq=1:N
  }
  
  
  
  
  ###########################
  #first compute the diganol
  ###########################
  for(i in newSeq)
  {
    loc=which(y1==i)
    ID=name1[loc]
    
    ID=intersect(ID,ID2)
    n=length(ID)
    within=0
    
    if (n==0)
      check=i
    else{
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
    
    
  }
  
  ##########################
  #compute the upper off-diaganol
  ##########################
  for(i in newSeq)
  {
    for(j in newSeq)
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
  
  
  
  
  return(list(W=With, E=EcoliWB))
}