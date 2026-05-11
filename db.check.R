#This function summarizes the DBSCAN clusterings obtainable by iteratively changing epsilon parameter.

#dat is a scaled matrix/data.frame
#eps is a value for Epsilon parameter
#minp is a value for MinPoints parameter
#prec measures the variations between epsilon parameters in the different clusterings tried

db.check=function(dat,eps,minp,prec=0.005){
  dat=as.matrix(dat)
  if (is.na(eps)){
    eps=0
  }
  if (minp>nrow(dat)){
    minp=3
  }#MinPoints must be less than n
  seq1=seq(eps,eps+10,by=.2)
  ok=F
  i=1
  while (ok==F){
    db=dbscan(dat,seq1[i],minp)
    if (any(db$cluster==0)==T){
      i=i+1
    }else{
      max=seq1[i]
      ok=T
    }
  } #Maximum value of epsilon, when no noise is found due to the high value of epsilon, the maximum is determined 
  seq2=seq(0,max,by=prec)
  res=matrix(ncol=9,nrow=length(seq2))#Results matrix
  res=as.data.frame(res)
  colnames(res)=c("Eps","MinPoints","Noise","Clusters","Dim. 1", "Dim. 2","Dim. 3", "Dim. 4","Dim. 5")
  for (j in 1:length(seq2)){
    db=dbscan(dat,seq2[j],minp)
    res[j,1]=db$eps
    res[j,2]=db$minPts
    ncl=max(db$cluster)
    res[j,4]=ncl
    a=hist(db$cluster,breaks=-1:ncl,plot=F)$counts
    res[j,3]=a[1]
    if (ncl>=5){
      a=sort(a[2:length(a)])
      res[j,5:9]=a[1:5]
    }else{
      if (ncl==0){
        res[j,5:9]=NA
      }else{
        elim=5-ncl
        res[j,c((10-elim):9)]=NA
        res[j,5:(4+ncl)]=a[2:length(a)]
      }
    }
  }
  return(res)
}
