#Finds mismatches between clustering and ground truth

#k is the clustering size
#methods is a vector of distance methods used for hierarchical clustering
#linkages is a vector of linkage methods used for hierarchical clustering
#eps is a vector of Epsilon parameters used for DBSCAN
#mins is a vector of MinPoints parameters used for DBSCAN
#means is the agg1 dataset
#logical, used to select a ground truth between Diet/Binocularity and Aquatic (both with k=3), Aquatic->TRUE

mismatch=function (k,methods,linkages,eps,mins,means,a){
  library(Jacc)
  library(dbscan)
  library(dplyr)
  
  err=data.frame(Name=c("km.e.1","km.e.2","km.e.3","km.e.4","km.e.5","km.e.6","km.r.1","km.r.2","km.r.3","km.r.4","km.r.5","km.r.6","km.p.1","km.p.2","km.p.3","km.p.4","km.p.5","km.p.6","km.c.1","km.c.2","km.c.3","km.c.4","km.c.5","km.c.6","km.m.1","km.m.2","km.m.3","km.m.4","km.m.5","km.m.6","hc.e.1","hc.e.2","hc.e.3","hc.e.4","hc.e.5","hc.e.6","hc.r.1","hc.r.2","hc.r.3","hc.r.4","hc.r.5","hc.r.6","hc.p.1","hc.p.2","hc.p.3","hc.p.4","hc.p.5","hc.p.6","hc.c.1","hc.c.2","hc.c.3","hc.c.4","hc.c.5","hc.c.6","hc.m.1","hc.m.2","hc.m.3","hc.m.4","hc.m.5","hc.m.6","db.e.1","db.e.2","db.e.3","db.e.4","db.e.5","db.e.6","db.r.1","db.r.2","db.r.3","db.r.4","db.r.5","db.r.6","db.p.1","db.p.2","db.p.3","db.p.4","db.p.5","db.p.6","db.c.1","db.c.2","db.c.3","db.c.4","db.c.5","db.c.6","db.m.1","db.m.2","db.m.3","db.m.4","db.m.5","db.m.6"),Shape=c(rep("Ellipsoid",6),rep("Round",6),rep("Pyramidal",6),rep("Complex",6),rep("Microglia",6),rep("Ellipsoid",6),rep("Round",6),rep("Pyramidal",6),rep("Complex",6),rep("Microglia",6),rep("Ellipsoid",6),rep("Round",6),rep("Pyramidal",6),rep("Complex",6),rep("Microglia",6)),Layer=c(rep(1:6,5),rep(1:6,5),rep(1:6,5)),Type=c(rep("km",30),rep("hc",30),rep("db",30)))
  err$Shape=as.factor(err$Shape)
  err$Layer=as.factor(err$Layer)
  err$Type=as.factor(err$Type)
  errs=c()#Setting up the output data.frame and errs
  
  #K-Means
  for (i in 1:30){
    s=subset(means,means$Shape==err$Shape[i] & means$Layer==err$Layer[i])
    if (k==2){
      par=s$Cetacean
    }
    if (k==3){
      if (a==T){
        par=s$Aquatic
      }
      else{
        par=s$Diet_Binocularity
      }
    }
    if (k==5){
      par=s$Order
    }
    cl=kmeans(scale(s[,10:19]),centers = k,nstart = 100)#Clustering
    ja=Jacc(cl$cluster,par)#Cluster matching based on F1-Score maximization 
    mod=ja$Results.Bijective$Matches
    for (j in 1:nrow(s)){
      nm=cl$cluster[j]#Calculated cluster
      if (par[j]!=mod[nm]){#If mismatch
        im=s$Image[j]
        cat=mod[nm]#Wrong partition
        v=paste(as.character(im),"-",as.character(cat))#Mismatch name
        if (length(errs)==0){#First mismatch
          errs=c(errs,v)
          err=data.frame(err,c(rep(0,i-1),1,rep(NA,90-i)))#Mismatch recording
          colnames(err)=c(colnames(err)[1:ncol(err)-1],v)
        }else{
          w=which(errs==v)
          if (length(w)==1){#New type of mismatch?
            err[i,w+4]=1
          }else{
            errs=c(errs,v)
            err=data.frame(err,c(rep(0,i-1),1,rep(NA,90-i)))
            colnames(err)=c(colnames(err)[1:ncol(err)-1],v)
          }
        }
      }
    }
    w=which(is.na(err[i,])==TRUE)
    err[i,w]=0
  }
  for (i in 1:30){#Hierarchical clustering
    s=subset(means,means$Shape==err$Shape[i] & means$Layer==err$Layer[i])
    if (k==2){
      par=s$Cetacean
      met=methods[1+(i-1)*4]
      lin=linkages[1+(i-1)*4]#Specification of distance and linkage methods
    }
    if (k==3){
      if (a==T){
        par=s$Aquatic
        met=methods[2+(i-1)*4]
        lin=linkages[2+(i-1)*4]
      }
      else{
        par=s$Diet_Binocularity
        met=methods[3+(i-1)*4]
        lin=linkages[3+(i-1)*4]
      }
    }
    if (k==5){
      par=s$Order
      met=methods[i*4]
      lin=linkages[i*4]
    }
    cl=hclust(dist(scale(s[,10:19]),method = met),method = lin)#Clustering
    ja=Jacc(cutree(cl,k),par)#Cluster matching based on F1-Score maximization 
    mod=ja$Results.Bijective$Matches
    for (j in 1:nrow(s)){
      nm=cutree(cl,k)[j]#Calculated cluster
      if (par[j]!=mod[nm]){#If mismatch
        im=s$Image[j]
        cat=mod[nm]#Wrong partition
        v=paste(as.character(im),"-",as.character(cat))#Mismatch name
        if (length(errs)==0){#First mismatch
          errs=c(errs,v)
          err=data.frame(err,c(rep(0,i+29),1,rep(NA,60-i)))#Mismatch recording
          colnames(err)=c(colnames(err)[1:ncol(err)-1],v)
        }else{
          w=which(errs==v)
          if (length(w)==1){#New type of mismatch?
            err[i+30,w+4]=1
          }else{
            errs=c(errs,v)
            err=data.frame(err,c(rep(0,i+29),1,rep(NA,60-i)))
            colnames(err)=c(colnames(err)[1:ncol(err)-1],v)
          }
        }
      }
    }
    w=which(is.na(err[i+30,])==TRUE)
    err[i+30,w]=0
  }
  nn=which(is.na(eps)==T)#Identification of invalid DBSCAN clusterings
  vect=c(1:30)
  if (k==2){
    b=1
  }
  if (k==3){
    if (a==T){
      b=2
    }else{
      b=3
    }
  }
  if (k==5){
    b=4
  }
  nn.2=c()
  for (i in 1:30){
    i2=b+(i-1)*4
    if (i2 %in% nn){
      nn.2=c(nn.2,i)
    }
  }
  if (length(nn.2)>0){
    vect=vect[-nn.2]#Dropping vect values corresponding to invalid clusterings
  }
  for (i in vect){
    s=subset(means,means$Shape==err$Shape[i] & means$Layer==err$Layer[i])
    if (k==2){
      par=s$Cetacean
      ep=eps[1+(i-1)*4]
      min=mins[1+(i-1)*4]#Specification of Epsilon and MinPoints parameters
    }
    if (k==3){
      if (a==T){
        par=s$Aquatic
        ep=eps[2+(i-1)*4]
        min=mins[2+(i-1)*4]
      }
      else{
        par=s$Diet_Binocularity
        ep=eps[3+(i-1)*4]
        min=mins[3+(i-1)*4]
      }
    }
    if (k==5){
      par=s$Order
      ep=eps[i*4]
      min=mins[i*4]
    }
    cl=dbscan(scale(s[,10:19]),ep,min)#Clustering
    if (length(which(cl$cluster==0))==0){
      ja=Jacc(cl$cluster,par)  
    }else{
      ja=Jacc(cl$cluster,par,unassigned = c(0))
    }#Cluster matching based on F1-Score maximization 
    mod=ja$Results.Bijective$Matches
    for (j in 1:nrow(s)){
      nm=cl$cluster[j]#Calculated cluster
      if (nm==0){#If noise
        im=s$Image[j]
        cat="Noise"
      }else{#If mismatch
        if (par[j]!=mod[nm]){
          im=s$Image[j]
          cat=mod[nm]
        }
      }
      v=paste(as.character(im),"-",as.character(cat))#Mismatch name
      if (length(errs)==0){#First mismatch
        errs=c(errs,v)
        err=data.frame(err,c(rep(0,i+59),1,rep(NA,30-i)))#Mismatch recording
        colnames(err)=c(colnames(err)[1:ncol(err)-1],v)
      }else{
        w=which(errs==v)
        if (length(w)==1){#New type of mismatch?
          err[i+60,w+4]=1
        }else{
          errs=c(errs,v)
          err=data.frame(err,c(rep(0,i+59),1,rep(NA,30-i)))
          colnames(err)=c(colnames(err)[1:ncol(err)-1],v)
        }
      }
    }
    w=which(is.na(err[i+60,])==TRUE)
    err[i+60,w]=0
  }
  nn.2=nn.2+60
  err[nn.2,5:ncol(err)]=NA#Treatment of invalid clusterings
  w=dplyr::select(err,contains("Noise"))
  w.c=dplyr::select(err,!contains("Noise"))
  err=cbind(w.c,w)
  ww=c(colMeans(na.omit(w.c[,5:ncol(w.c)])),colMeans(na.omit(w[which(err$Type=="db"),])))#Expected occurrence of specific mismatches and noise istances
  names(ww)=c(colnames(w.c[,5:ncol(w.c)]),colnames(w))
  nom=colnames(err)
  for (h in 1:(ncol(err)-4)){#Marking mismatches by frequency
    if (ww[h]>.8){
      tag="***"
    }
    if (ww[h]>.6 & ww[h]<=.8){
      tag="**"
    }
    if (ww[h]>.4 & ww[h]<=.6){
      tag="*"
    }
    if (ww[h]>.2 & ww[h]<=.4){
      tag=".."
    }
    if (ww[h]<=.2){
      tag="."
    }
    nom[h+4]=paste(nom[h+4],tag)  
  }
  colnames(err)=nom
  ee=list(Errors=err,Weight=ww)
  return(ee)
}
