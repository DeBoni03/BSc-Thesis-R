db.ground.truth=function(data,k,ground_truth){
  library(cluster)
  library(clusterCrit)
  library(ggfortify)
  source("db.estimate.R")
  source("db.check.R")
  n=nrow(data)
  s1=data[,1:9]
  s1=as.integer(s1[ground_truth][[1]])#Ground truth
  s2=scale(data[,10:19])#Scaled data
  d=dist(s2)
  pc=prcomp(s2,rank.=2)#Principal components for Biplot
  ris=db.estimate(d,2)#Epsilon and MinPoints selection
  if (is.na(ris$Epsilon.2)){
    ris$Epsilon.2=0
  }
  db=dbscan(s2,ris$Epsilon.2,ris$MinPoints.2)#DBSCAN clustering
  if (max(db$cluster)!=k){#If optimal parameters don't bring the desired k, parameters are recalculated 
    ok=F
    i=0
    while (ok==F & i<=5){
      ris=db.estimate(d,2,w=F,show=F,p=10,auto=T,which=i)#Changes in the 'which' parameter lead to changes in Epsilon and MinPoints
      if (is.na(ris$Epsilon.2)){
        ris$Epsilon.2=0
      }
      db=dbscan(s2,ris$Epsilon.2,ris$MinPoints.2)#DBSCAN clustering
      check=db.check(s2,ris$Epsilon.2,ris$MinPoints.2)#Checking for different configurations of Epsilon given different configurations of MinPoints
      if (any(check$Clusters==k)){#Are there clusterings with the desired k?
        check=subset(check,check$Clusters==k)
        db.r=dbscan(s2,check$Eps[1],check$MinPoints[1])#DBSCAN clustering for whether better clusterings are not found
        if (any(check$Noise<n/5) | (any(check$Noise<n/5)==F & n<=30)){#Are there low-noise clusterings?
          if (any(check$Noise<n/5) ){
            check=subset(check,check$Noise<n/5)
            db.r=dbscan(s2,check$Eps[1],check$MinPoints[1])#DBSCAN clustering for whether better clusterings are not found
            a=c()
            for (h in 1:nrow(check)){
              sor=sort(na.omit(t(check[h,5:9])),decreasing = T)
              a[h]=(sor[1]+check$Noise[h])/n
            }
            if (any(a<.8)){#Are there balanced clusterings?
              v=c()
              for (j in 1:ncol(t(check[,5:9]))){
                v[j]=var(na.omit(t(check[,5:9])[,j]))
              }
              v=(v-min(v))/(max(v)-min(v))#Normalized variances of top 5 clusters dimensions
              noi=(check$Noise-min(check$Noise))/(max(check$Noise)-min(check$Noise))#Normalized noise frequencies
              if (any(is.na(v)==F) & any(is.na(noi)==F)){
                ind=v/4+noi#Low-noise/low clusters dimensions' heterogeneity index: low values are preferred
                m=min(ind)
                w=which(ind==m)[1]
              } else{
                w=1
              }
              db=dbscan(s2,check[w,]$Eps,check[w,]$MinPoints)#DBSCAN clustering
              ok=T
            }
          }
        }
      }
      if (ok==F){
        i=i+1#Index incrementation: new MinPoints configurations are tried
      }
    }
  }else{
    ok=T#Cycle exit
  }
  if (ok==T | exists("db.r")){
    if (exists("db.r")){
      db=db.r
    }
    print(autoplot(pc,col=(db$cluster+c(rep(1,n)))))
    print(autoplot(pc,col=as.numeric(s1)))#Biplots
    if (any(db$cluster==0)){#Validation indexes calculations: different cases due to eventual presence of noise
      w=which(db$cluster==0)
      db$cluster[w]=9
      t=table(db$cluster[-w],s1[-w])
      if (dim(t)[1]==dim(t)[2]){
        f=function(tab) {
          r=tryCatch({
            fisher.test(tab)
          }, error = function(e) {
            message("Error. Simulating p-value.")
            return(fisher.test(tab, simulate.p.value = TRUE))
          })
          return(r)
        }
        p=f(t)$p.value
        cs=chisq.test(t)$statistic
        chin=cs/(nrow(s2[-w,])*(k-1))
      }else{
        p=NA
        chin=NA
      }
      cvi=cluster.stats(d,clustering=db$cluster,alt.clustering = s1,noisecluster=T)
      crit=c(k,cvi$avg.silwidth,length(which(silhouette(db$cluster[-w],dist(s2[-w,]))[,3]<0))/nrow(s2[-w,]),cvi$ch,cvi$dunn2,cvi$vi/log(nrow(s2[-w,])),cvi$corrected.rand,chin,p) #Validation indexes
    }
    else{
      t=table(db$cluster,s1)
      f=function(tab) {
        r=tryCatch({
          fisher.test(tab)
        }, error = function(e) {
          message("Error. Simulating p-value.")
          return(fisher.test(tab, simulate.p.value = TRUE))
        })
        return(r)
      }
      p=f(t)$p.value
      cs=chisq.test(t)$statistic
      chin=cs/(nrow(s2)*(k-1))
      cvi=cluster.stats(d,clustering=db$cluster,alt.clustering = s1,noisecluster=F)
      crit=c(k,cvi$avg.silwidth,length(which(silhouette(db$cluster,d[,3]<0)))/n,cvi$ch,cvi$dunn2,cvi$vi/log(n),cvi$corrected.rand,chin,p) #Validation indexes
    }
    names(crit)=c("K","Avg. Sil.","Neg.Sil.","CH","Dunn*","VI","ARI","Chi","P-value")
    l=list("DBSCAN"=db,"CVI"=crit)
  }
  else{#No valid DBSCAN clustering has been found
    crit=c(rep(NA,9))
    l=list("CVI"=crit)
  }
  return(l)
}
