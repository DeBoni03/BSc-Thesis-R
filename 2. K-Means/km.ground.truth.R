km.ground.truth=function(data,k,ground_truth){
  library(cluster)
  library(clusterCrit)
  library(ggfortify)
  s1=data[,1:9]
  s1=as.integer(s1[ground_truth][[1]])#Ground truth
  s2=scale(data[,10:19])#Scaled data
  d=dist(s2)
  pc=prcomp(s2,rank.=2)#Principal components for Biplot
  km=kmeans(s2,centers=k,nstart=100)#Clustering
  print(autoplot(pc,col=km$cluster))
  print(autoplot(pc,col=as.numeric(s1)))#Biplots
  t=table(km$cluster,s1)
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
  cvi=cluster.stats(d,clustering=km$cluster,alt.clustering = s1)
  crit=c(k,cvi$avg.silwidth,length(which(silhouette(km$cluster,d)[,3]<0))/nrow(s2),cvi$ch,cvi$dunn2,cvi$vi/log(nrow(s2)),cvi$corrected.rand,chin,p) #Validation indexes
  names(crit)=c("K","Avg. Sil.","Neg.Sil.","CH","Dunn*","VI","ARI","Chi","P-value")
  l=list("K-Means"=km,"CVI"=crit)
  return(l)
}
