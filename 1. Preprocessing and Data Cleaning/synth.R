synth=function(df,context){
  library(MASS)
  nam=colnames(df)
  groups=split(df,df[context])
  synth_list=lapply(groups,function(group){
    n=nrow(group)
    logn.col=c(5:8,10,14:15,28)
    poi.cols=c(16:27)
    norm.cols=c(30)
    norm.int.cols=c(31:32)
    beta.cols=c(12:13)
    log.data=log(group[logn.col]+1e-16)
    mu=colMeans(log.data)
    sigma=cov(log.data)
    synth.log=as.data.frame(mvrnorm(n=n,mu=mu,Sigma = sigma))
    synth.log=exp(synth.log)
    while (length(which(synth.log[,4]>synth.log[,3]))!=0){
      w=which(synth.log[,4]>synth.log[,3])
      synth.log[w,3:4]=synth.log[w,c(4,3)]
    }
    synth.poi=apply(group[poi.cols],2,function(col){
      n=length(col)
      new=rpois(n,mean(c(mean(col),var(col))))
      return(new)
    })
    while(length(which(sum(synth.poi[,1:6])>sum(synth.poi[,7:12])))!=0){
      w=which(sum(synth.poi[,1:6])>sum(synth.poi[,7:12]))
      synth.poi[w,]=synth.poi[w,c(7:12,1:6)]
    }
    no=as.vector(group[norm.cols])[[1]]
    synth.nor=rnorm(n,mean(no),sd(no))
    while (length(which(synth.nor<0 | synth.nor>255))!=0){
      w=which(synth.nor<0)
      synth.nor[w]=rnorm(length(w),mean(synth.nor),sd(synth.nor))
      w=which(synth.nor>255)
      synth.nor[w]=rnorm(length(w),mean(synth.nor),sd(synth.nor))
    }
    synth.nor=as.data.frame.vector(synth.nor)
    synth.nor.int=group[norm.int.cols]
    new=round(mvrnorm(n,colMeans(synth.nor.int),cov(synth.nor.int)))
    int.1=as.vector(new[,1])
    int.2=as.vector(new[,2])
    while (length(which(int.1>synth.nor))!=0){
      w=which(int.1>synth.nor)
      int.1[w]=round(rnorm(length(w),mean(int.1),sd(int.1)))
      while (length(which(int.1<0 | int.1>255))!=0){
        w=which(int.1<0)
        int.1[w]=rnorm(length(w),mean(int.1),sd(int.1))
        w=which(int.1>255)
        int.1[w]=rnorm(length(w),mean(int.1),sd(int.1))
      }
    }
    while (length(which(int.2<synth.nor))!=0){
      w=which(int.2<synth.nor)
      int.2[w]=round(rnorm(length(w),mean(int.2),sd(int.2)))
      while (length(which(int.2<0 | int.2>255))!=0){
        w=which(int.2<0)
        int.2[w]=rnorm(length(w),mean(int.2),sd(int.2))
        w=which(int.2>255)
        int.2[w]=rnorm(length(w),mean(int.2),sd(int.2))
      }
    }
    synth.nor.int=as.data.frame(cbind(int.1,int.2))
    synth.beta=apply(group[beta.cols],2,function(col){
      col=pmin(pmax(col, 1e-16), 1 - 1e-16)
      fit=fitdistr(col, "beta", start = list(shape1 = 2, shape2 = 2))
      new =rbeta(n, fit$estimate["shape1"], fit$estimate["shape2"])
      return(new)
    })
    area=synth.beta[,1]*synth.log[,5]
    eccentricity=sqrt(1-(synth.log[,4]/synth.log[,3])^2)
    eq=sqrt(area*4/pi)
    synth=cbind(group[,1:3],area,synth.log[,1:4],eccentricity,synth.log[,5],eq,synth.beta,synth.log[,6:7],synth.poi,synth.log[,8],group[,29],synth.nor,synth.nor.int)
    colnames(synth)=nam
    synth[context]=group[1,context]
    return(synth)
  })
  do.call(rbind,synth_list)
}
