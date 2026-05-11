#Epsilon and Min Points selection according to Starczewski et al. 2020

#di is a dist object
#k is a starting parameter, 2 is recommended
#p is the dimensionality of the dataset
#w is a logical parameter to unlock some selection methods for Min Points
#show is a logical parameter to show/hide plots
#auto is a logical parameter: should the selction method for MinPoints be added in the function's call?
#which is the selection method for MinPoints
db.estimate=function(di,k=2,p,w=T,show=T,auto=F,which=NA){
  library(dbscan)
  library(signal)
  n=attr(di,"Size")#Sample size extraction from di
  if (k>n){#MinPoints must be less than n 
    ris=list(Epsilon.1=NA,MinPoints.1=NA,Epsilon.2=NA,MinPoints.2=3)
    return(ris)
  }
  else{
    if (show==T){
      kNNdistplot(di,k)
    }
    d=c()
    for(i in 1:k){
      d=c(d,kNNdist(di,i))
    }
    d=sort(d)#Sorting k-dist across k=1,...,k
    if (show==T){
      plot(d,type="l")
    }
    v.start=n*(k-1)
    v.stop=n*k#Limits of the interval in which is expected to find the knee point
    x1=v.start
    x2=v.stop
    y1=d[x1]
    y2=d[x2]#Values corresponding to v.start and v.stop
    A1=(y2-y1)/(x2-x1)
    B1=y1-A1*x1
    if ((is.finite(B1)[1]==TRUE) & (is.finite(A1)[1]==TRUE) & show==T){
      abline(B1,A1,col="green")
    }#Line passing through p1=(x1,y1) and p2
    A2=-A1
    B2=A1*(x1+x2)+B1
    if ((is.finite(B2)[1]==TRUE) & (is.finite(A2)[1]==TRUE) & show==T){
      abline(B2,A2,col="blue")
    }#Line intercepting halfway between p1 and p2 with opposite slope if compared with previous line. Intercepts with d also, in p3
    xt=c(1:(n*k))
    r.2=xt*A2+c(rep(B2,n*k))
    x3=which.min(abs(d-r.2))
    y3=d[x3]#Determination of p3
    if (n>=500){
      l=31
    }
    if (n>=300 & n<500){
      l=21
    }
    if (n>=100 & n<300){
      l=11
    }
    if (n>=30 & n<100){
      l=7
    }
    if (n<30){
      l=5
    }
    der=sgolayfilt(d,p=3,l,1) #Estimate of a derivative for d
    m=der[x3]
    A3=m
    B3=y3-m*x3
    if ((is.finite(B3)[1]==TRUE) & (is.finite(A3)[1]==TRUE) & show==T){
      abline(B3,A3,col="red")
    }#Line tangent to d in p3
    r.3=A3*xt[x3:x2]+c(rep(B3,(x2-x3+1)))
    M=d[x3:x2]-r.3#Differences between d and the third line between x3 and x2
    ya=mean(M)
    xa=which.min(abs(c(rep(ya,length(M)))-M))#pa point determination
    eps=d[x3+xa]#Epsilon estimate
    if (attr(di,"method")=="euclidean"){
      d1=abs(x2-x3)+abs(y2-y3)
      d2=abs(x1-x3)+abs(y1-y3)
      dp=d1/d2
    } else if (attr(di,"method")=="manhattan"){
      d1=sqrt((x2-x3)^2+(y2-y3)^2)
      d2=sqrt((x1-x3)^2+(y1-y3)^2)
      dp=d1/d2
    }
    if (dp>=4){
      b=(xa+length(M))/2
      eps=d[x3+as.integer(b)]
    }#Estimate adjustment with bias
    if (k!=2){
      ris=list(Epsilon.1=eps)
      return(ris)
    } else{
      if (w==F){#Selection method for MinPoints
        if (auto==T){
          input=which
        }
        else{
          input=readline("Which selection method for MinPoints (MP)? (0: Estimate, 1: MP=3, 2: MP=4, 3: MP=5, 4: MP=2*dim(X), 5: MP=log(n) )")
        }
        if (input==1){
          minp=3
        } else if (input==2){
          minp=4
        } else if (input==3){
          minp=5
        } else if (input==4){
          minp=2*p
        } else if (input==5){
          minp=as.integer(log(attr(di,"Size")))
        } else{
          minp=max(as.integer(dp-.5),3)
        }
        
      }else{
        minp=max(as.integer(dp-.5),3)
      }
      db=db.estimate(di,minp,p)#Recursive call to db.estimate
      ris=list(Epsilon.1=as.numeric(eps),MinPoints.1=as.numeric(k),Epsilon.2=as.numeric(db$Epsilon.1),MinPoints.2=as.numeric(minp))#Epsilon 2 and MinPoints 2 are the correct estimates for Epsilon and MinPoints
      return(ris)
    }
  }
}
