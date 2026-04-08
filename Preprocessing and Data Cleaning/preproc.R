#PARAMETERS
#data: observations' data.frame.
#Image: vector listing all of the Image's modalities.
#Animal: vector listing all of the animal species corresponding to each Image modality.
#ID: vector listing all of the subjects corresponding to each Image modality.
#microglia: should the modality "Microglia" be encluded in Shape? (Cellls with MajorAxisLength in [5,7) micron -> Microglia).
#cut: should problematic observations be deleted?
#threshold: vector listing the minimum and maximum thersholds to apply on the data. If only one threshold should be considered, mark the other one as NA.
#Example: type 1 -> threshold=c(.05,.95) or c(.05,NA) or c(NA,.95); type=2 -> threshold=c(34,200) or c(34,NA) or c(NA,200).
#type=3 doesn't require the specification of the threshold parameter.
#by: column of data on which the threshold(s) is going to be applied.
#type: 1=thresholding using quantiles, 2=thresholding using specific variable's values, 3=5 micron threshold on MajorAxisLength.
#tol= tolerance between calculated values and registered values.

preproc=function(data,Image,Animal,ID,microglia=F,cut=F,threshold=NA,by=NA,type=NA,tol=.01){
  #INPUT CONTROL
  data=as.data.frame(data)
  ima=as.vector(Image)
  ani=as.vector(Animal)
  id=as.vector(ID)
  mic=as.logical(microglia)
  cut=as.logical(cut)
  thr=as.vector(threshold)
  by=as.character(by)
  type=as.integer(type)
  tol=as.numeric(tol)
  if (is.na(type)==F){
    if ((type %in% c(1,2,3))==F){
      warning("type must be equal to 1, 2 o 3!")
    }
  }
  
  #InvAR, Convex Circularity & Shape
  InvAR=data$MinorAxisLength/data$MajorAxisLength
  ConvexCircularity=(4*pi*data$ConvexArea)/(data$ConvexPerimeter^2)
  Shape=(InvAR>=.75)+(ConvexCircularity<=.85)*2
  if(microglia==T){
    w=which(data$MajorAxisLength<7 & data$MajorAxisLength>=5)
    Shape[w]=4
    Shape=factor(Shape)
    levels(Shape)=c("Ellipsoid","Round","Pyramidal","Complex","Microglia")
  }else{
    Shape=factor(Shape)
    levels(Shape)=c("Ellipsoid","Round","Pyramidal","Complex")
  }
  data=cbind(data,InvAR,ConvexCircularity,Shape)
  
  #ID,Animal
  ID=factor(data$Image)
  levels(ID)=id
  Animal=factor(data$Image)
  levels(Animal)=ani
  data=cbind(data,ID,Animal)
  
  #Ngb 50 e Ngb 100
  Ngb_50=data$Dir_50_1+data$Dir_50_2+data$Dir_50_3+data$Dir_50_4+data$Dir_50_5+data$Dir_50_6
  Ngb_100=data$Dir_100_1+data$Dir_100_2+data$Dir_100_3+data$Dir_100_4+data$Dir_100_5+data$Dir_100_6
  data=cbind(data,Ngb_50,Ngb_100)
  
  #VARIABLES CONTROL
  problems=data.frame(q=NA,var=NA,prob=NA,out=NA)#Incoherences report
  which_obs=list()#List of the problematic observations segregated by cause
  i=0 
  data_cln=data
  
  #Image, Layer, Ngb 50, Ngb 100, Min Intensity, Max Intensity must be integers 
  a=any(round(data$Image)!=data$Image)#Are there problematic observatios?
  if (a==T){
    w=which(round(data$Image)!=data$Image)#Which are they?
    i=i+1
    problems[i,1:3]=c(length(w),"Image","Not an integer")#Updating report
    if (cut==T){
      if (any(round(data_cln$Image)!=data_cln$Image)){
        w1=which(round(data_cln$Image)!=data_cln$Image)
        data_cln=data_cln[-w1,]#Deleting problematic observations
        problems$out[i]=length(w1)
      }
    }
    if (cut==F){
      which_obs[[i]]=w#No cut
    }
  }
  a=any(round(data$Layer)!=data$Layer)
  if (a==T){
    w=which(round(data$Layer)!=data$Layer)
    i=i+1
    problems[i,1:3]=c(length(w),"Layer","Not an integer")
    if (cut==T){
      if (any(round(data_cln$Layer)!=data_cln$Layer)){
        w1=which(round(data_cln$Layer)!=data_cln$Layer)
        data_cln=data_cln[-w1,]
        problems$out[i]=length(w1)
      }
    }
    if (cut==F){
      which_obs[[i]]=w
    }
  }
  a=any(round(data$Ngb_50)!=data$Ngb_50)
  if (a==T){
    w=which(round(data$Ngb_50)!=data$Ngb_50)
    i=i+1
    problems[i,1:3]=c(length(w),"Ngb_50","Not an integer")
    if (cut==T){
      if (any(round(data_cln$Ngb_50)!=data_cln$Ngb_50)){
        w1=which(round(data_cln$Ngb_50)!=data_cln$Ngb_50)
        data_cln=data_cln[-w1,]
        problems$out[i]=length(w1)
      }
    }
    if (cut==F){
      which_obs[[i]]=w
    }
  }
  a=any(round(data$Ngb_100)!=data$Ngb_100)
  if (a==T){
    w=which(round(data$Ngb_100)!=data$Ngb_100)
    i=i+1
    problems[i,1:3]=c(length(w),"Ngb_100","Not an integer")
    if (cut==T){
      if (any(round(data_cln$Ngb_100)!=data_cln$Ngb_100)){
        w1=which(round(data_cln$Ngb_100)!=data_cln$Ngb_100)
        data_cln=data_cln[-w1,]
        problems$out[i]=length(w1)
      }  
    }
    if (cut==F){
      which_obs[[i]]=w
    }
  }
  a=any(round(data$Min_Intensity)!=data$Min_Intensity)
  if (a==T){
    w=which(round(data$Min_Intensity)!=data$Min_Intensity)
    i=i+1
    problems[i,1:3]=c(length(w),"Min_Intensity","Not an integer")
    if (cut==T){
      if (any(round(data_cln$Min_Intensity)!=data_cln$Min_Intensity)){
        w1=which(round(data_cln$Min_Intensity)!=data_cln$Min_Intensity)
        data_cln=data_cln[-w1,]
        problems$out[i]=length(w1)
      }    
    }
    if (cut==F){
      which_obs[[i]]=w
    }
  }
  a=any(round(data$Max_Intensity)!=data$Max_Intensity)
  if (a==T){
    w=which(round(data$Max_Intensity)!=data$Max_Intensity)
    i=i+1
    problems[i,1:3]=c(length(w),"Max_Intensity","Not an integer")
    if (cut==T){
      if (any(round(data_cln$Max_Intensity)!=data_cln$Max_Intensity)){
        w1=which(round(data_cln$Max_Intensity)!=data_cln$Max_Intensity)
        data_cln=data_cln[-w1,]
        problems$out[i]=length(w1)
      }
    }
    if (cut==F){
      which_obs[[i]]=w
    }
  }
  
  #Layer must be an integer between 1 and 6
  a=any((data$Layer %in% c(1:6))==F)
  if (a==T){
    w=which((data$Layer %in% c(1:6))==F)
    i=i+1
    problems[i,1:3]=c(length(w),"Layer","Not an actual Layer")
    if (cut==T){
      if (any((data_cln$Layer %in% c(1:6))==F)){
        w1=which((data_cln$Layer %in% c(1:6))==F)
        data_cln=data_cln[-w1,]
        problems$out[i]=length(w1)
      }
    }
    if (cut==F){
      which_obs[[i]]=w
    }
  }
  
  #Area, MajorAxisLength, MinorAxisLength, ConvexArea, EquivDiameter, Perimeter e ConvexPerimeter must be >0
  a=any(data$Area<=0)
  if (a==T){
    w=which(data$Area<=0)
    i=i+1
    problems[i,1:3]=c(length(w),"Area","Less than or equal to zero")
    if (cut==T){
      if (any(data_cln$Area<=0)){
        w1=which(data_cln$Area<=0)
        data_cln=data_cln[-w1,]
        problems$out[i]=length(w1)
      }
    }
    if (cut==F){
      which_obs[[i]]=w
    }
  }
  a=any(data$MajorAxisLength<=0)
  if (a==T){
    w=which(data$MajorAxisLength<=0)
    i=i+1
    problems[i,1:3]=c(length(w),"MajorAxisLength","Less than or equal to zero")
    if (cut==T){
      if (any(data_cln$MajorAxisLength<=0)){
        w1=which(data_cln$MajorAxisLength<=0)
        data_cln=data_cln[-w1,]
        problems$out[i]=length(w1)
      }
    }
    if (cut==F){
      which_obs[[i]]=w
    }
  }
  a=any(data$MinorAxisLength<=0)
  if (a==T){
    w=which(data$MinorAxisLength<=0)
    i=i+1
    problems[i,1:3]=c(length(w),"MinorAxisLength","Less than or equal to zero")
    if (cut==T){
      if (any(data_cln$MinorAxisLength<=0)){
        w1=which(data_cln$MinorAxisLength<=0)
        data_cln=data_cln[-w1,]
        problems$out[i]=length(w1)
      }
    }
    if (cut==F){
      which_obs[[i]]=w
    }
  }
  a=any(data$ConvexArea<=0)
  if (a==T){
    w=which(data$ConvexArea<=0)
    i=i+1
    problems[i,1:3]=c(length(w),"ConvexArea","Less than or equal to zero")
    if (cut==T){
      if (any(data_cln$ConvexArea<=0)){
        w1=which(data_cln$ConvexArea<=0)
        data_cln=data_cln[-w1,]
        problems$out[i]=length(w1)
      }
    }
    if (cut==F){
      which_obs[[i]]=w
    }
  }
  a=any(data$EquivDiameter<=0)
  if (a==T){
    w=which(data$EquivDiameter<=0)
    i=i+1
    problems[i,1:3]=c(length(w),"EquivDiameter","Less than or equal to zero")
    if (cut==T){
      if (any(data_cln$EquivDiameter<=0)){
        w1=which(data_cln$EquivDiameter<=0)
        data_cln=data_cln[-w1,]
        problems$out[i]=length(w1)
      }
    }
    if (cut==F){
      which_obs[[i]]=w
    }
  }
  a=any(data$Perimeter<=0)
  if (a==T){
    w=which(data$Perimeter<=0)
    i=i+1
    problems[i,1:3]=c(length(w),"Perimeter","Less than or equal to zero")
    if (cut==T){
      if (any(data_cln$Perimeter<=0)){
        w1=which(data_cln$Perimeter<=0)
        data_cln=data_cln[-w1,]
        problems$out[i]=length(w1)
      }
    }
    if (cut==F){
      which_obs[[i]]=w
    }
  }
  a=any(data$ConvexPerimeter<=0)
  if (a==T){
    w=which(data$ConvexPerimeter<=0)
    i=i+1
    problems[i,1:3]=c(length(w),"ConvexPerimeter","Less than or equal to zero")
    if (cut==T){
      if (any(data_cln$ConvexPerimeter<=0)){
        w1=which(data_cln$ConvexPerimeter<=0)
        data_cln=data_cln[-w1,]
        problems$out[i]=length(w1)
      }
    }
    if (cut==F){
      which_obs[[i]]=w
    }
  }
  
  #MajorAxisLength must be more than or equal to MinorAxisLength
  #Area must be less than or equal to ConvexArea
  #Ngb 50 must be less than or equal to Ngb 100
  #Min/Mean/Max Intensity must be increasing
  #Solidity must be more than or equal to Extent
  a=any(data$MajorAxisLength<data$MinorAxisLength)
  if (a==T){
    w=which(data$MajorAxisLength<data$MinorAxisLength)
    i=i+1
    problems[i,1:3]=c(length(w),"MajorAxisLength","Less than MinorAxisLength")
    if (cut==T){
      if (any(data_cln$MajorAxisLength<data_cln$MinorAxisLength)){
        w1=which(data_cln$MajorAxisLength<data_cln$MinorAxisLength)
        data_cln=data_cln[-w1,]
        problems$out[i]=length(w1)
      }
    }
    if (cut==F){
      which_obs[[i]]=w
    }
  }
  a=any(data$ConvexArea<data$Area)
  if (a==T){
    w=which(data$ConvexArea<data$Area)
    i=i+1
    problems[i,1:3]=c(length(w),"ConvexArea","Less than Area")
    if (cut==T){
      if (any(data_cln$ConvexArea<data_cln$Area)){
        w1=which(data_cln$ConvexArea<data_cln$Area)
        data_cln=data_cln[-w1,]
        problems$out[i]=length(w1)
      }
    }
    if (cut==F){
      which_obs[[i]]=w
    }
  }
  a=any(data$Ngb_100<data$Ngb_50)
  if (a==T){
    w=which(data$Ngb_100<data$Ngb_50)
    i=i+1
    problems[i,1:3]=c(length(w),"Ngb_100","Less than Ngb_50")
    if (cut==T){
      if (any(data_cln$Ngb_100<data_cln$Ngb_50)){
        w1=which(data_cln$Ngb_100<data_cln$Ngb_50)
        data_cln=data_cln[-w1,]
        problems$out[i]=length(w1)
      }
    }
    if (cut==F){
      which_obs[[i]]=w
    }
  }
  a=any(data$Mean_Intensity<data$Min_Intensity)
  if (a==T){
    w=which(data$Mean_Intensity<data$Min_Intensity)
    i=i+1
    problems[i,1:3]=c(length(w),"Mean_Intensity","Less than Min_Intensity")
    if (cut==T){
      if (any(data_cln$Mean_Intensity<data_cln$Min_Intensity)){
        w1=which(data_cln$Mean_Intensity<data_cln$Min_Intensity)
        data_cln=data_cln[-w1,]
        problems$out[i]=length(w1)
      }
    }
    if (cut==F){
      which_obs[[i]]=w
    }
  }
  a=any(data$Max_Intensity<data$Min_Intensity)
  if (a==T){
    w=which(data$Max_Intensity<data$Min_Intensity)
    i=i+1
    problems[i,1:3]=c(length(w),"Max_Intensity","Less than Min_Intensity")
    if (cut==T){
      if (any(data_cln$Max_Intensity<data_cln$Min_Intensity)){
        w1=which(data_cln$Max_Intensity<data_cln$Min_Intensity)
        data_cln=data_cln[-w1,]
        problems$out[i]=length(w1)
      }
    }
    if (cut==F){
      which_obs[[i]]=w
    }
  }
  a=any(data$Max_Intensity<data$Mean_Intensity)
  if (a==T){
    w=which(data$Max_Intensity<data$Mean_Intensity)
    i=i+1
    problems[i,1:3]=c(length(w),"Max_Intensity","Less than Mean_Intensity")
    if (cut==T){
      if (any(data_cln$Max_Intensity<data_cln$Mean_Intensity)){
        w1=which(data_cln$Max_Intensity<data_cln$Mean_Intensity)
        data_cln=data_cln[-w1,]
        problems$out[i]=length(w1)
      }
    }
    if (cut==F){
      which_obs[[i]]=w
    }
  }
  a=any(data$Solidity<data$Extent)
  if (a==T){
    w=which(data$Solidity<data$Extent)
    i=i+1
    problems[i,1:3]=c(length(w),"Solidity","Less than Extent")
    if (cut==T){
      if (any(data_cln$Solidity<data_cln$Extent)){
        w1=which(data_cln$Solidity<data_cln$Extent)
        data_cln=data_cln[-w1,]
        problems$out[i]=length(w1)
      }
    }
    if (cut==F){
      which_obs[[i]]=w
    }
  }
  
  #Eccentricity ed Extent devono stare tra zero ed uno
  a=any(data$Eccentricity<0 | data$Eccentricity>1)
  if (a==T){
    w=which(data$Eccentricity<0 | data$Eccentricity>1)
    i=i+1
    problems[i,1:3]=c(length(w),"Eccentricity","Out of domain")
    if (cut==T){
      if (any(data_cln$Eccentricity<0 | data_cln$Eccentricity>1)){
        w1=which(data_cln$Eccentricity<0 | data_cln$Eccentricity>1)
        data_cln=data_cln[-w1,]
        problems$out[i]=length(w1)
      }
    }
    if (cut==F){
      which_obs[[i]]=w
    }
  }
  a=any(data$Extent<0 | data$Extent>1)
  if (a==T){
    w=which(data$Extent<0 | data$Extent>1)
    i=i+1
    problems[i,1:3]=c(length(w),"Extent","Out of domain")
    if (cut==T){
      if (any(data_cln$Extent<0 | data_cln$Extent>1)){
        w1=which(data_cln$Extent<0 | data_cln$Extent>1)
        data_cln=data_cln[-w1,]
        problems$out[i]=length(w1)
      }
    }
    if (cut==F){
      which_obs[[i]]=w
    }
  }
  a=any(data$ConvexCircularity<0 | data$ConvexCircularity>1)
  if (a==T){
    w=which(data$ConvexCircularity<0 | data$ConvexCircularity>1)
    i=i+1
    problems[i,1:3]=c(length(w),"ConvexCircularity","Out of domain")
    if (cut==T){
      if (any(data_cln$ConvexCircularity<0 | data_cln$ConvexCircularity>1)){
        w1=which(data_cln$ConvexCircularity<0 | data_cln$ConvexCircularity>1)
        data_cln=data_cln[-w1,]
        problems$out[i]=length(w1)
      }
    }
    if (cut==F){
      which_obs[[i]]=w
    }
  }
  
  #Eccentricity, EquivDiameter e Solidity must approximately follow their own formulas
  a=any(data$Eccentricity<sqrt(1-(data$MinorAxisLength/data$MajorAxisLength)^2)-tol | data$Eccentricity>sqrt(1-(data$MinorAxisLength/data$MajorAxisLength)^2)+tol)
  if (is.na(a)==T){
    a=F
  }
  if (a==T){
    w=which(data$Eccentricity<sqrt(1-(data$MinorAxisLength/data$MajorAxisLength)^2)-tol | data$Eccentricity>sqrt(1-(data$MinorAxisLength/data$MajorAxisLength)^2)+tol)
    i=i+1
    problems[i,1:3]=c(length(w),"Eccentricity","Inconsistency between formula and observed value")
    if (cut==T){
      if (any(data_cln$Eccentricity<sqrt(1-(data_cln$MinorAxisLength/data_cln$MajorAxisLength)^2)-tol | data_cln$Eccentricity>sqrt(1-(data_cln$MinorAxisLength/data_cln$MajorAxisLength)^2)+tol)==T){
        w1=which(data_cln$Eccentricity<sqrt(1-(data_cln$MinorAxisLength/data_cln$MajorAxisLength)^2)-tol | data_cln$Eccentricity>sqrt(1-(data_cln$MinorAxisLength/data_cln$MajorAxisLength)^2)+tol)
        data_cln=data_cln[-w1,]
        problems$out[i]=length(w1)
      }
    }
    if (cut==F){
      which_obs[[i]]=w
    }
  }
  a=any(data$Area<pi*(data$EquivDiameter/2)^2-tol | data$Area>pi*(data$EquivDiameter/2)^2+tol)
  if (is.na(a)==T){
    a=F
  }
  if (a==T){
    w=which(data$Area<pi*(data$EquivDiameter/2)^2-tol | data$Area>pi*(data$EquivDiameter/2)^2+tol)
    i=i+1
    problems[i,1:3]=c(length(w),"EquivDiameter","Inconsistency between formula and observed value")
    if (cut==T){
      if (any(data_cln$Area<pi*(data_cln$EquivDiameter/2)^2-tol | data_cln$Area>pi*(data_cln$EquivDiameter/2)^2+tol)==T){
        w1=which(data_cln$Area<pi*(data_cln$EquivDiameter/2)^2-tol | data_cln$Area>pi*(data_cln$EquivDiameter/2)^2+tol)
        data_cln=data_cln[-w1,]
        problems$out[i]=length(w1)
      }
    }
    if (cut==F){
      which_obs[[i]]=w
    }
  }
  if (is.na(a)==T){
    a=F
  }
  a=any(data$Solidity<data$Area/data$ConvexArea-tol | data$Solidity>data$Area/data$ConvexArea+tol)
  if (a==T){
    w=which(data$Solidity<data$Area/data$ConvexArea-tol | data$Solidity>data$Area/data$ConvexArea+tol)
    i=i+1
    problems[i,1:3]=c(length(w),"Solidity","Inconsistency between formula and observed value")
    if (cut==T){
      if (any(data_cln$Solidity<data_cln$Area/data_cln$ConvexArea-tol | data_cln$Solidity>data_cln$Area/data_cln$ConvexArea+tol)==T){
        w1=which(data_cln$Solidity<data_cln$Area/data_cln$ConvexArea-tol | data_cln$Solidity>data_cln$Area/data_cln$ConvexArea+tol)
        data_cln=data_cln[-w1,]
        problems$out[i]=length(w1)
      }
    }
    if (cut==F){
      which_obs[[i]]=w
    }
  }
  
  #THRESHOLDING
  if (any(is.na(threshold)==F)){
    if (type==1){
      if (is.na(thr[1])==F){
        a=any(data[[by]]<quantile(data_cln[[by]],thr[1]))
        if (a==T){
          w=which(data[[by]]<quantile(data_cln[[by]],thr[1]))
          i=i+1
          problems[i,1:3]=c(length(w),by,"Smaller than threshold")
          if (any(data_cln[[by]]<quantile(data_cln[[by]],thr[1]))==T){
            w1=which(data_cln[[by]]<quantile(data_cln[[by]],thr[1]))
            data_cln=data_cln[-w1,]
            problems$out[i]=length(w1)
          }
        }
      }
      if (is.na(thr[2])==F){
        a=any(data[[by]]>quantile(data_cln[[by]],thr[2]))
        if (a==T){
          w=which(data[[by]]>quantile(data_cln[[by]],thr[2]))
          i=i+1
          problems[i,1:3]=c(length(w),by,"Bigger than threshold")
          if (any(data_cln[[by]]>quantile(data_cln[[by]],thr[2]))==T){
            w1=which(data_cln[[by]]>quantile(data_cln[[by]],thr[2]))
            data_cln=data_cln[-w1,]
            problems$out[i]=length(w1)
          }
        }
      }
    }
    if (type==2){
      if (is.na(thr[1])==F){
        a=any(data[[by]]<thr[1])
        if (a==T){
          w=which(data[[by]]<thr[1])
          i=i+1
          problems[i,1:3]=c(length(w),by,"Smaller than threshold")
          if (any(data_cln[[by]]<thr[1])==T){
            w1=which(data_cln[[by]]<thr[1])
            data_cln=data_cln[-w1,]
            problems$out[i]=length(w1)
          }
        }
      }
      if (is.na(thr[2])==F){
        a=any(data[[by]]>thr[2])
        if (a==T){
          w=which(data[[by]]>thr[2])
          i=i+1
          problems[i,1:3]=c(length(w),by,"Bigger than threshold")
          if (any(data_cln[[by]]>thr[2])==T){
            w1=which(data_cln[[by]]>thr[2])
            data_cln=data_cln[-w1,]
            problems$out[i]=length(w1) 
          }
        }
      }
    }
  }
  if (is.na(type)==F){
    if (type==3){
      a=any(data$MajorAxisLength<5)
      if (a==T){
        w=which(data$MajorAxisLength<5)
        i=i+1
        problems[i,1:3]=c(length(w),"MajorAxisLength","Smaller than threshold")
        if (any(data_cln$MajorAxisLength<5)==T){
          w1=which(data_cln$MajorAxisLength<5)
          data_cln=data_cln[-w1,]
          problems$out[i]=length(w1)
        }
      }
    }
  }
  
  Data=data.frame(Shape=data_cln$Shape,Layer=data_cln$Layer,Animal=data_cln$Animal,ID=data_cln$ID,Image=data_cln$Image,Centroid_1=data_cln$Centroid_1,Centroid_2=data_cln$Centroid_2,Area=data_cln$Area,ConvexArea=data_cln$ConvexArea,MajorAxisLength=data_cln$MajorAxisLength,MinorAxisLength=data_cln$MinorAxisLength,Perimeter=data_cln$Perimeter,ConvexPerimeter=data_cln$ConvexPerimeter,EquivDiameter=data_cln$EquivDiameter,Eccentricity=data_cln$Eccentricity,Solidity=data_cln$Solidity,Extent=data_cln$Extent,InvAR=data_cln$InvAR,ConvexCircularity=data_cln$ConvexCircularity,Ngb_50=data_cln$Ngb_50,Ngb_100=data_cln$Ngb_100,Min_Intensity=data_cln$Min_Intensity,Mean_Intensity=data_cln$Mean_Intensity,Max_Intensity=data_cln$Max_Intensity,Subject=data_cln$Subject,Distance_from_outer_layer=data_cln$Distance_from_outer_layer,Cell_class=data_cln$Cell_class)
  Warnings=problems
  colnames(Warnings)=c("N","Variable","Problem","Out")
  Which_Obs=which_obs
  if (cut==T){
    r=list(Data=Data,Warnings=Warnings)
    return(r)
  }
  else{
    r=list(Data=Data,Warnings=Warnings,Which_Obs=Which_Obs)
    return(r)
  }
}
