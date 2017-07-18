#-------------------------------------------------------------------------------------------# 
# Metrics 
# Calculates SCOR and richness/diversity measures 
# from the sampled plankton world
# 
# Written by Lee Hsiang Liow 2015
# with contributions from Bjarte Hannisdal
#
# This file is provided as part of the electronic supplementary material to
# Hannisdal B, Haaga KA, Reitan T, Diego D, and Liow LH, 2017,
# Common species link global ecosystems to climate change: dynamical evidence
# in the planktonic fossil record. Proc. R. Soc. B 284: 2010722
# http://dx.doi.org/10.1098/rspb.2017.0722
#-------------------------------------------------------------------------------------------#

rawmetrics<- function(pres, t){
  # this function returns the total number of sampled species and range-through counts
  
  print("Calculating raw richness and range-through")
  #ptm <- proc.time()[3]
  
  rawsp.n= NULL
  rangethru=NULL
  sp.occ=vector("list", length(t))
  for (i in 1:t){
    rawsp.n[i]=length(unique(pres$pres.data[[i]]$s))# number of raw species sampled
    temp <- unique(pres$pres.data[[i]]$s) # unique species list in each time interval for RT
    sp.occ[[i]]<-cbind(rep(i, length(temp)),  temp)
  }  

  all.sp.occ=as.data.frame(do.call(rbind, sp.occ))
  colnames(all.sp.occ)=c("time", "sp")
  sp=unique(all.sp.occ$s)
  sp=as.data.frame(sp)
  
  for ( i in 1:(dim(sp)[1]) ){
    sp$minage[i]=min(all.sp.occ$time[which(all.sp.occ$s==sp$sp[i])])
    sp$maxage[i]=max(all.sp.occ$time[which(all.sp.occ$s==sp$sp[i])])
  }
  for (i in 1:t){
    rangethru[i]=length(which(sp$minage<=i & sp$maxage>=i))
  }

  output<-list(rawsp=rawsp.n, RT=rangethru)
  
  #etm=round(proc.time()[3] - ptm)
  #print(paste(as.character(etm), " seconds (", as.character(round(etm/60)), " mins.)"))
  print(" ")
  
  return(output) 
}

#-------------------------------------------------------------------------------------------#

SCOR<- function(pres, t){
  # this function returns the SCOR estimates and std. errors as well as the number of cells sampled
  
  print("Calculating SCOR")
  #ptm <- proc.time()[3]
  
  cell.n=NULL
  sl.n=NULL
  sesl.n=NULL
  sl2.n=NULL
  sesl2.n=NULL

  for (i in 1:t){
    cell.n[i]=length(unique(pres$pres.data[[i]]$cell))# number of cells sampled
    x.n <- with(pres$pres.data[[i]], tapply(cell,s, function(i) length(unique(i))))# number of sites observed occupied for given species
    p.n <- x.n/cell.n[i]
    sl.n[i] <-  sum(-log(1-p.n)) # this is the relative density
    sesl.n[i] <- sqrt(sum((p.n/(1-p.n))/(cell.n[i])))      # this is SE

  }
  output<-list( scor=sl.n, scorE=sesl.n)
  
  #etm=round(proc.time()[3] - ptm)
  #print(paste(as.character(etm), " seconds (", as.character(round(etm/60)), " mins.)"))
  print(" ")
  
  return(output) 
}

#-------------------------------------------------------------------------------------------#

RARE<- function(pres, t, repeats){
  # this function returns the classical rarefaction richness estimates
  
  print("Estimating richness using classical rarefaction")
  ptm <- proc.time()[3]
  
  samples=NULL
  for ( i in 1:t){
    temp <- unique(pres$pres.data[[i]])
    samples[i]=dim(temp)[1]
  }
  lower=min(samples) # this is the number of unique species-cell in the worst time interval

  CR.sp=matrix(NA, nrow=t,ncol=repeats )
  CR.est=NULL

    for(j in 1:repeats){
        for ( i in 1:t){
          temp <- unique(pres$pres.data[[i]])
          random = sample(seq(1:dim(temp)[1]), lower, replace=F)
          # see http://onlinelibrary.wiley.com/enhanced/doi/10.1046/j.1461-0248.2001.00230.x/
          subsample=pres$pres.data[[i]][random,]
          CR.sp[i, j]=length(unique(subsample$s))
        }
    }

  CR.est=apply(CR.sp, 1, mean)  
  output<-list( CR.data=CR.sp, CR.est=CR.est)
  
  etm=round(proc.time()[3] - ptm)
  print(paste(as.character(etm), " seconds (", as.character(round(etm/60)), " mins.)"))
  print(" ")
  
  return(output)
}

#-------------------------------------------------------------------------------------------#

O2W<- function(pres, t, repeats){
  # this function returns the occurrences-squared-weighted richness estimates
  
  print("Estimating richness using O2W subsampling")
  ptm <- proc.time()[3]
  
  O2W.sp=matrix(NA, nrow=t,ncol=repeats )
  O2W.est=NULL
  
  samples=NULL
  for ( i in 1:t){
    temp <- unique(pres$pres.data[[i]])
    counts.in.list=with(temp, tapply(s, cell, length))
    samples[i]=sum(counts.in.list^2)
  }
  lower=min(samples) # this is the squared number to be reached
  
  for ( j in 1:repeats){
    for ( i in 1:t){
     data.temp <- unique(pres$pres.data[[i]]$cell) # unique samples (cells)
      occ=0
      sp=NA
    
        while(occ<lower){ # while still below the squared occurrence target
        random = sample(data.temp, 1, replace=F) # cell being sampled
        data.temp=setdiff(data.temp,random) # take that cell away and rename the data
        subsample=pres$pres.data[[i]][which(pres$pres.data[[i]]$cell==random),]# the sampled cell with its species
        sp=append(sp,unique(subsample$s)) # compress the list into PBDB "list"
        occ=occ+dim(subsample)[1]^2 # compute the occurrence squared             
        }
     
    sps=unique(sp)
    O2W.sp[i, j]=length(sps)-1
    }
  }
  O2W.est=apply(O2W.sp, 1, mean)
  output<-list( O2W.data=O2W.sp, O2W.est=O2W.est)
  
  etm=round(proc.time()[3] - ptm)
  print(paste(as.character(etm), " seconds (", as.character(round(etm/60)), " mins.)"))
  print(" ")
  
  return(output)
}
