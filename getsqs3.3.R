#-------------------------------------------------------------------------------------------# 
# getsqs3.3: prepare sampled poseidon data and run Alroy's sqs3.3 R script
# use abundance = F for "occurrences" like in PBDB (and most applications)
# use abundance = T for abundances as recommended by Alroy
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

getsqs3.3=function(pres, abundance){
  
  if (abundance==T){      
    print("Estimating richness with SQS version 3.3 - using abundance")
  }else{
    print("Estimating richness with SQS version 3.3 - using occurrence")}
  ptm <- proc.time()[3]
  
  SQSres=NULL # prep result vector for SQS using collections
  rawrows=NULL # result vector for rows of data
  specimens=NULL # direct output from sqs3.3
  H=NULL # Shannon's H
  
  # parameter settings for SQS
  quo=0.6 # quorum level (0.6 is "typical")
  dominant="exclude" # Alroy's recommendation
  trials=Ntrial # Alroy's recommendation is 1000
  ignore.singletons=F 
  t=105
  
  cell.n=NULL
  x.n=vector("list", length(t))   
  for ( i in 1:t){
    cell.n[i]=length(unique(pres$pres.data[[i]]$cell))# number of cells sampled
    if (abundance==F){
      x.n[[i]] <- with(pres$pres.data[[i]], tapply(cell,s, function(i) length(unique(i))))
    }else{
      x.n[[i]] <- with(pres$pres.data[[i]], tapply(cell,s, function(i) length(i)))
    }
  }
  
  for (i in 1:t){
    if (length(x.n[[i]])>0){ # sometimes there are no species in the list e.g. in the case of increasingly poor preservation back in time
      x=x.n[[i]][which(x.n[[i]]>0)] # this is to remove the NA's that have been created by removing rare species
      temp=sqs(ab=x, q=quo, trials=trials, ignore.singletons=ignore.singletons, dominant=dominant) 
      SQSres[i]=temp[3]    # temp [3] is sqs 
      rawrows[i]=temp[16]
      specimens[i]=temp[12]
      H[i]=temp[9]
    }else{
      SQSres[i]=NA 
      rawrows[i]=temp[16] # same as specimens, just to check
      specimens[i]=temp[12]
      H[i]=temp[9]
    } 
  } 
  output<-list(SQSres=SQSres, rawrows=rawrows, specimens=specimens, H=H)
  # possible output from SQS 3.3:
  # [1] raw richness 
  # [2] Good's u 
  # [3] subsampled richness 
  # [4] subsampled u 
  # [5] Chao 1 
  # [6] subsampled Chao 1 
  # [7] k
  # [8] Fisher's alpha 
  # [9] Shannon's H 
  # [10] Hurlbert's PIE 
  # [11] dominance 
  # [12] specimens 
  # [13] singletons 
  # [14] doubletons 
  # [15] specimens drawn 
  # [16] datarows 
  
  etm=round(proc.time()[3] - ptm)
  print(paste(as.character(etm), " seconds (", as.character(round(etm/60)), " mins.)"))
  print(" ")
  
  return(output) 
}
