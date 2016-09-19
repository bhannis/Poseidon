#-------------------------------------------------------------------------------------------# 
# convert2pbdb: 
# format sampled poseidon data for use with Alroy's sqs4.3 perl script
#
# input:
# pres is sampled poseidon data from sampling() function
# use unique = T to make species unique in each "collection" (= our cell)
# ave is rate parameter of gamma distribution used to assign collections to "references"
#
# output:
# a file "myDownload.csv" mimicking a PBDB download
# a file "FR2_bins" mimicking a time scale
#
# Modified from version convert2pbdb12.06.2015.R
# Written by Lee Hsiang Liow 2015
# (with contributions from Bjarte Hannisdal)
#
# This file is provided as auxiliary supplementary material to
# Hannisdal B, Haaga KA, Reitan T, Diego D, and Liow LH, 
# "Common Species Link Global Ecosystems to Climate Change".
#-------------------------------------------------------------------------------------------#

convert2pbdb=function(pres, t, unique, ave)
{             
  print("Converting data to PBDB format for SQS 4.3")
  ptm <- proc.time()[3]
  
  if(unique==T){
  pbdbdata=unique(pres$pres.data[[1]])# make unique species in each "collection" (== our cell)
  }else{pbdbdata=pres$pres.data[[1]]}# or not

  pbdbdata$ref=0 # make new column of reference numbers
  temp <- unique(pres$pres.data[[1]]$cell) # cells that are sampled
  k=1 # ticker
  no=1
  while (length(temp)>no){                
  no=ceiling(rgamma(1, 1,ave))# generates a cluster of cells to combine into a reference
  if(length(temp)<no) break() # stop if cells lefts are less than draw then go straight to assigning last numbers
  oneref=sample(temp, no) # samples the cells that will be combined
  temp=setdiff(temp, oneref) # take off the ones already done
  pbdbdata$ref[which(pbdbdata$cell %in% oneref)]=k # assign the reference number
  k=k+1
  }
  pbdbdata$ref[which(pbdbdata$cell %in% temp)]=k# remaining ones
  o=order(pbdbdata$cell) # Alroy's data is organized according to collections, so this has to be ordered!
  pbdbdata=pbdbdata[o,]


  for (i in 2:t){  
    if(unique==T){
      pbdb=unique(pres$pres.data[[i]])# make unique species in each "collection" (== our cell)
    }else{pbdb=pres$pres.data[[i]]}# or not
    
    pbdb$ref=0
    temp <- unique(pres$pres.data[[i]]$cell)
    no=1  
    while (length(temp)>no){
    k=k+1
    no=ceiling(rgamma(1, 1,ave)) #generates a cluster of cells to combine into a reference
    if(length(temp)<no) break()
    oneref=sample(temp, no) # samples the cells that will be combined
    temp=setdiff(temp, oneref)
    pbdb$ref[which(pbdb$cell %in% oneref)]=k # assign the reference number 
    }
    pbdb$ref[which(pbdb$cell %in% temp)]=k
    pbdb$cell=pbdb$cell+i*1e03 # make collection numbers different in each time interval to please sqs4.3

    o=order(pbdb$cell)
    pbdb=pbdb[o,]
    
    pbdbdata=rbind(pbdbdata, pbdb)
  }
                
  colnames(pbdbdata)=c("occurrence.genus_name","collection_no", "FR2_bin", "collection.reference_no")
                
  pbdbdata$source_database="poseidon"
  pbdbdata$collection.authorizer="dummy"
  pbdbdata$occurrence.genus_reso="dummy"
  pbdbdata$original.genus_reso="dummy"
  pbdbdata$original.genus_name="dummy"
  pbdbdata=pbdbdata[,c(2,5,6,7,1,8,9,4,3)]
  rownames(pbdbdata)=NULL
  
  write.csv(pbdbdata,file="myDownload.csv",row.names=F, quote=F) # mimicking a PBDB download
  write(matrix(1:t,ncol=1), file="FR2_bins", ncolumns=1) # mimicking a time scale file
  
  etm=round(proc.time()[3] - ptm)
  print(paste(as.character(etm), " seconds (", as.character(round(etm/60)), " mins.)"))
  print(" ")
}
