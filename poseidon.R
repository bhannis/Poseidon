#-------------------------------------------------------------------------------------------# 
# Poseidon
# 
# Modified from version poseidon12.06.2015.R
# Written by Lee Hsiang Liow 2015
# (with contributions from Bjarte Hannisdal)
#
# This file is provided as auxiliary supplementary material to
# Hannisdal B, Haaga KA, Reitan T, Diego D, and Liow LH, 
# "Common Species Link Global Ecosystems to Climate Change".
#-------------------------------------------------------------------------------------------#

poseidon=function(Nt, mu, sig, cell, max_sp, min_sp, max_indiv, min_indiv){

print("Generating plankton world")
ptm <- proc.time()[3]

#-------------------------------------------------------------------------------------------#
# Species richness changes through time

x=sin(seq(0,6.5*pi, pi/16)) # sinusoidal
x=(x-min(x))/(max(x)-min(x)) # scale to [0 1]
t=seq(1,length(x),1) # time intervals
sp=x*(max_sp-min_sp)+min_sp
sp=round(sp)     # sp. richness

#  Count loss or gain of species  
spLG=(sp[-1]-sp[1:(length(sp)-1)]) #  

gain=NA # gains from t=2 through t=t
loss=NA  # losses from t=2 through t=t       
for (i in 1:(length(t)-1))
{
  if (spLG[i]>0)
  {
    gain[i]=spLG[i]          
    loss[i]=0 
  } else {
    loss[i]=-spLG[i]
    gain[i]=0
  }
}

# Make matrix of species
Species.list = vector("list", length(t))     

# establish the time slice
for (j in 1:1){
  spid=seq(1:sp[1])     
  splist=as.data.frame(spid)
  splist$timbin=1
  Species.list[[j]] = splist  
}

# subsequent species

for (j in 2:length(t)){
  splist=Species.list[[j-1]]      # splist is the previous cohort
  
  if (loss[j-1]<gain[j-1]){ # this is the gain loop
        ins= seq(max(splist$spid)+1, max(splist$spid)+gain[j-1], 1)  # species IDs to insert
    
        tempA=splist$spid
        spid=c(tempA, ins)
        splist=as.data.frame(spid)
        splist$timebin=j
    
    }else{ # this is the loss loop
        del=sample(splist$spid,loss[j-1], replace=F)    # species IDs to delete
    
        tempA=splist$spid
        spid=setdiff(tempA, del)
        splist=as.data.frame(spid)
        splist$timebin=j
  }
  
  Species.list[[j]] = splist
}

#-------------------------------------------------------------------------------------------#
# Distribute individuals among species according to lognormal RAD

series=sin((seq(3.6*pi+.5,14*pi, pi/27))) # sinusoidal
series=(series-min(series))/(max(series)-min(series)) # scale to [0 1]
indiv=series*(max_indiv-min_indiv)+min_indiv
indiv=indiv[1:Nt]  # varying abundance through time

k.n=NA
Raw.list = vector("list", length(t)) # each row is cell id, and species ID and time interval
PJ=NULL

# pass mu and sig on log scale
plognorm=function(x,sig) plnorm(x, log(mu),log(sig),log = F)
lognorm=function(x,sig) dlnorm(x, log(mu),log(sig),log = F) 

for (j in 1:length(t))   # for each time interval
{
  splist=Species.list[[j]]
  n=round(sp[j]) # the number of species in system
  
  top= 0.99 # area under the curve to reach, this is a decision that we have to make, 0.95 is a good compromise
  a=0 # set the cummulative sum to zero
  z=0 # z is the x axis we need to look at
  
  while (a< top){ # keep adding xaxis until we reach the predetermined curve area
    z=z+1 # add for every loop to xaxis we want to look at
    pp=plognorm(seq(0,z),sig[j])
    a=tail(pp, n=1) # sum up the probability
  }
  
  pp=plognorm(seq(0,z),sig[j]) # this is the part of the curve to be used.probability
  aa=lognorm(seq(0,z),sig[j]) # this is the part of the curve to be used.density
  b=seq(0,z)
  
  # extraction
  # n is z, lengthout is species
  y=seq(1,z, length.out=n) # extract for those species that are in the system from poseidon
  bb=lognorm(y,sig[j])
  A=indiv[j]/sum(bb)
  counts=round(A*bb) 
  
  stick=function(species, counts){
    s=rep(species,counts)
    return(s)}
  
  s=stick(splist$spid[1:length(splist$spid)], ceiling(counts[1:length(splist$spid)]))
  s=as.data.frame(s)
  
  s$cell=sample(seq(1,cell,1),dim(s)[1],replace=T) # put in cells or sites (not all cells equally likely)
  s$time=j
  k.n[j]=length(unique(s$cell))
  
  Raw.list[[j]]=s
  
  # calculate Pielou's J evenness on true data
  newlist=Raw.list[[j]]
  nn=tapply(newlist$s, newlist$s, length)
  freq = nn / sum(nn)
  PJ[j] = (-1 * sum(freq * log(freq)))/log(n)
}

#-------------------------------------------------------------------------------------------#

output<-list(Jumbled.list=Raw.list,t=t,indiv=indiv, sp=sp, k.n=k.n, PJ=PJ)

etm=round(proc.time()[3] - ptm)
print(paste("Elapsed time: ", as.character(etm), " seconds (", as.character(round(etm/60)), " mins.)"))
print(" ")

return(output) 
}