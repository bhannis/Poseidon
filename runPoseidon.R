#-------------------------------------------------------------------------------------------#
# Set up and run Poseidon simulations
# Written by Lee Hsiang Liow 2015
# (with contributions from Bjarte Hannisdal)
#
# In order to run shareholder quorum subsampling
# the user needs to obtain Alroy's SQS scripts (see below)
#
# This file is provided as auxiliary supplementary material to
# Hannisdal B, Haaga KA, Reitan T, Diego D, and Liow LH, 
# "Common Species Link Global Ecosystems to Climate Change".
#-------------------------------------------------------------------------------------------#

# set the poseidon scaling parameters
Ncell=1000 # number of spatial cells (sites)
max_sp=400 # maximum richness (number of species)
min_sp=200 # minimum richness
maxabund=3e04 # maximum abundance (number of individuals)
minabund=1e04 # minimum abundance
Nt=105 # number of time steps

# set the shape of the rank-abundance distribution (RAD)
mu=3 # locality parameter
sigmin=4.0 # 
sigmax=4.0 # if equal to sigmin, then constant RAD shape
sig=runif(Nt,sigmin,sigmax) # shape parameter, max variability when sigmin=2 and sigmax=6

# set the spatial sampling parameters
sitemin=0.10 # proportion of cells sampled (1 = all cells sampled)
sitemax=0.40 # sampling at t=end. If > sitemin, then linear increase through time
sitevar=0.0 # stddev of gaussian noise on siteprop (= spatial sampling variability, max=0.4)
siteprop=seq(sitemin,sitemax,length=Nt)+rnorm(Nt, 0, sitevar) 
siteprop[which(siteprop<0.02)]=0.02 # impose a minimum of 20 sites
siteprop[which(siteprop>1)]=1 # impose a maximum of Ncell sites

# set proportion of species randomly lost (e.g. dissolution, or biostrat whim)
randrm=F # toggle on or off
randmin=0.5 # 
randmax=0.5 # if equal to randmin, then constant proportion lost
randprop=runif(Nt,randmin,randmax) # max variability when randmin=0.3 and randmax=0.7

# set the number of trials/iterations for all subsampling methods
Ntrial=100 # Note: this must be set separately in SQS 4.3 perl script!

#-------------------------------------------------------------------------------------------#

# The poseidon() function generates the "plankton world"
source("poseidon.R") 
realworld1=poseidon(Nt,mu,sig,Ncell,max_sp,min_sp,maxabund,minabund) #
 
# The sampling() function generates sampled "fossil" world from realworld1
source("sampling.R")
pres=sampling(Jumbled.list=realworld1$Jumbled.list, 
                  t=length(realworld1$t),  
                  siteprop,
                  randrm,
                  randprop) 

#-------------------------------------------------------------------------------------------#
# The metrics.R source file includes routines for calculating
# SCOR, raw sampled richness, range-through (in "rawmetrics"),
# classic rarefaction ("rare"), and occurrences-squared-weighted (O2W)
# Note that for rarefaction and O2W we follow conventions and convert
# abundance data to lists of taxon occurrence in "collections" (=cells)
source("metrics.R")
SCOR=SCOR(pres,t=Nt)
raw=rawmetrics(pres,t=Nt)
rare=RARE(pres,Nt,Ntrial)
O2W=O2W(pres,Nt,Ntrial)

#-------------------------------------------------------------------------------------------#
# Shareholder Quorum Subsampling
# In order to run these, the user needs to obtain Alroy's SQS scripts
# sqs3.3.R (http://bio.mq.edu.au/~jalroy/SQS.html)
# and quorumSubsample4-3.pl (by contacting Alroy directly)

run_sqs=F # 
if(run_sqs==T){
  source("sqs3.3.R") # Alroy's code
  
  source("getsqs3.3.R")# Our code to call sqs3.3 after some data formatting
  alroyres=getsqs3.3(pres, abundance=T) 
  # Note: abundance=T here corresponds to unique=F in convert2pbdb() below
  
  # convert data to Paleobiology Database format for sqs4.3 perl script
  source("convert2pbdb.R")
  convert2pbdb(pres,Nt, unique=T, 1) 
  
  print("Estimating richness using SQS version 4.3")
  ptm <- proc.time()[3]
  
  # The user has to edit the perl script to
  # set the number of trials to match Ntrial above
  # set the quorum level to 0.6
  # use the "fr2" time scale
  system(sprintf("perl quorumSubsample4-3.pl")) # Alroy's code
  sqsres=read.table("quorumSubsample4.txt", header=T, sep="\t")
  
  etm=round(proc.time()[3] - ptm)
  print(paste(as.character(etm), " seconds (", as.character(round(etm/60)), " mins.)"))  
  print(" ")
  
}else{
  # these are just plot dummies in absence of SQS scripts
  alroyres=list(SQSres=NA)
  sqsres=list(Subsampled.diversity=NA)
}
#-------------------------------------------------------------------------------------------#
# bundle some of the output for easy plotting

dat=data.frame(realworld1$t, realworld1$indiv, realworld1$sp,  realworld1$PJ, 
               SCOR$scor, raw$rawsp, raw$RT, rare$CR.est, O2W$O2W.est,
               alroyres$SQSres, sqsres$Subsampled.diversity,
               siteprop,sig) # 

colnames(dat)=c("t", "ab", "sp", "PJ",
                "scor","raw", "rt", "cr", "o2w",
                "sqs3", "sqs4",
                "siteprop","sig")

nrmd<-function(x){
  # normalize data to zero mean and unit std. dev.
  nx=(x-mean(x,na.rm=T))/sd(x,na.rm=T)
  return(nx)
}

# A few quick plots comparing some of the richness estimators
# to true richness, and comparing SCOR to true abundance:

par(mfrow=c(2,2))
yl=c(-3,3) # y-axis limits for normalized data

plot(dat$t, nrmd(dat$sp), ylim=yl, col="black", type="l", lwd=2, 
     xlab="time bin", ylab="normalized")
lines(dat$t, nrmd(dat$raw), lty=1, col="red", lwd=2)
lines(dat$t, nrmd(dat$rt), lty=1, col="green", lwd=2)
title("raw richness, range-through")

plot(dat$t, nrmd(dat$sp), ylim=yl, col="black", type="l", lwd=2, 
     xlab="time bin", ylab="normalized")
lines(dat$t, nrmd(dat$cr), lty=1, col="red", lwd=2)
lines(dat$t, nrmd(dat$o2w), lty=1, col="green", lwd=2)
title("rarefaction, O2W")

plot(dat$t, nrmd(dat$sp), ylim=yl, col="black", type="l", lwd=2, 
     xlab="time bin", ylab="normalized")
lines(dat$t, nrmd(dat$sqs3), lty=1, col="red", lwd=2)
lines(dat$t, nrmd(dat$sqs4), lty=1, col="green", lwd=2)
title("SQS 3.3, SQS 4.3")

plot(dat$t, nrmd(dat$ab), ylim=yl, type="l", col="black", lwd=2, 
     xlab="time bin", ylab="normalized")
SCOR$scor[which(SCOR$scor==Inf)]=NA # strip away Infs
lines(nrmd(SCOR$scor), lty=1, col="red", lwd=2)
title("SCOR")

