#-------------------------------------------------------------------------------------------# 
# Sampling: subjects Poseidon plankton world to sampling or "preservation"
# 
# input: 
# Jumbled.list from poseidon()
# siteprop is the proportion of sites to sample (vector of length t)
# randrm is to an option to say whether species will be removed randomly
# randprop is the proportion of random sp. to remove (vector of length t)
# 
# output:
# k.n, the number of cells sampled 
# scor.input[[i]], the no. cells in which a sp. is sampled, where i is time
# pres.data which is the sampled Jumbled.list 
# 
# Modified from version preservation12.06.2015.R
# Written by Lee Hsiang Liow 2015
# (with contributions from Bjarte Hannisdal)
#
# This file is provided as auxiliary supplementary material to
# Hannisdal B, Haaga KA, Reitan T, Diego D, and Liow LH, 
# "Common Species Link Global Ecosystems to Climate Change".
#-------------------------------------------------------------------------------------------#

sampling<- function(Jumbled.list, t, siteprop, randrm, randprop){

print("Preserving and sampling fossil plankton world")
#ptm <- proc.time()[3]
  
Data.list= vector("list", t)  # prepare results list
k.n=NULL
if (randrm==T){print("Random removal of species")}

for (i in 1:t){
        
        newlist=Jumbled.list[[i]]
        cells=unique(newlist$cell)
        keep=sample(cells,(siteprop[i])*length(cells), replace=F) # sample a proportion of cells
        Data.list[[i]]=newlist[which(newlist$cell %in% keep),]
        k.n[i]=length(unique(Data.list[[i]]$cell))  
          
        
        if (randrm==T){
          newlist=Data.list[[i]]
          sp=unique(newlist$s)
          keep=sample(sp, (1-randprop[i])*length(sp), replace=F) # sample a proportion of species
          Data.list[[i]]=newlist[which(newlist$s %in% keep),]
          k.n[i]=length(unique(Data.list[[i]]$cell))  
        }
      }
      

output<-list(pres.data=Data.list)

#etm=round(proc.time()[3] - ptm)
#print(paste("Elapsed time: ", as.character(etm), " seconds (", as.character(round(etm/60)), " mins.)"))
print(" ")

return(output) 

}