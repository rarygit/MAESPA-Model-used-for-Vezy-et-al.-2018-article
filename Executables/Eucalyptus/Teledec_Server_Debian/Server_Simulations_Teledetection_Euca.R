#####################################################################################################
####-------------------------------Rstudio Linux server script ----------------------------------####
# Aim: Script to launch multiple simulations in parallel into a server using RSTUDIO.
# Author: Remi Vezy
# Date: 19/01/2017
####################################################################################################
# Before running the script, you have to enter several arguments:
# 1- number of CPU cores to leave free from any simulation:
nbCoresToLeave= 5
# 2- the name of maespa in use:
V_MAESPA= "maespa"
# 3- the path of the simulations (path to plots 1-7):
Sim.Path= "/donnees/goplus/users/rvezy/Simulations/2-Simulations_FullPeriod2010/V2.Iter20/Sim 1"


# Loading packages or install them:
library(parallel) # parallel comes with base R, no need to download
pckg= c("iterators","snow","foreach","doSNOW","data.table")
# librariesPath= "/home/rvezy/R_routines/R_library"
if(any(lapply(pckg, require, character.only = TRUE)==F)){
  install.packages(pckg, Ncpus= length(pckg))  
}


Make.MAESPA.executable= function(Sim.Path,V_MAESPA){
  # Aim: give permissions to execute maespa iteratively within each plot
  # Author: Rémi Vezy
  # Date: 19/01/2017
  Plots= list.files(Sim.Path)
  for(i in 1:length(Plots)){
    setwd(file.path(Sim.Path,Plots[i]))
    system(paste0("chmod u+x ",V_MAESPA))
  }
}

RunMAESPA.Parallel.whole= function(Sim.Path,nbCoresToLeave,V_MAESPA){
  # Aim: Running MAESPA in true parallel computing
  # Author: Rémi Vezy
  # Date: 19/01/2017
  Plots= list.files(Sim.Path)
  
  Simuls= vector(mode = "list", length = length(Plots))
  for(i in 1:length(Plots)){
    filesnames= list.files(file.path(Sim.Path,Plots[i]), recursive = T)
    Nsim= as.numeric(gsub("_confile.dat", "",filesnames[grep("confile", filesnames)], perl = T))
    Simuls[[i]]= data.frame(Plot= rep(Plots[i], length(Nsim)), Simul= Nsim)
  }
  Simuls2= as.data.frame(rbindlist(Simuls))
  
  Projectwd= getwd()
  NbCores= detectCores()-nbCoresToLeave # Leave "nbCoresToLeave" free processors
  cl<-makeCluster(NbCores)
  registerDoSNOW(cl)
  foreach(i= 1:nrow(Simuls2)) %dopar% {
    setwd(file.path(Sim.Path,Simuls2[i,1]))
    system2(paste0("./",V_MAESPA), input= as.character(Simuls2[i,2]), wait=T, invisible= F,
            stdout= paste(Simuls2[i,2],"_Console_Outputs.txt", sep=""))
  }
  stopCluster(cl)
  setwd(Projectwd)
}

Make.MAESPA.executable(Sim.Path,V_MAESPA)
RunMAESPA.Parallel.whole(Sim.Path,nbCoresToLeave,V_MAESPA)