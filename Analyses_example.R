# R-code: Analyses r-scripts
# manuscrpt titled "Information arms race explains plant-herbivore chemical communication in ecological communities”
# Authors: Zu P*, del-Val E, Boege K, Schuman M, Stevenson P, Alejandro Zaldivar-Riveron, Saavedra S. 
# correspondence to:  Pengjuan.zu@gmail.com
# cite this paper: 


#rm(list = ls())

# load library
library(mgcv) 
library(ggplot2)
library(reshape2)
library(tidyverse)
source("toolbox.r") # load functions related to mutual information

## Set work directory where you have these files
setwd()

# proposed information arms-race mechanism "A(min)-P(max)", and three alternatives
Mech_vectors <- c("A(min)-P(max)", "A(min)-P(min)","A(max)-P(min)", "A(max)-P(max)")

######################################
### settings for simulation number ###
######################################
## set a small number to play with the code
# the bigger number, the longer time for running

# set for evolutionary time steps
N_sim <- 100 ## N_sim <- 1000  # in our analyses 
# set for calculating cumulative mutual information PV
Ns = 100 ## Ns = 5000 in our analyses
# set number of random sampling round, cumulative mutual informationfor AV 
N_test <- 30
# choose echanism for simulation
mech <- 1 # for A(min)-P(max) mechanism. can choose 2,3,4 for other mechanism

###################################################
## load field observation matrix AP_obs, PV_obs ###
###################################################
AP_obs <- read.csv("AP_obs.csv", header = TRUE, as.is = TRUE, row.names = 1)
PV_obs <- read.csv("PV_obs.csv", header = TRUE, as.is = TRUE, row.names = 1)

nA <- nrow(AP_obs) # number of herbivores
nP <- ncol(AP_obs) # number of plants
nV <- ncol(PV_obs) # number of VOCs

## generate matrix for simulation starting
# starting matrices do not affect simulation. 
# Can also use more specialised/generalized matrices for starting
# Make sure colSum > 0. i.e., each plant interacts with at least one herbivore; each VOC is emited from at least one plant 
# here we use random matrix for illustration 
AP <- matrix(rbinom(nA*nP,1,0.5),nA,nP)
PV <- matrix(rbinom(nP*nV,1,0.5),nP,nV)
colSums(AP) # Make sure all colSum > 0
colSums(PV) # Make sure all colSum > 0
# element mutation circles
M1 <- 0.2*nP*nV # 20% of PV matrix. However, results remain the same if setup other %. eg, 5%
M2 <- 0.2*nA*nP # 20% of AP matrix. However, results remain the same if setup other % 

########################################
### Analyses1: Simulation process #####
#######################################
# create plot
plot(-1, xlim = c(0,N_sim), ylim = c(0,1), ylab = Mech_vectors[mech] )
# calculate and plot entropy and fitness based on field observations
{
  AV_obs <- as.matrix(AP_obs) %*% as.matrix(PV_obs)
  # function "H_A_VOC" in toolbox.r to calculate entropy when feeding matrix
  H_V   <- H_A_VOC(PV_obs)[["Hn_V"]]
  H_P_V <- H_A_VOC(PV_obs)[["Hn_S_V"]]
  H_V_P <- H_A_VOC(PV_obs)[["Hn_V_S"]]
  H_A_V <- H_A_VOC(AV_obs)[["Hn_S_V"]]
  H_V_A <- H_A_VOC(AV_obs)[["Hn_V_S"]]
  H_A_P <- H_A_VOC(AP_obs)[["Hn_S_V"]]
  H_P_A <- H_A_VOC(AP_obs)[["Hn_V_S"]]
  
  # fitness function. function details in toolbox.r
  E_plant <- F_plt(mech)
  E_animal<- F_anm(mech)
  
  # write to E, A_AP, A_VOC
  E_obs <- c(0, E_plant, E_animal, H_V, H_P_V, H_V_P, H_A_V, H_V_A, H_A_P, H_P_A)
  names(E_obs) <- c("N","E_plant", "E_animal", "H_V", "H_P_V", "H_V_P", "H_A_V", "H_V_A", "H_A_P", "H_P_A")
  
  ## observed values as lines
  # conditional entropy for PV, AV and AP
  abline(h = H_P_V, lty = "dotted", pch = 8) 
  abline(h = H_A_V, lty = "solid", pch = 8)
  abline(h = H_A_P, lty = "dashed", pch = 8) 
  # fitness for animal and plant
  abline(h = E_animal, col = "blue") 
  #  abline(h = E_plant, col = "green3") # the same as H(A|V) in mech = 1
}
## save values from field observation
#saveRDS(E_obs, paste("E_obs",year,nA, nP, nV,".RDATA", sep = "_"))

# Simulated values
# create list to store matrix and fitness file 
A_PV <- list()  # matrix to store PV after each simulation
A_AP  <- list() # matrix to store AP after each simulation
E = matrix(NA, nrow = N_sim +1, ncol = 10) # dataset to store variables after each simulation
colnames(E) <- c("N", "E_plant", "E_animal", "H_V", "H_P_V", "H_V_P", "H_A_V", "H_V_A", "H_A_P", "H_P_A")

## first measurements based on starting simuation matrices
n = 1
{
  AV <- AP %*% PV
  # function in MutInf_functions.r to calculate entropy when feeding matrix
  H_V   <- H_A_VOC(PV)[["Hn_V"]]
  H_P_V <- H_A_VOC(PV)[["Hn_S_V"]]
  H_V_P <- H_A_VOC(PV)[["Hn_V_S"]]
  H_A_V <- H_A_VOC(AV)[["Hn_S_V"]]
  H_V_A <- H_A_VOC(AV)[["Hn_V_S"]]
  H_A_P <- H_A_VOC(AP)[["Hn_S_V"]]
  H_P_A <- H_A_VOC(AP)[["Hn_V_S"]]
  
  # fitness function
  E_plant <- F_plt(mech) 
  E_animal<- F_anm(mech) 
  
  # plot initial values
  points((2*n),  H_P_V, col = "red", pch = 24) 
  points((2*n),  H_A_V, col = "red", pch=  21)
  points((2*n),  H_A_P, col = "red", pch=  22) 
  points((2*n),  E_animal, col = "blue", pch = 4) 

  #write to E, A_AP, A_PV
  E[n,] <- c(n, E_plant, E_animal, H_V, 
             H_P_V, H_V_P, H_A_V, H_V_A, H_A_P, H_P_A)
  A_AP[[n]] <- AP
  A_PV[[n]]<- PV
}


### herbivores and plants play in turn to maximize their fitness
for (n in 1:(N_sim/2))
{
  print(n)
  ## Herbivore circles: herbivore maximize its fitness
  for (m in 1:M2) {
    AP_new <- random.sample_1element1(AP) #function details see toolbox.r
    AV <- AP_new %*% PV
    H_V   <- H_A_VOC(PV)[["Hn_V"]]
    H_P_V <- H_A_VOC(PV)[["Hn_S_V"]]
    H_V_P <- H_A_VOC(PV)[["Hn_V_S"]]
    H_A_V <- H_A_VOC(AV)[["Hn_S_V"]]
    H_V_A <- H_A_VOC(AV)[["Hn_V_S"]]
    H_A_P <- H_A_VOC(AP_new)[["Hn_S_V"]]  #AP_new
    H_P_A <- H_A_VOC(AP_new)[["Hn_V_S"]]  #AP_new
    
    # fitness function
    E_plant_new <- F_plt(mech)
    E_animal_new <- F_anm(mech)
    
    # update based on optimization mechanism
    if(E_animal_new > E_animal)
    {E_animal <- E_animal_new
    E_plant  <- E_plant_new
    AP <- AP_new 
    }
  }
  
  # use the matrix for calculation
  AV <- AP %*% PV
  H_V   <- H_A_VOC(PV)[["Hn_V"]]
  H_P_V <- H_A_VOC(PV)[["Hn_S_V"]]
  H_V_P <- H_A_VOC(PV)[["Hn_V_S"]]
  H_A_V <- H_A_VOC(AV)[["Hn_S_V"]]
  H_V_A <- H_A_VOC(AV)[["Hn_V_S"]]
  H_A_P <- H_A_VOC(AP)[["Hn_S_V"]]
  H_P_A <- H_A_VOC(AP)[["Hn_V_S"]]
  
  # fitness function
  E_plant <- F_plt(mech)
  E_animal<- F_anm(mech)
  
  #write to E, A_AP, A_PV
  E[(2*n),] <- c((2*n), E_plant, E_animal, H_V, H_P_V, H_V_P, H_A_V, H_V_A, H_A_P, H_P_A)
  A_AP[[(2*n)]] <- AP
  A_PV[[(2*n)]]<- PV
  
  # plot simulated values
  points((2*n),  H_P_V, col = "red", pch = 24) 
  points((2*n),  H_A_V, col = "red", pch=  21)
  points((2*n),  H_A_P, col = "red", pch=  22) 
  points((2*n),  E_animal, col = "blue", pch = 4) 

  
  ## Plant circles: plant maximize its fitness
  for (m in 1:M1) {
    PV_new <- random.sample_1element1(PV)
    AV <- AP %*% PV_new
    H_V   <- H_A_VOC(PV_new)[["Hn_V"]]    #update to PV_new
    H_P_V <- H_A_VOC(PV_new)[["Hn_S_V"]]
    H_V_P <- H_A_VOC(PV_new)[["Hn_V_S"]]
    H_A_V <- H_A_VOC(AV)[["Hn_S_V"]]
    H_V_A <- H_A_VOC(AV)[["Hn_V_S"]]
    H_A_P <- H_A_VOC(AP)[["Hn_S_V"]]
    H_P_A <- H_A_VOC(AP)[["Hn_V_S"]]
    
    # fitness function
    E_plant_new <- F_plt(mech)
    E_animal_new <- F_anm(mech)
    
    if(E_plant_new > E_plant)
    { E_plant <- E_plant_new
    E_animal<- E_animal_new
    PV <- PV_new
    }
  }
  
  # use the matrix for calculation
  AV <- AP %*% PV
  H_V   <- H_A_VOC(PV)[["Hn_V"]]
  H_P_V <- H_A_VOC(PV)[["Hn_S_V"]]
  H_V_P <- H_A_VOC(PV)[["Hn_V_S"]]
  H_A_V <- H_A_VOC(AV)[["Hn_S_V"]]
  H_V_A <- H_A_VOC(AV)[["Hn_V_S"]]
  H_A_P <- H_A_VOC(AP)[["Hn_S_V"]]
  H_P_A <- H_A_VOC(AP)[["Hn_V_S"]]
  
  # fitness function
  E_plant <- F_plt(mech)
  E_animal<- F_anm(mech)
  
  #write to E, A_AP, A_PV
  E[(2*n+1),] <- c((2*n+1), E_plant, E_animal, H_V, H_P_V, H_V_P, H_A_V, H_V_A, H_A_P, H_P_A)
  A_AP[[(2*n+1)]] <- AP
  A_PV[[(2*n+1)]]<- PV
  
  # plot simulated values
  points((2*n),  H_P_V, col = "red", pch = 24) 
  points((2*n),  H_A_V, col = "red", pch=  21)
  points((2*n),  H_A_P, col = "red", pch=  22) 
  points((2*n),  E_animal, col = "blue", pch = 4) 
}

## save output 
#saveRDS(A_PV, paste("M_PV", Mech_vectors[mech], nA, nP, nV, "test.RDATA",sep = "-"))
#saveRDS(A_AP, paste("M_AP", Mech_vectors[mech], nA, nP, nV, "test.RDATA",sep = "-"))
#saveRDS(E, paste("E", Mech_vectors[mech], nA, nP, nV,"test.RDATA",sep = "-"))



##########################################################################
##### Analyses2: cumulative mutual information along number of VOCs  #####
##########################################################################

########################################
## for plant coding: PV matrix
# will need to run it for observed matrix and for PV-matrix resulted from the 4 mechanism 
# we use observed matrix for illustration
PV <- PV_obs
AP <- AP_obs

I_PV_Ls_r = matrix(NA, nrow = Ns, ncol = nV)
for(j in 1:nV) {
  print(paste("Library size = ", j, "simulation number =", Ns))
  index_list <- replicate(n= Ns, sample(x = 1:ncol(PV), size = j, replace = FALSE), simplify = FALSE)
  #dim(I) = c(Ns,2)
  for (i in 1:length(index_list)) {
    P_vi <- Pro_A_Lib(PV[, index_list[[i]]])[[1]] # function see toolbox.r
    IVi  <- Pro_A_Lib(PV[, index_list[[i]]])[[2]]
    # plant species level
    I_PV_Ls_r[i,j] <- mut.inf_pr(P_vi, IVi)
  }
}

I_PV <- cbind.data.frame(I_PV_Ls_r, "type" = "PlantSp", "Mech" = "obs") # change here if it is other mechanism
#saveRDS(I_PV, paste("I_PV", Mech_vectors[mech],".RDATA",sep = "-"))

## plot 
I_plant <- melt(I_PV)
colnames(I_plant)[3]  <- "VOC_number"
colnames(I_plant)[4]  <- "Mut.Inf"
I_agg <- I_plant %>% 
  group_by(Mech, VOC_number) %>% 
  dplyr::summarise(Max.Mut.Inf = max(Mut.Inf)) # max for each library
I_agg$VOC_number <- as.numeric(I_agg$VOC_number) 

# when the I(Lsize = n+1) < I(Lsize=n), give the previous value
# This is to correct for sampling effort
I_agg_max <- NULL
for (m in 1: 5) {
  I_agg1 <- I_agg
  for (i in 2:nV) {
    if(I_agg1$Max.Mut.Inf[i] < I_agg1$Max.Mut.Inf[i-1]) {I_agg1$Max.Mut.Inf[i] <- I_agg1$Max.Mut.Inf[i-1]}
  }
  I_agg_max <- rbind.data.frame(I_agg_max, I_agg1)
}

ggplot(I_agg_max, aes(y = Max.Mut.Inf, x = VOC_number, color = Mech, group = Mech)) +
  geom_line(size = 1, alpha = 0.6) +
  scale_x_continuous(breaks=seq(0, nV, by = 5), limits = c(0,nV)) +
  ggtitle("Plant coding")


########################################
## for herbivore decoding AV matrix
### function
M_rdm <- function(r, c) {
  M <- matrix(runif((r*c), 0, 1), nrow = r, ncol = c)
  return(M)
}

M_rdm_all <- replicate(N_test, M_rdm(nA, nV))

### assign the matrices for test
AV <- as.matrix(AP) %*% as.matrix(PV)
Pla <- rowSums(AP)
AV_sc <- AV/Pla

I_AV_result <- NULL
for (nr in 1:N_test) {
  print(nr)
  Mr <- M_rdm_all[[nr]]
  AV <- (AV_sc > Mr)*1 # if the value bigger than random then assign 1, otherwise, assign 0
  I_AV_Ls_r = matrix(NA, nrow = Ns, ncol = nV)
  for(j in 1:nV) {
    print(paste("Library size = ", j, "simulation number =", Ns))
    index_list <- replicate(n= Ns, sample(x = 1:ncol(AV), size = j, replace = FALSE), simplify = FALSE)
    for (i in 1:length(index_list)) {
      P_vi <- Pro_A_Lib(AV[, index_list[[i]]])[[1]]
      IVi  <- Pro_A_Lib(AV[, index_list[[i]]])[[2]]
      I_AV_Ls_r[i,j] <- mut.inf_pr(P_vi, IVi)
    }
  }
  I_AV <- cbind.data.frame(I_AV_Ls_r, "Type" = "AV", "Mech" = "obs", "SimN" = nr) # # change here if it is other mechanism
  I_AV_melt <- melt(I_AV, id.vars = c("Type", "Mech", "SimN"))
  I_AV_result <- rbind.data.frame(I_AV_result, I_AV_melt)
}

#saveRDS(I_AV_result, paste("I_AV", Mech_vectors[mech],".RDATA",sep = "-"))

I2 <- I_AV_result
colnames(I2)[4]  <- "VOC_number"
colnames(I2)[5]  <- "Mut.Inf"
I2$VOC_number <- as.numeric(I2$VOC_number) 

#####
# max for each round of simulation and for each VOC_number
I2_agg1 <- I2 %>% 
  group_by(Mech ,SimN, VOC_number) %>% 
  dplyr::summarise(Max.I.Sim = max(Mut.Inf)) #%>%
#### due to the non-normal distribution of I_agg_max, take the median value for final analysis
I2_agg <- I2_agg1 %>% 
  group_by(Mech, VOC_number)  %>% 
  dplyr::summarise(Mean.Mut.Inf = median(Max.I.Sim))

########### with max value if VOC_number n +1 
I2_agg_max <- NULL
for (i in 2:31 ) {
  I2_agg$Mean.Mut.Inf[i] <- max(I2_agg$Mean.Mut.Inf[i-1], I2_agg$Mean.Mut.Inf[i] )
}
I2_agg_max <- rbind.data.frame(I2_agg_max, I2_agg)

#plot
ggplot(I2_agg_max, aes(y = Mean.Mut.Inf, x = VOC_number, color = Mech, group = Mech)) +
  geom_line(size = 1, alpha = 0.6) +
  scale_x_continuous(breaks=seq(0, nV, by = 5), limits = c(0, nV)) +
  ggtitle("Herbivore decoding")
