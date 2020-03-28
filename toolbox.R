# R-code: Functions for the main r-scripts
# manuscrpt titled"Information arms race explains plant-herbivore chemical communication in ecological communities”
# Authors: Zu P*, Del-Val E, Boege K, Schuman M, Stevenson P, Alejandro Zaldivar-Riveron, Saavedra S. 

library(mgcv)
############################################
## Function of fitness F_plant, F_animal ###
############################################
# for different simulation mechanisms 
# i.e., Mech_vectors <- c("A(min)-P(max)", "A(min)-P(min)","A(max)-P(min)", "A(max)-P(max)")
# we have different fitness functions accordingly to adjust max and min 
# for animal
F_anm <- function(mech){
  E1 <- (1-H_V_A) # n=1 for A_min & P_max
  E2 <- (1-H_V_A) # n=2 for A_min & P_min
  E3 <- (H_V_A)   # n=3 for A_max & P_min
  E4 <- (H_V_A)   # n=4 for A_max & P_max
  E <- list(E1,E2,E3, E4)
  return(E[[mech]])
}

# for plant
F_plt <- function(mech){
  E1 <- (H_A_V)   # n=1 for A_min & P_max
  E2 <- (1-H_A_V) # n=2 for A_min & P_min
  E3 <- (1-H_A_V) # n=3 for A_max & P_min
  E4 <- (H_A_V)   # n=4 for A_max & P_max
  E <- list(E1,E2,E3, E4)
  return(E[[mech]])
}


############################################
##   Function of entropy calculation   ###
############################################
## Calculate entropy: 
# P_vector is a probability vector 
entropy <- function(P_vector){
  entropy <- 0
  for (i in 1:length(P_vector)) {
    if(P_vector[i] != 0) {  #zero check in class
      entropy0 <- -(P_vector[i]*log(P_vector[i], base = length(P_vector))) # for binary data set 
    }else{
      entropy0 <- 0
    }
    entropy <- entropy + entropy0
  }
  return(entropy)  #dont forget return!
}

# for average individual VOC when feeding a matrix
H_A_VOC <- function(x1) {
  x <- x1[, colSums(x1)!=0] # remove columns that only contains 0
  NI <- nrow(x) # number of species
  NJ <- ncol(x) # number of VOCs
  if ( is.null(NJ) == TRUE) {
    Hn_V <- 1
    Hn_V_S <- 1
    P_s_v <- x/sum(x)
    Hn_S_V <- entropy(P_s_v)
  }
  if (is.null(NJ) == FALSE) {
    # Average mut.ind VOC
    { P_v = c()
    P_s_v = matrix(NA, NI, NJ)
    Hn_S_vj = c()
    P_v_s = matrix(NA, NI, NJ)
    Hn_V_si = c()
    for (j in 1: NJ) {
      P_v[j] <- colSums(x)[j]/sum(x)
      for (i in 1: NI) {
        P_s_v[i, j] <- x[i, j]/colSums(x)[j]
        P_v_s[i, j] <- x[i, j]/rowSums(x)[i]
      }
    }
    }
    Hn_S_vj <- apply(P_s_v, 2, entropy)  # scale matrix so that colsum is 1
    Hn_S_V <- sum(P_v*Hn_S_vj)
    Hn_V_si <- apply(P_v_s, 1, entropy)
    Hn_V_S  <- sum((1/NI)*Hn_V_si)
    Hn_V <- entropy(P_v)
  }
  results <- list(Hn_S_V, Hn_V_S, Hn_V)
  names(results) <- c("Hn_S_V", "Hn_V_S", "Hn_V")
  return(results)
}


##############################################################
##   Function of entropy calculation with multiple columns ###
##############################################################
# whole VOC library (lib)
Pro_A_Lib <- function(x) { 
  # need: library(mgcv) This library can calculate identical rows (or col) and give a index for unique rows
  L_unique <- uniquecombs(x)
  index <- attr(L_unique, "index")
  NI_L <- length(index) #library size, as rows in probability matrix
  NJ_L <- length(unique(index)) # unique library number, as columns in probability matrix
  P_s_l = matrix(NA, nrow = NI_L, ncol=NJ_L)
  P_l = c()
  for(j in 1:NJ_L) {
    count_L <- length(index[which(index==unique(index)[j])]) # number of cases that share the same library
    for(i in 1:NI_L) {
      if(index[i] == unique(index)[j]){
        P_s_l[i,j] <- 1/count_L
      }else{P_s_l[i,j] <- 0}
    }
    P_l[j] <- count_L/NI_L
  }
  results <- list(P_l, P_s_l)
  names(results) <- c("P_l", "P_s_l")
  return(results)
}

#######################################
##   calculating mutual information ###
#######################################
## fomular for mutual information between species (S) and VOC (V):  
# mut.inf = Hm(S) - Hn(S|V)
# where Hn(S|V) = sum(P_vj * Hn(S|vj))
# P_si, P_vj are probability for s and v respectively. 
# P_si_vj is probability matrix with i rows and j columns. P(si|vj): given vj, the probability of observing si

mut.inf_pr <- function(P_vj, P_si_vj) {  # one more parameter can be added depending on Hm_S <- entropy(P_si)
  Hn_S_vj <- apply(scale.matrix(P_si_vj), 2, entropy)  # scale matrix so that colsum is 1. Function see below
  Hn_S_V <- sum(P_vj*Hn_S_vj)
  mut.inf <- 1 - Hn_S_V
  return(mut.inf)
}

# standardize matrix so that colSum =1
scale.matrix <- function(x) {
  A1 = matrix(NA, nrow = nrow(x), ncol = ncol(x))
  for (j in 1: ncol(x)) {
    for (i in 1: nrow(x)) {
      if(colSums(x)[j] != 0) {A1[i,j] <- x[i, j]/colSums(x)[j]}
      else{A1[i,j] <- 0}
    }
  }
  return(A1)
}


############################################
## Function of random mutation in matrix ###
############################################
# change one element of the matrix (0, 1) flip
random.sample_1element1 <- function(x) {
  n = 1
  repeat {
    print(paste( "Sample",n))
    l <- sample(length(x), 1) #sample(length(HP), length(HP))
    x_new <- x
    x_new[l] <- 1- x[l]
    l_r <- Length_unique_A(x_new)[[1]]
    l_c <- Length_unique_A(x_new)[[2]]
    n <- n +1
    print(paste("l_row",l_r))
    print(paste("l_col", l_c))
    # exit if the condition is met
    if (min(rowSums(x_new)) >0 & min(colSums(x_new)) >0 )  break # to do check
  }
  return(x_new)
}


### check unique col and row patterns to make sure no redundant lines
# whole VOC library (lib)
Length_unique_A <- function(x) { 
  # need: library(mgcv) This library can calculate identical rows (or col) and give a index for unique rows
  L_uniqueR <- uniquecombs(x)
  indexR <- attr(L_uniqueR, "index")
  #  NI_L <- length(index) #library size, as rows in probability matrix
  NJ_R <- length(unique(indexR)) # unique library number, as columns in probability matrix
  
  L_uniqueC <- uniquecombs(t(x))
  indexC <- attr(L_uniqueC, "index")
  #  NI_L <- length(index) #library size, as rows in probability matrix
  NJ_C <- length(unique(indexC)) # unique library number, as columns in probability matrix
  
  results <- list(NJ_R, NJ_C)
  names(results) <- c("NJ_R", "NJ_C")
  return(results)
}


#####################
### Functions End ###
##################### 

