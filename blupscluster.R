## This will be a simulation for approaches to estimating posterior means.


### TODO: Make a second version of this with varying cluster sizes and number of clusters
### TODO: Fix biased sampling based on positive vs negative (clusters need to be sampled differently not within clusters)
library(lme4)
library(svylme)
library(tidyverse)
#rcloud.execute.asset("HelpFunctions.R")
source("HelpFunctions.R")

N_vec <- c(50, 100, 200)
cluster_sizes <- c(20, 50, 100)
n_per <- c(2, 10, 20) ##
n_dfs <- 5 ## When we want to start generating different datasets we can run add a loop
n_iter <- 25


### Maybe declared later
N <- 100
cluster_size <- 50
n <- 5
### This will be for different size clusters. 

######  This is temporary
#pairwise_indx <- sapply(1:N, function(cluster_id){
#  local_data <- filter(df, g == cluster_id)
#  ## g is the cluster id, t is within cluster
#  pairs <- as.vector(combn(local_data$t, 2))
#  local_data$id[pairs] 
#  ## If the within cluster ids are not 1:cluster_size, this needs work
#})
#pairwise_vals <- get_pairwise_values(df, pairwise_indx, fit=fit, N=N)
#######

set.seed(10112)
df<-data.frame(x=rnorm(N*cluster_size),
               g=rep(1:N,each=cluster_size), ## g is the cluster id
               t=rep(1:cluster_size,N), ## t is an id within cluster
               id=1:(N*cluster_size))

df$u<-with(df, rnorm(N)[g])
df$e<-rnorm(N*cluster_size,s=2)
df$y<-with(df, x+u+e)

prob_cluster <- exp(abs(df$u/2)-2.2)[df$t==1]  
prob_outcome <- exp(-1.5 +(df$u>0)*0.5)
M_vec <- round(runif(N, n-0.5, cluster_size+0.5)) 

df$w_cluster <- 1/rep(prob_cluster, each=cluster_size) ## Full sample weights
## Level 2 weights
df$w1 <- cluster_size/n
df$w3 <- rep(exp(1/(M_vec/60+0.5)), each=cluster_size)
df$w4 <- (1/prob_outcome)

##### These are the approximate targets
fit <- lmer(y~x+(1|g), data=df, REML=FALSE)
BLUPS <- coef(fit)$g[,1] # - coef(summary(fit))[1,1]
tau <- 1# 2*fit@theta ## NO
# How to get random effects variance
qi<-sapply(fit@cnms,length)    
L<-as.matrix(Matrix::bdiag(lapply(qi,function(i) matrix(1,i,i))))
######(need indicator for where thetas go in the matrix)
ThInd<-which((L==1) & lower.tri(L,diag=TRUE))
#### getAnywhere("VarCorr.merMod")
#### mkVarCorr

#######

s2 <- fit@devcomp$cmp["sigmaML"]^2
Xhi_inv <- solve(s2 * 
                   matrix(c(1+tau^2, tau^2, tau^2, 1+tau^2), nrow=2, ncol=2)) 
L_mat <- chol(Xhi_inv)

pairwise_indx <- get_pairwise_indx(df, N, 1:N)
pairwise_BLUPS <- get_pairwise_values(df, pairwise_indx, fit, N)

## These can help get the true targets. e.g. variance will be larger
#test_design <- survey::svydesign(id=~g+id, data=df, weights=~wtest+wtest2)
#full_pairs <- svy2lme(y~x+(1|g), design=test_design)
#### ## We will cheat and use estimates from fit instead of full_pairs

output <- map_dfr(1:n_iter, function(iter){
  ## oversample extreme `u` to bias random-intercept variance
  ind_cluster <- rbinom(N, 1, prob_cluster)==1
  N1 <- sum(ind_cluster)
  
  m_vec <- round(runif(N, 1.5, M_vec/2+0.5))
  #### We need the weights for the 2 other sampling schemes.
  df$w2 <- rep(M_vec/m_vec, each=cluster_size)
  
  ind_obs1 <- 1:n ## sample n from each cluster
  ## Different sampling schemes
  ind_obs2 <- map(m_vec, seq) ## deterministic
  ## random based on cluster size
  ind_obs3 <- map(rbinom(N , M_vec, exp(-1/(M_vec/60+0.5))), function(M){seq(max(M, 2))})
  ## Based on outcome
  ind_obs4 <- map(seq(1,N), function(M){
    local_obs <- 0
    while(length(local_obs) < 3){
      local_obs<-which(rbinom(M_vec[M], 1, prob_outcome[df$g==M])==1)
    }
    return(local_obs)
  })## potentially infinite loop, lazysolution
  
  #### Create 4 sampled clusters.
  sdf1 <- subset(df, (g %in% (1:N)[ind_cluster]) & (t %in% ind_obs1))
  sdf2 <- get_subset_data(df, N, ind_cluster, ind_obs2)
  sdf3 <- get_subset_data(df, N, ind_cluster, ind_obs3)
  sdf4 <- get_subset_data(df, N, ind_cluster, ind_obs4)
  
  ######### Next the pairwise ordering is decided
  pairwise_indx1 <- get_pairwise_indx(sdf1, N, ind_cluster) 
  ## deepdive here to see if clusters are correctly order
  pairwise_indx2 <- get_pairwise_indx(sdf2, N, ind_cluster)
  pairwise_indx3 <- get_pairwise_indx(sdf3, N, ind_cluster)
  pairwise_indx4 <- get_pairwise_indx(sdf4, N, ind_cluster)
  
  ### TODO: Testing
  b1 <- get_predictions_1(df, pairwise_indx1, df$w1,
                          L_mat, tau[1], N1, M_vec = rep(cluster_size, each=N)) 
  b2 <- get_predictions_1(df, pairwise_indx2, df$w2,
                          L_mat, tau[1], N1, M_vec = M_vec) 
  b3 <- get_predictions_1(df, pairwise_indx3, df$w3,
                          L_mat, tau[1], N1, M_vec = M_vec) 
  b4 <- get_predictions_1(df, pairwise_indx4, df$w4, L_mat, tau[1], N1, M_vec = M_vec) 
  
  b5 <- get_predictions_2(df, pairwise_indx1, df$w1,
                          filter(sdf1, t==1)$w_cluster, tau=tau, N1=N1)[1:N1]
  b6 <- get_predictions_2(df, pairwise_indx2, df$w2,
                          filter(sdf2, t==1)$w_cluster, tau=tau, N1=N1)[1:N1]
  b7 <- get_predictions_2(df, pairwise_indx3, df$w3,
                          slice(group_by(sdf3, g), 1)$w_cluster, tau=tau, N1=N1)[1:N1]
  b8 <- get_predictions_2(df, pairwise_indx4, df$w4, slice(group_by(sdf4, g), 1)$w_cluster, tau=tau, N1=N1)[1:N1]
  
  b9 <- coef(lmer(y~x+(1|g), data=sdf1, REML=FALSE))$g[,1]
  
  b11 <- get_predictions_1(df, pairwise_indx1, df$w1,
                           L_mat, tau[1], N1, M_vec = rep(cluster_size, each=N), cluster_weights=FALSE)
  b12 <- get_predictions_1(df, pairwise_indx2, df$w2,
                           L_mat, tau[1], N1, M_vec = M_vec, cluster_weights=FALSE)
  b13 <- get_predictions_1(df, pairwise_indx3, df$w3,
                           L_mat, tau[1], N1, M_vec = M_vec, cluster_weights=FALSE) 
  b14 <- get_predictions_1(df, pairwise_indx4, df$w4,
                           L_mat, tau[1], N1, M_vec = M_vec, cluster_weights=FALSE) 
  
  data.frame("iter" = iter, "cluster" = unique(sdf1$g),
             "Cluster_SRS" = b1, "Cluster_SRSP" = b2,
             "Cluster_Size" = b3, "Cluster_Latent" = b4,
             "mi"=n, "Mi"=cluster_size, "N_sample" = nrow(sdf1),
             "V2_SRS" = b5, "V2_SRSP" = b6,
             "V2_Size" = b7, "V2_Latent" = b8, "Model"=b9,
             "NoC_SRS" = b11, "NoC_SRSP" = b12, "NoC_Size"=b13,
             "NoC_Latent"=b14)
})

write_csv(output, "Data/blups.csv")

