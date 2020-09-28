library(lme4)
library(tidyverse)
library(svylme)

N <- 1000
n <- 10
N_vec <- c(1000, 10000, 100000)
n_per <- c(2, 10, 25) 
cluster_size <- 20
### Setup
set.seed(10012)
### svylme example
df<-data.frame(x=rnorm(N*cluster_size),
               g=rep(1:N,each=cluster_size), ## g is the cluster id
               t=rep(1:cluster_size,N), ## t is an id within cluster
               id=1:(N*cluster_size))

df$u<-with(df, rnorm(N)[g])

df$y<-with(df, x+u+rnorm(N,s=2))

## oversample extreme `u` to bias random-intercept variance
pg <- exp(abs(df$u/2)-2.2)[df$t==1]  

in1<-rbinom(N,1,pg)==1
in2<-rep(1:5, length(in1)) ## sample 5 from each cluster

sdf<-subset(df, (g %in% (1:N)[in1]) & (t %in% in2)) 

p1<-rep(pg[in1],each=5)
N2<-rep(cluster_size, nrow(sdf))

## Population values
lme4::lmer(y~x+(1|g), data=df, REML=FALSE)

## Naive estimator: higher intercept variance
lme4::lmer(y~x+(1|g), data=sdf, REML=FALSE)

##pairwise estimator
sdf$w1<-1/p1
sdf$w2<-cluster_size/5
df$wtest <- 1 ## Doesn't change the estimates of beta because it is equal prob sampling
df$wtest2 <- 1

design<-survey::svydesign(id=~g+id, data=sdf, weights=~w1+w2)
pair <- svy2lme(y~x+(1|g), design=design)

###### ###### ###### ###### ###### ###### ###### ######
###### ###### ###### ###### ###### ###### ###### ###### 
######   CHECKING PERFORMANCE ##### ##### #####  ###### 
###### ###### ###### ###### ###### ###### ###### ###### 
###### ###### ###### ###### ###### ###### ###### ###### 

fit <- lmer(y~x+(1|g), data=df, REML=FALSE)

#### We may later try to get brute force estimator
#### For now try to get pairwise likelihood estimator
#### Then I compare to svy2lme
test_design <- survey::svydesign(id=~g+id, data=df, weights=~wtest+wtest2)
full_pairs <- svy2lme(y~x+(1|g), design=test_design)

pairwise_indx <- sapply(1:N, function(cluster_id){
  local_data <- filter(df, g == cluster_id)
  ## g is the cluster id, t is within cluster
  pairs <- as.vector(combn(local_data$t, 2))
  local_data$id[pairs] ## If the within cluster ids are not 1:cluster_size, this needs work
}) ## If ID are not 1:N*M then this needs work

##### Brute force all the terms
Yp <- df$y[pairwise_indx]
Xp <- model.matrix(~x, df)[pairwise_indx, ]
Zp <- model.matrix(~1, df)[pairwise_indx] ## can add z here
tau <- sqrt(full_pairs$L)
Xhi_inv <- solve(full_pairs$s2[1] * matrix(c(1+tau^2, tau^2, tau^2, 1+tau^2), nrow=2, ncol=2)) 
G <- diag(tau[1]^2, nrow=N*1)
### NEED TO GET G MATRIX,  ^^ CAN BE CORRECTED TO USE FULL DATA, e.g. lmer or full_pairs

L_mat <- chol(Xhi_inv)
XL_mat <- matrix(sapply(1:(nrow(Xp)/2), function(ind){
  crossprod(Xp[c(2*ind - 1, 2*ind), ], L_mat)
}), ncol=2, byrow = TRUE)
LY_p <-  matrix(sapply(1:(nrow(Xp)/2), function(ind){
  crossprod(Yp[c(2*ind - 1, 2*ind)], L_mat)
}), ncol=1, byrow=TRUE)
##### END OF DEFINITIONS. Missing G.
## This isn't exact but it's pretty close. 
beta_P <- solve(crossprod(XL_mat)) %*% crossprod(XL_mat, LY_p)
###### 

BLUPS <- coef(fit)$g[,1] # - coef(summary(fit))[1,1]
## We may have to subtract the fixed intercept: probably

resids <- Yp - Xp%*%beta_P
ZL_vec <-  matrix(sapply(1:(nrow(Xp)/2), function(ind){ ## fix when it works
  crossprod(Zp[c(2*ind - 1, 2*ind)], L_mat) ## adjust when it's a matrix
}), ncol=1, byrow = TRUE) ## change 1 for dimension z

ZL_mat <- matrix(0, nrow=nrow(Xp), ncol=N) ## Need to extend for matrices
cluster_indx <- df$g[pairwise_indx]
for(i in 1:nrow(Xp)){
  ZL_mat[i, cluster_indx[i]] <- ZL_vec[i]
}

Lresi_p <- matrix(sapply(1:(nrow(Xp)/2), function(ind){
  crossprod(resids[c(2*ind - 1, 2*ind)], L_mat)
}), ncol=1, byrow=TRUE)

## Everything works.
bis  <- G %*% crossprod(ZL_mat, Lresi_p) / 19 ## This is the number of times each observation is included

ZL_mat <- matrix(ZL_vec[1:380], ncol=1)
Gstar <- as.matrix(G[1,1])

subset_indx <- which(pairwise_indx[1:380] %in% c(1,2,3,4,5))
b1  <- Gstar %*% crossprod(ZL_mat, Lresi_p[1:380]) / 19
bis_un  <- Gstar %*% crossprod(ZL_mat[subset_indx], Lresi_p[subset_indx])*4 / 19
bis_small  <- Gstar %*% crossprod(ZL_mat[subset_indx], Lresi_p[subset_indx]) / 19

