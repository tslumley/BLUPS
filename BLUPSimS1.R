## This will be a simulation for approaches to estimating posterior means.
library(lme4)
library(svylme)
library(purrr)

N <- 200 ## number of clusters
cluser_sizes <- c(20, 50, 100)
n_per <- c(2, 10, 20) ##
n_dfs <- 5 ## When we want to start generating different datasets we can run add a loop
n_iter <- 1000

set.seed(10012)
df<-data.frame(x=rnorm(N*cluster_size),
               g=rep(1:N,each=cluster_size), ## g is the cluster id
               t=rep(1:cluster_size,N), ## t is an id within cluster
               id=1:(N*cluster_size))

df$u<-with(df, rnorm(N)[g])
df$y<-with(df, x+u+rnorm(N,s=2))
df$wtest <- 1 ## Full sample weights
df$wtest2 <- 1

##### Takes a long time, these are the targets
fit <- lmer(y~x+(1|g), data=df, REML=FALSE)
test_design <- survey::svydesign(id=~g+id, data=df, weights=~wtest+wtest2)
full_pairs <- svy2lme(y~x+(1|g), design=test_design)
#### ## Can cheat and use estimates from fit instead of full_pairs


cluster_size <- 50
n <- 2
M_vec <- round(runif(N, n-0.5, cluster_size+0.5))
m_vec <- round(runif(N, 1.5, cluster_size/2+0.5))
### This will be for different size clusters. 

map_dfr(n_iter, function(iter){
  ## oversample extreme `u` to bias random-intercept variance
  prob_cluster <-exp(abs(df$u/2)-2.2)[df$t==1]  
  
  ind_cluster <- rbinom(N,1, prob_cluster)==1
  ind_obs1 <- rep(1:n, length(ind_cluster)) ## sample n from each cluster
  
  ## Different sampling scheme
  ind_obs2 <- unlist(map(m_vec[ind_cluster], seq))
  
  sdf1 <- subset(df, (g %in% (1:N)[ind_cluster]) & (t %in% ind_obs1))
  sdf2 <- subset(df, (g %in% (1:N)[ind_cluster]) & (t %in% ind_obs2)) 
  
  p1 <-rep(prob_cluster[ind_cluster], each=5)
  N_1 <-rep(cluster_size, nrow(sdf1))
  N_2 <-rep(cluster_size, nrow(sdf2))
  
})


