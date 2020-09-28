# New file PairwiseBLUPs.R
# New file SimPairwiseBLUPS

### TODO
### 1) ideally show that as n -> infty that bias -> 0 and order
### 2) root n (bhat - blup) distribution for various n and  vairance estimates
library(lme4)
library(svylme)
library(purrr)
library(dplyr)
library(ggplot2)

#rcloud.execute.asset("HelpFunctions.R")
source("HelpFunctions.R")

N <- 10 ## number of clusters
n_iter <- 100
#cluster_size <- 20
cluster_vec <- c(50)#c(10, 25, 50, 100)#, 500)#, 1000)

##### Takes a long time, these would give better results
#test_design <- survey::svydesign(id=~g+id, data=df, weights=~wtest+wtest2)
#full_pairs <- svy2lme(y~x+(1|g), design=test_design)
#### ## We will cheat and use estimates from fit instead of full_pairs

set.seed(10012)
output2 <- map_dfr(cluster_vec, function(cluster_size){
  map_dfr(1:n_iter, function(iter){
    df<-data.frame(x=rnorm(N*cluster_size),
                   g=rep(1:N,each=cluster_size), ## g is the cluster id
                   t=rep(1:cluster_size,N), ## t is an id within cluster
                   id=1:(N*cluster_size))
    
    df$u<-with(df, rnorm(N)[g]) ### tau = 1
    df$e<-rnorm(N*cluster_size,s=2) ## sigma =2
    df$y<-with(df, x+u+e)
    
    ### Model fit
    fit <- lmer(y~x+(1|g), data=df, REML=FALSE)
    BLUPS <- coef(fit)$g[,1] # - coef(summary(fit))[1,1]
    ###
    
    pairwise_indx <- get_pairwise_indx(df, N, 1)
    pairwise_BLUPS <- get_pairwise_values(df, pairwise_indx, fit, 1)
    #        pairwise_indx <- get_pairwise_indx(df, N, 1:N)
    #        pairwise_BLUPS <- get_predictions_2(df, pairwise_indx, rep(1, cluster_size), rep(1, N), tau=1, N1=N)[1]
    
    tibble("n" = cluster_size, "iter"= iter,
           "b1" = BLUPS[1], "btilde" = pairwise_BLUPS[1]*2/cluster_size)
    
  })})


write_csv(output2, "Data/pairwise.csv")
# output2 %>%
#   ggplot() +
#   geom_histogram(aes(x=btilde - b1), binwidth=0.5) + 
#   theme_bw() + xlab("Difference")
# 
# output %>%
#   filter(n==10) %>%
#   ggplot() +
#   geom_histogram(aes(x=btilde - b1), binwidth=0.5) + 
#   theme_bw() + xlab("Difference")
# ggsave("~/BLUPS_HIST.pdf")
# 
# output2 %>%
#   filter(n==50) %>%
#   ggplot() +
#   geom_point(aes(x=btilde , y=b1), alpha=0.5) + 
#   theme_bw() + xlab("Pairwise") + ylab("Model") 
# ggsave("~/BLUPS_LINES.pdf")
# 
# plot
# 
# 
# rcloud.download.file("~/BLUPS_HIST.pdf")
# rcloud.download.file("~/BLUPS_LINES.pdf")
