# New file HelpFunctions.R
get_predictions_1 <- function(sdf, pairwise_indx, weights, L_mat, tau, N1, M_vec, cluster_weights=TRUE){
  ## sdf = data frame
  ## N1 = the number of clusters
  ## M_vec = a vector denoting the cluster sizes
  ## tau = random effect sqrt(variance)
  Yp <- sdf$y[pairwise_indx]
  Xp <- model.matrix(~x, sdf)[pairwise_indx, ]
  Zp <- model.matrix(~1, sdf)[pairwise_indx] ## can add z here
  G <- diag(tau, nrow=N1*1) ## sqrt of tau, 
  W_mat <- diag(sqrt(weights[pairwise_indx]*sdf$w_cluster[pairwise_indx]))
  if(!cluster_weights){W_mat <- diag(sqrt(weights[pairwise_indx]))}
  
  cluster_indx <- sdf$g[pairwise_indx]
  
  XL_mat <- matrix(
    sapply(1:(nrow(Xp)/2), function(ind){
      crossprod(Xp[c(2*ind - 1, 2*ind), ], L_mat)
    }),
    ncol=2, byrow = TRUE)
  LY_p <-  matrix(
    sapply(1:(nrow(Xp)/2), function(ind){
      crossprod(Yp[c(2*ind - 1, 2*ind)], L_mat)
    }),
    ncol=1, byrow=TRUE)
  ZL_vec <-  matrix(
    sapply(1:(nrow(Xp)/2), function(ind){ ## fix when it works
      crossprod(Zp[c(2*ind - 1, 2*ind)], 
                W_mat[c(2*ind - 1, 2*ind),c(2*ind - 1, 2*ind)]%*%L_mat) 
      ## adjust when it's a matrix
    }),
    ncol=1, byrow = TRUE) ## change 1 for dimension z
  
  ZL_mat <- matrix(0, nrow=nrow(Xp), ncol=N1*1) ## Need to extend for matrices
  cluster_uni <- unique(cluster_indx)
  for(i in 1:nrow(Xp)){ ## 
    ZL_mat[i, which(cluster_indx[i]==cluster_uni)] <- 
      ZL_vec[i,1]/(M_vec[cluster_indx[i]]-1)
  } ## really shitty code using rank. need to improve
  ## I put the m-1 here, as it is slighty easier.
  
  
  ## This isn't exact but it's pretty close. Can use lme4 output 
  beta_P <- solve(crossprod(XL_mat)) %*% crossprod(XL_mat, LY_p)
  #### This may be wrong
  
  ###### 
  resids <- Yp - Xp%*%beta_P
  Lresi_p <- matrix(
    sapply(1:(nrow(Xp)/2), function(ind){
      crossprod(resids[c(2*ind - 1, 2*ind)],
                W_mat[c(2*ind - 1, 2*ind),c(2*ind - 1, 2*ind)]%*%L_mat)
    }),
    ncol=1, byrow=TRUE)
  ########
  
  bis  <- G %*% crossprod(ZL_mat, Lresi_p)
  ## This is the number of times each observation is included
  
  
}

## We will solve the weighted penalized problem using linear regression
get_predictions_2 <- function(sdf, pairwise_indx, weights1, weights2, L_mat, tau, N1, M_vec){
  
  Yp <- sdf$y[pairwise_indx]
  Xp <- model.matrix(~x, sdf)[pairwise_indx, ]
  Z_vec <- model.matrix(~1, sdf)[pairwise_indx] ## can add z here
  G <- diag(tau, nrow=N1*1) ## sqrt of tau
  W1 <- diag(weights1[pairwise_indx]) ## weights for the residuals
  W2 <- diag(weights2) ## weights for the predictions
  
  cluster_indx <- sdf$g[pairwise_indx]  
  cluster_uni <- unique(cluster_indx)
  
  Z_mat <- matrix(0, nrow=nrow(Xp), ncol=N1*1) ## Need to extend for matrices
  for(i in 1:nrow(Xp)){ ## 
    Z_mat[i, which(cluster_indx[i]==cluster_uni)] <- Z_vec[i]#[i,1]
  } ## really shitty code using rank. need to improve
  
  
  ### COMBINE THE PIECES!!
  Y_vec1 <- G %*% t(Z_mat) %*% W1 %*% Yp  
  ### G^T * Z^T * W_1 * Y
  Y_vec2 <- t(Xp) %*% W1 %*% Yp
  ## X^T W_1 Y
  
  Mat_11 <- G %*% t(Z_mat) %*% W1 %*% Z_mat %*% G + W2
  Mat_12 <- G %*% t(Z_mat) %*% W1 %*% Xp
  Mat_21 <- t(Mat_12)
  Mat_22 <- t(Xp) %*% W1 %*% Xp
  Big_mat <- rbind(cbind(Mat_11, Mat_12), cbind(Mat_21, Mat_22)) 
  
  solve(Big_mat) %*% rbind(Y_vec1, Y_vec2) ## returns(mu, beta)
}

get_pairwise_values <- function(df, pairwise_indx, fit, N, tau=1){
  ### Pairwise fit
  Yp <- df$y[pairwise_indx]
  Xp <- model.matrix(~x, df)[pairwise_indx, ]
  Zp <- model.matrix(~1, df)[pairwise_indx] ## can add z here
  #tau <- 2*fit@theta
  s2 <- fit@devcomp$cmp["sigmaML"]^2
  Xhi_inv <- solve(s2 * 
                     matrix(c(1+tau^2, tau^2, tau^2, 1+tau^2),
                            nrow=2, ncol=2)) 
  G <- diag(tau[1], nrow=N*1)
  
  cluster_indx <- df$g[pairwise_indx]
  
  L_mat <- chol(Xhi_inv)
  
  #beta_P <- solve(crossprod(XL_mat)) %*% crossprod(XL_mat, LY_p)
  beta_P <- summary(fit)$coef[,1]
  ###### 
  resids <- Yp - Xp%*%beta_P
  ZL_vec <-  matrix(sapply(1:(nrow(Xp)/2), function(ind){ 
    ## fix when it works
    crossprod(Zp[c(2*ind - 1, 2*ind)], L_mat) ## adjust when it's a matrix
  }), ncol=1, byrow = TRUE) ## change 1 for dimension z
  
  ZL_mat <- matrix(0, nrow=nrow(Xp), ncol=N*1) ## Need to extend for matrices
  
  cluster_uni <- unique(cluster_indx)
  
  for(i in 1:nrow(Xp)){ ## 
    ZL_mat[i, which(cluster_indx[i]==cluster_uni)] <- 
      ZL_vec[i,1]/(length(filter(df, g==cluster_indx[i])$t)-1)
  }
  
  Lresi_p <- matrix(sapply(1:(nrow(Xp)/2), function(ind){
    crossprod(resids[c(2*ind - 1, 2*ind)], L_mat)
  }), ncol=1, byrow=TRUE)
  
  b1  <- G %*% crossprod(ZL_mat, Lresi_p)    
  b1   
}

get_pairwise_indx <- function(sdf, N, ind_cluster){
  unlist(lapply((1:N)[ind_cluster], function(cluster_id){
    local_data <- filter(sdf, g == cluster_id)
    #pairs <- as.vector(combn(local_data$t, 2))
    pairs <- as.vector(combn(1:length(local_data$t), 2))
    local_data$id[pairs]
  }))
}

get_subset_data <- function(df, N, ind_cluster, ind_obs){
  map_dfr((1:N)[ind_cluster], function(cluster_id){
    filter(df, (g == cluster_id) & (t %in% ind_obs[[cluster_id]])) 
  })
}
