## Bayesian ordered probit regression (BOPR)

## necessary packages
library(truncnorm)  # for truncated normal dstn
library(progress)  # for progress bar
library(mvtnorm)  # for multivariate normal dstn

## function parameters
# formula: formula for BOPR
# train_df
# test_df
# N: total number of iterations for the MCMC algorithm
# burn.in: number of iterations for the burn-in period of the MCMC algorithm
# verbose: if TRUE, it will print the progress bar for the MCMC algorithm

bopr = function(formula, train_df, test_df, N, burn.in, verbose = TRUE) {
  
  X_train = model.matrix(Y ~ .,data=train_df)
  row.names(X_train) = NULL
  X_test = model.matrix(Y ~ .,data=test_df)
  row.names(X_test) = NULL
  Y_train = model.response(model.frame(formula, train_df))
  
  stopifnot("the response variable must be a factor" = is.factor(Y_train))
  stopifnot("the response variable must have at least 3 ordered categories" = length(levels(Y_train)) >= 3)
  Y_train = as.numeric(train_df$Y)
  
  ## Initialize output
  
  J = length(unique(Y_train))  # total number of ordered categories
  us_p = matrix(nrow = N, ncol = J - 1)  # there are J - 1 thresholds
  theta_p = matrix(nrow=N, ncol=ncol(X_train))
  
  ## Initial values for Gibbs sampling
  
  theta_p[1,] = rep(1, ncol(X_train))
  us_p[1,] = 0:(J - 2)  # -1 because there are J - 1 thresholds and -1 because it starts from 0
  us_p[,1] = 0  # necessary restriction for us1
  lower_us = c(-Inf, us_p[1,])
  lower = lower_us[Y_train]
  upper_us = c(us_p[1,], Inf)
  upper = upper_us[Y_train]
  Z = rtruncnorm(n = length(Y_train), a = lower, b = upper, mean = X_train %*% theta_p[1,], sd = 1)
  
  
  ## MCMC
  
  pb = progress_bar$new(
    format = "  MCMC [:bar] :percent eta: :eta",
    total = N-1, clear = FALSE, width= 60)
  
  for(i in 2:N) {
    if(verbose) pb$tick()
    
    ## update threshold parameters
    for(j in 2:ncol(us_p)){
      if(j == ncol(us_p)){
        us_p[i,j] = runif(1, min = max(Z[which(Y_train==j)], us_p[i,j-1]), max = min(Z[which(Y_train==j+1)], Inf))
      } else{
        us_p[i,j] = runif(1, min = max(Z[which(Y_train==j)], us_p[i,j-1]), max = min(Z[which(Y_train==j+1)], us_p[i-1,j+1]))
      }
    }
    
    ## sample Z
    lower_us = c(-Inf, us_p[i,])
    lower = lower_us[Y_train]
    upper_us = c(us_p[i,], Inf)
    upper = upper_us[Y_train]
    Z = rtruncnorm(n = length(Y_train), a = lower, b = upper, mean = X_train %*% theta_p[i-1,], sd = 1)
    
    ## update theta
    theta_hat = solve(t(X_train) %*% X_train) %*% (t(X_train) %*% Z)
    theta_sigma = solve(t(X_train) %*% X_train)
    theta_p[i,] = rmvnorm(n=1, mean=theta_hat, sigma=theta_sigma)
    
  }
  
  us_p = us_p[(burn.in+1):N,]
  theta_p = theta_p[(burn.in+1):N,]
  
  ## Initialize a matrix to store the predicted values for each MCMC sample
  fx_train = matrix(nrow = (N - burn.in), ncol = length(Y_train))
  fx_test = matrix(nrow = (N - burn.in), ncol = nrow(test_df))
  
  ## Calculate the predicted values f(x) for the i-th sample
  for (i in 1:(N - burn.in)) {  # Loop through each MCMC sample
    fx_train[i,] = X_train %*% theta_p[i,]
    fx_test[i,] = X_test %*% theta_p[i,]
  }
  
  ## Predicting probabilities for each ordered category
  
  p_train = list()  # for train
  for(j in 1:J){
    if(j == 1){p_train[[j]] = pnorm(us_p[,j] - fx_train)}
    else if(j == J){p_train[[j]] = 1 - pnorm(us_p[,j-1] - fx_train)}
    else{p_train[[j]] = pnorm(us_p[,j] - fx_train) - pnorm(us_p[,j-1] - fx_train)}
  }
  train_probs = colMeans(p_train[[1]])  # predicted probabilities
  for(j in 2:J){
    train_probs = cbind(train_probs, colMeans(p_train[[j]]))
  }
  colnames(train_probs) = c(1:J)
  train_preds = factor(apply(train_probs, 1, which.max), levels=c(1:J))  # predicted classes which have the maximum prob. for each obs.
  
  p_test = list()  # for test
  for(j in 1:J){
    if(j == 1){p_test[[j]] = pnorm(us_p[,j] - fx_test)}
    else if(j == J){p_test[[j]] = 1 - pnorm(us_p[,j-1] - fx_test)}
    else{p_test[[j]] = pnorm(us_p[,j] - fx_test) - pnorm(us_p[,j-1] - fx_test)}
  }
  test_probs = colMeans(p_test[[1]])  # predicted probabilities
  for(j in 2:J){
    test_probs = cbind(test_probs, colMeans(p_test[[j]]))
  }
  colnames(test_probs) = c(1:J)
  test_preds = factor(apply(test_probs, 1, which.max), levels=c(1:J))  # predicted classes which have the maximum prob. for each obs.
  
  out = list(us_p = us_p,
             theta_p = theta_p,
             p_train = p_train, 
             p_test = p_test, 
             fx_train = fx_train,
             fx_test = fx_test,
             train_probs = train_probs,
             train_preds = train_preds,
             test_probs = test_probs,
             test_preds = test_preds,
             formula = formula)
  
  class(out) = "bopr"
  return(out)
  
}

## Outputs:
# us_p: MCMC samples (after burn-in) from the posterior distributions of the threshold parameters
# theta_p: MCMC samples from the posterior distributions of the regression coefficients
# p_train: MCMC samples of the probabilities for each ordered category for the train set
# p_test: MCMC samples of the probabilities for each ordered category for the test set
# fx_train: MCMC samples of the fitted f(x) for the train set
# fx_test: MCMC samples of the fitted f(x) for the test set
# train_probs: Column wise averages of p_train (probability predictions)
# train_preds: Class predictions for the train set
# test_probs: Column wise averages of p_test (probability predictions)
# test_preds: Class predictions for the test set
# formula: formula used for BOPR

