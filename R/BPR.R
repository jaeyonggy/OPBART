## Bayesian probit regression

## necessary packages
library(truncnorm)  # for truncated normal dstn
library(progress)  # for progress bar
library(mvtnorm)  # for multivariate normal dstn

## function parameters
# formula: formula for BPR
# train_df
# test_df
# N: total number of iterations for the MCMC algorithm
# burn.in: number of iterations for the burn-in period of the MCMC algorithm
# verbose: if TRUE, it will print the progress bar for the MCMC algorithm

bpr = function(formula, train_df, test_df, N, burn.in, verbose = TRUE) {
  
  X_train = model.matrix(Y ~ .,data=train_df)
  row.names(X_train) = NULL
  X_test = model.matrix(Y ~ .,data=test_df)
  row.names(X_test) = NULL
  Y_train = model.response(model.frame(formula, train_df))
  Y_test  = model.response(model.frame(formula, test_df))
  
  stopifnot("the response variable must be a factor" = is.factor(Y_train))
  stopifnot("the response variable must be a factor" = is.factor(Y_test))
  stopifnot("the response variable must be a binary categorical variable" = length(levels(Y_train)) == 2)
  Y_train = as.numeric(train_df$Y) - 1
  Y_test = as.numeric(test_df$Y) - 1
  
  ## Initialize output
  
  theta_p = matrix(nrow=N, ncol=ncol(X_train))
  
  ## Initial values for Gibbs sampling
  
  theta_p[1,] = rep(1, ncol(X_train))
  lower = ifelse(Y_train == 0, -Inf, 0)
  upper = ifelse(Y_train == 0, 0, Inf)
  
  ## MCMC
  
  pb = progress_bar$new(
    format = "  MCMC [:bar] :percent eta: :eta",
    total = N-1, clear = FALSE, width= 60)
  
  for(i in 2:N) {
    if(verbose) pb$tick()
    
    ## sample Z
    Z = rtruncnorm(n = length(Y_train), a = lower, b = upper, mean = X_train %*% theta_p[i-1,], sd = 1)
    
    ## update theta
    theta_hat = solve(t(X_train) %*% X_train) %*% (t(X_train) %*% Z)
    theta_sigma = solve(t(X_train) %*% X_train)
    theta_p[i,] = rmvnorm(n=1, mean=theta_hat, sigma=theta_sigma)
    
  }
  
  theta_p = theta_p[(burn.in+1):N,]
  
  ## Initialize a matrix to store the predicted values for each MCMC sample
  fx_train = matrix(nrow = (N - burn.in), ncol = length(Y_train))
  fx_test = matrix(nrow = (N - burn.in), ncol = length(Y_test))
  
  ## Calculate the predicted values f(x) for the i-th sample
  for (i in 1:(N - burn.in)) {  # Loop through each MCMC sample
    fx_train[i,] = X_train %*% theta_p[i,]
    fx_test[i,] = X_test %*% theta_p[i,]
  }
  
  ## Predicting probabilities for each ordered category
  
  p_train = pnorm(fx_train)  # MCMC samples of P(y=1)
  train_probs = 1 - colMeans(p_train)  # P(y=0)
  train_probs = cbind(train_probs, colMeans(p_train))  # P(y=1)
  colnames(train_probs) = c(0, 1)
  train_preds = factor(apply(train_probs, 1, which.max) - 1, levels=c(0, 1))
  p_test = pnorm(fx_test)
  test_probs = 1 - colMeans(p_test)
  test_probs = cbind(test_probs, colMeans(p_test))
  colnames(test_probs) = c(0, 1)
  test_preds = factor(apply(test_probs, 1, which.max) - 1, levels=c(0, 1))
  
  out = list(theta_p = theta_p,
             p_train = p_train, 
             p_test = p_test, 
             fx_train = fx_train,
             fx_test = fx_test,
             train_probs = train_probs,
             train_preds = train_preds,
             test_probs = test_probs,
             test_preds = test_preds,
             formula = formula)
  
  class(out) = "bpr"
  return(out)
  
}

## Outputs:
# theta_p: MCMC samples from the posterior distributions of the regression coefficients
# p_train: MCMC samples of the probabilities for P(y=1) for the train set
# p_test: MCMC samples of the probabilities for P(y=1) for the test set
# fx_train: MCMC samples of the fitted f(x) for the train set
# fx_test: MCMC samples of the fitted f(x) for the test set
# train_probs: Column wise averages of p_train (probability predictions)
# train_preds: Class predictions for the train set
# test_probs: Column wise averages of p_test (probability predictions)
# test_preds: Class predictions for the test set
# formula: formula used for BPR

