## semiparametric Probit BART (semi-PBART) regression
## - Built under SoftBART model

## necessary packages
library(SoftBart)  # for BART
library(truncnorm)  # for truncated normal dstn
library(progress)  # for progress bar
library(mvtnorm)  # for multivariate normal dstn

## necessary preprocessing function for BART (from SoftBart implementation)
dummy_assign = function(dummy) {
  terms = attr(dummy$terms, "term.labels")
  group = list()
  j     = 0
  for(k in terms) {
    if(k %in% dummy$facVars) {
      group[[k]] = rep(j, length(dummy$lvls[[k]]))
    } else {
      group[[k]] = j
    }
    j = j + 1
  }
  return(do.call(c, group))
}

## WARNING: marginal effects can be calculated for continuous covariates only at the moment

## function parameters
# formula: formula for the nonparametric part of semi-PBART
# linear_formula: a formula for the semiparametric part of semi-PBART
# one_unit: you can specify the one unit to a desired unit such as 0.1
# train_df
# test_df
# num_tree: number of trees to train with for BART algorithm
# k: hyperparameter for BART
# hypers: a set of hyperparameters that can be given by function Hypers() from SoftBart package
# opts: a set of optimizers that can be given by the function Opts() from SoftBart package
# remove_intercept: if TRUE, it will remove the intercept term in the linear_formula
# verbose: if TRUE, it will print the progress bar for the MCMC algorithm

smpbart = function(formula, linear_formula, train_df, test_df, one_unit = 1, num_tree = 20,
                   k = 1, hypers = NULL, opts = NULL, remove_intercept = TRUE, verbose = TRUE) {
  
  ## Get design matricies and groups for categorical
  
  dv = dummyVars(formula, train_df)
  terms = attr(dv$terms, "term.labels")
  group = dummy_assign(dv)
  suppressWarnings({
    X_train = predict(dv, train_df)
    X_test  = predict(dv, test_df)
  })
  Y_train = model.response(model.frame(formula, train_df))
  
  stopifnot("the response variable must be a factor" = is.factor(Y_train))
  stopifnot("the response variable must be a binary categorical variable" = length(levels(Y_train)) == 2)
  Y_train = as.numeric(Y_train) - 1
  
  W_train = model.matrix(linear_formula, train_df)
  W_test = model.matrix(linear_formula, test_df)
  
  remove_intercept =TRUE
  if(remove_intercept) {
    idx = which(colnames(W_train) != "(Intercept)")
    W_train = W_train[,idx]
    W_test = W_test[,idx]    
  }
  
  W_train = as.matrix(W_train)  # when there's only one covariate it stops being a matrix so we force it to be a matrix
  W_test = as.matrix(W_test)
  
  ## Set up hypers
  
  if(is.null(hypers)) {
    hypers = Hypers(X = X_train, Y = Y_train)
  }
  hypers$sigma_mu = 3 / k / sqrt(num_tree)
  hypers$sigma = 1
  hypers$sigma_hat = 1
  hypers$num_tree = num_tree
  hypers$group = group
  
  ## Set up opts
  
  if(is.null(opts)) {
    opts = Opts()
  }
  opts$update_sigma = FALSE
  opts$num_print = 2147483647
  
  ## Normalize!
  
  make_01_norm = function(x) {
    a = min(x)
    b = max(x)
    return(function(y) (y - a) / (b - a))
  }
  
  ecdfs   = list()
  for(i in 1:ncol(X_train)) {
    ecdfs[[i]] = ecdf(X_train[,i])
    if(length(unique(X_train[,i])) == 1) ecdfs[[i]] = identity
    if(length(unique(X_train[,i])) == 2) ecdfs[[i]] = make_01_norm(X_train[,i])
  }
  for(i in 1:ncol(X_train)) {
    X_train[,i] = ecdfs[[i]](X_train[,i])
    X_test[,i] = ecdfs[[i]](X_test[,i])
  }
  
  ## Resetting the indexes of the covariate matrix
  row.names(X_train) = NULL
  row.names(X_test) = NULL
  row.names(W_train) = NULL
  row.names(W_test) = NULL
  
  ## Make forest ----
  
  smpbart_forest = MakeForest(hypers, opts, FALSE)
  
  ## Initialize output
  
  N = opts$num_burn + opts$num_save
  burn.in = opts$num_burn
  
  fx_train  = matrix(NA, nrow = N, ncol = length(Y_train))
  fx_test   = matrix(NA, nrow = N, ncol = nrow(test_df))
  sigma_mu  = numeric(N)
  varcounts = matrix(NA, nrow = N, ncol = length(terms))
  theta_p = matrix(nrow=N, ncol=ncol(W_train))
  
  ## Initialize values for Gibbs sampling
  
  fx_train[1,] = as.numeric(smpbart_forest$do_predict(X_train))  # initial value
  theta_p[1,] = rep(1, ncol(W_train))
  lower = ifelse(Y_train == 0, -Inf, 0)
  upper = ifelse(Y_train == 0, 0, Inf)
  
  ## MCMC
  
  pb = progress_bar$new(
    format = "  MCMC [:bar] :percent eta: :eta",
    total = N-1, clear = FALSE, width= 60)
  
  for(i in 2:N) {
    if(verbose) pb$tick()
    
    ## sample Z
    Z = rtruncnorm(n = length(Y_train), a = lower, b = upper, mean = fx_train[i-1,] + W_train %*% theta_p[i-1,], sd = 1)
    
    ## update theta
    Z_star = Z - fx_train[i-1,]
    theta_hat = solve(t(W_train) %*% W_train) %*% (t(W_train) %*% Z_star)
    theta_sigma = solve(t(W_train) %*% W_train)
    theta_p[i,] = rmvnorm(n=1, mean=theta_hat, sigma=theta_sigma)
    
    ## update f(x)
    Z_tilde = Z - W_train %*% theta_p[i,]
    fx_train[i,] = smpbart_forest$do_gibbs(X_train, Z_tilde, X_train, 1)
    
    ## saving after burn-in period
    if(i > burn.in){
      fx_test[i,] = smpbart_forest$do_predict(X_test)
      sigma_mu[i] = smpbart_forest$get_sigma_mu()
      varcounts[i,] = smpbart_forest$get_counts()
    }
  }
  
  
  ## MCMC samples after burn-in period
  
  fx_train  = fx_train[(burn.in+1):N,]
  fx_test   = fx_test[(burn.in+1):N,]
  sigma_mu  = sigma_mu[(burn.in+1):N]
  varcounts = varcounts[(burn.in+1):N,]
  theta_p = theta_p[(burn.in+1):N,]
  theta_p = as.matrix(theta_p)  # if there is only one covariate it becomes a vector so we force it to be a matrix
  
  ## Initialize a matrix to store the h(w) for each MCMC sample
  hw_train = matrix(nrow = (N - burn.in), ncol = length(Y_train))
  hw_test = matrix(nrow = (N - burn.in), ncol = nrow(test_df))
  
  ## Calculate the predicted values h(w) for the ith sample
  for (i in 1:(N - burn.in)) {  # Loop through each MCMC sample
    hw_train[i,] = W_train %*% theta_p[i,]
    hw_test[i,] = W_test %*% theta_p[i,]
  }
  
  ## Predicting probabilities
  
  p_train = pnorm(fx_train + hw_train)  # MCMC samples of P(y=1)
  train_probs = 1 - colMeans(p_train)  # P(y=0)
  train_probs = cbind(train_probs, colMeans(p_train))  # P(y=1)
  colnames(train_probs) = c(0, 1)
  train_preds = factor(apply(train_probs, 1, which.max) - 1, levels=c(0, 1))
  p_test = pnorm(fx_test + hw_test)
  test_probs = 1 - colMeans(p_test)
  test_probs = cbind(test_probs, colMeans(p_test))
  colnames(test_probs) = c(0, 1)
  test_preds = factor(apply(test_probs, 1, which.max) - 1, levels=c(0, 1))
  
  ## Calculating marginal effects (ME) for each w
  # for train set
  fx_train_mean = rowMeans(fx_train)  # when calculating MEs all remaining variables assume their respective average values
  hw_train_mean = rowMeans(hw_train)
  
  ME_train_list = list()  # list to store the ME matrix for each theta
  
  for(w in 1:ncol(theta_p)){
    ME_train = -pnorm(fx_train_mean + hw_train_mean)*(theta_p[,w]*one_unit)
    ME_train = cbind(ME_train, pnorm(fx_train_mean + hw_train_mean)*(theta_p[,w]*one_unit))
    ME_train_list[[w]] = ME_train
  }
  
  # for test set
  fx_test_mean = rowMeans(fx_test)  # when calculating MEs all remaining variables assume their respective average values
  hw_test_mean = rowMeans(hw_test)
  
  ME_test_list = list()  # list to store the ME matrix for each theta
  
  for(w in 1:ncol(theta_p)){
    ME_test = -pnorm(fx_test_mean + hw_test_mean)*(theta_p[,w]*one_unit)
    ME_test = cbind(ME_test, pnorm(fx_test_mean + hw_test_mean)*(theta_p[,w]*one_unit))
    ME_test_list[[w]] = ME_test
  }
  
  
  ## Outputs
  
  colnames(varcounts) = terms
  
  out = list(sigma_mu = sigma_mu, 
             var_counts = varcounts, 
             theta_p = theta_p,
             p_train = p_train, 
             p_test = p_test, 
             fx_train = fx_train,
             fx_test = fx_test,
             hw_train = hw_train,
             hw_test = hw_test,
             train_probs = train_probs,
             train_preds = train_preds,
             test_probs = test_probs,
             test_preds = test_preds,
             ME_train_list = ME_train_list,
             ME_test_list = ME_test_list,
             formula = formula,
             linear_formula = linear_formula,
             ecdfs = ecdfs,
             opts = opts,
             forest = smpbart_forest,
             dv = dv)
  
  class(out) = "smpbart"
  return(out)
  
}

## Outputs:
# sigma_mu: refer to SoftBart package
# var_counts: refer to SoftBart package
# theta_p: MCMC samples from the posterior distributions of the regression coefficients
# p_train: MCMC samples of the probabilities for P(y=1) for the train set
# p_test: MCMC samples of the probabilities for P(y=1) for the test set
# fx_train: MCMC samples of the fitted f(x) for the train set
# fx_test: MCMC samples of the fitted f(x) for the test set
# hw_train: MCMC samples of the fitted h(w) for the train set
# hw_test: MCMC samples of the fitted h(w) for the test set
# train_probs: Column wise averages of p_train (probability predictions)
# train_preds: Class predictions for the train set
# test_probs: Column wise averages of p_test (probability predictions)
# test_preds: Class predictions for the test set
# ME_train_list: Marginal effects of w covariates from train set
# ME_test_list: Marginal effects of w covariates from test set
# formula: formula used for semi-PBART
# linear_formula: linear_formula used for semi-PBART
# ecdfs: refer to SoftBart package
# opts: opts used for semi-PBART
# forest: refer to SoftBart package
# dv: refer to SoftBart package
