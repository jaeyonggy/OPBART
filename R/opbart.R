## Ordered probit BART regression
## - Built under SoftBART model

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


## function parameters
# formula: formula for OPBART
# train_df
# test_df
# num_tree: number of trees to train with for BART algorithm
# k: hyperparameter for BART
# hypers: a set of hyperparameters that can be given by function Hypers() from SoftBart package
# opts: a set of optimizers that can be given by the function Opts() from SoftBart package
# verbose: if TRUE, it will print the progress bar for the MCMC algorithm

opbart = function(formula, train_df, test_df, num_tree = 20,
                              k = 1, hypers = NULL, opts = NULL, verbose = TRUE) {
  
  ## Get design matricies and groups for categorical
  
  dv = dummyVars(formula, train_df)
  terms = attr(dv$terms, "term.labels")
  group = dummy_assign(dv)
  suppressWarnings({
    X_train = predict(dv, train_df)
    X_test  = predict(dv, test_df)
  })
  Y_train = model.response(model.frame(formula, train_df))
  Y_test  = model.response(model.frame(formula, test_df))
  
  stopifnot("the response variable must be a factor" = is.factor(Y_train))
  stopifnot("the response variable must be a factor" = is.factor(Y_test))
  stopifnot("the response variable must have at least 3 ordered categories" = length(levels(Y_train)) >= 3)
  Y_train = as.numeric(Y_train)  # need -1 if the factor data starts from 0
  Y_test  = as.numeric(Y_test)  # need -1 if the factor data starts from 0
  
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
  
  ## Make forest ----
  
  opbart_forest = MakeForest(hypers, opts, FALSE)
  
  ## Initialize output
  
  N = opts$num_burn + opts$num_save
  burn.in = opts$num_burn
  
  fx_train  = matrix(NA, nrow = N, ncol = length(Y_train))
  fx_test   = matrix(NA, nrow = N, ncol = length(Y_test))
  sigma_mu  = numeric(N)
  varcounts = matrix(NA, nrow = N, ncol = length(terms))
  J = length(unique(Y_train))  # total number of ordered categories
  us_p = matrix(nrow = N, ncol = J - 1)  # there are J - 1 thresholds
  
  
  ## Initialize values for Gibbs sampling
  
  fx_train[1,] = as.numeric(opbart_forest$do_predict(X_train))  # initial value
  us_p[1,] = 0:(J - 2)  # -1 because there are J - 1 thresholds and -1 because it starts from 0
  us_p[,1] = 0  # necessary restriction for us1
  lower_us = c(-Inf, us_p[1,])
  lower = lower_us[Y_train]
  upper_us = c(us_p[1,], Inf)
  upper = upper_us[Y_train]
  Z = rtruncnorm(n = length(Y_train), a = lower, b = upper, mean = fx_train[1,], sd = 1)
  
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
    Z = rtruncnorm(n = length(Y_train), a = lower, b = upper, mean = fx_train[i-1,], sd = 1)
    
    ## update mu
    fx_train[i,] = opbart_forest$do_gibbs(X_train, Z, X_train, 1)
    
    ## saving after burn-in period
    if(i > burn.in){
      fx_test[i,] = opbart_forest$do_predict(X_test)
      sigma_mu[i]   = opbart_forest$get_sigma_mu()
      varcounts[i,] = opbart_forest$get_counts()
    }
  }
  
  
  ## MCMC samples after burn-in period
  
  fx_train  = fx_train[(burn.in+1):N,]
  fx_test   = fx_test[(burn.in+1):N,]
  sigma_mu  = sigma_mu[(burn.in+1):N]
  varcounts = varcounts[(burn.in+1):N,]
  us_p = us_p[(burn.in+1):N,]
  
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
  
  ## Outputs
  
  colnames(varcounts) = terms
  
  out = list(sigma_mu = sigma_mu, 
             var_counts = varcounts, 
             us_p = us_p,
             p_train = p_train, 
             p_test = p_test, 
             fx_train = fx_train,
             fx_test = fx_test,
             train_probs = train_probs,
             train_preds = train_preds,
             test_probs = test_probs,
             test_preds = test_preds,
             formula = formula,
             ecdfs = ecdfs,
             opts = opts,
             forest = opbart_forest,
             dv = dv)
  
  class(out) = "opbart"
  return(out)
  
}

## Outputs:
# sigma_mu: refer to SoftBart package
# var_counts: refer to SoftBart package
# us_p: MCMC samples (after burn-in) from the posterior distributions of the threshold parameters
# p_train: MCMC samples of the probabilities for each ordered category for the train set
# p_test: MCMC samples of the probabilities for each ordered category for the test set
# fx_train: MCMC samples of the fitted f(x) for the train set
# fx_test: MCMC samples of the fitted f(x) for the test set
# train_probs: Column wise averages of p_train (probability predictions)
# train_preds: Class predictions for the train set
# test_probs: Column wise averages of p_test (probability predictions)
# test_preds: Class predictions for the test set
# formula: formula used for OPBART
# ecdfs: refer to SoftBart package
# opts: opts used for OPBART
# forest: refer to SoftBart package
# dv: refer to SoftBart package



