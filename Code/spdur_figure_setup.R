spdur_roc_plot <- function(ModelResults, modelname){
  require(plotROC)
  require(pROC)
  require(dplyr)
  library(purrr)
  library(tibble)
  library(magrittr)
  Y_pred =lapply(ModelResults, function(x) FUN = predict(x, type = "response"))
  
  Y_obs = lapply(ModelResults, function(x) FUN = x$Y[,"fail"])
  
  roc = Map(function(x, y) roc(x,y), Y_obs, Y_pred)
  
  names(roc) <- modelname
  
  roc_df = lapply(roc, function(x) FUN= data.frame(
    plotx = x$specificities,
    ploty = rev(x$sensitivities),
    name = paste("AUC =",
                 sprintf("%.3f",x$auc)))) %>%
    map_df(., rbind, .id="modelname")
  return(roc_df)
}

# set the hazard function
hazard <- function(ti, lambda, cure, alpha, out, dist) {
  
  # minimum value for P, to avoid divide by 0 errors
  p_min <- 1e-16
  
  lambda <- as.vector(lambda)
  cure   <- as.vector(cure)
  
  ht <- vector("numeric", length(ti))
  if (dist=="weibull") {
    
    st       <- exp(-(lambda * ti)^alpha)
    cure.t   <- cure / pmax(p_min, (st + cure * (1 - st))) # pmax to avoid dividing it by 0
    atrisk.t <- 1 - cure.t
    ft       <- lambda * alpha * (lambda * ti)^(alpha-1) * exp(-(lambda * ti)^alpha)
    
    ht       <- atrisk.t * ft / pmax(p_min, (cure.t + atrisk.t * st))
    
  } else if (dist=="loglog") {
    
    st       <- 1/(1+(lambda * ti)^alpha)
    cure.t   <- cure / pmax(p_min, (st + cure * (1 - st))) # pmax to avoid dividing it by 0
    atrisk.t <- 1 - cure.t
    ft       <- (lambda * alpha * (lambda * ti)^(alpha-1)) / ((1 + (lambda * ti)^alpha)^2)  
    
    ht       <- atrisk.t * ft / pmax(p_min, (cure.t + atrisk.t * st))
    
  } else {
    
    stop(paste0("Unrecognized distribution: ", dist))
    
  }
  
  ht
}

atrisk <- function(ti, lambda, cure, alpha, out, dist) {
  
  # minimum value for P, to avoid divide by 0 errors
  p_min <- 1e-16
  
  lambda <- as.vector(lambda)
  cure   <- as.vector(cure)
  
  ht <- vector("numeric", length(ti))
  if (dist=="weibull") {
    
    st       <- exp(-(lambda * ti)^alpha)
    cure.t   <- cure / pmax(p_min, (st + cure * (1 - st))) # pmax to avoid dividing it by 0
    atrisk.t <- 1 - cure.t
    ft       <- lambda * alpha * (lambda * ti)^(alpha-1) * exp(-(lambda * ti)^alpha)
    
    ht       <- atrisk.t * ft / pmax(p_min, (cure.t + atrisk.t * st))
    
  } else if (dist=="loglog") {
    
    st       <- 1/(1+(lambda * ti)^alpha)
    cure.t   <- cure / pmax(p_min, (st + cure * (1 - st))) # pmax to avoid dividing it by 0
    atrisk.t <- 1 - cure.t
    ft       <- (lambda * alpha * (lambda * ti)^(alpha-1)) / ((1 + (lambda * ti)^alpha)^2)  
    
    ht       <- atrisk.t * ft / pmax(p_min, (cure.t + atrisk.t * st))
    
  } else {
    
    stop(paste0("Unrecognized distribution: ", dist))
    
  }
  
  atrisk.t
}


### cured t
cure <- function(ti, lambda, cure, alpha, out, dist) {
  
  # minimum value for P, to avoid divide by 0 errors
  p_min <- 1e-16
  
  lambda <- as.vector(lambda)
  cure   <- as.vector(cure)
  
  ht <- vector("numeric", length(ti))
  if (dist=="weibull") {
    
    st       <- exp(-(lambda * ti)^alpha)
    cure.t   <- cure / pmax(p_min, (st + cure * (1 - st))) # pmax to avoid dividing it by 0
    atrisk.t <- 1 - cure.t
    ft       <- lambda * alpha * (lambda * ti)^(alpha-1) * exp(-(lambda * ti)^alpha)
    
    ht       <- atrisk.t * ft / pmax(p_min, (cure.t + atrisk.t * st))
    
  } else if (dist=="loglog") {
    
    st       <- 1/(1+(lambda * ti)^alpha)
    cure.t   <- cure / pmax(p_min, (st + cure * (1 - st))) # pmax to avoid dividing it by 0
    atrisk.t <- 1 - cure.t
    ft       <- (lambda * alpha * (lambda * ti)^(alpha-1)) / ((1 + (lambda * ti)^alpha)^2)  
    
    ht       <- atrisk.t * ft / pmax(p_min, (cure.t + atrisk.t * st))
    
  } else {
    
    stop(paste0("Unrecognized distribution: ", dist))
    
  }
  
  cure.t
}


st <- function(ti, lambda, cure, alpha, out, dist) {
  
  # minimum value for P, to avoid divide by 0 errors
  p_min <- 1e-16
  
  lambda <- as.vector(lambda)
  cure   <- as.vector(cure)
  
  ht <- vector("numeric", length(ti))
  if (dist=="weibull") {
    
    st       <- exp(-(lambda * ti)^alpha)
    cure.t   <- cure / pmax(p_min, (st + cure * (1 - st))) # pmax to avoid dividing it by 0
    atrisk.t <- 1 - cure.t
    ft       <- lambda * alpha * (lambda * ti)^(alpha-1) * exp(-(lambda * ti)^alpha)
    
    ht       <- atrisk.t * ft / pmax(p_min, (cure.t + atrisk.t * st))
    
  } else if (dist=="loglog") {
    
    st       <- 1/(1+(lambda * ti)^alpha)
    cure.t   <- cure / pmax(p_min, (st + cure * (1 - st))) # pmax to avoid dividing it by 0
    atrisk.t <- 1 - cure.t
    ft       <- (lambda * alpha * (lambda * ti)^(alpha-1)) / ((1 + (lambda * ti)^alpha)^2)  
    
    ht       <- atrisk.t * ft / pmax(p_min, (cure.t + atrisk.t * st))
    
  } else {
    
    stop(paste0("Unrecognized distribution: ", dist))
    
  }
  
  st
}

hazard_plot_df <- function(x, t = NULL, ci=TRUE, n=1000, xvals=NULL, zvals=NULL,
                           xpred_name, xpred_value, scenario,...) {
  
  # Set t vector if needed to 1.2 * max observed duration; lower limit is 1
  if (is.null(t)) {
    max_t <- round(max(x$Y[x$Y[, "last"]==1, "duration"]) * 1.2)
    t     <- seq(1, max_t, length.out=100)
  } 
  
  # Extract covariate matrices
  dur.dat  <- x$mf.dur
  risk.dat <- x$mf.risk 
  X <- model.matrix(attr(x$mf.dur, 'terms'), data=x$mf.dur)
  Z <- model.matrix(attr(x$mf.risk, 'terms'), data=x$mf.risk)
  
  # Extract coefficient point estimates
  beta  <- coef(x, model = "duration")
  gamma <- coef(x, model = "risk")
  alpha <- coef(x, model = "distr")
  alpha <- exp(-alpha)
  beta_vcv  <- vcov(x, "duration")
  gamma_vcv <- vcov(x, "risk")
  alpha_vcv <- vcov(x, "distr")
  
  if (is.null(xvals)) {
    X[, xpred_name] <- xpred_value
    X_vals <- apply(X, 2, mean)		
  } else if (!length(xvals)==ncol(X) && length(xvals)-ncol(X)==-1) {
    stop("Incorrect length for xvals, did you forget 1 for intercept term?")			
  } else if (!length(xvals)==ncol(X)) {
    stop("Incorrect length for xvals")
  } else {
    X_vals <- xvals
  }
  
  if (is.null(zvals)) {
    Z[, xpred_name] <- xpred_value
    Z_vals <- apply(Z, 2, mean)		
  } else if (!length(zvals)==ncol(Z) && length(zvals)-ncol(Z)==-1) {
    stop("Incorrect length for zvals, did you forget 1 for intercept term?")			
  } else if (!length(zvals)==ncol(Z)) {
    stop("Incorrect length for zvals")
  } else {
    Z_vals <- zvals
  }
  
  # Calculate hazard using point estimates only
  lambda <- exp(-X_vals %*% beta)
  cure   <- 1 - plogis(Z_vals %*% gamma)
  
  ht <- hazard(ti = t, lambda = lambda, cure = cure, alpha = alpha,
               out = NULL, dist = x$distr)
  
  if (ci==TRUE) {
    set.seed(1234)
    Coef_smpl <- MASS::mvrnorm(n = n, mu = coef(x, "full"), Sigma = vcov(x, "full"))
    
    b_idx <- 1:x$n.terms$duration
    g_idx <- (max(b_idx) + 1):(max(b_idx) + x$n.terms$risk)
    a_idx <- (max(g_idx) + 1):ncol(Coef_smpl)
    
    Beta  <- Coef_smpl[, b_idx]
    Gamma <- Coef_smpl[, g_idx]
    A     <- Coef_smpl[, a_idx]
    Alpha <- exp(-A)
    
    lambda <- exp(-tcrossprod(X_vals, Beta))
    cure   <- 1 - plogis(tcrossprod(Z_vals, Gamma))
    
    sims <- matrix(nrow = length(t), ncol = n)
    hmat <- matrix(nrow = length(t), ncol = 3)
    
    for (i in 1:n) {
      sims[, i] <- hazard(ti = t, lambda = lambda[i], cure = cure[i], 
                          alpha = Alpha[i], out = NULL, dist = x$distr)
    }
    
    hmat[, 1] <- ht
    hmat[, 2] <- apply(sims, 1, quantile, probs = 0.05)
    hmat[, 3] <- apply(sims, 1, quantile, probs = 0.95)
    
    ## use ggplot
    sur_df <- data.frame(Time = t,
                         surv = hmat[, 1],
                         lo95 = hmat[, 2],
                         up95 = hmat[, 3],
                         scenario = scenario)
    return(sur_df)
  } else {
    # Plot without CIs
    
    sur_df <- data.frame(Time = t, 
                         surv=ht,
                         scenario = scenario)
    return(sur_df)
  }
  
  invisible(NULL)
}


### plot predicted
hazardpred_plot_df <- function(x, t = NULL, ci=TRUE, n=1000, xvals=NULL, zvals=NULL,...) {
  
  # Set t vector if needed to 1.2 * max observed duration; lower limit is 1
  if (is.null(t)) {
    max_t <- round(max(x$Y[x$Y[, "last"]==1, "duration"]) * 1.2)
    t     <- seq(1, max_t, length.out=100)
  } 
  
  # Extract covariate matrices
  dur.dat  <- x$mf.dur
  risk.dat <- x$mf.risk 
  X <- model.matrix(attr(x$mf.dur, 'terms'), data=x$mf.dur)
  Z <- model.matrix(attr(x$mf.risk, 'terms'), data=x$mf.risk)
  
  # Extract coefficient point estimates
  beta  <- coef(x, model = "duration")
  gamma <- coef(x, model = "risk")
  alpha <- coef(x, model = "distr")
  alpha <- exp(-alpha)
  beta_vcv  <- vcov(x, "duration")
  gamma_vcv <- vcov(x, "risk")
  alpha_vcv <- vcov(x, "distr")
  
  if (is.null(xvals)) {
    X_vals <- apply(X, 2, mean)		
  } else if (!length(xvals)==ncol(X) && length(xvals)-ncol(X)==-1) {
    stop("Incorrect length for xvals, did you forget 1 for intercept term?")			
  } else if (!length(xvals)==ncol(X)) {
    stop("Incorrect length for xvals")
  } else {
    X_vals <- xvals
  }
  
  if (is.null(zvals)) {
    Z_vals <- apply(Z, 2, mean)		
  } else if (!length(zvals)==ncol(Z) && length(zvals)-ncol(Z)==-1) {
    stop("Incorrect length for zvals, did you forget 1 for intercept term?")			
  } else if (!length(zvals)==ncol(Z)) {
    stop("Incorrect length for zvals")
  } else {
    Z_vals <- zvals
  }
  
  # Calculate hazard using point estimates only
  lambda <- exp(-X_vals %*% beta)
  cure   <- 1 - plogis(Z_vals %*% gamma)
  
  ht <- hazard(ti = t, lambda = lambda, cure = cure, alpha = alpha,
               out = NULL, dist = x$distr)
  
  if (ci==TRUE) {
    set.seed(1234)
    Coef_smpl <- MASS::mvrnorm(n = n, mu = coef(x, "full"), Sigma = vcov(x, "full"))
    
    b_idx <- 1:x$n.terms$duration
    g_idx <- (max(b_idx) + 1):(max(b_idx) + x$n.terms$risk)
    a_idx <- (max(g_idx) + 1):ncol(Coef_smpl)
    
    Beta  <- Coef_smpl[, b_idx]
    Gamma <- Coef_smpl[, g_idx]
    A     <- Coef_smpl[, a_idx]
    Alpha <- exp(-A)
    
    lambda <- exp(-tcrossprod(X_vals, Beta))
    cure   <- 1 - plogis(tcrossprod(Z_vals, Gamma))
    
    sims <- matrix(nrow = length(t), ncol = n)
    hmat <- matrix(nrow = length(t), ncol = 3)
    
    for (i in 1:n) {
      sims[, i] <- hazard(ti = t, lambda = lambda[i], cure = cure[i], 
                          alpha = Alpha[i], out = NULL, dist = x$distr)
    }
    
    hmat[, 1] <- ht
    hmat[, 2] <- apply(sims, 1, quantile, probs = 0.05)
    hmat[, 3] <- apply(sims, 1, quantile, probs = 0.95)
    
    ## use ggplot
    sur_df <- data.frame(Time = t,
                         surv = hmat[, 1],
                         lo95 = hmat[, 2],
                         up95 = hmat[, 3])
    return(sur_df)
  } else {
    # Plot without CIs
    
    sur_df <- data.frame(Time = t, 
                         surv=ht)
    return(sur_df)
  }
  
  invisible(NULL)
}

## conditional risk
risk_plot_df <- function(x, t = NULL, ci=TRUE, n=1000, xvals=NULL, zvals=NULL,
                           xpred_name, xpred_value, scenario,...) {
  
  # Set t vector if needed to 1.2 * max observed duration; lower limit is 1
  if (is.null(t)) {
    max_t <- round(max(x$Y[x$Y[, "last"]==1, "duration"]) * 1.2)
    t     <- seq(1, max_t, length.out=100)
  } 
  
  # Extract covariate matrices
  dur.dat  <- x$mf.dur
  risk.dat <- x$mf.risk 
  X <- model.matrix(attr(x$mf.dur, 'terms'), data=x$mf.dur)
  Z <- model.matrix(attr(x$mf.risk, 'terms'), data=x$mf.risk)
  
  # Extract coefficient point estimates
  beta  <- coef(x, model = "duration")
  gamma <- coef(x, model = "risk")
  alpha <- coef(x, model = "distr")
  alpha <- exp(-alpha)
  beta_vcv  <- vcov(x, "duration")
  gamma_vcv <- vcov(x, "risk")
  alpha_vcv <- vcov(x, "distr")
  
  if (is.null(xvals)) {
    X[, xpred_name] <- xpred_value
    X_vals <- apply(X, 2, mean)		
  } else if (!length(xvals)==ncol(X) && length(xvals)-ncol(X)==-1) {
    stop("Incorrect length for xvals, did you forget 1 for intercept term?")			
  } else if (!length(xvals)==ncol(X)) {
    stop("Incorrect length for xvals")
  } else {
    X_vals <- xvals
  }
  
  if (is.null(zvals)) {
    Z[, xpred_name] <- xpred_value
    Z_vals <- apply(Z, 2, mean)		
  } else if (!length(zvals)==ncol(Z) && length(zvals)-ncol(Z)==-1) {
    stop("Incorrect length for zvals, did you forget 1 for intercept term?")			
  } else if (!length(zvals)==ncol(Z)) {
    stop("Incorrect length for zvals")
  } else {
    Z_vals <- zvals
  }
  
  # Calculate hazard using point estimates only
  lambda <- exp(-X_vals %*% beta)
  cure   <- 1 - plogis(Z_vals %*% gamma)
  
  atrisk.t <- atrisk(ti = t, lambda = lambda, cure = cure, alpha = alpha,
               out = NULL, dist = x$distr)
  
  if (ci==TRUE) {
    set.seed(1234)
    Coef_smpl <- MASS::mvrnorm(n = n, mu = coef(x, "full"), Sigma = vcov(x, "full"))
    
    b_idx <- 1:x$n.terms$duration
    g_idx <- (max(b_idx) + 1):(max(b_idx) + x$n.terms$risk)
    a_idx <- (max(g_idx) + 1):ncol(Coef_smpl)
    
    Beta  <- Coef_smpl[, b_idx]
    Gamma <- Coef_smpl[, g_idx]
    A     <- Coef_smpl[, a_idx]
    Alpha <- exp(-A)
    
    lambda <- exp(-tcrossprod(X_vals, Beta))
    cure   <- 1 - plogis(tcrossprod(Z_vals, Gamma))
    
    sims <- matrix(nrow = length(t), ncol = n)
    hmat <- matrix(nrow = length(t), ncol = 3)
    
    for (i in 1:n) {
      sims[, i] <- atrisk(ti = t, lambda = lambda[i], cure = cure[i], 
                          alpha = Alpha[i], out = NULL, dist = x$distr)
    }
    
    hmat[, 1] <- atrisk.t
    hmat[, 2] <- apply(sims, 1, quantile, probs = 0.05)
    hmat[, 3] <- apply(sims, 1, quantile, probs = 0.95)
    
    ## use ggplot
    sur_df <- data.frame(Time = t,
                         atrisk = hmat[, 1],
                         lo95 = hmat[, 2],
                         up95 = hmat[, 3],
                         scenario = scenario)
    return(sur_df)
  } else {
    # Plot without CIs
    
    sur_df <- data.frame(Time = t, 
                         atrisk=atrisk.t,
                         scenario = scenario)
    return(sur_df)
  }
  
  invisible(NULL)
}

## surivial estimates
st_plot_df <- function(x, t = NULL, ci=TRUE, n=1000, xvals=NULL, zvals=NULL,
                         xpred_name, xpred_value, scenario,...) {
  
  # Set t vector if needed to 1.2 * max observed duration; lower limit is 1
  if (is.null(t)) {
    max_t <- round(max(x$Y[x$Y[, "last"]==1, "duration"]) * 1.2)
    t     <- seq(1, max_t, length.out=100)
  } 
  
  # Extract covariate matrices
  dur.dat  <- x$mf.dur
  risk.dat <- x$mf.risk 
  X <- model.matrix(attr(x$mf.dur, 'terms'), data=x$mf.dur)
  Z <- model.matrix(attr(x$mf.risk, 'terms'), data=x$mf.risk)
  
  # Extract coefficient point estimates
  beta  <- coef(x, model = "duration")
  gamma <- coef(x, model = "risk")
  alpha <- coef(x, model = "distr")
  alpha <- exp(-alpha)
  beta_vcv  <- vcov(x, "duration")
  gamma_vcv <- vcov(x, "risk")
  alpha_vcv <- vcov(x, "distr")
  
  if (is.null(xvals)) {
    X[, xpred_name] <- xpred_value
    X_vals <- apply(X, 2, mean)		
  } else if (!length(xvals)==ncol(X) && length(xvals)-ncol(X)==-1) {
    stop("Incorrect length for xvals, did you forget 1 for intercept term?")			
  } else if (!length(xvals)==ncol(X)) {
    stop("Incorrect length for xvals")
  } else {
    X_vals <- xvals
  }
  
  if (is.null(zvals)) {
    Z[, xpred_name] <- xpred_value
    Z_vals <- apply(Z, 2, mean)		
  } else if (!length(zvals)==ncol(Z) && length(zvals)-ncol(Z)==-1) {
    stop("Incorrect length for zvals, did you forget 1 for intercept term?")			
  } else if (!length(zvals)==ncol(Z)) {
    stop("Incorrect length for zvals")
  } else {
    Z_vals <- zvals
  }
  
  # Calculate hazard using point estimates only
  lambda <- exp(-X_vals %*% beta)
  cure   <- 1 - plogis(Z_vals %*% gamma)
  
  st.t <- st(ti = t, lambda = lambda, cure = cure, alpha = alpha,
                     out = NULL, dist = x$distr)
  
  if (ci==TRUE) {
    set.seed(1234)
    Coef_smpl <- MASS::mvrnorm(n = n, mu = coef(x, "full"), Sigma = vcov(x, "full"))
    
    b_idx <- 1:x$n.terms$duration
    g_idx <- (max(b_idx) + 1):(max(b_idx) + x$n.terms$risk)
    a_idx <- (max(g_idx) + 1):ncol(Coef_smpl)
    
    Beta  <- Coef_smpl[, b_idx]
    Gamma <- Coef_smpl[, g_idx]
    A     <- Coef_smpl[, a_idx]
    Alpha <- exp(-A)
    
    lambda <- exp(-tcrossprod(X_vals, Beta))
    cure   <- 1 - plogis(tcrossprod(Z_vals, Gamma))
    
    sims <- matrix(nrow = length(t), ncol = n)
    hmat <- matrix(nrow = length(t), ncol = 3)
    
    for (i in 1:n) {
      sims[, i] <- st(ti = t, lambda = lambda[i], cure = cure[i], 
                          alpha = Alpha[i], out = NULL, dist = x$distr)
    }
    
    hmat[, 1] <- st.t
    hmat[, 2] <- apply(sims, 1, quantile, probs = 0.05)
    hmat[, 3] <- apply(sims, 1, quantile, probs = 0.95)
    
    ## use ggplot
    sur_df <- data.frame(Time = t,
                         st = hmat[, 1],
                         lo95 = hmat[, 2],
                         up95 = hmat[, 3],
                         scenario = scenario)
    return(sur_df)
  } else {
    # Plot without CIs
    
    sur_df <- data.frame(Time = t, 
                         st=atrisk.t,
                         scenario = scenario)
    return(sur_df)
  }
  
  invisible(NULL)
}

# cured (t)
## surivial estimates
cure_plot_df <- function(x, t = NULL, ci=TRUE, n=1000, xvals=NULL, zvals=NULL,
                       xpred_name, xpred_value, scenario,...) {
  
  # Set t vector if needed to 1.2 * max observed duration; lower limit is 1
  if (is.null(t)) {
    max_t <- round(max(x$Y[x$Y[, "last"]==1, "duration"]) * 1.2)
    t     <- seq(1, max_t, length.out=100)
  } 
  
  # Extract covariate matrices
  dur.dat  <- x$mf.dur
  risk.dat <- x$mf.risk 
  X <- model.matrix(attr(x$mf.dur, 'terms'), data=x$mf.dur)
  Z <- model.matrix(attr(x$mf.risk, 'terms'), data=x$mf.risk)
  
  # Extract coefficient point estimates
  beta  <- coef(x, model = "duration")
  gamma <- coef(x, model = "risk")
  alpha <- coef(x, model = "distr")
  alpha <- exp(-alpha)
  beta_vcv  <- vcov(x, "duration")
  gamma_vcv <- vcov(x, "risk")
  alpha_vcv <- vcov(x, "distr")
  
  if (is.null(xvals)) {
    X[, xpred_name] <- xpred_value
    X_vals <- apply(X, 2, mean)		
  } else if (!length(xvals)==ncol(X) && length(xvals)-ncol(X)==-1) {
    stop("Incorrect length for xvals, did you forget 1 for intercept term?")			
  } else if (!length(xvals)==ncol(X)) {
    stop("Incorrect length for xvals")
  } else {
    X_vals <- xvals
  }
  
  if (is.null(zvals)) {
    Z[, xpred_name] <- xpred_value
    Z_vals <- apply(Z, 2, mean)		
  } else if (!length(zvals)==ncol(Z) && length(zvals)-ncol(Z)==-1) {
    stop("Incorrect length for zvals, did you forget 1 for intercept term?")			
  } else if (!length(zvals)==ncol(Z)) {
    stop("Incorrect length for zvals")
  } else {
    Z_vals <- zvals
  }
  
  # Calculate hazard using point estimates only
  lambda <- exp(-X_vals %*% beta)
  cure   <- 1 - plogis(Z_vals %*% gamma)
  
  cure.t <- cure(ti = t, lambda = lambda, cure = cure, alpha = alpha,
             out = NULL, dist = x$distr)
  
  if (ci==TRUE) {
    set.seed(1234)
    Coef_smpl <- MASS::mvrnorm(n = n, mu = coef(x, "full"), Sigma = vcov(x, "full"))
    
    b_idx <- 1:x$n.terms$duration
    g_idx <- (max(b_idx) + 1):(max(b_idx) + x$n.terms$risk)
    a_idx <- (max(g_idx) + 1):ncol(Coef_smpl)
    
    Beta  <- Coef_smpl[, b_idx]
    Gamma <- Coef_smpl[, g_idx]
    A     <- Coef_smpl[, a_idx]
    Alpha <- exp(-A)
    
    lambda <- exp(-tcrossprod(X_vals, Beta))
    cure   <- 1 - plogis(tcrossprod(Z_vals, Gamma))
    
    sims <- matrix(nrow = length(t), ncol = n)
    hmat <- matrix(nrow = length(t), ncol = 3)
    
    for (i in 1:n) {
      sims[, i] <- cure(ti = t, lambda = lambda[i], cure = cure[i], 
                      alpha = Alpha[i], out = NULL, dist = x$distr)
    }
    
    hmat[, 1] <- cure.t
    hmat[, 2] <- apply(sims, 1, quantile, probs = 0.05)
    hmat[, 3] <- apply(sims, 1, quantile, probs = 0.95)
    
    ## use ggplot
    sur_df <- data.frame(Time = t,
                         cure = hmat[, 1],
                         lo95 = hmat[, 2],
                         up95 = hmat[, 3],
                         scenario = scenario)
    return(sur_df)
  } else {
    # Plot without CIs
    
    sur_df <- data.frame(Time = t, 
                         cure=atrisk.t,
                         scenario = scenario)
    return(sur_df)
  }
  
  invisible(NULL)
}


glm_roc_plot <- function(ModelResults, modelname){
  require(plotROC)
  require(pROC)
  require(dplyr)
  library(purrr)
  library(tibble)
  library(magrittr)
  Y_pred =lapply(ModelResults, function(x) FUN = predict(x, type = "response"))
  
  Y_obs = lapply(ModelResults, function(x) FUN = x$y)
  
  roc = Map(function(x, y) roc(x,y), Y_obs, Y_pred)
  
  names(roc) <- modelname
  
  roc_df = lapply(roc, function(x) FUN= data.frame(
    plotx = x$specificities,
    ploty = rev(x$sensitivities),
    name = paste("AUC =",
                 sprintf("%.3f",x$auc)))) %>%
    map_df(., rbind, .id="modelname")
  return(roc_df)
}


spdur_pr_plot <- function(ModelResults, modelname){
  require(plotROC)
  require(pROC)
  require(dplyr)
  library(purrr)
  library(tibble)
  library(magrittr)
  Y_pred =lapply(ModelResults, function(x) FUN = predict(x, type = "response"))
  Y_obs = lapply(ModelResults, function(x) FUN = x$Y[,"fail"])
  
  pr_df = Map(function(x, y)  prediction(x, y), Y_pred, Y_obs)
  pr_df <- lapply(pr_df,function(x) FUN = performance(x, "prec", "rec"))
  
  names(pr_df) <- modelname
  
  pr_df = lapply(pr_df, function(x) FUN= data.frame(
    rec = x@x.values[[1]],
    prec = x@y.values[[1]]))%>%
    map_df(., rbind, .id="modelname")
  pr_df[1,3] <- 0
  return(pr_df)
}


glm_pr_plot <- function(ModelResults, modelname){
  require(plotROC)
  require(pROC)
  require(dplyr)
  library(purrr)
  library(tibble)
  library(magrittr)
  Y_pred =lapply(ModelResults, function(x) FUN = predict(x, type = "response"))
  Y_obs = lapply(ModelResults, function(x) FUN = x$y)
  
  pr_df = Map(function(x, y)  prediction(x, y), Y_pred, Y_obs)
  pr_df <- lapply(pr_df,function(x) FUN = performance(x, "prec", "rec"))
  
  names(pr_df) <- modelname
  
  pr_df = lapply(pr_df, function(x) FUN= data.frame(
    rec = x@x.values[[1]],
    prec = x@y.values[[1]]))%>%
    map_df(., rbind, .id="modelname")
  pr_df[1,3] <- 0
  return(pr_df)
}



## a function to extract coefficient and build ggplot
plot_spdm_coef <- function(x, subvar = FALSE, variables = NULL){
  library(spduration)
  library(ggplot2)
  library(dplyr)
  library(ggthemes)
  library(cowplot)

  
  
  beta  <- coef(x, model = "duration")
  gamma <- coef(x, model = "risk")
  alpha <- coef(x, model = "distr")
  
  beta_se  <- sqrt(diag(vcov(x, "duration")))
  gamma_se <- sqrt(diag(vcov(x, "risk")))
  alpha_se <- sqrt(diag(vcov(x, "distr")))
  
  beta_df <- data.frame(beta = beta,
                        vars = names(beta),
                        se = beta_se,
                        type = "beta") %>%
    dplyr::mutate(ad = .5/se) %>%
    dplyr::mutate(beta = beta*ad,
                  se = se*ad) %>%
    dplyr::mutate(lo95 = beta - 1.96*se,
                  hi95 = beta +1.96*se) %>%
    dplyr::mutate(color_sig = ifelse(lo95 < 0 & hi95 > 0, "#ef8a62", "#2c7bb6"))
  
  beta_df <- beta_df %>%
    dplyr::mutate(vars = as.character(vars)) %>%
    dplyr::arrange(vars)
  
  if(subvar == TRUE){
    beta_df <- beta_df[grep(paste(variables, collapse = "|"), beta_df$vars),]
    #re-order the variables: if factor, the order of factor levels is used; if character, an alphabetical order ist used
    beta_df$vars <- factor(beta_df$vars, levels = variables)
  } 
  
  
  gamma_df <- data.frame(beta = gamma,
                         vars = names(gamma),
                         se = gamma_se,
                         type = "gamma")%>%
    dplyr::mutate(ad = .5/se) %>%
    dplyr::mutate(beta = beta*ad,
                  se = se*ad) %>%
    dplyr::mutate(lo95 = beta - 1.96*se,
                  hi95 = beta +1.96*se) #%>%
   # dplyr::mutate(color_sig = ifelse(lo95 < 0 & hi95 > 0, "#ef8a62", "#2c7bb6"))
  
  gamma_df <- gamma_df %>%
    dplyr::mutate(vars = as.character(vars)) %>%
    dplyr::arrange(vars)
  
  if(subvar == TRUE){
    gamma_df <- gamma_df[grep(paste(variables, collapse = "|"), gamma_df$vars),]
    #re-order the variables: if factor, the order of factor levels is used; if character, an alphabetical order ist used
    gamma_df$vars <- factor(gamma_df$vars, levels = variables)
  } 
  
  
  alpha_df <- data.frame(beta = alpha,
                         vars = names(alpha),
                         se = alpha_se,
                         type = "alpha")%>%
    dplyr::mutate(ad = .5/se) %>%
    dplyr::mutate(beta = beta*ad,
                  se = se*ad) %>%
    dplyr::mutate(lo95 = beta - 1.96*se,
                  hi95 = beta +1.96*se) #%>%
    #dplyr::mutate(color_sig = ifelse(lo95 < 0 & hi95 > 0, "#ef8a62", "#2c7bb6"))
   
  
  coef_df <- list(beta_df, gamma_df, alpha_df)
  plot_obj <- list()
  
  for (ii in seq_along(coef_df)){
    plot_obj[[ii]] <- ggplot(coef_df[[ii]], aes(x=vars, y=beta, color = "red"))+ #, color=color_sig)) +
      geom_hline(yintercept = 0, colour = gray(1/2), lty = 2) + 
      geom_point(size=3) + 
      geom_linerange(aes(ymin=lo95, ymax=hi95),alpha = 1, size = 1.5) + 
      coord_flip()  +
      theme_bw() +
      xlab('') + ylab('') +
      theme(legend.position="none",
            legend.title=element_blank(),
            axis.text = element_text(size=14),
            text = element_text(size=14),
            plot.title = element_text(hjust = .5, size = 16, face = "bold"))
  }
  return(plot_obj)
}
