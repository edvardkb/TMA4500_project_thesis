library(RTMB)
library(Matrix)

set.seed(1)


smooth_min <- function(x, U, eps = 1e-6) {
  0.5 * (x + U - sqrt((x - U)^2 + eps^2))
}

smooth_max <- function(x, L, eps = 1e-6) {
  0.5 * (x + L + sqrt((x - L)^2 + eps^2))
}

smooth_clamp <- function(x, L, U, eps = 1e-6) {
  smooth_min(smooth_max(x, L, eps), U, eps)
}


safe_exp <- function(x) exp(smooth_clamp(x, -100, 100))




SV_Model_SHASH_inference <- function(data, fixed_SV_params = list(), burnin_time = 50, norm_latent_variables = TRUE){
  # Computes MLE of the marginal likelihood function of the SV model with SHASH noise.
  # Args:
  #   data: numeric list, complete time series data (y_1:n) which the SV model are to be fitted with.
  #   normal_transformed_latent_variables: bool, when true the joint likelihood is f(y_1:n, h_1, z_2:n) with
  #                                        z_2:n being standard normal distributed. When false it is f(y_1:n, h_1:n).
  #   fixed_SV_params: numeric list, contains parameters of the distribution of the volatility that is to be fixed.
  # Returns:
  #   summary: list, contains the fixed parameter estimates.
  
  
  # ---------- Define helper functions --------------
  
  E_SHASH <- function(sigma, epsilon, delta){
    # Helper function used to center eta.
    # Computes expected value of a SHASH distributed r.v. given the parameters sigma, epsilon, delta.
    P <- function(q){
      exp(1 / 4) / sqrt(8 * pi) * (besselK(1 / 4, (q + 1) / 2) + besselK(1 / 4, (q - 1) / 2))
    }
    
    E_SHASH <- sigma * sinh(epsilon / delta) * P(1 / delta)
    
    return(E_SHASH)
  }
  
  dSHASH <- function(x, sigma, epsilon, delta, log = FALSE) {
    # Helper function used for the distribution of latent variables
    # when treated as SHASH distributed.
    # Returns the pdf (log(pdf) if log = TRUE) of SHASH given the parameters. 
    z  <- asinh(x / sigma)
    w  <- delta * z - epsilon
    s  <- sinh(w)
    c  <- cosh(w)
    logdens <-  - 0.5*log(2*pi) - 0.5*log1p((x / sigma)^2) + log(delta) + log(c) - 0.5*s^2 - log(sigma)
    if (log) logdens else exp(logdens)
  }
 
  Eta <- function(Z, sigma, epsilon, delta){
    # Transforms a standard normal distributed r.v. Z into the noise term
    # of the stochastic volatility process. Noise term is assumed to follow
    # a centered SHASH distribution, i.e. E(Eta(Z|sigma, epsilon, delta)) = 0.
    # Args:
    #  Z: numeric, realization of a standard normal r.v.
    #  sigma: numeric, scale parameter for SHASH.
    #  epsilon: numeric, skew parameter for SHASH.
    #  delta: numeric, tail weight parameter of SHASH.
    # Returns:
    #  numeric, realization of a centered SASH distributed r.v.
    
    SHASH <- sigma * sinh((asinh(Z) + epsilon) / delta)
    
    return(SHASH - E_SHASH(sigma, epsilon, delta))
  }
  
  
  # --------- Create the joint nll which is to be minimized -------------
  f_tilde <- function(par){
    # nll of the joint pdf. Used for RTMB.
    
    # Unpack and transform AR par
    mu    <- par$mu
    alpha <- par$alpha
    phi   <- tanh(par$psi)        # |phi|<1
    u     <- par$u                # random effects vector (length n)
    
    
    # Unpack and transform Vol par
    sigma <- exp(par$log_sigma)
    epsilon <- par$epsilon
    delta <- exp(par$log_delta)
    
    
    # y is define as the observed data
    n <- length(data)
    m <- length(u) - n
    y <- data
    y <- OBS(y)
    
    
    if(norm_latent_variables){
      nll <- - sum(dnorm(u[1:n+m], log = TRUE))
      
      # Add the nll of y
      
      
      h <- alpha
      if (m>0){
        for (i in 1:m){
          h <- alpha + phi * (h - alpha) + Eta(u[i], sigma, epsilon, delta)
        }
      }
      
      
      for (i in (m+1):(m+n)){
        h <- alpha + phi * (h - alpha) + Eta(u[i], sigma, epsilon, delta)
        nll <- nll - dnorm(y[i-m], mean = mu, sd = safe_exp(0.5*h), log = TRUE)
      }
      
    }
    else{
      nll <- - dSHASH(u[1] - alpha + E_SHASH(sigma, epsilon, delta), sigma, epsilon, delta, log = TRUE)
      
      # Observations: y_t | h_t ~ N(mu, exp(h_t))   (dropping constants)
      nll <- nll + 0.5 * sum( u[(m+1):(m+n)] + (y - mu)^2 * safe_exp(-u[(m+1):(m+n)]) )
      
      # Transitions: h_{t+1} | h_t
      for (i in 1:(n+m-1)){
        translational_transform <- u[i+1] - alpha + E_SHASH(sigma, epsilon, delta) - phi * (u[i] - alpha)
        nll <- nll - dSHASH(translational_transform, sigma, epsilon, delta, log = TRUE)
      }
    }
    
    
    
    
    return(nll)
  }
  
  
  # -------- Find initial guesses for the parameters used in joint nll ---------
  
  # Initiate a list of the latent variables.
  u <- numeric(length(data) + burnin_time)
  
  
  # Transforming AR and Vol par to an unconstrained format.
  transformed_AR_par <- list(
    mu         = mean(data),
    alpha      = log(var(data) + 1e-12),
    psi        = atanh(0.5),
    u          = u  # initialized as a numeric list of zeros, i.e. the expected value of z_(-m+1):n
  )
  transformed_SV_par <- list(
    log_sigma = 0,  #Harder to initialize based on observed data. 
    epsilon = 0,
    log_delta = 0
  )
  
  # Combining them into one list.
  par0 <- c(transformed_AR_par, transformed_SV_par)
  
  
  # --------- If specified, state which parameters are to be fixed -----------------
  if(length(fixed_SV_params) != 0){
    
    map <- list()
    
    if("mu" %in% names(fixed_SV_params)){
      par0$mu <- fixed_SV_params$mu
      map$mu <- factor(NA)
    }
    if("alpha" %in% names(fixed_SV_params)){
      par0$alpha <- fixed_SV_params$alpha
      map$alpha <- factor(NA)
    }
    if("phi" %in% names(fixed_SV_params)){
      par0$psi <- atan(fixed_SV_params$alpha)
      map$psi <- factor(NA)
    }
    if("epsilon" %in% names(fixed_SV_params)){
      par0$epsilon <- fixed_SV_params$epsilon
      map$epsilon <- factor(NA)
    }
    if("sigma" %in% names(fixed_SV_params)){
      par0$log_sigma <- log(fixed_SV_params$sigma)
      map$log_sigma <- factor(NA)
    }
    if("delta" %in% names(fixed_SV_params)){
      par0$log_delta <- log(fixed_SV_params$delta)
      map$log_delta <- factor(NA)
    }
    obj <- MakeADFun(
      func     = f_tilde,      # f is def as the nll. The closure captures y as data
      parameters = par0,       # Starting parameters
      map = map,
      random   = "u",
      random.start = expression(last.par[random]),
      inner.control=list(maxit=500, smartsearch=TRUE, trace=1),
      silent = FALSE
    )
  }
  else{
    obj <- MakeADFun(
      func     = f_tilde,      # f is def as the nll. The closure captures y as data
      parameters = par0,       # Starting parameters
      random   = "u",
      random.start = expression(last.par.best[random] + 0.2 * rnorm(length(last.par.best[random]))),
      inner.control = list(
        smartsearch = TRUE,
        maxit = 2000,
        trace = 1,
        silent = FALSE,
        tol = 1e-10,
        grad.tol = 1e-10,
        step.tol = 1e-10,
        tol10 = 0,     # IMPORTANT: don't early-exit on "no improvement"
        ustep = 0.1,
        u0 = 1e-2,
        power = 0.5,
        super = FALSE
      )
    )
  }
  
  
  # ------- Wrapper for catching bad parameters yielding NaN objective or NaN gradient -------
  bad <- new.env(parent = emptyenv())
  
  fn_wrap <- function(par) {
    val <- obj$fn(par)
    #H <- obj$env$spHess(par, random = TRUE)
    #cat("any nonfinite in H:", any(!is.finite(as.matrix(H))), "\n")
    if (!is.finite(val)) {
      bad$which <- "objective"
      bad$par   <- par
      
      bad$val   <- val
      bad$full  <- obj$env$last.par  # includes random effects at the last inner state
      
      cat("NaN eval in", bad$which, "\n")
      # 1) Joint nll at the current full state:
      f_joint <- obj$env$f(bad$full)
      cat("Joint nll: ", f_joint, "\n")
      
      Huu <- obj$env$spHess(par = bad$full, random = TRUE)
      Huu <- (Huu + Matrix::t(Huu)) * 0.5   # enforce symmetry numerically
      
      chol_ok <- try(Matrix::Cholesky(Huu, LDL = FALSE), silent = TRUE)
      if (inherits(chol_ok, "try-error")) {
        cat("Huu is NOT SPD at this state (Laplace/Cholesky will fail)\n")
      } else {
        # logdet(H) from Cholesky (what Laplace wants)
        ldet_chol <- 2 * sum(log(diag(as(chol_ok, "Matrix"))))
        cat("Cholesky logdet(Huu):", ldet_chol, "\n")
      }
      
      if (!is.null(Huu@x) && any(!is.finite(Huu@x))) {
        print("Huu contains NaN/Inf (stored entries)")
      }
      else{
        n <- nrow(Huu)
        k <- min(6L, n - 1L)
        
        Afun <- function(x, args = NULL) {
          drop(Huu %*% x)   # drop -> plain numeric vector
        }
        
        ev <- RSpectra::eigs_sym(
          A = Afun,
          n = n,
          k = k,
          which = "SA",
          opts = list(retvec = FALSE)
        )$values
        
        cat("6 lowest eigenvalues in Huu:", sort(ev), "\n")
        
        
      }
    }
    val
  }
  
  gr_wrap <- function(par) {
    print("Gradient was called")
    g <- obj$gr(par)
    
    if (any(!is.finite(g))) {
      bad$which <- "gradient"
      bad$par   <- par
      bad$g     <- g
      bad$full  <- obj$env$last.par
      
    }
    
    g
  }
  
  
  # ------- Begin outer optimization routine -----------
  
  
  opt <- nlminb(start = obj$par, objective = fn_wrap, gradient = gr_wrap,
                control = list(trace = 0))
  
  
  rep <- sdreport(obj)
  
  print(rep)
  print(opt$message)
  
  if (length(bad$full) != 0){
    cat("NaN eval in", bad$which, "\n")
    # 1) Joint nll at the current full state:
    f_joint <- obj$env$f(bad$full)
    cat("Joint nll: ", f_joint, "\n")
    
    # 2) Hessian wrt random effects + logdet:
    Huu <- obj$env$spHess(par = bad$full, random = TRUE)
    ldet <- as.numeric(determinant(Huu, logarithm = TRUE)$modulus)
    
    cat("Log det hessian:", ldet, "\n")
    
    # 3) Inner gradient size:
    g_all <- obj$env$f(bad$full, order = 1)[4:103]  #Broken
    cat("Maximum inner gradient component:", max(g_all), "\n")
    
    saveRDS(bad$full, file = "estimated_params.rds")
  }
  
  
  summary <- rep$par.fixed
  
  
  # ------------  Undo the transformations ---------------------
  if("epsilon" %in% names(summary) && nn_epsilon){
    summary["epsilon"] <- sqrt(summary["epsilon"])
    names(summary)[names(summary) == "epsilon"] <- "epsilon"
  }
  if("psi" %in% names(summary)){
    summary["psi"] <- tanh(summary["psi"])
    names(summary)[names(summary) == "psi"] <- "phi"
  }
  if("log_sigma" %in% names(summary)){
    summary["log_sigma"] <- exp(summary["log_sigma"])
    names(summary)[names(summary) == "log_sigma"] <- "sigma"
  }
  if("log_delta" %in% names(summary)){
    summary["log_delta"] <- exp(summary["log_delta"])
    names(summary)[names(summary) == "log_delta"] <- "delta"
  }
  
  return(summary)
}

generateData <- function(n, AR_par, SHASH_par, random = TRUE){
  
  SHASH_mean_sd <- function(sigma, epsilon, delta){
    
    P <- function(q){
      exp(1 / 4) / sqrt(8 * pi) * (besselK(1 / 4, (q + 1) / 2) + besselK(1 / 4, (q - 1) / 2))
    }
    
    EX <- sinh(epsilon / delta) * P(1 / delta)
    
    EX2 <- 1 / 2 * (cosh(2 * epsilon / delta) * P(2 / delta) - 1)
    
    return(list(mean = sigma * EX, sd = sigma * sqrt(EX2 - EX^2)))
  }
  
  rSHASH <- function(n, sigma = 1, epsilon = 0, delta = 1){
    
    Z <- rnorm(n)
    shash <- sigma * sinh((asinh(Z) + epsilon) / delta)
    
    return(shash)
  }
  
  # Unpack AR parameters.
  alpha <- AR_par$alpha
  phi <- AR_par$phi
  mu <- AR_par$mu
  
  # Unpack SV parameters.
  sigma <- SHASH_par$sigma
  epsilon <- SHASH_par$epsilon
  delta <- SHASH_par$delta
  
  
  epsilon_tilde <- rnorm(n)
  
  
  SHASH_moments <- SHASH_mean_sd(sigma, epsilon, delta)
  eta <- rSHASH(n, sigma, epsilon, delta) - SHASH_moments$mean
  
  h <- numeric(n)
  y <- numeric(n)
  
  
  if (random){ #We make the assumption that the infinite weighted sum is normal distributed
    h[1] <- rnorm(1, alpha, SHASH_moments$sd/sqrt(1-phi^2))
  }
  else{
    h[1] <-  alpha
  }
  
  y[1] <- mu + exp(h[1] / 2) * epsilon_tilde[1]
  
  for (t in 1:(n-1)) {
    h[t+1] <- alpha + phi * (h[t] - alpha) +  eta[t+1]
    y[t+1] <- mu + exp(h[t+1] / 2) * epsilon_tilde[t+1]
  }
  
  return(list(y = y, h = h))
}

# ---------------- Defining true parameters of the model and creating synthetic data -----------

# True AR(1) parameters
true_AR_params <- list(
  alpha = 2,
  phi = 0.95,
  mu = 0
)
# True SV parameters
true_SV_params <- list(
  sigma = 1,
  epsilon = 0,
  delta = 1
)


# Create synthetic data
n <- 100
data <- generateData(n, true_AR_params, true_SV_params, random = FALSE)
saveRDS(data, file = "synthetic_data.rds")
y <- data$y
h <- data$h

plot(y)

# Possible estimator for alpha: log(var(y))
# Possible estimator for mu: E(y)

# --------------------------- Running RTMB --------------------------
fixed_SV_params <- list(epsilon = true_SV_params$epsilon, delta = true_SV_params$delta) #List of fixed parameter values. Remove element if it should be treated
#as a free variable
summary <- SV_Model_SHASH_inference(y, fixed_SV_params = fixed_SV_params, norm_latent_variables = FALSE)

print(summary)



# delta = 0.4 resterende standard cond. (mu fixed) -> ustabil
# delta = 0.3 resterende standard cond. (mu, alpha fixed) -> ustabil

