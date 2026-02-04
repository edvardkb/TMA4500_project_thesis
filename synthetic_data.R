set.seed(1)

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
  #eta <- rgamma(n, shape = 4, rate = 1) - 4
  #eta <- rsn(n, alpha  = 0.2)
  
  
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


# ---- Parameters used for the synthetic log vol plot  ----

# True AR(1) parameters
true_AR_params <- list(
  alpha = -4.4,
  phi = 0.975,
  mu = 0
)
# True SV parameters
true_SV_params <- list(
  sigma = 0.001,
  epsilon = 0.150,
  delta = 0.330
)
