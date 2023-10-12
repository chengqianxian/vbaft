
# This is an example for simulation of a sample size, 300
# First let us create three NULL matrices for saving results
replicate.total <- 500
aft.size.300.cen.0 <- matrix(0, nrow = 8, ncol = replicate.total)
aft.size.300.cen.10 <- matrix(0, nrow = 8, ncol = replicate.total)
aft.size.300.cen.30 <- matrix(0, nrow = 8, ncol = replicate.total)

# using while loop to get replicate.total converged cases
i <- 1
seed.con <- numeric(replicate.total)
j <- 1

start.time <- Sys.time()
while (i <= replicate.total) {
  # generate true survival times
  set.seed(j)
  x1 <- rnorm(300, 1, 0.2)
  x2 <- rbinom(300, 1, 0.5)
  z <- rlogis(300)
  beta0 <- 0.5
  beta1 <- 0.2
  beta2 <- 0.8
  b <- 0.8
  Y <- beta0 + beta1 * x1 + beta2 * x2 + b * z 
  T <- exp(Y)
  
  # generate censoring times
  set.seed(j)
  cen.time.10 <- runif(300, 0, 48)
  cen.time.30 <- runif(300, 0, 17)
  
  # obtain observed time
  T.10 <- pmin(T, cen.time.10)
  T.30 <- pmin(T, cen.time.30)
  
  # obtain censoring indicator
  delta <- rep(1, 300)
  delta.10 <- ifelse(T == T.10, 1, 0)
  delta.30 <- ifelse(T == T.30, 1, 0)
  
  # create X matrix
  X <- matrix(c(rep(1, 300), x1, x2), nrow = 300)
  
  # priors, use non-informative priors
  mu_0 <- c(0, 0, 0)
  v_0 <-  0.1
  alpha_0 <- 11
  omega_0 <- 10
  
  # running models
  vb.aft.model.cen.0 <- vb_aft_final(log(T), X, delta, mu_0, v_0, alpha_0, 
                                     omega_0, max_iteration = 100, threshold = 0.01)
  if(vb.aft.model.cen.0$iteration_used > 100) j = j + 1
  else{
    vb.aft.model.cen.10 <- vb_aft_final(log(T.10), X, delta.10, mu_0, v_0, alpha_0, 
                                       omega_0, max_iteration = 100, threshold = 0.01)
    if(vb.aft.model.cen.10$iteration_used > 100) j = j + 1
    else{
      vb.aft.model.cen.30 <- vb_aft_final(log(T.30), X, delta.30, mu_0, v_0, alpha_0, 
                                          omega_0, max_iteration = 100, threshold = 0.01)
      if(vb.aft.model.cen.10$iteration_used <= 100){
        aft.size.300.cen.0[, i] <- c(vb.aft.model.cen.0$alpha, 
                                     vb.aft.model.cen.0$omega, 
                                     vb.aft.model.cen.0$mu_vector,
                                     diag(vb.aft.model.cen.0$Sigma))
        aft.size.300.cen.10[, i] <- c(vb.aft.model.cen.10$alpha, 
                                      vb.aft.model.cen.10$omega, 
                                      vb.aft.model.cen.10$mu_vector,
                                      diag(vb.aft.model.cen.10$Sigma))
        aft.size.300.cen.30[, i] <- c(vb.aft.model.cen.30$alpha, 
                                      vb.aft.model.cen.30$omega, 
                                      vb.aft.model.cen.30$mu_vector,
                                      diag(vb.aft.model.cen.30$Sigma))
        seed.con[i] <- j
        j = j + 1
        i = i + 1
      }
      else j = j + 1
    }
  }
}
end.time  <- Sys.time()
end.time - start.time


write.csv(aft.size.300.cen.0, "aft.size.300.cen.0.weak.csv")
write.csv(aft.size.300.cen.10, "aft.size.300.cen.15.weak.csv")
write.csv(aft.size.300.cen.30, "aft.size.300.cen.30.weak.csv")

