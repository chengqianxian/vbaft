vb.aft.get.elbo <- function(Y, X, delta, mu_0, v_0, alpha_0, omega_0, curr_mu_vector, curr_b, alpha, omega, Sigma){
  res <- 0
  for (i in 1:nrow(X)) {
    C_i <- ifelse((Y[i] - sum(X[i, ] * curr_mu_vector)) / curr_b <= -5, 0,
                  ifelse((Y[i] - sum(X[i, ] * curr_mu_vector)) / curr_b <= -1.701 & 
                           (Y[i] - sum(X[i, ] * curr_mu_vector)) / curr_b > -5, 0.0426,
                         ifelse((Y[i] - sum(X[i, ] * curr_mu_vector)) / curr_b <= 0 & 
                                  (Y[i] - sum(X[i, ] * curr_mu_vector)) / curr_b > -1.701, 0.3052,
                                ifelse((Y[i] - sum(X[i, ] * curr_mu_vector)) / curr_b <= 1.702 & 
                                         (Y[i] - sum(X[i, ] * curr_mu_vector)) / curr_b > 0, 0.6950,
                                       ifelse((Y[i] - sum(X[i, ] * curr_mu_vector)) / curr_b <= 5 & 
                                                (Y[i] - sum(X[i, ] * curr_mu_vector)) / curr_b > 1.702, 0.9574, 1)))))
    res_i <- (delta[i] - C_i * (1 + delta[i])) * (Y[i] - sum(X[i, ] * curr_mu_vector))
    res <- res + res_i
  }
  r <- sum(delta)
  expect_inverse_b <- alpha / omega
  expect_log_b <- log(omega) - digamma(alpha)
  expect_log_likelihood <- -r * expect_log_b + expect_inverse_b * res # the first term in the ELBO
  
  diff_beta <- (- v_0 * (sum(diag(Sigma)) + sum((curr_mu_vector - mu_0) * (curr_mu_vector - mu_0))) + log(det(Sigma))) / 2 # the second term in the ELBO
  
  diff_b <- (alpha - alpha_0) * expect_log_b + (omega - omega_0) * expect_inverse_b - alpha * log(omega) # the third term
  elbo <- expect_log_likelihood + diff_beta + diff_b
  return(elbo)
}
