vb_aft_final <- function(Y, X, delta, mu_0, v_0, alpha_0, omega_0, max_iteration = 100, threshold = 0.0001){
  n <- nrow(X)
  p <- ncol(X)
  
  alpha <- get_alpha_star(alpha_0, delta) # fixed always
  omega <- omega_0 # initialization
  mu_vector <- mu_0 # initialization
  
  converged <- FALSE
  iteration <- 0
  curr_mu_vector <- mu_0
  curr_b <- omega / (alpha - 1)
  curr_elbo <- 0
  while (converged == FALSE & iteration <= max_iteration) {
    iteration <- iteration + 1
    Sigma <- get_Sigma_star(curr_b, curr_mu_vector, alpha, omega, delta, X, v_0, Y)
    mu_vector <- get_mu_star(curr_b, curr_mu_vector, alpha, omega, delta, X, v_0, Sigma, mu_0, Y)
    omega <- get_omega_star(curr_b, omega_0, delta, Y, X, mu_vector)
    elbo <- vb.aft.get.elbo(Y, X, delta, mu_0, v_0, alpha_0, omega_0, curr_mu_vector, curr_b, alpha, omega, Sigma)
    converged_1 <- ifelse(abs(elbo - curr_elbo) <= threshold, TRUE, FALSE)
    converged_2 <- ifelse(sum(abs(mu_vector - curr_mu_vector)) <= threshold, TRUE, FALSE)
    converged <- ifelse(converged_1 == FALSE & converged_2 == FALSE, FALSE, TRUE)
    curr_elbo <- elbo
    curr_mu_vector <- mu_vector
    curr_b <- omega / (alpha - 1)
  }
  if(iteration > max_iteration & converged == FALSE) warning("The algorithm can not converge and the max iteration has achieved")
  return(list(
    converged_ELBO = curr_elbo,
    iteration_used = iteration,
    # a = a,
    alpha = alpha,
    Sigma = Sigma,
    # G = G,
    # m_vector = m_vector,
    # b = b,
    omega = omega,
    mu_vector = mu_vector))
}
