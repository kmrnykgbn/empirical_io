log_prodction <-  function(l, k , omega, eta, beta_0, beta_l, beta_k) {
  ln_y = beta_0 + beta_l*l + beta_k*k +omega + eta 
  return(ln_y)
}

log_labor_choice <- function(k, wage, omega, beta_0, beta_l, beta_k, sigma_eta) { 
  n = length(k)
  eta = rnorm(n, 0, sigma_eta)
  ln_l = log( ((1/wage) * beta_l *exp(beta_0+omega) * exp(k)^beta_k)^(1/(1-beta_l)) )
  return(ln_l)
}

log_labor_choice_error <- function(k, wage, omega, beta_0, beta_l, beta_k, iota, sigma_eta) { 
  n = length(k)
  eta = rnorm(n, 0, sigma_eta)
  ln_l = log( ((1/wage) * beta_l *exp(beta_0+omega+iota) * exp(k)^beta_k)^(1/(1-beta_l)) )
  return(ln_l)
}

investment_choice <- function(k, omega, gamma, delta){
  I = (delta + gamma*omega)*k
  return(I)
}
