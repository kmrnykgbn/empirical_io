compute_indirect_utility <- function(df, beta, sigma, mu, omega){
  K = length(beta)
  nu = as.matrix(df[c("v_x_1", "v_x_2", "v_x_3")])
  x = as.matrix(df[c("x_1", "x_2", "x_3")])
  alpha=-exp(mu + omega*df$v_p)
  u = alpha*df$p
  for(k in 1:K){
    beta_t = beta[k] + sigma[k]*nu[,k]
    add_u = beta_t * x[,k]
    u = u + add_u
  }
  
  return(u)
}

compute_choice <- function(df, e_t, beta, sigma, mu, omega) {
  # compute utility
  df <- df %>%
    mutate(u = compute_indirect_utility(df, beta, sigma, mu, omega))
  
  # add q
  df <- df %>%
    mutate(e = e_t, u_e = u + e) %>%
    group_by(t, i) %>%
    mutate(q = if_else(u_e == max(u_e), 1, 0)) %>%
    ungroup()
  
  return(df)
}

compute_share <- function(X, M, V, e, beta, sigma, mu, omega) {
  # merge
  df <- left_join(V, M, by="t")
  df <- left_join(df, X, by="j") %>%
    select(t, i, j, v_x_1, v_x_2, v_x_3, v_p, x_1, everything())
  
  # compute_choice
  df_choice <- compute_choice(df, e, beta, sigma, mu, omega)
  
  N <- length(df_choice$i %>% unique())
  df_share <- df_choice %>%
    group_by(t, j) %>%
    mutate(s=sum(q)/N) %>%
    ungroup()
  df_share <- df_share %>%
    group_by(t, i) %>%
    mutate(y=log(s/s[1])) %>%
    ungroup() %>%
    replace_na(list(y = 0))

  # summarize
  df_share <- df_share %>%
    select(i, t, j, x_1, x_2, x_3, xi, c, p, q, s, y) %>%
    group_by(t, j) %>%
    select(t, j, x_1, x_2, x_3, xi, c, p, q, s, y) %>%
    summarise(t=mean(t), j=mean(j), x_1=mean(x_1), x_2=mean(x_2), x_3=mean(x_3),
              xi=mean(xi), c=mean(c), p=mean(p), q=sum(q),
              s=mean(s), y=mean(y)) %>%
    group_by()
   
  return(df_share)
}

NLLS_objective <- function(theta, df_share, X, M, V_mcmc, e_mcmc) {
  df_share_mcmc <- compute_share(X, M, V_mcmc, e_mcmc,
                                 theta[1:3], theta[4:6], theta[7], theta[8])
  MSE <- (df_share - df_share_mcmc) %>%
    summarize(mean(s^2))
  
  return(as.numeric(MSE))
}

