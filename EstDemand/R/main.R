library(tidyverse)
library(evd)
library(doParallel)

source("R/functions.R")

#set seed
set.seed(1)

J<-10
K<-3
T<-100
N<-500
L<-500
beta <- rnorm(K)
beta[1]<-4
beta
sigma<-abs(rnorm(K)); sigma
mu <- 0.5
omega<-1


price_xi<-1
prop_jt<-0.6
sd_x<-0.5
sd_c <- 0.05
sd_p<-0.05

# create X_data
X <- tibble(j=1:10, x_1=1,
            x_2=rnorm(J, 0, sd_x),
            x_3=rnorm(J, 0, sd_x))
X <- rbind(0, X)
X
# create M_data
rlnorm(T*J)
M <- expand_grid(t = 1:T, j = 1:J) %>%
  mutate(xi = 0,
         c = rlnorm(T*J,0,sd_c),
         p = rlnorm(T*J,0,sd_p) + c) 

M <- M %>% 
  group_by(t) %>%
  dplyr::sample_frac(prop_jt) %>%
  ungroup()

add_M <- cbind(1:T, as.data.frame(matrix(0, T, 4)))
names(add_M) <- names(M)
M <- rbind(add_M, M) %>%
  select(j, t, xi, c, p) %>%
  arrange(t, j)


# create V_data
V <- expand_grid(t=1:T, i=1:N) %>%
  mutate(v_x_1=rnorm(N*T), 
         v_x_2=rnorm(N*T),
         v_x_3=rnorm(N*T),
         v_p=rnorm(N*T)) %>%
  select(i, t, v_x_1, v_x_2, v_x_3, v_p)
V

# merge
df <- left_join(V, M, by="t")
df <- left_join(df, X, by="j") %>%
  select(t, i, j, v_x_1, v_x_2, v_x_3, v_p, x_1, everything())
df

# gumbel
e <- rgumbel(nrow(df))
head(e)


# compute utility
u = compute_indirect_utility(df, beta, sigma, mu, omega)
head(u)

# compute_choice
df_choice <- compute_choice(df, e, beta, sigma, mu, omega)
df_choice
summary(df_choice)

# compute_share
df_share <- compute_share(X, M, V, e, beta, sigma, mu, omega)
df_share
summary(df_share)

reg_res = lm(y ~ -1 +x_1 + x_2+ x_3 + p, data=df_share)
summary(reg_res)

# create V_mcmc
V_mcmc <- expand_grid(t=1:T, i=1:N) %>%
  mutate(v_x_1=rnorm(N*T), 
         v_x_2=rnorm(N*T),
         v_x_3=rnorm(N*T),
         v_p=rnorm(N*T)) %>%
  select(i, t, v_x_1, v_x_2, v_x_3, v_p)
V_mcmc

# e_mcmc
e_mcmc <- rgumbel(nrow(df))
head(e_mcmc)

# compute mcmc share 
df_share_mcmc <- compute_share(X, M, V_mcmc, e_mcmc, beta, sigma, mu, omega)
df_share_mcmc

# set parameter
theta <- c(beta, sigma, mu, omega)
theta_name <- c('beta1', 'beta2', 'beta3', 'sigma1', 'sigma2', 'sigma3',
                'mu', 'omega')

# mcmc
scale_vec <- seq(0.5, 1.5, 0.1)

s_time <- Sys.time()
registerDoParallel()
estimate <- foreach(i=1:length(theta), .combine = c) %do% {

  # result
  scaled_param = rep(0, length=length(scale_vec))
  res = rep(0, length=length(scale_vec))
  
  for (j in 1:length(scale_vec)) {
    theta_cp <- theta
    param <- theta[i]*scale_vec[j]
    theta_cp[i] <- param
    MSE <- NLLS_objective(theta_cp, df_share, X, M, V_mcmc, e_mcmc)
    scaled_param[j] <- param
    res[j] <- MSE
  }
  res_df <- tibble(param=scaled_param, objective=res)
  
  #plot
  g <- ggplot(res_df, aes(x=param, y=objective)) + geom_point() + xlab(theta_name[i])
  plot(g)
  ggsave(filename=paste("./tablefigures/",theta_name[i],".png", sep = ""),
         width=5,height=4,dpi=100)
  

}
e_time <- Sys.time()
print(e_time-s_time)

# optimize
res_optim = optim(par=theta, fn=NLLS_objective,df_share=df_share,
                  X=X, M=M, V_mcmc=V_mcmc, e_mcmc=e_mcmc,
                  method='Nelder-Mead')
res_df = tibble(true=theta, estimate=res_optim$par)
res_df
