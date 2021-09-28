library(tidyverse)
library(tseries)
library(psych)
library(np)
library(gmm)

source("R/functions.R")
sourceCpp("src/function.cpp")

# parameters
beta_0 = 1.0
beta_l = 0.2
beta_k = 0.7
alpha = 0.7
sigma_eta = 0.2
sigma_nu = 0.5
sigma_w = 0.1
delta = 0.05

# create data frame
n = 1000
t = 10 
set.seed(1)

df <- tibble(j=1:n,
             t=1,
             k=rnorm(n,1,0.5),
             omega=2*alpha*rnorm(n,0,sigma_nu),
             wage=0.5,
             iota=rnorm(n, 0, 0.05),
             l=log_labor_choice(k, wage, omega, beta_0, beta_l, beta_k, sigma_eta),
             l_error=log_labor_choice_error(k, 0.5, omega, beta_0,
                                            beta_l, beta_k, iota, sigma_eta),
             I=investment_choice(exp(k), omega, 0.1, delta),
             eta=rnorm(n,0,sigma_eta),
             y=log_prodction(l, k, omega, eta, beta_0, beta_l, beta_k),
             y_error=log_prodction(l_error, k, omega, eta, beta_0, beta_l, beta_k))

df
df_T <- df
for (time in 2:t){
  omega_b = df_T[which(df_T[["t"]] == time-1), ]$omega
  k_b = df_T[which(df_T[["t"]] == time-1), "k"]$k
  I_b = df_T[which(df_T[["t"]] == time-1), "I"]$I
  add_df <- tibble(j = 1:n,
                   t = time,
                   k = log((1-delta)*exp(k_b)+I_b),
                   omega = alpha * omega_b  + rnorm(n, 0, sigma_nu) , 
                   wage = 0.5,
                   iota = rnorm(n, 0, 0.05),
                   l = log_labor_choice(k, wage, omega, beta_0, beta_l, beta_k, sigma_eta),
                   l_error = log_labor_choice_error(k, 0.5, omega, beta_0,
                                                  beta_l, beta_k, iota, sigma_eta),
                   I = investment_choice(exp(k), omega, 0.1, delta),
                   eta = rnorm(n,0,sigma_eta),
                   y = log_prodction(l, k, omega, eta, beta_0, beta_l, beta_k),
                   y_error = log_prodction(l_error, k, omega, eta, beta_0, beta_l, beta_k))
 df_T <- bind_rows(df_T, add_df)
}

df_T
describe(df_T)[c("mean", "sd", "min", "max")]


# regression
reg_m = lm(y_error~l_error+k, df_T)
summary(reg_m)

df_within = df_T %>%
  select(t, j, y_error, l_error, k) %>%
  group_by(j) %>%
  arrange(j, t) %>%
  mutate(dy_error = y_error - lag(y_error, default=0),
         dk = k - lag(k, default=0),
         dl_error = l_error - lag(l_error, default=0)) %>%
  ungroup()

#within
within_reg = lm(dy_error~-1 + dl_error +dk, df_within)
summary(within_reg)

#partial regression(l_error_version)
#bw <- npplregbw(formula=y~1+l_error+k+I|k+I, data=df_T)
#bw <- npplregbw(formula=y~1+l_error|k+I, data=df_T, itmax=100, nmulti=1)
pl_error <- npplreg(bws = bw, txdat = as.data.frame(df_T[, c("l_error", "k", "I")]), tydat = df_T$y, tzdat = as.data.frame(df_T[, c("k", "I")]))
summary(pl_error)

#partial regression(l_version)
#bw <- npplregbw(formula=y~1+l|k+I, data=df_T, itmax=5, nmulti=1) # collaboration of bandwidth to reduce calc time
#pl <- npplreg(bws = bw, txdat = as.data.frame(df_T[, "l"]), tydat = df_T$y, tzdat = as.data.frame(df_T[, c("k", "I")]))
#pl <- lm(y ~ 1 + l + poly(k, I, degree = 10), data = df_T)
#summary(pl)

# plot y_error and fitted
df_y_error = tibble(fitted=fitted(pl_error),
                    actual=df_T$y_error)
g <- ggplot(df_y_error, aes(x=fitted, y=actual)) + geom_point()
plot(g)

# second_step estimate
beta_l_hat <- coef(pl_error) 
phi_1_t <- fitted(pl_error) -beta_l_hat*df$l_error
lag_df <-  df_T %>% 
  mutate(phi_1_t = phi_1_t) %>%
  group_by(j) %>% 
  arrange(j) %>%
  mutate(lag_phi = lag(phi_1_t, default=0)) %>%
  ungroup()


df_T_1st <- expand_grid(i=1:n, t=1:t) %>%
  mutate(y_error_tilde = df_T$y_error-beta_l_hat*df_T$l_error,
         phi_t_1=lag_df$lag_phi) 
df_T_1st

# get moment
obj_val <- objective_OP_2nd(alpha, beta_0, beta_k, df_T, df_T_1st)
obj_val

# draw_obj_func_graph
param_seq <- seq(0, 1, 0.1)
res_alpha = rep(0, length=length(param_seq))
res_beta0 = rep(0, length=length(param_seq))
res_betak = rep(0, length=length(param_seq))
for (i in 1:length(param_seq)) {
  res_alpha[i] = objective_OP_2nd(param_seq[i], beta_0, beta_k, df_T, df_T_1st)
  res_beta0[i] = objective_OP_2nd(alpha, param_seq[i], beta_k, df_T, df_T_1st)
  res_betak[i] = objective_OP_2nd(alpha, beta_0, param_seq[i], df_T, df_T_1st)
}

g <- ggplot(tibble(alpha=param_seq, Objective=res_alpha),
            aes(x=alpha, y=Objective)) + geom_point()
plot(g)

g <- ggplot(tibble(beta0=param_seq, Objective=res_beta0),
            aes(x=beta0, y=Objective)) + geom_point()
plot(g)

g <- ggplot(tibble(betak=param_seq, Objective=res_betak),
            aes(x=betak, y=Objective)) + geom_point()
plot(g)

