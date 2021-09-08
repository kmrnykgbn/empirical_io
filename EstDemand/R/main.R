library(tidyverse)
library(tseries)
library(psych)
library(np)
library(gmm)

source("R/functions.R")

# parameters
beta_0 = 1
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

df <- tibble(i=1:n,
             t=1,
             k=rnorm(n,1,0.5),
             omega=1.401*rnorm(n,0,sigma_nu),
             wage=0.5,
             iota=rnorm(n, 0, 0.05),
             l=log_labor_choice(k, 0.5, omega, beta_0, beta_l, beta_k, sigma_eta),
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
  add_df <- tibble(i=1:n,
                   t=time,
                   k=(1-delta) * k_b + I_b,
                   omega=alpha * omega_b  + rnorm(n, 0, sigma_nu) ,
                   wage=0.5,
                   iota=rnorm(n, 0, 0.05),
                   l=log_labor_choice(k, 0.5, omega, beta_0, beta_l, beta_k, sigma_eta),
                   l_error=log_labor_choice_error(k, 0.5, omega, beta_0,
                                                  beta_l, beta_k, iota, sigma_eta),
                   I=investment_choice(k, omega, 0.1, delta),
                   eta=rnorm(n,0,sigma_eta),
                   y=log_prodction(l, k, omega, eta, beta_0, beta_l, beta_k),
                   y_error=log_prodction(l_error, k, omega, eta, beta_0, beta_l, beta_k))
 df_T <- bind_rows(df_T, add_df)
}

df_T
describe(df_T)


# regression
reg_m = lm(y_error~l_error+k, df_T)
summary(reg_m)

df_T = df_T %>%
  mutate(diff_y = lag(y, default=0),
         diff_k = lag(k, default=0),
         diff_l = lag(l, default=0))

#within
within_reg = lm(diff_y~-1 + diff_l + diff_k, df_T)
summary(within_reg)

#partial regression(l_error_version)
bw <- npplregbw(formula=y~1+l_error|k+I, data=df_T, itmax=5, nmulti=1)
pl_error <- npplreg(bws = bw, txdat = as.data.frame(df_T[, "l_error"]), tydat = df_T$y, tzdat = as.data.frame(df_T[, c("k", "I")]))
#pl_error <- lm(y ~ 1 + l_error + poly(k, I, degree = 2), data = df_T)
summary(pl_error)

#partial regression(l_version)
bw <- npplregbw(formula=y~1+l|k+I, data=df_T, itmax=5, nmulti=1) # collaboration of bandwidth to reduce calc time
pl <- npplreg(bws = bw, txdat = as.data.frame(df_T[, "l"]), tydat = df_T$y, tzdat = as.data.frame(df_T[, c("k", "I")]))
#pl <- lm(y ~ 1 + l + poly(k, I, degree = 10), data = df_T)
summary(pl)

#df_1st
df_1st <- expand.grid(j=1:1000, t=1:10) %>%
  mutate(y_error_tilde)
