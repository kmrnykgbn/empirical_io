library(tidyverse)

set.seed(1)

J <- 10
K <- 3
T <- 100
N <- 500
L <- 500

beta <- rnorm(K)
beta[1] <- 4
sigma <- abs(rnorm(K))
mu <- 0.5
omega <- 1

price_xi <- 1
sd_x <- 2
sd_xi <- 0.5
sd_c <- 0.05
ssd_p <- 0.05

j <- 1:J
X <- tibble(
  j = j,
  x_1 = 1,
  x_2 = rnorm(J, mean = 0, sd = sd_x),
  x_3 = rnorm(J, mean = 0, sd = sd_x)
)
X <- rbind(0, X)

t <- 1:T
M <- expand_grid(t, j)
M <- mutate(M, 
            xi =rnorm(T*J, mean = 0, sd = sd_xi),
            c = rlnorm(T*J, mean = 0, sd = sd_c),
            p = 0)
M <- M[, c(2, 1, 3, 4, 5)]
M <- M %>%
  group_by(t) %>%
  dplyr::sample_n(as.integer(rdunif(1, 1, J))) %>%
  ungroup()
outside <- tibble(j = 0, 
                  t = 1:T,
                  xi = 0,
                  c = 0,
                  p=0)
M <- rbind(M, outside)
M <- M[order(M$t, M$j), ]

i <- 1:N
V <- expand_grid(t, i)
V <- V[, c(2,1)]
V <- mutate(V,
            v_x_1 = rnorm(N*T),
            v_x_2 = rnorm(N*T),
            v_x_3 = rnorm(N*T),
            v_p = rnorm(N*T))

compute_indirect_utility <- function(df, beta, sigma, mu, omega){
  data <- mutate(df, 
                 beta_1 = beta[1] + sigma[1]*v_x_1,
                 beta_2 = beta[2] + sigma[2]*v_x_2,
                 beta_3 = beta[3] + sigma[3]*v_x_3,
                 alpha = -exp(mu+omega*v_p),
                 u = beta_1*x_1 + beta_2*x_2 + beta_3*x_3 + alpha*p + xi)
  return(data$u)
}

compute_choice_smooth <- function(df, beta, sigma, mu, omega){
  data <- mutate(df, u = compute_indirect_utility(df, beta, sigma, mu, omega))
  data_ <- mutate(data, exp_u = exp(u))
  data_ <- data_ %>%
    group_by(i, t) %>%
    mutate(q_sum = sum(exp_u)) %>%
    ungroup()
  data_ <- mutate(data_, q=exp_u/(q_sum))
  data <- mutate(data, q=data_$q)
  return(data)
}

compute_share_smooth <- function(data){
  df <- data %>%
    group_by(t, j) %>%
    mutate(s=sum(q)/N) %>%
    ungroup()
  return(df)
}

compute_derivative_share_smooth <- function(X, M, V, beta, sigma, mu, omega){
  df_ <- left_join(M, X, by="j")
  df <- left_join(V, df_, by="t")
  df <- df[, c(2,1,7,3,4,5,6,11,12,13,8,9,10)]
  
  u <- compute_indirect_utility(df, beta, sigma, mu, omega)
  df <- mutate(df, u=u)
  
  df <- df %>%
    group_by(t, i) %>%
    mutate(sigma = exp(u)/sum(exp(u))) %>%
    ungroup()
  
  df <- mutate(df, alpha = -exp(mu + omega*v_p))
  #df <- mutate(df, alpha_sigma = alpha*sigma)
  
  T <- length(unique(df$t))
  N <- length(unique(df$i)) 
  
  result <- list()
  
  for(l in 1:T){
    dd <- filter(df, t == l & j > 0)
    num_market <- length(unique(dd$j))
    
    market <- unique(dd$j)
    
    sigma <- data.frame()
    
    for(m in 1:num_market){
      sigma <- rbind(sigma, t(dd$sigma[dd$j==market[m]]))
    }
    
    sigma <- as.matrix(sigma)
    alpha <- as.matrix(dd$alpha[dd$j==market[1]])
    kata <- matrix(1, num_market, 1)
    alpha_mat <- kata %*% t(alpha)
    
    alpha_sigma <- sigma * alpha_mat
    ass_ <- alpha_sigma %*% t(sigma)
    ass <- -1/N * ass_
    
    as_one <- matrix(1, N, 1)
    
    as_ <- (alpha_sigma %*% as_one) * 1/N
    
    if(nrow(as_)>1){
      as <- diag(as.vector(as_))
    }
    else{
      as <- as_
    }
    
    res <- as + ass
    
    result <- c(result, list(res))
  }
  
  return(result)
}

derivative_share_smooth <- compute_derivative_share_smooth(X, M, V, beta, sigma, mu, omega)

Delta <- list()
for(l in 1:T){
  dd <- filter(df, t==l & j>0)
  num_market <- length(unique(dd$j))
  
  mat <- diag(1, nrow=num_market, ncol=num_market)
  Delta <- c(Delta, list(mat))
}



lambda <- 1e-6
p <- M[M$j>0, "p"]
logp <- log(rep(1, dim(p)[1]))






