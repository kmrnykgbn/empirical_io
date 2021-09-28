library(tidyverse)
library(evd)

generate_data <- function (n,k, beta){
  df <- expand.grid(i=1:n, k=1:k)
  
  # add x
  df <- df %>%
    mutate(x=ifelse(k==1, 0, 1)) %>%
    dplyr::arrange(i)
  
  # add epsilon
  set.seed(1);
  df <- df %>%
    mutate(e=rgumbel(n*k))
  
  #plot
  plot(y=df$e, x=1:(n*k),
       type = "l", ylab = "e", xlab = "time")
  
  # add lattent
  df <- df %>%
    mutate(latent=beta*x+e,
           y=)
  
  # add y 
  df <- df %>%
    group_by(i) %>%
    mutate(y=ifelse(k[which(latent==max(latent))]==k, 1, 0))

  return(df)
}


loglikelihood_A1 <- function(beta, df) {
  df <- df %>%
    mutate(p=exp(beta*x)/(exp(beta*x[k==1]) + exp(beta*x[k==2]))) %>%
    filter(k==1)
  lf <- 0.0
  for (i in 1:1000) {
    add_lf <- df$y[[i]]*log(df$p[[i]]) + (1 - df$y[[i]])* log(1-df$p[[i]])
    lf <- lf + add_lf
  }
  return(lf)
} 

