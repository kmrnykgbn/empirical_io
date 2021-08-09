source("R/functions.R")

#generate data
df <- generate_data(1000, 2, 0.2)

# estmate MLE
beta <- seq(0, 1, 0.1)
i <- 1
Loglikelihood <- rep(0, length=length(beta))
for (b in beta) {
  Loglikelihood[i] <- loglikelihood_A1(b, df)  
  i <- i + 1
}
#plot the log likelihood
ggplot2::qplot(beta, Loglikelihood)

# maximize the loglikelihood
MLE = optim(par=runif(1,min=-1,max=+1), fn=loglikelihood_A1,df=df, method='Brent', lower = -1,
      upper = 1, control = list(fnscale = -1))
MLE

