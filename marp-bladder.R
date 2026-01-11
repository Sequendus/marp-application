
library(marp)
library(tidyverse)

intervals <- bladder1 |>
  filter(id == 46) |>
  arrange(start) |>
  mutate(interval = stop - start) |>
  pull(interval)

# dat <- rgamma(20,3,0.01)


dat <- intervals
dat

set.seed(2026)

# set parameters
m <- 80 # number of iterations for MLE optimization
t <- as.numeric(quantile(dat, probs = c(0.3, 0.5, 0.7)))
# time intervals
B <- 300 # number of bootstraps
BB <- 60 # number of double-bootstrapps
alpha <- 0.05 # confidence level
y <- mean(dat) 
# model_gen <- 2 # specifying the data generating model (if known)

# step one: fitting differnt renewal models
res1 <- marp::poisson_rp(dat,t,y)
res2 <- marp::gamma_rp(dat,t,m,y)
res3 <- marp::loglogis_rp(dat,t,m,y)
res4 <- marp::weibull_rp(dat,t,m,y)
res5 <- marp::lognorm_rp(dat,t,y)
res6 <- marp::bpt_rp(dat,t,m,y)

# step two: model selection and obtain model-averaged estimates
res <- marp::marp(dat,t,m,y)
res

# step three: construct different confidence intervals (including model-averaged CIs)
ci <- marp::marp_confint(dat,m,t,B,BB,alpha,y, 6)
ci

