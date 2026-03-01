library(marp)
library(stringr)
library(lubridate)
library(dplyr)
library(edfReader)
library(tidyverse)

# Run this in a different R session (so completely different environment) from marp-seizure.R

data_fuller <- read.table("data/FULLER2.DAT", skip = 25)

dat <- data_fuller$V1

set.seed(2026)


hist(dat, prob = TRUE)

# set parameters
m <- 100 # number of iterations for MLE optimization
t <- quantile(dat, probs = seq(0.1, 0.9, length.out = 7))
B <- 20 # number of bootstraps use 1000 for paper
BB <- 20 # number of double-bootstrapps use 200 for paper
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

ci
res


# == Placeholder histogram graph==

pdf("graphs/multipanel_figure.pdf")
par(mfrow = c(2, 2))


hist(dat,
     probability = TRUE,
     xlab = "Time to failure",
     col = 'white',
     border = 'black',
     main = '(a)')

hist(dat,
     probability = TRUE,
     xlab = "Time to failure",
     col = 'white',
     border = 'black',
     main = '(b)')

hist(dat,
     probability = TRUE,
     xlab = "Time to failure",
     col = 'white',
     border = 'black',
     main = '(c)')

hist(dat,
     probability = TRUE,
     xlab = "Time to failure",
     col = 'white',
     border = 'black',
     main = '(d)')

dev.off()




# --- Empirical log-hazard vs estimates graph---

pdf("graphs/loghaz_plot_fuller.pdf")

# 1) KDE for f(t); ensure coverage across t to avoid NA from extrapolation
d <- density(dat, from = min(t), to = max(t))
# Interpolate, with rule=2 to allow linear extrapolation at the ends if needed
f_t <- approx(d$x, d$y, xout = t, rule = 2)$y

# 2) Survival S(t) via ECDF: S(t) = P(T > t) = 1 - F(t)
Fhat <- ecdf(dat)
# use right-continuity convention via Fhat(t) then 1 - Fhat(t)
S_t <- 1 - Fhat(t)

# 3) Protect against 0 or negative values before log
eps <- .Machine$double.eps
f_t_pos <- pmax(f_t, eps)
S_t_pos <- pmax(S_t, eps)

# 4) Compute log-hazard directly to avoid intermediate 'haz_emp'
# log h(t) = log f(t) - log S(t)
loghaz_emp <- log(f_t_pos) - log(S_t_pos)

# --- Plot ---
y_all <- c(res$haz_best, res$haz_aic, loghaz_emp)
y_rng <- range(y_all, na.rm = TRUE)

plot(t, res$haz_best,
     type = "l",
     lwd = 3,
     ylim = y_rng,
     xlab = "Time",
     ylab = "Log hazard",
     main = "Log-Hazard Comparison")

lines(t, res$haz_aic, lwd = 3, lty = 2)
lines(t, loghaz_emp, lwd = 3, lty = 3)

legend("topleft",
       legend = c("Best model", "AIC-averaged", "Empirical"),
       lty = c(1, 2, 3),
       lwd = 3,
       bty = "n")

dev.off()


length(t)
length(lower)
length(upper)






# ------Confidence intervals graph---------------

pdf("graphs/CI_plot_fuller.pdf")

# Extract series
haz   <- res$haz_aic
lower <- ci$student_CI$haz_lower_ma
upper <- ci$student_CI$haz_upper_ma


# Define a semi-transparent gray for the CI band
ci_col <- rgb(0, 0, 0, 0.15)

# Plot an empty canvas first to avoid drawing the line twice
plot(t, haz,
     type = "n",
     lwd  = 3,
     ylim = range(lower, upper, haz, finite = TRUE),
     xlab = "Time",
     ylab = "Log hazard",
     main = "Model-averaged Log-Hazard with 95% CI")

# Add shaded CI band (upper on the left path, lower on the return path)
polygon(c(t, rev(t)),
        c(upper, rev(lower)),
        col = ci_col,
        border = NA)

# Draw the line on top
lines(t, haz, lwd = 3)

# Legend that visually matches the line and the shaded band
legend("topleft",
       legend = c("AIC-averaged", "95% CI"),
       lty    = c(1, NA),
       lwd    = c(3, NA),
       col    = c("black", NA),
       fill   = c(NA, ci_col),
       border = NA,
       bty    = "n")

dev.off()