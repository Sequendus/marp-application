library(marp)
library(stringr)
library(lubridate)
library(dplyr)
library(edfReader)
library(tidyverse)
library(grid)

# ==== Ignore this unused code for subject 24 ==== 
# source("chb24-pipeline.R")
# edf_dir <- "data/seizure-data/chb24"   # directory containing chb24_*.edf
# seizure_df <- parse_seizure_text_with_edf(txt, edf_dir, tz = "UTC")



# ==== Pre-process data ====

# get function to pre-process data
source("pipelines/general-seizure-pipeline.R")

# read seizure data from subject 15
txt <- readLines("data/seizure-data/chb15-summary.txt")

# Run the parser to get data into suitable format (vector) for input into marp functions
seizure_df <- parse_seizure_text(txt)

# Visualise data
ggplot(seizure_df, aes(x = global_start_sec)) +
  geom_rug(
    sides = "b",
    color = "black",
    size = 0.3,                
    length = unit(10, "pt")     
  ) +
  labs(
    title = "Seizure Onsets (Global Time)",
    x = "Global Time (seconds since first file)",
    y = NULL
  ) +
  theme_minimal()


# Optional: histogram of inter event times (times between events, not global times)
# before removing clusters
global_secs1 <- seizure_df$global_start_sec
dat1 <- diff(global_secs1)
summary(dat1)
hist(dat1)


# Seizures often occur in short bursts (clusters), meaning the events within
# these bursts are not independent of each other.
# However, renewal process models require the inter-event times to be i.i.d.
# To approximate independence, we remove closely spaced events by collapsing
# any seizures that occur within a minimum gap threshold (here: 1800 sec).\
df_renewal <- collapse_clusters(seizure_df, min_gap_sec = 1800)


# Visualise dat after removing clusters
ggplot(df_renewal, aes(x = global_start_sec, y = 0)) +
  geom_point(size = 3, color = "steelblue") +
  geom_rug(sides = "b") +
  labs(
    title = "Seizure Onsets (Global Time)",
    x = "Global Time (seconds since first file)",
    y = NULL
  ) +
  theme_minimal()


# Extract inter event times (times between events, not global times)
# and visualise inter event times in a histogram
global_secs <- df_renewal$global_start_sec
dat <- diff(global_secs)
summary(dat)
hist(dat, breaks=8)



# ==== Start analysis ====
set.seed(2026)


# set parameters
m <- 80 # number of iterations for MLE optimization
t <- quantile(dat, probs = seq(0.1, 0.9, length.out = 4))
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
# NOTE: this might take a while, so feel free to skip it
# as the graphs can be done without this
ci <- marp::marp_confint(dat,m,t,B,BB,alpha,y, 6)
ci

summary(dat)
hist(dat)


# ==== graphs for paper ====

pdf("loghaz_plot_seizure.pdf")

# --- Empirical log-hazard ---

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

