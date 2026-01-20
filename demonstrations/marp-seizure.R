library(marp)
library(stringr)
library(lubridate)
library(dplyr)
library(edfReader)
library(tidyverse)


source("general-seizure-pipeline.R")
# source("chb24-pipeline.R")


# read summary text
txt <- readLines("data/seizure-data/chb15-summary.txt")

# edf_dir <- "data/seizure-data/chb24"   # directory containing chb24_*.edf


# Run the parser
seizure_df <- parse_seizure_text(txt)
# seizure_df <- parse_seizure_text_with_edf(txt, edf_dir, tz = "UTC")


df_renewal <- collapse_clusters(seizure_df, min_gap_sec = 1800)


# Numeric seconds since first file:
global_secs1 <- seizure_df$global_start_sec
dat1 <- diff(global_secs1)
summary(dat1)
hist(dat1)

global_secs <- df_renewal$global_start_sec
dat <- diff(global_secs)
summary(dat)
hist(dat, breaks=8)

ggplot(seizure_df, aes(x = global_start_sec, y = 0)) +
  geom_point(size = 3, color = "steelblue") +
  geom_rug(sides = "b") +
  labs(
    title = "Seizure Onsets (Global Time)",
    x = "Global Time (seconds since first file)",
    y = NULL
  ) +
  theme_minimal()


ggplot(df_renewal, aes(x = global_start_sec, y = 0)) +
  geom_point(size = 3, color = "steelblue") +
  geom_rug(sides = "b") +
  labs(
    title = "Seizure Onsets (Global Time)",
    x = "Global Time (seconds since first file)",
    y = NULL
  ) +
  theme_minimal()

set.seed(2026)


# set parameters
m <- 80 # number of iterations for MLE optimization
t <- as.numeric(quantile(dat, probs = c(0.25, 0.5, 0.75))) # time intervals
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

summary(dat)

