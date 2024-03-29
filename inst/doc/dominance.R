## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
num.datasets <- 5

## -----------------------------------------------------------------------------
library("metagam")

## -----------------------------------------------------------------------------
## simulate datasets
set.seed(123)
datasets <- lapply(
  seq_len(num.datasets), 
  function(x) mgcv::gamSim(scale = x, verbose = FALSE))

## -----------------------------------------------------------------------------
df <- datasets[[1]]
df[df$x2<0.2,] <- NA
datasets[[1]] <- df

## -----------------------------------------------------------------------------
df <- datasets[[2]]
df[df$x2 > 0.8, ] <- NA
datasets[[2]] <- df

## -----------------------------------------------------------------------------
## fit a generalized additive model to each dataset separately
models <- lapply(datasets, function(dat){
  ## Full fit using mgcv
  gamfit <- mgcv::gam(y ~ s(x0, bs = "cr") + s(x1, bs = "cr") + s(x2, bs = "cr"), data = dat)
  ## Extract the necessary components for performing a meta-analysis
  ## This removes all subject-specific data
  strip_rawdata(gamfit)
})

## -----------------------------------------------------------------------------
names(models) <- c("A", "B", "C", "D", "E")

## -----------------------------------------------------------------------------
meta_analysis <- metagam(models, grid_size = 50, terms = "s(x2)")

## -----------------------------------------------------------------------------
plot_dominance(meta_analysis)

