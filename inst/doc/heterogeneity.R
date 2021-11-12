## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library("metagam")

## -----------------------------------------------------------------------------
library("mgcv")
set.seed(1233)
shifts <- c(0, .5, 1, 0, -1)
datasets <- lapply(shifts, function(x) {
  ## Simulate data
  dat <- gamSim(scale = .1, verbose = FALSE)
  ## Add a shift
  dat$y <- dat$y + x * dat$x2^2
  ## Return data
  dat
})

## -----------------------------------------------------------------------------
models <- lapply(datasets, function(dat){
  b <- gam(y ~ s(x2, bs = "cr"), data = dat)
  strip_rawdata(b)  
})

## -----------------------------------------------------------------------------
meta_analysis <- metagam(models, type = "response")

## -----------------------------------------------------------------------------
plot(meta_analysis)

## -----------------------------------------------------------------------------
plot(meta_analysis) + 
  ggplot2::scale_colour_brewer(palette = "Set1")

## -----------------------------------------------------------------------------
plot_heterogeneity(meta_analysis)

## -----------------------------------------------------------------------------
plot_heterogeneity(meta_analysis, type = "p")

