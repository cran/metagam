## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- message=FALSE-----------------------------------------------------------
library(metagam)
library(mgcv)
## simulate three datasets
set.seed(123)
datasets <- lapply(1:3, function(x) gamSim(scale = 3, verbose = FALSE))

## -----------------------------------------------------------------------------
## Data location 1
fit1 <- gam(y ~ s(x0, k = 8, bs = "cr") + s(x1, bs = "cr") + s(x2, bs = "cr"), 
            data = datasets[[1]])

## Data location 2, use P-splines for the first and third term
fit2 <- gam(y ~ s(x0, bs = "ps") + s(x1, k = 20, bs = "cr") + s(x2, bs = "bs"), 
            data = datasets[[2]])

## Data location 3, use maximum likelihood for smoothing
fit3 <- gam(y ~ s(x0, bs = "cr") + s(x1, bs = "cr") + s(x2, bs = "cr"), 
            data = datasets[[3]], method = "ML")

## -----------------------------------------------------------------------------
## Data location 1
fit_no_raw1 <- strip_rawdata(fit1)

## Data location 2
fit_no_raw2 <- strip_rawdata(fit2)

## Data location 3
fit_no_raw3 <- strip_rawdata(fit3)

## -----------------------------------------------------------------------------
summary(fit_no_raw1)

## -----------------------------------------------------------------------------
models <- list(cohort1 = fit_no_raw1, 
               cohort2 = fit_no_raw2, 
               cohort3 = fit_no_raw3)

## -----------------------------------------------------------------------------
metafit <- metagam(models, terms = "s(x0)")
summary(metafit)

## ----out.width="60%", fig.align="center"--------------------------------------
plot(metafit)

