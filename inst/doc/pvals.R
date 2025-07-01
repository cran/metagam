## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 4
)

## ----setup--------------------------------------------------------------------
library("metagam")
library("mgcv")
library("metafor")

## ----echo=FALSE---------------------------------------------------------------
simulate_data <- function(){
  dat1 <- data.frame(x = runif(100))
  dat1$y <- dat1$x + rnorm(100, sd = .7)
  
  dat2 <- dat1
  dat2$y <- -dat2$x + .1 * dat2$x^2 + rnorm(100, sd = .7)
  
  list(dat1 = dat1, dat2 = dat2)
}
set.seed(1)
dd <- simulate_data()


for(i in seq_along(dd)){
  plot(dd[[i]]$x, dd[[i]]$y, 
       xlab = "x", ylab = "y", main = paste("Dataset", i))
}




## -----------------------------------------------------------------------------
mod1 <- gam(y ~ s(x, bs = "cr"), data = dd$dat1)
summary(mod1)

mod2 <- gam(y ~ s(x, bs = "cr"), data = dd$dat2)
summary(mod2)

## -----------------------------------------------------------------------------
models <- list(strip_rawdata(mod1), strip_rawdata(mod2))
metafit <- metagam(models)

## -----------------------------------------------------------------------------
plot(metafit, ci = "pointwise", legend = FALSE)

## -----------------------------------------------------------------------------
metafit$pvals

## ----eval=FALSE---------------------------------------------------------------
# library("metap")
# allmetap(p = unlist(lapply(metafit$pvals, function(x) as.data.frame(x)[, "p-value"])),
#          method = "all")

## ----echo=FALSE---------------------------------------------------------------
structure(list(p = list(logitp = 5.61819601706602e-07, maximump = 1.25713383795263e-06, 
    meanp = NA_real_, meanz = 6.75835070011165e-08, minimump = 7.22516917847225e-06, 
    sumlog = 8.23242533258099e-08, sump = 6.32623953067482e-07, 
    sumz = structure(4.8107525303013e-08, .Dim = c(1L, 1L))), 
    valid = list(logitp = 2L, maximump = 2L, meanp = 2L, meanz = 2L, 
        minimump = 2L, sumlog = 2L, sump = 2L, sumz = 2L), eponym = c("", 
    "", "", "", "Tippett", "Fisher", "Edgington", "Stouffer")), row.names = c("logitp", 
"maximump", "meanp", "meanz", "minimump", "sumlog", "sump", "sumz"
), class = c("allmetap", "data.frame"))

## -----------------------------------------------------------------------------
mod1 <- lm(y ~ x, data = dd$dat1)
mod2 <- lm(y ~ x, data = dd$dat2)

estimates <- c(coef(mod1)[["x"]], coef(mod2)[["x"]])
sampling_variances <- c(vcov(mod1)[["x", "x"]], vcov(mod2)[["x", "x"]])
rma(estimates, vi = sampling_variances)

## -----------------------------------------------------------------------------
set.seed(123)
dat <- gamSim(verbose = FALSE)
mod <- gam(y ~ s(x0, bs = "cr"), data = dat)
plot(mod)

## -----------------------------------------------------------------------------
newdat <- with(dat, data.frame(x0 = seq(min(x0), max(x0), length = 200)))
masd <- getmasd(mod, newdat = newdat, nsim = 1000, term = "s(x0)")
(crit <- quantile(masd, prob = .95, type = 8))

## -----------------------------------------------------------------------------
fit <- predict(mod, newdata = newdat, se.fit = TRUE)
dat <- data.frame(
  x0 = newdat$x0,
  pred = fit$fit,
  ci_pt_lb = fit$fit + qnorm(.025) * fit$se.fit,
  ci_pt_ub = fit$fit + qnorm(.975) * fit$se.fit,
  ci_sim_lb = fit$fit - crit * fit$se.fit,
  ci_sim_ub = fit$fit + crit * fit$se.fit
)

plot(dat$x0, dat$pred, type = "l",
     ylim = range(dat$ci_sim_lb, dat$ci_sim_ub), 
     xlab = "x0", ylab = "s(x0)")
polygon(
  x = c(rev(dat$x0), dat$x0), y = c(rev(dat$ci_sim_ub), dat$ci_sim_lb),
  col = "gray80", border = NA
)
polygon(
  x = c(rev(dat$x0), dat$x0), y = c(rev(dat$ci_pt_ub), dat$ci_pt_lb),
  col = "gray60", border = NA
)
lines(dat$x0, dat$pred)

legend(x = "bottom", legend = c("Pointwise", "Simultaneous"),
       fill = c("gray60", "gray80"))


## -----------------------------------------------------------------------------
set.seed(124)
datasets <- lapply(1:5, function(x) gamSim(scale = 5, verbose = FALSE))

models <- lapply(datasets, function(dat){
  model <- gam(y ~ s(x2, bs = "cr"), data = dat)
  strip_rawdata(model)
})
names(models) <- paste("Model", letters[1:5])

meta_analysis <- metagam(models, terms = "s(x2)", grid_size = 100,
                         nsim = 10000, ci_alpha = .05)

## ----fig.height=6, fig.width=6------------------------------------------------
plot(meta_analysis, ci = "both", legend = TRUE)

## -----------------------------------------------------------------------------
summary(meta_analysis)

## ----eval=FALSE, echo=FALSE---------------------------------------------------
# library(parallel)
# cl <- makeCluster(10)
# pvals <- parLapply(
#   cl = cl, X = 1:100, fun = function(x){
#     models <- lapply(1:5, function(i){
#       dat <- data.frame(x = runif(100), y = rnorm(100))
#       mod <- mgcv::gam(y ~ s(x, bs = "cr"), data = dat)
#       metagam::strip_rawdata(mod)
#       })
#     fit <- metagam::metagam(models, nsim = 1000, ci_alpha = .05)
#     fit$simulation_results$`s(x)`$pval
#   }
# )
# stopCluster(cl)
# png("figures/quantile_plot.png")
# plot(qunif(seq(from = 0, to = 1, length.out = 100)),
#      sort(as.numeric(unlist(pvals))),
#      xlab = "Theoretical quantile",
#      ylab = "Data quantile")
# abline(0, 1)
# dev.off()

## ----eval=FALSE---------------------------------------------------------------
# library(parallel)
# cl <- makeCluster(10)
# pvals <- parLapply(
#   cl = cl, X = 1:100, fun = function(x){
#     models <- lapply(1:5, function(i){
#       dat <- data.frame(x = runif(100), y = rnorm(100))
#       mod <- mgcv::gam(y ~ s(x, bs = "cr"), data = dat)
#       metagam::strip_rawdata(mod)
#       })
#     fit <- metagam::metagam(models, nsim = 1000, ci_alpha = .05)
#     fit$simulation_results$`s(x)`$pval
#   }
# )
# stopCluster(cl)
# 
# plot(qunif(seq(from = 0, to = 1, length.out = 100)),
#      sort(as.numeric(unlist(pvals))),
#      xlab = "Theoretical quantile",
#      ylab = "Data quantile")
# abline(0, 1)

