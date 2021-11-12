## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library("metagam")
library("mgcv")
library("metafor")
library("ggplot2")
theme_set(theme_bw())

## ---- echo=FALSE--------------------------------------------------------------
simulate_data <- function(){
  dat1 <- data.frame(x = runif(100))
  dat1$y <- dat1$x + rnorm(100, sd = .7)
  
  dat2 <- dat1
  dat2$y <- -dat2$x + .1 * dat2$x^2 + rnorm(100, sd = .7)
  
  list(dat1 = dat1, dat2 = dat2)
}
set.seed(1)
dd <- simulate_data()

ggplot(dd$dat1, aes(x = x, y = y)) + 
  geom_point() + 
  ggtitle("Dataset 1")

ggplot(dd$dat2, aes(x = x, y = y)) + 
  geom_point() + 
  ggtitle("Dataset 2")


## -----------------------------------------------------------------------------
mod1 <- gam(y ~ s(x, bs = "cr"), data = dd$dat1)
summary(mod1)

mod2 <- gam(y ~ s(x, bs = "cr"), data = dd$dat2)
summary(mod2)

## -----------------------------------------------------------------------------
models <- list(strip_rawdata(mod1), strip_rawdata(mod2))
metafit <- metagam(models)

## -----------------------------------------------------------------------------
ggplot(metafit$meta_estimates, aes(x = x, y = estimate, ymin = ci.lb, ymax = ci.ub)) + 
  geom_line() + 
  geom_ribbon(alpha = .3) + 
  geom_hline(yintercept = 0, linetype = "dashed")

## -----------------------------------------------------------------------------
metafit$pvals

## ---- eval=FALSE--------------------------------------------------------------
#  library("metap")
#  allmetap(p = metafit$pvals$`p-value`, method = "all")

## ---- echo=FALSE--------------------------------------------------------------
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
grid <- with(dat, data.frame(x0 = seq(min(x0), max(x0), length = 200)))
masd <- getmasd(mod, grid, 10000, terms = "s(x0)")
(crit <- quantile(masd, prob = .95, type = 8))

## -----------------------------------------------------------------------------
fit <- predict(mod, newdata = grid, se.fit = TRUE)
dat <- data.frame(
  x0 = grid$x0,
  pred = fit$fit,
  ci_pt_lb = fit$fit + qnorm(.025) * fit$se.fit,
  ci_pt_ub = fit$fit + qnorm(.975) * fit$se.fit,
  ci_sim_lb = fit$fit - crit * fit$se.fit,
  ci_sim_ub = fit$fit + crit * fit$se.fit
)

ggplot(dat, aes(x = x0, y = pred)) + 
  geom_line() + 
  geom_ribbon(aes(ymin = ci_pt_lb, ymax = ci_pt_ub), alpha = .3) + 
  geom_ribbon(aes(ymin = ci_sim_lb, ymax = ci_sim_ub), alpha = .3)

## -----------------------------------------------------------------------------
set.seed(124)
datasets <- lapply(1:5, function(x) gamSim(scale = 5, verbose = FALSE))

models <- lapply(datasets, function(dat){
  model <- gam(y ~ s(x2, bs = "cr"), data = dat)
  strip_rawdata(model)
})

meta_analysis <- metagam(models, terms = "s(x2)", grid_size = 100,
                         nsim = 10000, ci_alpha = .05)

## -----------------------------------------------------------------------------
sim_conf <- meta_analysis$sim_ci
ptw_conf <- meta_analysis$meta_estimates
ggplot(sim_conf, aes(x = x2)) + 
  geom_line(aes(y = pred)) + 
  geom_ribbon(aes(ymin = ci_sim_lb, ymax = ci_sim_ub), alpha = .3) + 
  geom_ribbon(data = ptw_conf, aes(ymin = ci.lb, ymax = ci.ub), alpha = .3)

## -----------------------------------------------------------------------------
summary(meta_analysis)

## ---- eval=FALSE, echo=FALSE--------------------------------------------------
#  library(parallel)
#  cl <- makeCluster(10)
#  pvals <- parLapply(
#    cl = cl, X = 1:100, fun = function(x){
#      models <- lapply(1:5, function(i){
#        dat <- data.frame(x = runif(100), y = rnorm(100))
#        mod <- mgcv::gam(y ~ s(x, bs = "cr"), data = dat)
#        metagam::strip_rawdata(mod)
#        })
#      metagam::metagam(models, nsim = 1000, ci_alpha = .05)$meta_pval
#    }
#  )
#  stopCluster(cl)
#  png("figures/quantile_plot.png")
#  plot(qunif(seq(from = 0, to = 1, length.out = 100)),
#       sort(as.numeric(unlist(pvals))),
#       xlab = "Theoretical quantile",
#       ylab = "Data quantile")
#  dev.off()

## ---- eval=FALSE--------------------------------------------------------------
#  library(parallel)
#  cl <- makeCluster(10)
#  pvals <- parLapply(
#    cl = cl, X = 1:100, fun = function(x){
#      models <- lapply(1:5, function(i){
#        dat <- data.frame(x = runif(100), y = rnorm(100))
#        mod <- mgcv::gam(y ~ s(x, bs = "cr"), data = dat)
#        metagam::strip_rawdata(mod)
#        })
#      metagam::metagam(models, nsim = 1000, ci_alpha = .05)$meta_pval
#    }
#  )
#  stopCluster(cl)
#  
#  plot(qunif(seq(from = 0, to = 1, length.out = 100)),
#       sort(as.numeric(unlist(pvals))),
#       xlab = "Theoretical quantile",
#       ylab = "Data quantile")

