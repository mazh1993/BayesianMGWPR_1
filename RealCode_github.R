##### Bayesian MGWPR ####

set.seed(1007)
library(INLA)
library(readxl)

data_use <- read_xlsx(path = "../georgia_2022.xlsx", sheet = "data_use", col_names = TRUE)
load("../GA_map.Rdata")
data_use_0 <- data_use[, c(1:3,5:8)]
data_use_0$County <- as.character(data_use_0$County)
names(data_use_0) <- c("county", "deaths", "PM2.5", 
                       "commute-driving","uninsured", "inequity", "unemployment")

library(ggplot2)
GAmap <- fortify(mapdata)

library(dplyr)
long1 <- aggregate(GAmap$long, by = list(subregion = GAmap$subregion), mean)
lat1 <- aggregate(GAmap$lat, by = list(subregion = GAmap$subregion), mean)
GAdata <- data.frame("SR_ID" = c(1:159), "Latitude" = lat1$x, "Longitud" = long1$x)
GAdata <- cbind(GAdata, data_use_0)


library(INLA)
inla.rgeneric.beta.model <-
  function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const", "log.prior", "quit"),
           theta = NULL) {
    envir = parent.env(environment())
  
    interpret.theta <- function() {
      return(list(prec = exp(theta[1L]),
                  rho = 1/(1+exp(-theta[2L])),
                  mu0 = theta[3L]))
    }
    graph <- function() {
      return (inla.as.sparse(matrix(1, n, n)))
    }
    Q <- function() {
      require(Matrix)
      param <- interpret.theta()
      cov.mat <- 1/param$prec * exp(-D/(param$rho))  
      prec.mat <- solve(cov.mat)
      return(inla.as.sparse(prec.mat))
    }
    mu <- function() {
      return (rep(interpret.theta()$mu, n))
    }
    log.norm.const <- function() {
      return (numeric())
    }
    log.prior <- function() {
      param = interpret.theta()
      res <- dgamma(param$prec, 0.01, 0.01, log = T) + log(param$prec) +
        log(1) + log(param$rho) + log(1 - param$rho) +
        dnorm(param$mu0, 0, 1, log = T)
      return(res)
    }
    initial <- function() {
      return(c(0, 0, 0))
    }
    quit <- function() {
      return(invisible())
    }
    
    if (!length(theta) || is.null(theta)) theta = initial()
    
    val <- do.call(match.arg(cmd), args = list())
    return(val)
  }


data.inla <- GAdata
data.inla$y <- as.vector(data.inla$deaths)
data.inla$x0 <- rep(1, length(data.inla$y))
data.inla$x1 <- as.vector(scale(data.inla$PM2.5))
data.inla$x2 <- as.vector(scale(data.inla$`commute-driving`))
data.inla$x3 <- as.vector(scale(data.inla$uninsured))


Dist0 <- as.matrix(dist(data.frame(data.inla$Longitud, data.inla$Latitude)))
Dmax <- max(Dist0); Dmax
Dist01 <- Dist0/Dmax
D <- Dist01
NN <- nrow(D); NN

beta.model <- inla.rgeneric.define(inla.rgeneric.beta.model,
                                   D = D,
                                   n = NN)

f.beta <- y ~ -1 + f(idx1, x1, model = beta.model, n = NN) + 
  f(idx2, x2, model = beta.model, n = NN) + f(idx3, x3, model = beta.model, n = NN) +
  f(idx0, x0, model = beta.model, n = NN)

m.beta <- inla(f.beta,
               safe = FALSE, 
               data = data.frame(data.inla, idx1 = 1:NN, idx2 = 1:NN, idx3 = 1:NN,
                                 idx0 = 1:NN),
               family = "poisson", 
               control.predictor = list(compute = T),
               control.compute = list(cpo = TRUE, waic = TRUE, dic = TRUE),
               verbose = T)


# parameter estimates
beta1.mean <- m.beta$summary.random$idx1$mean
beta2.mean <- m.beta$summary.random$idx2$mean
beta3.mean <- m.beta$summary.random$idx3$mean
beta0.mean <- m.beta$summary.random$idx0$mean


c(mean(beta1.mean), mean(beta2.mean), mean(beta3.mean), mean(beta0.mean))
c(sd(beta1.mean), sd(beta2.mean), sd(beta3.mean), sd(beta0.mean))

beta1.sd <- m.beta$summary.random$idx1$sd
beta2.sd <- m.beta$summary.random$idx2$sd
beta3.sd <- m.beta$summary.random$idx3$sd
beta0.sd <- m.beta$summary.random$idx0$sd
c(mean(beta1.sd), mean(beta2.sd), mean(beta3.sd), mean(beta0.sd))

beta1.q1 <- m.beta$summary.random$idx1$`0.025quant`
beta2.q1 <- m.beta$summary.random$idx2$`0.025quant`
beta3.q1 <- m.beta$summary.random$idx3$`0.025quant`
beta0.q1 <- m.beta$summary.random$idx0$`0.025quant`
c(mean(beta1.q1), mean(beta2.q1), mean(beta3.q1),  mean(beta0.q1))


beta1.q3 <- m.beta$summary.random$idx1$`0.975quant`
beta2.q3 <- m.beta$summary.random$idx2$`0.975quant`
beta3.q3 <- m.beta$summary.random$idx3$`0.975quant`
beta0.q3 <- m.beta$summary.random$idx0$`0.975quant`
c(mean(beta1.q3), mean(beta2.q3), mean(beta3.q3), mean(beta0.q3))


