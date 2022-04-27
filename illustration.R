rm(list=ls())
setwd("C:/Users/JINHU/Documents/ProjectCode/CFO-lateonset")
source("CFO_tox_utils.R")

### 1. Example for CFO design
phi <- 0.3
cns <- c(3, 6, 3)
cys <- c(0, 1, 1)
alp.prior <- phi
bet.prior <- 1-phi

gam1 <- optim.gamma.fn(cns[1], cns[2], phi, "L", alp.prior, bet.prior)$gamma;gam1
gam2 <- optim.gamma.fn(cns[2], cns[3], phi, "R", alp.prior, bet.prior)$gamma;gam2
OR.v1 <- OR.values(phi, cys[1], cns[1], cys[2], cns[2], alp.prior, bet.prior, type="L");OR.v1
OR.v2 <- OR.values(phi, cys[2], cns[2], cys[3], cns[3], alp.prior, bet.prior, type="R");OR.v2


### 2. Imputation example

tau <- 90
t1 <- 10
t2 <- 80
phi <- 0.3
tilde.p <- (1+phi/2)/(3+1)

impute.fn <- function(t, tau, tilde.p){
    num <- tilde.p * (1-t/tau) 
    den <- tilde.p * (1-t/tau)  + (1-tilde.p)
    num/den
}

y1.hat <- impute.fn(t1, tau, tilde.p);y1.hat
y2.hat <- impute.fn(t2, tau, tilde.p);y2.hat