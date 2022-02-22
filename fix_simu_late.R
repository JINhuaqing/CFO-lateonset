#setwd("C:/Users/JINHU/Documents/ProjectCode/CFO-lateonset")
setwd("/home/r5user5/MyResearch/CFO-lateonset")
library(magrittr)
library(parallel)
source("Lateonset_utils.R")
source("utilities.R")
source("./fixScs2.R")


nsimu <- 5000
seeds <- 1:nsimu
tau <- 3
accrual <- 6
tite.dist <- 2
accrual.dist <- 1
init.dose=1
add.args=list(alp.prior=target, bet.prior=1-target, CV=0.95, suspend=F, crmCI.CV=0.80, impute.method=2)
add.args.frac <- list(alp.prior=target, bet.prior=1-target, CV=0.95, suspend=F, crmCI.CV=0.80, impute.method=1)
#   design: the phase I design, 1: CFO, 2: TITE-CRM, 3:TITE-BOIN

idx <- 1
for (idx in 1:8){
p.true <- p.trues[[idx]]

run.fn <- function(i){
    set.seed(seeds[i]) #10
    print(i)
    cfo.res <- Simu.Fn(target, p.true, tau, cohortsize, ncohort, 
                     accrual, tite.dist, accrual.dist, design=1, add.args=add.args)
    fcfo.res <- Simu.Fn(target, p.true, tau, cohortsize, ncohort, 
                     accrual, tite.dist, accrual.dist, design=1, add.args=add.args.frac)
    crm.res <- Simu.Fn(target, p.true, tau, cohortsize, ncohort, 
                     accrual, tite.dist, accrual.dist, design=2, add.args=add.args)
    boin.res <- Simu.Fn(target, p.true, tau, cohortsize, ncohort, 
                     accrual, tite.dist, accrual.dist, design=3, add.args=add.args)
    ress <- list(
                 cfo=cfo.res,
                 crm = crm.res, 
                 fcfo=fcfo.res,
                 boin = boin.res, 
                 paras=list(p.true=p.true, 
                             add.args=add.args,
                             target=target,
                             ncohort=ncohort,
                             cohortsize=cohortsize)
        )
    ress
    
}


results <- mclapply(1:nsimu, run.fn, mc.cores=75)
file.name <- paste0("./phaseI-late/results/", "Simu", nsimu, "_fix2_",  idx, ".RData")
save(results, file=file.name)

crm.ress <- lapply(1:nsimu, function(i)results[[i]]$crm)
boin.ress <- lapply(1:nsimu, function(i)results[[i]]$boin)
cfo.ress <- lapply(1:nsimu, function(i)results[[i]]$cfo)
fcfo.ress <- lapply(1:nsimu, function(i)results[[i]]$fcfo)
sum.all <- list(
                CFO = phase1.post.fn(cfo.ress),
                fCFO = phase1.post.fn(fcfo.ress),
                BOIN = phase1.post.fn(boin.ress),
                CRM = phase1.post.fn(crm.ress)
                )
print(p.true)
print(phase.I.pretty.tb(sum.all))
}

