setwd("/home/r5user5/MyResearch/CFO-lateonset")
library(magrittr)
library(parallel)
source("utilities.R")
source("./phaseI/crm_utils.R")
source("./phaseI/boin_utils.R")
source("CFO_tox_utils.R")


target <- 0.20
ncohort <- 12
cohortsize <- 3

add.args <- list(alp.prior=target, bet.prior=1-target)

# target = 0.2
# 0.05: mu = 0.35
# 0.07: mu = 0.49
# 0.10: mu = 0.65
# 0.15: mu = 0.85

mu <- 0.21
run.fn <- function(k){
    set.seed(seeds[k])
    print(k)
    p.true.all <- gen.rand.doses(6, target, mu1=mu, mu2=mu)
    p.true <- p.true.all$p.true
    tmtd <- p.true.all$mtd.level

    cfo.res <- CFO.simu.fn(target, p.true, ncohort=ncohort, cohortsize=cohortsize, add.args=add.args)
    crm.res <- crm.simu.fn(target=target, p.true=p.true, cohortsize=cohortsize, ncohort=ncohort)
    boin.res <- boin.simu.fn(target=target, p.true=p.true, ncohort=ncohort, cohortsize)

    ress <- list(
                 cfo=cfo.res,
                 crm = crm.res, 
                 boin = boin.res, 
                 paras=list(p.true=p.true, 
                         mtd=tmtd, 
                         add.args=add.args,
                         target=target,
                         ncohort=ncohort,
                         cohortsize=cohortsize)
                 )
    ress
}


nsimu <- 100
seeds <- 1:nsimu
file.name <- paste0("./phaseI/results/", "Simu_", nsimu, "_random_05", ".RData")
results <- mclapply(1:nsimu, run.fn, mc.cores=75)
post.process.random(results)
save(results, file=file.name)

