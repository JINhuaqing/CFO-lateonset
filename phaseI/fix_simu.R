setwd("/home/r5user5/MyResearch/CFO-lateonset")
library(magrittr)
library(parallel)

source("utilities.R")
source("./phaseI/crm_utils.R")
source("./phaseI/boin_utils.R")
source("CFO_tox_utils.R")
source("./fixScs3.R")

init.level <- 1
nsimu <- 5000

add.args <- list(alp.prior=target, bet.prior=1-target)

for (idx in 1:8){
p.true <- p.trues[[idx]]
tmtd <- MTD.level(target, p.true)


run.fn <- function(i){
    set.seed(seeds[i]) #10
    print(i)
    cfo.res <- CFO.simu.fn(target, p.true, ncohort=ncohort, cohortsize=cohortsize,
                                init.level=init.level,  add.args=add.args)
    crm.res <- crm.simu.fn(target=target, p.true=p.true, 
                              init.level=init.level, cohortsize=cohortsize, ncohort=ncohort)
    boin.res <- boin.simu.fn(target=target, p.true=p.true, ncohort=ncohort, 
                                init.level=init.level, cohortsize=cohortsize)
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

seeds <- 1:nsimu

results <- mclapply(1:nsimu, run.fn, mc.cores=80)
file.name <- paste0("./phaseI/results/", "Simu_", nsimu, "_fix3_",  idx, ".RData")
save(results, file=file.name)

crm.ress <- lapply(1:nsimu, function(i)results[[i]]$crm)
boin.ress <- lapply(1:nsimu, function(i)results[[i]]$boin)
cfo.ress <- lapply(1:nsimu, function(i)results[[i]]$cfo)
sum.all <- list(
                CFO = phase1.post.fn(cfo.ress),
                BOIN = phase1.post.fn(boin.ress),
                CRM = phase1.post.fn(crm.ress)
                )
print(tmtd)
print(p.true)
phase.I.pretty.tb(sum.all)
}
