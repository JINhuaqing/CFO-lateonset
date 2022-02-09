setwd("../")
#setwd("C:/Users/Dell/Documents/ProjectCode/phaseI/Rcode")
library(magrittr)
library(parallel)

source("utilities.R")
source("./phaseI/crm_utils.R")
source("./phaseI/boin_utils.R")
source("ORM_utils.R")


target <- 0.33
ncohort <- 10
cohortsize <- 3
init.level <- 1
nsimu <- 5000

add.args <- list(alp.prior=target, bet.prior=1-target)
p.trues <- list()
p.trues[[1]] <- c(0.33, 0.45, 0.58, 0.70, 0.80) 
p.trues[[2]] <- c(0.18, 0.33, 0.52, 0.60, 0.70) 
p.trues[[3]] <- c(0.12, 0.20, 0.33, 0.40, 0.50) 
p.trues[[4]] <- c(0.01, 0.02, 0.03, 0.33, 0.50) 
p.trues[[5]] <- c(0.00, 0.00, 0.05, 0.10, 0.33) 

idx <- 4
p.true <- p.trues[[idx]]
tmtd <- MTD.level(target, p.true)


run.fn <- function(i){
    set.seed(seeds[i]) #10
    print(i)
    orm.res <- ORM.simu.fn(target, p.true, ncohort=ncohort, cohortsize=cohortsize,
                                init.level=init.level,  add.args=add.args)
    crm.res <- crm.simu.fn(target=target, p.true=p.true, 
                              init.level=init.level, cohortsize=cohortsize, ncohort=ncohort)
    boin.res <- boin.simu.fn(target=target, p.true=p.true, ncohort=ncohort, 
                                init.level=init.level, cohortsize=cohortsize)
    ress <- list(
                 orm=orm.res,
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

results <- mclapply(1:nsimu, run.fn, mc.cores=75)
file.name <- paste0("../results/JRSSC-R/", "Hard_MTDSimu_", nsimu, "fix_",  idx, ".RData")
save(results, file=file.name)

crm.ress <- lapply(1:nsimu, function(i)results[[i]]$crm)
boin.ress <- lapply(1:nsimu, function(i)results[[i]]$boin)
orm.ress <- lapply(1:nsimu, function(i)results[[i]]$orm)
sum.all <- list(
                BOIN = phase1.post.fn(boin.ress),
                ORM = phase1.post.fn(orm.ress),
                CRM = phase1.post.fn(crm.ress)
                )
print(tmtd)
phase.I.pretty.tb(sum.all)

