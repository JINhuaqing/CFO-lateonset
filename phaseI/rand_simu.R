setwd("/home/r5user5/MyResearch/CFO-lateonset")
library(magrittr)
library(parallel)
source("utilities.R")
source("./phaseI/crm_utils.R")
source("./phaseI/boin_utils.R")
source("CFO_tox_utils.R")


target <- 0.30
#mus <- c(0.37, 0.50, 0.66, 0.86)
mus <- c(0.26, 0.41, 0.56, 0.74)
ncohort <- 12
cohortsize <- 3
ndose <- 7

add.args <- list(alp.prior=target, bet.prior=1-target)
diff <- 0.05
diffs <- c(0.05, 0.07, 0.10, 0.15)


for (diff in diffs){
mu <- mus[which(diffs==diff)]

# target = 0.2, Dose Level 6
# 0.05: mu = 0.35
# 0.07: mu = 0.49
# 0.10: mu = 0.65
# 0.15: mu = 0.85
# target = 0.2, Dose Level 7
# 0.05: mu = 0.37
# 0.07: mu = 0.50
# 0.10: mu = 0.66
# 0.15: mu = 0.86
# target = 0.3, Dose Level 7
# 0.05: mu = 0.26
# 0.07: mu = 0.41
# 0.10: mu = 0.56
# 0.15: mu = 0.74

run.fn <- function(k){
    set.seed(seeds[k])
    print(k)
    p.true.all <- gen.rand.doses(ndose, target, mu1=mu, mu2=mu)
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


nsimu <- 5000
seeds <- 1:nsimu
file.name <- paste0("./phaseI/results/", "Simu_", nsimu, "_phi_", 100*target,  "_random_", 100*diff,  ".RData")
print(file.name)
results <- mclapply(1:nsimu, run.fn, mc.cores=75)
post.process.random(results)
save(results, file=file.name)
}
