setwd("/home/r5user5/MyResearch/CFO-lateonset")
library(magrittr)
library(parallel)
source("utilities.R")
source("Lateonset_utils.R")


target <- 0.20
mus <- c(0.37, 0.50, 0.66, 0.86)
#mus <- c(0.26, 0.41, 0.56, 0.74)
ncohort <- 12
cohortsize <- 3
ndose <- 7
tau <- 3
accrual <- 6
tite.dist <- 2
accrual.dist <- 1
init.dose <- 1
add.args <- list(alp.prior=target, bet.prior=1-target, CV=0.95, suspend=F, crmCI.CV=0.80, impute.method=2)
add.args.frac <- list(alp.prior=target, bet.prior=1-target, CV=0.95, suspend=F, crmCI.CV=0.80, impute.method=1)

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
                 fcfo=fcfo.res,
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
file.name <- paste0("./phaseI-late/results/", "Simu_", nsimu, "_phi_", 100*target,  "_random_", 100*diff,  ".RData")
print(file.name)
results <- mclapply(1:nsimu, run.fn, mc.cores=75)
post.process.random(results)
save(results, file=file.name)
}
