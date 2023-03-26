# This file is to check the robustness of selection of the prior (\alpha, \beta) in imputation (Mar 26, 2023)
setwd("/home/MyResearch/CFOs/")
library(magrittr)
library(parallel)
source("utilities.R")
source("Lateonset_utils.R")
source("CFO_tox_utils.R")


target <- 0.30
nsimu <- 5000
#mus <- c(0.37, 0.50, 0.66, 0.86)
mus <- c(0.26, 0.41, 0.56, 0.74)
ncohort <- 12
cohortsize <- 3
ndose <- 7
tau <- 3
accrual <- 6
tite.dist <- 2
accrual.dist <- 1
init.dose <- 1
add.args1 <- list(alp.prior=target, bet.prior=1-target, CV=0.95, suspend=F, crmCI.CV=0.80, impute.method=2)
add.args2 <- list(alp.prior=1, bet.prior=1, CV=0.95, suspend=F, crmCI.CV=0.80, impute.method=2)
add.args3 <- list(alp.prior=0.5, bet.prior=0.5, CV=0.95, suspend=F, crmCI.CV=0.80, impute.method=2)
add.args4 <- list(alp.prior=0.5*target, bet.prior=1-0.5*target, CV=0.95, suspend=F, crmCI.CV=0.80, impute.method=2)

diff <- 0.05
diffs <- c(0.05, 0.07, 0.10, 0.15)


#for (diff in diffs){
mu <- mus[which(diffs==diff)]

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

    cfo.res1 <- Simu.Fn(target, p.true, tau, cohortsize, ncohort, 
                     accrual, tite.dist, accrual.dist, design=1, add.args=add.args1)
    cfo.res2 <- Simu.Fn(target, p.true, tau, cohortsize, ncohort, 
                     accrual, tite.dist, accrual.dist, design=1, add.args=add.args2)
    cfo.res3 <- Simu.Fn(target, p.true, tau, cohortsize, ncohort, 
                     accrual, tite.dist, accrual.dist, design=1, add.args=add.args3)
    cfo.res4 <- Simu.Fn(target, p.true, tau, cohortsize, ncohort, 
                     accrual, tite.dist, accrual.dist, design=1, add.args=add.args4)

    ress <- list(
                 cfo1=cfo.res1,
                 cfo2=cfo.res2,
                 cfo3=cfo.res3,
                 cfo4=cfo.res4,
                 paras=list(p.true=p.true, 
                         mtd=tmtd, 
                         target=target,
                         ncohort=ncohort,
                         cohortsize=cohortsize)
                 )
    ress
}


seeds <- 1:nsimu
file.name <- paste0("./PharStatR1/results/", "Simu_robustalpbet", nsimu, "_phi_", 100*target,  "_random_", 100*diff,  ".RData")
print(file.name)
results <- mclapply(1:nsimu, run.fn, mc.cores=20)
print(post.process.random(results))
save(results, file=file.name)
#}
