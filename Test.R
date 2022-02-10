#setwd("C:/Users/JINHU/Documents/ProjectCode/CFO-lateonset")
setwd("/home/r5user5/MyResearch/CFO-lateonset")
source("Lateonset_utils.R")
source("utilities.R")


phi <- 0.3
ps <- c(0.2, 0.3, 0.4, 0.5, 0.6)
tau <- 3
cohortsize <- 3
ncohort <- 12
accrual <- 6
tite.dist <- 2
accrual.dist <- 0
init.dose=1
add.args=list(CV=0.9, suspend=TRUE)
design <- 1
res <- Simu.Fn(phi, ps, tau, cohortsize, ncohort, 
                     accrual, tite.dist, accrual.dist, design=design, add.args=add.args)

ress <- list()
for (i in 1:100){
res <- Simu.Fn(phi, ps, tau, cohortsize, ncohort, 
                     accrual, tite.dist, accrual.dist, design=3, add.args=add.args)
ress[[i]] <- res
}
length(ress)
phase1.post.fn(ress)

