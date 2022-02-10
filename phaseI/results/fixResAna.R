setwd("C:/Users/JINHU/Documents/ProjectCode/CFO-lateonset")
source("utilities.R")
source("fixScs3.R")

idx <- 1
nsimu <- 5000
for (idx in 1:8){
file.name <- paste0("./phaseI/results/", "Simu_", 5000, "_fix3_",  idx, ".RData")
load(file.name)

crm.ress <- lapply(1:nsimu, function(i)results[[i]]$crm)
boin.ress <- lapply(1:nsimu, function(i)results[[i]]$boin)
cfo.ress <- lapply(1:nsimu, function(i)results[[i]]$cfo)
sum.all <- list(
                CFO = phase1.post.fn(cfo.ress),
                BOIN = phase1.post.fn(boin.ress),
                CRM = phase1.post.fn(crm.ress)
                )
print(p.trues[[idx]])
#print(phase.I.pretty.tb(sum.all))
    
}
res <- phase.I.pretty.tb(sum.all)

latex.out.fn(res, prefix="test")
