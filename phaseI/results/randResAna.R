setwd("C:/Users/JINHU/Documents/ProjectCode/CFO-lateonset")
source("utilities.R")

nsimu = 5000
target <- 0.2
diff <- 0.05
res.ls <- list()
flag <- 1
for (diff in c(0.05, 0.07, 0.1, 0.15)){
file.name <- paste0("./phaseI/results/", "Simu_", nsimu, "_phi_", 100*target,  "_random_", 100*diff,  ".RData")
load(file.name)
res.ls[[flag]] <-  post.process.random(results)
flag <- flag + 1
}


Ress <- list()
for (nam in names(res.ls[[1]])){
    
    cRes <- do.call(rbind, lapply(res.ls, function(x){x[[nam]]}))
    colnames(cRes) <- row.names(res.ls[[1]])
    Ress[[nam]] <- cRes
}

singlePlot <- function(cRes, main=""){
    plot(cRes[, 1]*100, type = "b", col=1, lwd=2, lty=1, pch=1, ylim=c(min(cRes)-0.02, max(cRes)+0.02)*100, 
         xlab="Prob diff around the target", main=main, 
         xaxt="n", ylab="Percentage (%)")
    lines(cRes[, 2]*100, type = "b", col=2, lwd=2, lty=2, pch=2)
    lines(cRes[, 3]*100, type = "b", col=3, lwd=2, lty=3, pch=3)
    axis(1, at=1:4, labels=c("0.05", "0.07", "0.10", "0.15"))
    #legend("topleft", toupper(colnames(cRes)), col=1:3, lty=1:3, pch=1:3, lwd=2)
}


fig.name <- paste0("./phaseI/results/", "Simu_", nsimu, "_phi_", 100*target,  "_random",  ".png")
png(filename=fig.name, unit="in", height=10, width=8, res=300)
par(mfrow=c(3, 2), oma = c(2,1,1,1))
singlePlot(Ress[[1]], main="MTD Selection")
singlePlot(Ress[[2]], main="MTD Allocation")
singlePlot(Ress[[3]], main="Overdose Selection")
singlePlot(Ress[[4]], main="Overdose Allocation")
singlePlot(Ress[[5]], main="Risk of High toxicity")
singlePlot(Ress[[6]], main="Average DLT rate")
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend("bottom", toupper(colnames(cRes)), col=1:3, lty=1:3, pch=1:3, lwd=2, 
       xpd=TRUE, horiz = TRUE, cex = 1, seg.len=1, bty = 'n')
par(mfrow=c(1, 1))
dev.off()