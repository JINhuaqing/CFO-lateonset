rm(list=ls())
setwd("/Users/hujin/ProjectCode/TITE-CFO/")
source("utilities.R")
source("fixScs3.R")

fix.type <- "fix3"

idx <- 1
p.trues[[idx]]
nsimu <- 5000
MTDs <- c()
sum.alls <- list()
for (idx in 1:8){
    file.name <- paste0("./PharStatR1/results/", "Simu_withCFO", 5000, "_", fix.type, "_",  idx, ".RData")
    load(file.name)
    
    crm.ress <- lapply(1:nsimu, function(i)results[[i]]$crm)
    boin.ress <- lapply(1:nsimu, function(i)results[[i]]$boin)
    cfo.ress <- lapply(1:nsimu, function(i)results[[i]]$cfo)
    ncfo.ress <- lapply(1:nsimu, function(i)results[[i]]$ncfo)
    sum.all <- list(
        CFO = phase1.post.fn(cfo.ress),
        nCFO = phase1.post.fn(ncfo.ress),
        BOIN = phase1.post.fn(boin.ress),
        CRM = phase1.post.fn(crm.ress)
    )
    MTDs <- c(MTDs, which.min(abs(p.trues[[idx]]-target)))
    sum.alls[[idx]] <- sum.all
    print(p.trues[[idx]])
    print(phase.I.pretty.tb(sum.all))
    
}



OneRes.fn <- function(sum.all, nams, MTD){
    if (missing(nams))
        nams <- names(sum.all)
    
    ndose <- length(sum.all[[1]]$Allocation)
    MTD.sels <- c()
    MTD.allos.num <- c()
    MTD.allos <- c()
    durations <- c()
    errStops <- c()
    overDose.sels <- c()
    overDose.allos.num <- c()
    overDose.allos <- c()
    DLT.per <- c()
    for (nam in nams){
        curRes <- sum.all[[nam]]
        MTD.sels <- c(MTD.sels, curRes$Selection[MTD])
        MTD.allos.num <- c(MTD.allos.num, curRes$Allocation[MTD])
        MTD.allos <- c(MTD.allos, 100*curRes$Allocation[MTD]/sum(curRes$Allocation))
        durations <- c(durations, curRes$tol.time)
        errStops <- c(errStops, curRes$errStop)
        if (MTD==ndose){
            overDose.sels <- c(overDose.sels, 0)
            overDose.allos.num <- c(overDose.allos.num, 0)
            overDose.allos <- c(overDose.allos, 0)
        }else{
            overDose.sels <- c(overDose.sels, sum(curRes$Selection[(MTD+1):ndose]))
            overDose.allos.num <- c(overDose.allos.num, sum(curRes$Allocation[(MTD+1):ndose]))
            overDose.allos <- c(overDose.allos, 100*sum(curRes$Allocation[(MTD+1):ndose])/sum(curRes$Allocation))
        }
        DLT.per <- c(DLT.per, curRes$DLTs)
    }
    
    res <- list(
        MTD.sels=MTD.sels, 
        MTD.allos.num=MTD.allos.num, 
        MTD.allos=MTD.allos,
        durations=durations, 
        errStops=errStops, 
        overDose.sels=overDose.sels, 
        overDose.allos.num=overDose.allos.num, 
        overDose.allos=overDose.allos,
        DLT.per=DLT.per, 
        nams=nams
    )
    res
}

allRes.Fn <- function(sum.alls, nams, MTDs){
    #sum.alls <- list(...)
    ress <- lapply(1:length(sum.alls), function(x)OneRes.fn(sum.alls[[x]], nams=nams, MTD=MTDs[x]))
    m.nams <- names(ress[[1]])
    all.res <- list()
    for (m.nam in m.nams){
        all.res[[m.nam]] <- do.call(rbind, lapply(ress, function(x)x[[m.nam]]))
    }
    all.res
}

plot.single.fix <- function(plot.res, m.nam, ylab="Percentage (%)", main=""){
    cur.plot.res <- plot.res[[m.nam]]
    nams <- plot.res[["nams"]][1, ]
    plot(cur.plot.res[, 1], type = "b", col=1, lwd=2, lty=1, pch=1, 
         ylim=c(min(cur.plot.res)-0.02, max(cur.plot.res)+0.02), 
         xlab="Scenarios", 
         xaxt="n", ylab=ylab, main=main)
    for (ix in 2:(dim(plot.res$nams)[2])){
        lines(cur.plot.res[, ix], type = "b", col=ix, lwd=ix, lty=ix, pch=ix)
    }
    #lines(cur.plot.res[, 3], type = "b", col=3, lwd=2, lty=3, pch=3)
    axis(1, at=1:dim(cur.plot.res)[1], labels=1:dim(cur.plot.res)[1])
}

nams <- c("nCFO", "CFO", "BOIN", "CRM")
plt.nams <- c("CFO", "TITE-CFO", "TITE-BOIN", "TITE-CRM")
plot.res <- allRes.Fn(sum.alls, nams=nams, MTDs=MTDs)
plot.res$durations[plot.res$durations==-1] = ncohort * 3 # !!!not if tau != 3


fig.name <- paste0("./PharStatR1/results/", "SimuLate_", nsimu, "_phi_", 100*target,  "_", fix.type,  ".png")
png(filename=fig.name, unit="in", height=7, width=8, res=300)
par(mfrow=c(2, 3), oma = c(2,1,1,1))
plot.single.fix(plot.res, m.nam="MTD.sels", main="(A) MTD Selection")
plot.single.fix(plot.res, m.nam="MTD.allos", main="(B) MTD Allocation")
plot.single.fix(plot.res, m.nam="overDose.sels", main="(C) Overdose Selection")
plot.single.fix(plot.res, m.nam="overDose.allos",  main="(D) Overdose Allocation")
plot.single.fix(plot.res, m.nam="DLT.per", main="(E) Average DLT Rate")
plot.single.fix(plot.res, m.nam="durations", main="(F) Average Trial Duration", ylab="Time (in months)")
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend("bottom", toupper(plt.nams), col=1:4, lty=1:4, pch=1:4, lwd=2, 
       xpd=TRUE, horiz = TRUE, cex = 1, seg.len=1)
par(mfrow=c(1, 1))
dev.off()
