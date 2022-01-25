library(magrittr)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(grid)

OBD.level <- function(phi, phiE, p.true, pE.true){
    if (p.true[1]>phi+0.1){
        OBD <- 99
        return(OBD)
    }
    MTD <- which.min(abs(phi - p.true))
    eff.idxs <- phiE > pE.true[1:MTD]
    if (sum(eff.idxs)==MTD){
        OBD <- 99
        return(OBD)
    }

    OBD <- which.max(pE.true[1:MTD])
    return(OBD)
}

MTD.level <- function(phi, p.true){
    if (p.true[1]>phi+0.1){
        MTD <- 99
        return(MTD)
    }
    MTD <- which.min(abs(phi - p.true))
    return(MTD)
}

phase.I.pretty.tb <- function(sum.all){
    # The errStops in fact is the none selection rate
    ndose <- length(sum.all[[1]]$Selection)
    tb <- NULL
    tox.nums <- c()
    errStops <- c()
    sub.nums <- c()
    DLT.rates <- c()
    for (res in sum.all){
        rw <- paste0(round(res$Selection, 1), "(", round(res$Allocation, 1), ")")
        tb <- rbind(tb, rw)
        tox.nums <- c(tox.nums, res$tol.toxs)
        sub.nums <- c(sub.nums, sum(res$tol.Subjs))
        errStops <- c(errStops, 100-sum(res$Selection))
        DLT.rates <- c(DLT.rates, res$DLTs)
    }
    
    tb.df <- as.data.frame(tb)
    rownames(tb.df) <- names(sum.all)
    names(tb.df) <- paste0("Level ", 1:ndose)
    
    tb.df["nTox"] <- round(tox.nums, 1)
    tb.df["nSub"] <- round(sub.nums, 1)
    tb.df["DLT.rate"] <- round(DLT.rates, 1)
    tb.df["NonSel.rate"] <- round(errStops, 1)
    tb.df
}

phase.I.II.pretty.tb <- function(sum.all){
    # The errStops in fact is the none selection rate
    ndose <- length(sum.all[[1]]$Selection)
    tb <- NULL
    tox.nums <- c()
    eff.nums <- c()
    errStops <- c()
    sub.nums <- c()
    for (res in sum.all){
        rw <- paste0(round(res$Selection, 1), "(", round(res$Allocation, 1), ")")
        tb <- rbind(tb, rw)
        tox.nums <- c(tox.nums, res$tol.toxs)
        eff.nums <- c(eff.nums, res$tol.effs)
        sub.nums <- c(sub.nums, sum(res$Allocation))
        errStops <- c(errStops, 100-sum(res$Selection))
    }
    
    tb.df <- as.data.frame(tb)
    rownames(tb.df) <- names(sum.all)
    names(tb.df) <- paste0("Level ", 1:ndose)
    
    tb.df["PerTox"] <- round(100*tox.nums/sub.nums, 1)
    tb.df["PerEff"] <- round(100*eff.nums/sub.nums, 1)
    tb.df["nSub"] <- round(sub.nums, 1)
    tb.df["NonSel.rate"] <- round(errStops, 1)
    tb.df
}

phase1.post.fn <- function(ress){
    numTrials <- length(ress)
    ndose <- length(ress[[1]]$dose.ns)
    
    Allo <- rep(0, ndose)
    Sel <- rep(0, ndose)
    toxs.cts <- rep(0, ndose)
    tol.Subjs <- 0
    nonErrStops <- 0
    for (res in ress){
        if (res$MTD != 99){
            nonErrStops <- nonErrStops + 1
            Sel[res$MTD] <- 1 + Sel[res$MTD]
        }
        Allo <- Allo + res$dose.ns
        toxs.cts <- res$DLT.ns + toxs.cts
        tol.Subjs <- tol.Subjs + sum(res$dose.ns)
    }
    
    sum.v <- list(Allocation=Allo, Selection=Sel*100,
    #sum.v <- list(Allocation=numTrials*100*Allo/tol.Subjs, Selection=Sel*100,
                  toxs.nums=toxs.cts,
                  DLTs = numTrials*100*sum(toxs.cts)/tol.Subjs,
                  tol.Subjs=tol.Subjs, errStop=100*(numTrials-nonErrStops),
                  tol.toxs=sum(toxs.cts))
    lapply(sum.v, function(i)i/numTrials)
}

phase12.post.fn <- function(ress){
    numTrials <- length(ress)
    ndose <- length(ress[[1]]$dose.ns)
    
    Allo <- rep(0, ndose)
    Sel <- rep(0, ndose)
    effs.cts <- rep(0, ndose)
    toxs.cts <- rep(0, ndose)
    tol.Subjs <- 0
    nonErrStops <- 0
    for (res in ress){
        if (res$OBD != 99){
            nonErrStops <- nonErrStops + 1
            Sel[res$OBD] <- 1 + Sel[res$OBD]
        }
        Allo <- Allo + res$dose.ns
        effs.cts <- res$eff.ns + effs.cts
        toxs.cts <- res$DLT.ns + toxs.cts
        tol.Subjs <- tol.Subjs + sum(res$dose.ns)
    }
    
    sum.v <- list(Allocation=Allo, Selection=Sel*100,
                  effs.nums=effs.cts, toxs.nums=toxs.cts,
                  tol.Subjs=tol.Subjs, errStop=100*(numTrials-nonErrStops),
                  tol.effs=sum(effs.cts), tol.toxs=sum(toxs.cts))
    lapply(sum.v, function(i)i/numTrials)
}

# return latex scr code for output
latex.out.fn <- function(res, prefix, tidx, suffix=NULL){
  MTDs <- res$MTDs.percent  
  MTDs <- format(round(MTDs*100, 1), nsmall=1)
  if (!missing(tidx)){
      MTDs[tidx] <- paste0("\\bf{", MTDs[tidx], "}")
  }
  av.dose <- res$av.dose
  tss <- res$t.dose
  tdlt <- res$t.DLT
  row1 <- paste("&", MTDs, collapse=" ") 
  row1 <- paste(prefix, "& MTD \\%", row1, "\\\\"); 
  
  row2 <- paste("&", format(round(av.dose, 1), nsmall=1), collapse=" ") 
  row2 <- paste("& ASN", row2, "&", format(round(tdlt, 1), nsmall=1), "&",
                format(round(tss, 1), nsmall=1), "\\\\")
  out <- paste(row1, "\n", row2, suffix)
  cat(out)
}

# Gnerate dose toxicity relations randomly
gen.rand.doses <- function(ndose, target, mu1=0.55, mu2=0.55, inv.fn=qnorm, fn=pnorm, MTD=NA){
    sigma0 <- 0.05
    sigma1 <- 0.35
    sigma2 <- 0.35
    raw.p.trues <- rep(0, ndose)
    if (is.na(MTD)){
        mtd.level <- sample.int(ndose, 1)
    }else{
        mtd.level <- MTD
    }
    
    raw.p.trues[mtd.level] <- rnorm(1, inv.fn(target), sigma0)
    
    if (mtd.level != 1){
        eps.minus <- rnorm(1, mu1, sigma1)
        if (raw.p.trues[mtd.level] > inv.fn(target)){
            dlt <- raw.p.trues[mtd.level] - inv.fn(2*target-fn(raw.p.trues[mtd.level]))
        }else{
            dlt <- 0 
        }
        raw.p.trues[mtd.level-1] <- raw.p.trues[mtd.level] - (dlt +  eps.minus**2)
    }
    
    if (mtd.level != ndose){
        eps.plus <- rnorm(1, mu2, sigma2)
        if (raw.p.trues[mtd.level] < inv.fn(target)){
            dlt <- inv.fn(2*target-fn(raw.p.trues[mtd.level])) - raw.p.trues[mtd.level]
        }else{
            dlt <- 0 
        }
        raw.p.trues[mtd.level+1] <- raw.p.trues[mtd.level] + (dlt +  eps.plus**2)
    }
    
    if ((mtd.level-2)>0){
        for (i in (mtd.level-2):1){
           eps.minus <- rnorm(1, mu1, sigma1)
           raw.p.trues[i]  <- raw.p.trues[i+1] - eps.minus**2
        }
    }
    if ((mtd.level+2)<=ndose){
        for (j in (mtd.level+2):ndose){
           eps.plus <- rnorm(1, mu2, sigma2)
           raw.p.trues[j] <- raw.p.trues[j-1] + eps.plus**2
        }
    }
    p.trues <- fn(raw.p.trues)
    list(p.trues=p.trues, mtd.level=mtd.level)
}



# compute the prob diff around the target level
prob.diff.fn <- function(res, target=0.3){
    p.trues <- res$p.trues
    mtd <- res$mtd.level
    ndose <- length(p.trues)
    
    diffs <- c()
    if (mtd!=1){
        diffs <- c(abs(p.trues[mtd-1]-target), diffs)
    }

    if (mtd!=ndose){
        diffs <- c(abs(p.trues[mtd+1]-target), diffs)
    }
    min(diffs)
}

# Compute the prob diff around the target level, separately
prob.diff.fn.sep <- function(res, target=0.3){
    p.trues <- res$p.trues
    mtd <- res$mtd.level
    ndose <- length(p.trues)
    
    diffsL <- NULL
    diffsU <- NULL
    if (mtd!=1){
        diffsL <- abs(p.trues[mtd-1]-target)
    }

    if (mtd!=ndose){
        diffsU <- abs(p.trues[mtd+1]-target)
    }
    list(L=diffsL, U=diffsU) 
}



#mu <- 0.69
#phi <- 0.33
#ress <- lapply(1:10000, function(i)gen.rand.doses(5, phi, mu1=mu, mu2=mu))
#sapply(ress, prob.diff.fn, target=phi) %>% mean
#tmp <- lapply(ress, prob.diff.fn.sep, target=0.3)
#sapply(tmp, function(i)i$L) %>% unlist %>% mean
#sapply(tmp, function(i)i$U) %>% unlist %>% mean




post.process.onemethod <- function(res, paras){
    tmtd <- paras$mtd
    target <- paras$target
    ndose <- length(res$p.true)
    rv <- rep(0, 7)
    rv[1] <- sum(res$MTD==tmtd)
    rv[2] <- res$dose.ns[tmtd]
    if (res$MTD == 99){
        rv[3] <- 0
    }else{
        rv[3] <- sum(res$MTD>tmtd)
    }
    if (tmtd==ndose){
        rv[4] <- 0
    }else{
        rv[4] <- sum(res$dose.ns[(tmtd+1):ndose])
    }
    rv[5] <- sum((sum(res$DLT.ns)/sum(res$dose.ns))>target)
    rv[6] <- sum(res$DLT.ns)
    rv[7] <- sum(res$dose.ns)
    names(rv) <- c(
        "MTD Sel", "MTD Allo", "Over Sel", 
        "Over Allo", "Risk of HT", "No DLT",
        "No Subject")
    rv
    
}
post.process.single <- function(result){
    nMethod <- length(result)-1
    paras <- result$paras
    tmtd <- paras$mtd
    methods <- names(result)[1:nMethod]
    res.v <- list()
    for (name in methods){
        res.v[[name]] <- post.process.onemethod(result[[name]], paras)
    }
    do.call(rbind, res.v)
}

# Function to handle the results under random scenario
post.process.random <- function(results){
    nsimu <- length(results)
    res.all <- 0
    for (result in results){
        res.all <- post.process.single(result) + res.all
    }
    res.all.df <- data.frame(res.all)
    final.res <- transmute(res.all.df, 
                           MTD.Sel=MTD.Sel/nsimu,
                           MTD.Allo=MTD.Allo/No.Subject,
                           Over.Sel=Over.Sel/nsimu,
                           Over.Allo=Over.Allo/No.Subject,
                           Risk.of.HT=Risk.of.HT/nsimu,
                           PerDLT=No.DLT/No.Subject
    )
    rownames(final.res) <- rownames(res.all.df)
    final.res
}


# function to plot the results under random scenarios
res.plot.fn <- function(res, filename, M.names, is.save=TRUE, angle=45){
  if (missing(filename)){
    is.save <- FALSE
  }
  myBarplot <- function(data, labels_vec, rot_angle, main, ylim=NULL) {
    nM <- dim(res)[1]
    plt <- barplot(data, col=2:(nM+1), xaxt="n", main=main, ylim=ylim)
    text(plt, par("usr")[3], labels = labels_vec, srt = rot_angle, adj = c(1.1,1.1), 
         xpd = TRUE, cex=0.8, font=1) 
  }
  if (missing(M.names)){
      M.names <- c("Butterfly-BB", "Butterfly-CRM", "Butterfly-BB-CRM", "CRM", "BOIN")
  }
  if (is.save){
    png(filename, width=360*2, height=240*2)
  }
  par(mfrow=c(2, 3))
  myBarplot(res$MTD.Sel, main="MTD selection %", labels_vec=M.names, rot_angle=angle)
  myBarplot(res$MTD.Allo, main="MTD Allocation %", labels_vec=M.names, rot_angle=angle)
  myBarplot(res$Over.Sel, main="Overdose selection %", labels_vec=M.names, rot_angle=angle)
  myBarplot(res$Over.Allo, main="Overdose Allocation %", labels_vec=M.names, rot_angle=angle)
  myBarplot(res$Risk.of.HT, main="Risk of high toxicity %", labels_vec=M.names, rot_angle=angle)
  myBarplot(res$PerDLT, main="Avergae DLT rate %", labels_vec=M.names, rot_angle=angle)
  par(mfrow=c(1, 1))
  if (is.save){
    dev.off()
  }
}

res.ggplot.fn <- function(res, filename, main, M.names, is.save=TRUE, angle=NULL){
    if (missing(filename)){
        is.save <- FALSE
    }
    myggBarplot <- function(data, angle, main, theylim) {
         g <- ggplot(mapping=aes(x=M.names, y=data, fill=M.names))
         g <- g + geom_bar(stat="identity") + guides(fill=FALSE) + ggtitle(main)
         g <- g + theme(axis.title = element_blank())
         if (!is.null(angle)){
             g <- g + theme(axis.text.x = element_text(angle=angle, vjust=1, hjust=1))
         }
         if (!missing(theylim)){
             g <- g + ylim(theylim)
         }
             
         g
    }
    if (missing(M.names)){
        M.names <- c("Butterfly-BB", "Butterfly-CRM", "Butterfly-BB-CRM", "CRM", "BOIN")
    }
    g1 <- myggBarplot(res$MTD.Sel, main="MTD selection %", angle=angle)
    g2 <- myggBarplot(res$MTD.Allo, main="MTD Allocation %", angle=angle)
    g3 <- myggBarplot(res$Over.Sel, main="Overdose selection %", angle=angle)
    g4 <- myggBarplot(res$Over.Allo, main="Overdose Allocation %", angle=angle)
    g5 <- myggBarplot(res$Risk.of.HT, main="Risk of high toxicity %", angle=angle)
    g6 <- myggBarplot(res$PerDLT, main="Avergae DLT rate %", angle=angle)
    gg <- grid.arrange(g1, g2, g3, g4, g5, g6, nrow=3, top=textGrob(main, gp=gpar(fontsize=20, fontface=2)))
    if (is.save){
        ggsave(filename, gg, width=8, height=12, unit="in")
    }
}


# Plot the randomly generated dose levels
gen.rand.doses.plot <- function(cases, phi){
    ncase <- length(cases)
    case <- cases[[1]]
    ndose <- length(case$p.trues)
    plot(1:ndose, case$p.trues, type="b", ylim=c(0, 1), ylab="Prob of toxicity", xlab="Dose levels", xaxt="n")
    points(case$mtd.level, case$p.trues[case$mtd.level], pch=19, col=1)
    axis(1, at=1:ndose, labels=1:ndose)
    for (i in 2:ncase){
        ccase <- cases[[i]]
        lines(1:ndose, ccase$p.trues, type="b", ylim=c(0, 1), lty=i)
        points(ccase$mtd.level, ccase$p.trues[ccase$mtd.level], pch=19, col=1)
        i <- i + 1
    }
    abline(h=phi, col=2, lwd=2)
    legend("topleft", legend=c("Target prob", "MTD level"), col=c(2, 1), pch=c(NA, 19), 
           lty=c(1, NA), lwd=c(2, 1))
    
}

# posterior probability of pj >= phi given data
post.prob.fn <- function(phi, y, n, alp.prior=0.1, bet.prior=0.1){
    alp <- alp.prior + y 
    bet <- bet.prior + n - y
    1 - pbeta(phi, alp, bet)
}


#plot phase I/II scenarios
phaseI_II.scenario.plot.fn <- function(p.true, pE.true, main="Plateau", OBD=NULL){
    ndose <- length(p.true)
    lvls <- 1:ndose
    cols <- c("black", "red")
    lwds <- c(2, 2)
    ltys <- c(1, 2)
    pchs <- c(15, 16)
    plot(1:ndose, p.true, ylab = "Probability", xlab="Dose level", ylim=c(0, 0.8), main=main, xaxt="n", 
         type="b", 
         lwd=lwds[1], col=cols[1], lty=ltys[1], pch=pchs[1])
    if (is.null(OBD)){
        axis(1, at=lvls, lvls)
    }else{
        axis(1, at=OBD, paste0(OBD, "*"), col.axis="red", font.axis=2)
        axis(1, at=lvls[-OBD], lvls[-OBD])
    }
    lines(1:ndose, pE.true, type="b", lwd=lwds[2], col=cols[2], lty=ltys[2], pch=pchs[2])
    #    legend("topleft", c("DLT rate", "Efficacy rate"), col=cols, lty=ltys, lwd=lwds, pch=pchs)
}


# Generate phase I/II scenarios
gen.rand.doses.plateau <- function(ndose, phi, psi, psi.U, mu1, mu2){
  # args:
  # ndose: number of dose levels
  # phi: DLT target rate
  # psi: the minimal acceptable efficacy rate 
  # psi.U: the upper bound of the efficacy rate
  # mu1, mu2: parameters for phase I trial
    k.OBD <- sample.int(ndose, 1)
    k.MTD <- sample(k.OBD:ndose, 1)
    
    ps <- gen.rand.doses(ndose, phi, mu1, mu2, MTD=k.MTD)
    
    q.OBD <- runif(1, psi, psi.U)
    if (k.OBD == 1){
        qs <- rep(q.OBD, ndose)
    }else if (k.OBD == ndose){
        qs.l <- sort(runif(ndose-1, 0, q.OBD))
        qs <- c(qs.l, q.OBD)
    }else {
        qs.l <- sort(runif(k.OBD-1, 0, q.OBD))
        qs.u <- rep(q.OBD, ndose-k.OBD)
        qs <- c(qs.l, q.OBD, qs.u)
    }
    
    res <- list(
        qs=qs,
        ps=ps$p.trues, 
        k.MTD=k.MTD, 
        k.OBD=k.OBD
    )
    res
}

# with xxx.U, generate sc such that has maximal utils in OBD for STEIN
gen.rand.doses.plateau.U<- function(ndose, phi, psi, psi.U, mu1, mu2){
    flag <- FALSE
    while (!flag) {
        sc <- gen.rand.doses.plateau(ndose, phi, psi, psi.U, mu1, mu2)
        u<-sc$qs-0.33*sc$ps-1.09*(sc$ps>phi) 	
        flag <- sc$k.OBD == which.max(u)
    }
    sc
}

gen.rand.doses.umbrella.U <- function(ndose, phi, psi, psi.U, mu1, mu2){
    flag <- FALSE
    while (!flag) {
        sc <- gen.rand.doses.umbrella(ndose, phi, psi, psi.U, mu1, mu2)
        u<-sc$qs-0.33*sc$ps-1.09*(sc$ps>phi) 	
        flag <- sc$k.OBD == which.max(u)
    }
    sc
}
gen.rand.doses.umbrella <- function(ndose, phi, psi, psi.U, mu1, mu2){
  # args:
  # ndose: number of dose levels
  # phi: DLT target rate
  # psi: the minimal acceptable efficacy rate 
  # psi.U: the upper bound of the efficacy rate
  # mu1, mu2: parameters for phase I trial
  
    k.OBD <- sample.int(ndose, 1)
    k.MTD <- sample(k.OBD:ndose, 1)
    
    ps <- gen.rand.doses(ndose, phi, mu1, mu2, MTD=k.MTD)
    
    q.OBD <- runif(1, psi, psi.U)
    if (k.OBD == 1){
        qs.u <- sort(runif(ndose-1, 0, q.OBD), decreasing=TRUE)
        qs <- c(q.OBD, qs.u)
    }else if (k.OBD == ndose){
        qs.l <- sort(runif(ndose-1, 0, q.OBD))
        qs <- c(qs.l, q.OBD)
    }else {
        qs.l <- sort(runif(k.OBD-1, 0, q.OBD))
        qs.u <- sort(runif(ndose-k.OBD, 0, q.OBD), decreasing=TRUE)
        qs <- c(qs.l, q.OBD, qs.u)
    }
    
    res <- list(
        qs=qs,
        ps=ps$p.trues, 
        k.MTD=k.MTD, 
        k.OBD=k.OBD
    )
    res
}


## Below, the three functions are used for processing results under random scenarios for phase I/II trials
phase12.post.process.single <- function(res){

    ndose <- length(res$p.true)
    
    #us<- res$pE.true-0.33*res$p.true-1.09*(res$p.true>res$target) 	
    #U.OBD <- which.max(us)

    rv <- rep(0, 9)
    rv[1] <- sum(res$OBD==res$k.OBD) # OBD sel
    rv[2] <- res$dose.ns[res$k.OBD] # OBD allo
    #rv[10] <- sum(res$OBD==U.OBD) # U.OBD sel
    #rv[11] <- res$dose.ns[U.OBD] # U.OBD allo
    if (res$OBD== 99){ # over sel
        rv[3] <- 0
    }else{
        rv[3] <- sum(res$OBD>res$k.MTD)
    }
    if (res$k.MTD==ndose){ # over allo
        rv[4] <- 0
    }else{
        rv[4] <- sum(res$dose.ns[(res$k.MTD+1):ndose])
    }
    # good bio dose
    GBDs <- 1:res$k.MTD
    GBDs[res$pE.true[1:res$k.MTD] > res$min.eff]
    rv[5] <- sum(res$OBD %in% GBDs) #GBD sel
    rv[6] <- sum(res$dose.ns[GBDs]) # GBD allo
    rv[7] <- sum(res$DLT.ns)
    rv[8] <- sum(res$eff.ns)
    rv[9] <- sum(res$dose.ns)
    names(rv) <- c(
            "OBD Sel", "OBD Allo", "Over Sel", 
            "Over Allo", "GBD Sel", "GBD Allo", 
            "No DLT", "No Eff", "No Subject")
        rv
}
    
phase12.post.process.random <- function(results){
    nsimu <- length(results)
    res.all <- 0
    for (result in results){
        res.all <- phase12.post.process.single(result) + res.all
    }
    res.all.df <- data.frame(t(res.all))
    final.res <- transmute(res.all.df, 
                           OBD.Sel=OBD.Sel/nsimu,
                           #UOBD.Sel=UOBD.Sel/nsimu,
                           OBD.Allo=OBD.Allo/No.Subject,
                           #UOBD.Allo=UOBD.Allo/No.Subject,
                           Over.Sel=Over.Sel/nsimu,
                           Over.Allo=Over.Allo/No.Subject,
                           GBD.Sel = GBD.Sel/nsimu, 
                           GBD.Allo = GBD.Allo/No.Subject, 
                           PerDLT=No.DLT/No.Subject,
                           PerEff=No.Eff/No.Subject,
                           AvgSub=No.Subject/nsimu
    )
    final.res
}

phase12.post.process.random.all <- function(..., names){
  #args:
  # ...: results for each method
  # names: name for each method
    all.DF <- c()
    for (ress in list(...)){
        all.DF <- rbind(all.DF, phase12.post.process.random(ress))
        
    }
    if (!missing(names))
        rownames(all.DF) <- names
    
    all.DF
}
