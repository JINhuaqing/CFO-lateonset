# this file is to do trials with suspending
# note I assume when suspending, patient accural is suspended.
setwd("/Users/hujin/ProjectCode/TITE-CFO")
source("./Lateonset_utils.R")

make.decision.CFO.fn <- function(phi, cys, cns, alp.prior, bet.prior, cover.doses, diag=FALSE){
    if (cover.doses[2] == 1){
        if (!diag){
            res <- 1
        }else{
            res <- list(
                gam2=NA,
                gam1=NA,
                OR.v2=NA,
                OR.v1=NA,
                cidx=1
            )
        }
    }else{
        if (is.na(cys[1]) & (cover.doses[3]==1)){
            if (!diag){
                res <- 2
            }else{
                res <- list(
                    gam2=NA,
                    gam1=NA,
                    OR.v2=NA,
                    OR.v1=NA,
                    cidx=2
                )
            }
        }else  if (is.na(cys[1]) & (!(cover.doses[3]==1))){
            gam2 <- optim.gamma.fn(cns[2], cns[3], phi, "R", alp.prior, bet.prior)$gamma 
            OR.v2 <- OR.values(phi, cys[2], cns[2], cys[3], cns[3], alp.prior, bet.prior, type="R")
            if (OR.v2>gam2){
                cidx <- 3
            }else{
                cidx <- 2
            }
            if (!diag){
                res <- cidx
            }else{
                res <- list(
                    gam2=gam2,
                    gam1=NA,
                    OR.v2=OR.v2,
                    OR.v1=NA,
                    cidx=cidx
                )
            }
        }else  if (is.na(cys[3]) | (cover.doses[3]==1)){
            gam1 <- optim.gamma.fn(cns[1], cns[2], phi, "L", alp.prior, bet.prior)$gamma 
            OR.v1 <- OR.values(phi, cys[1], cns[1], cys[2], cns[2], alp.prior, bet.prior, type="L")
            if (OR.v1>gam1){
                cidx <- 1
            }else{
                cidx <- 2
            }
            if (!diag){
                res <- cidx
            }else{
                res <- list(
                    gam2=NA,
                    gam1=gam1,
                    OR.v2=NA,
                    OR.v1=OR.v1,
                    cidx=cidx
                )
            }
            
        }else  if (!(is.na(cys[1]) | is.na(cys[3]) | cover.doses[3]==1)){
            gam1 <- optim.gamma.fn(cns[1], cns[2], phi, "L", alp.prior, bet.prior)$gamma 
            gam2 <- optim.gamma.fn(cns[2], cns[3], phi, "R", alp.prior, bet.prior)$gamma 
            OR.v1 <- OR.values(phi, cys[1], cns[1], cys[2], cns[2], alp.prior, bet.prior, type="L")
            OR.v2 <- OR.values(phi, cys[2], cns[2], cys[3], cns[3], alp.prior, bet.prior, type="R")
            v1 <- OR.v1 > gam1
            v2 <- OR.v2 > gam2
            if (v1 & !v2){
                cidx <- 1
            }else if (!v1 & v2){
                cidx <- 3
            }else{
                cidx <- 2
            }
            if (!diag){
                res <- cidx
            }else{
                res <- list(
                    gam2=gam2,
                    gam1=gam1,
                    OR.v2=OR.v2,
                    OR.v1=OR.v1,
                    cidx=cidx
                )
            }
        }
    }
    
    return(res)
}

# Propose the next dose level with CFO design
CFOlateonset.next.dose <- function(curDose, phi, tau, impute.method, 
                                   enter.times, dlt.times, current.t, doses, tover.doses, simu=FALSE, add.args=list()){
    #args:
    #   curDose: the current dose level
    #   phi: Target DLT rate
    #   tau: maximal observing win size
    #   impute.method: impute method: 1: fractional design, 2: TITE-BOIN way
    #   enter.times: enter.times of each subject
    #   dlt.times: enter.times of each subject, let 0 if no DLT, tau+1 if unobserved and obs win is not reached.
    #   current.t: The current time 
    #   doses: Dose level for each subject
    #   tover.doses: over dose index for each dose, 1 over-toxic, 0 safe
    #   simu: Whether simulation or not, if TRUE, return tover.dose also
    #return:
    #   dose: Recommend doso, 0 early stopping
    
    ndose <- length(tover.doses)
    ## Obtain effective results
    impute.res <-  TITEImpute2(enter.times, dlt.times, current.t, tau, doses, ndose, c(phi/2, 1-phi/2))
    y.raw <- impute.res$y.raw
    y.impute <- impute.res$y.impute
        
    cy <- sum(y.impute[doses==curDose])
    cn <- sum(doses==curDose)
    
    ### early stop
    add.args <- c(list(y=cy, n=cn), add.args)
    if (overdose.fn(phi, add.args)){
        tover.doses[curDose:ndose] <- 1
    }
    
    if (tover.doses[1] == 1){
        dose <- 0
        dose.chg.res <- NULL
    }else{
        ### the results for current 3 dose levels
        if (curDose==1){
            cys <- c(NA, cy, sum(y.impute[doses==2]))
            cns <- c(NA, cn, sum(doses==2))
            cover.doses <- c(NA, tover.doses[1:2])
        }else if (curDose==ndose){
            cys <- c(sum(y.impute[doses==(curDose-1)]), cy, NA)
            cns <- c(sum(doses==(curDose-1)), cn, NA)
            cover.doses <- c(tover.doses[(curDose-1):curDose], NA)
        }else {
            cys <- c(sum(y.impute[doses==(curDose-1)]), cy, sum(y.impute[doses==(curDose+1)]))
            cns <- c(sum(doses==(curDose-1)), cn, sum(doses==(curDose+1)))
            cover.doses <- tover.doses[(curDose-1):(curDose+1)]
            
        }
        
        dose.chg.res <- make.decision.CFO.fn(phi, cys, cns, add.args$alp.prior, add.args$bet.prior, cover.doses, diag=TRUE)
        dose.chg <-  dose.chg.res$cidx - 2
        dose <- dose.chg + curDose
        
    }
    
    if (simu){
        res <- list(dose=dose, tover.doses=tover.doses, cfo.res=dose.chg.res, 
                    y.impute=y.impute, y.raw=y.raw)
        return(res)
    }else{
        res <- list(dose=dose)
        return(res)
    }
}


Ex.Simu.Fn <- function(phi, ps, tau, cohortsize, ncohort,
                    accrual, tite.dist, accrual.dist,
                    init.dose=1, add.args=list()){
    #args:
    #   phi: Target DLT rate
    #   ps: True dlt rates for each dose level      
    #   tau: maximal observing win size
    #   cohortsize: cohort size
    #   ncohort: num of cohort
    #   accural: the accrual rate, i.e., the number of patients accrued in tau time 
    #   tite.dist: Time to event distibution, 1-uniform, 2-weibull, 3-log-log
    #   accrual.dist: accrual distribution  0: fixed, 1, uniform, 2: exponential,  
    #   design: the phase I design, 1: CFO, 2: TITE-CRM, 3:TITE-BOIN
    #   init.dose: Initial dose level, 1 by default
    #   add.args: additional parameters
    #       alp.prior, bet.prior: prior parameter for TITE
    #       CV: critical value of early stopping for CFO and TITE-BOIN
    #       crmCI.CV: coverage prob of CI for CRM
    #       ps.prior: prior skeleton for CRM
    #       impute.method: impute method: 1: fractional design, 2: TITE-BOIN way
    #       boinEps: The parameter set the boundary for BOIN, by default, boinEps = 0.4
    
    ndose <- length(ps)
    if (is.null(add.args$alp.prior)){
        add.args <- c(add.args, list(alp.prior=phi, bet.prior=1-phi))
    }
    
    if (is.null(add.args$CV)){
        add.args <- c(add.args, list(CV=0.95))
    }
    
    
    
    earlystop <- 0
    enter.times <- NULL # enter time of each subject
    dlt.times <- NULL # dlt time of each subject
    dlts <- NULL # dlt event for each subject
    doses <- NULL # dose level for each subject
    cfo.ress <- NULL
    tover.dosess <- NULL
    suspend.ts <- c()
    current.ts.b4suspend <- c()
    current.ts.aftersuspend <- c()
    y.imputes <- list()
    y.raws <- list()
    
    current.t<- 0
    curDose <- init.dose  #current dose level
    
    tover.doses <- rep(0, ndose)
    if (is.null(add.args$impute.method)){
         add.args$impute.method <- 2
    }

    
    for (i in 1:ncohort){
        curP <- ps[curDose]    
        
        if (accrual.dist==0){
            delta.times <- rep(0, cohortsize)
        }else if (accrual.dist == 1){
            delta.times <- cumsum(c(0, runif(cohortsize-1, 0,  2*tau/accrual)))
            
        }else if (accrual.dist == 2){
            delta.times <- cumsum(c(0, rexp(cohortsize-1, rate=accrual/tau)))
            
        }
        enter.times <- c(enter.times, current.t+delta.times)
        
        # obtain the results of the patients
        obscohort <- gen.tite(tite.dist, cohortsize, curP, alpha=0.5, tau=tau);
        dlt.times <- c(dlt.times, obscohort$t.tox);
        dlts <- c(dlts, obscohort$tox);
        doses <- c(doses, rep(curDose, cohortsize));
        
        
        
        # Move to next cohort 
        if (i != ncohort){
            if (accrual.dist==0){
                delta.time <- tau*cohortsize/accrual
            }else if (accrual.dist == 1){
                delta.time <- runif(1, 0, 2*tau/accrual)
                
            }else if (accrual.dist == 2){
                delta.time <- rexp(1, rate=accrual/tau)
                
            }
        }else{
            delta.time <- tau
        }
        current.t<- enter.times[length(enter.times)] + delta.time;
        current.ts.b4suspend <- c(current.ts.b4suspend, current.t)
        
        
        # CFO with suspend
        # we assume when suspending, the patient accrual is stopped.
        suspend.ts.tmp <- c()
        for (tmp.dose in c(curDose-1, curDose, curDose+1)){
            suspend.ts.tmp <- c(suspend.ts.tmp, suspDelta.Fn(dlt.times, enter.times, current.t, tmp.dose, doses))
        }
        suspend.t <- max(suspend.ts.tmp)
        current.t <-  suspend.t + current.t
        current.ts.aftersuspend <- c(current.ts.aftersuspend, current.t)
        suspend.ts <- c(suspend.ts, suspend.t)
        
        res <- CFOlateonset.next.dose(curDose, phi, tau, add.args$impute.method, enter.times, 
                                          dlt.times, current.t, doses, tover.doses, TRUE, add.args)
        tover.doses <- res$tover.doses
      
        
        cfo.ress[[i]] <- res$cfo.res
        tover.dosess <- rbind(tover.dosess, tover.doses)
        y.raws[[i]] <-  res$y.raw
        y.imputes[[i]] <-res$y.impute
        
        if (res$dose==0){
            earlystop <- 1
            break()
        }else{
            curDose <- res$dose
        }
        
    }
    
    tns <- NULL
    tys <- NULL
    assess.t <- enter.times + tau
    y.raw <- (dlt.times!=0)*1
    for (j in 1:ndose){
        tns <- c(tns, sum(doses==j))
        tys <- c(tys, sum(y.raw[doses==j]))
    }
    if (earlystop==0){
        MTD.res <- select.mtd(phi, tns, tys)
        MTD <- select.mtd(phi, tns, tys)$MTD
      
    }else{
        MTD <- 99
    }
    res <- list(MTD=MTD, 
                MTD.res=MTD.res,
                dose.ns=tns, DLT.ns=tys, ps=ps, target=phi, 
                doses=doses,
                total.time=current.t,
                cfo.ress=cfo.ress, 
                tover.dosess=tover.dosess,
                dlt.times=dlt.times,
                y.raws = y.raws,
                y.imputes = y.imputes,
                enter.times=enter.times, 
                suspend.ts =suspend.ts, 
                current.ts.b4suspend = current.ts.b4suspend, 
                current.ts.aftersuspend = current.ts.aftersuspend
                )
    res
}

target <- 0.3
ncohort <- 10
cohortsize <- 3
tau <- 3
accrual <- 6
tite.dist <- 2
accrual.dist <- 1
init.dose=1
add.args=list(alp.prior=target, bet.prior=1-target, CV=0.95, suspend=F, crmCI.CV=0.80, impute.method=2)
p.true <- c(0.1, 0.15, 0.3, 0.46, 0.6)

set.seed(0)
cfo.res <- Ex.Simu.Fn(target, p.true, tau, cohortsize, ncohort, 
                   accrual, tite.dist, accrual.dist, add.args=add.args)
cfo.res$doses
cfo.res$tover.dosess
cfo.res$MTD
cfo.res$y.imputes
cfo.res$y.raws
cfo.res$suspend.ts
cfo.res$dlt.times

all.res.time <- matrix(cfo.res$enter.times + cfo.res$dlt.times + 3*(cfo.res$dlt.times==0), ncol=3, byrow=T);all.res.time
DLTs = matrix(cfo.res$dlt.times, ncol=3, byrow=TRUE) !=0
enter.times <- matrix(cfo.res$enter.times, ncol=3, byrow = T)
trial.pro <- cbind(cfo.res$current.ts.b4suspend, 
      cfo.res$current.ts.aftersuspend,
      cfo.res$suspend.ts, 
      enter.times, 
      all.res.time, 
      DLTs, 
      cfo.res$doses[seq(1, ncohort*cohortsize, cohortsize)]
      )
colnames(trial.pro) <- c("Cur.t0", "Cur.t1", "Susp.t", 
                         "", "Ent.t", "", 
                         "", "Res.t", "", 
                        "", "DLT", "", "Doselevel")
trial.pro

cfo.res$enter.times[1:3]*30 
cfo.res$enter.times[4]*30
cfo.res$y.imputes[1]
cfo.res$cfo.ress[[1]]

cfo.res$enter.times[1:6]*30 + 90
cfo.res$enter.times[7]*30
cfo.res$y.imputes[2]
cfo.res$cfo.ress[[2]]

cfo.res$enter.times[1:9]*30 + 90
cfo.res$enter.times[10]*30
cfo.res$y.imputes[3]
cfo.res$cfo.ress[[3]]

cfo.res$enter.times[1:12]*30 + 90
cfo.res$enter.times[13]*30
cfo.res$y.imputes[4]
cfo.res$cfo.ress[[4]]
cfo.res$tover.dosess[4, ]

cfo.res$total.time * 30
cfo.res$MTD.res

# Plot the results
DLT.plot <- function(cfo.res, ptsIdx, delta=0.5){
    curDose <- cfo.res$doses[ptsIdx]
    enter.time <- cfo.res$enter.times[ptsIdx]*30
    t.dlt.time <- enter.time+(cfo.res$dlt.times[ptsIdx])*30
    lines(c(enter.time, enter.time, t.dlt.time), 
          c(curDose, curDose+delta, curDose+delta), 
          lty=2)
    #lines(c(enter.time, enter.time, enter.time, t.dlt.time), 
    #      c(curDose, curDose+delta, curDose+delta, curDose+delta), 
    #      lty=2)
    points(t.dlt.time, curDose+delta, pch=4)
    text(x=t.dlt.time, y=curDose+delta, labels=round(t.dlt.time, 0), pos=4)
}
suspend.plot <- function(cur.time1, cur.time2, cur.dose){
    lines(c(cur.time1, cur.time1), c(cur.dose, cur.dose+0.8), lty=4, lwd=2)
    lines(c(cur.time2, cur.time2), c(cur.dose, cur.dose+0.8), lty=4, lwd=2)
    text(x=(cur.time1+cur.time2)/2, y=cur.dose+0.4, labels="Suspend", pos=3)
    arrows(y0=cur.dose+0.4, y1=cur.dose+0.4, x0=cur.time1, x1=cur.time2, length=0.1, lwd=2)
    arrows(y1=cur.dose+0.4, y0=cur.dose+0.4, x1=cur.time1, x0=cur.time2, length=0.1, lwd=2)
}
jpeg("./PharStatR1/TrialExampleSusp.jpg", width=10, height=5, units="in", res=300)
y.out <- cfo.res$y.raws[[length(cfo.res$y.raws)]]
plot(cfo.res$enter.times[y.out==1]*30, cfo.res$doses[y.out==1], 
     ylim=c(0, 5), xlim=c(0, cfo.res$total.time*30), 
     ylab = "Dose level", xlab="Time (in days)", 
     pch=19, cex=1.2)
points(cfo.res$enter.times[y.out==0]*30, cfo.res$doses[y.out==0], 
     pch=1, cex=1.2)
dltIdxSeq <- which(y.out==1);dltIdxSeq
DLT.plot(cfo.res, 3, -0.5)
DLT.plot(cfo.res, 8, 0.75)
DLT.plot(cfo.res, 10, 0.5)
DLT.plot(cfo.res, 11, 0.25)
DLT.plot(cfo.res, 20, 0.75)
DLT.plot(cfo.res, 24, 0.5)
DLT.plot(cfo.res, 27, 0.25)

suspend.plot(cur.time1=cfo.res$enter.times[3] * 30,
             cur.time2=cfo.res$enter.times[4] * 30, 
             cur.dose=1.1)

suspend.plot(cur.time1=cfo.res$enter.times[6] * 30,
             cur.time2=cfo.res$enter.times[7] * 30, 
             cur.dose=2.1)

suspend.plot(cur.time1=cfo.res$enter.times[9] * 30,
             cur.time2=cfo.res$enter.times[10] * 30, 
             cur.dose=2.1)

legend("topright", c("No DLT", "DLT", "DLT time"), 
       pch=c(1, 19, 4), 
       lty=c(0, 0, 0))
dev.off()

