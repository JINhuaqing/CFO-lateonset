# dlt.time: time to event, if 0, no event
# delta.time, time interval of patient arrive
library(survival) # to obtain K-M curve
library(BOIN) #mtd.select
library(dfcrm) #TITE-CRM, getprior
source("CFO_tox_utils.R")

# Below 3 functions are to impute missing y 
#------------------------------------------------------------------------------------------
fracImpute <- function(enter.times, dlt.times, current.time, tau){ 
    
    #args:
    # enter.times: The enter times of the patients, a vector 
    # dlt.times: The DLT times of the patients, if no DLT, 0, a vector
    # current.time: Current time point: a value
    # tau: Observing window size
    
    #return:
    # ym: Imputed y for the patient with no DLT and follow-up time < tau
  
    assesstime = enter.times+tau;	
    dlt.times[dlt.times==0]= tau+1;
    yo = (dlt.times<=tau)*(assesstime<=current.time)+(dlt.times<=(current.time-enter.times))*(current.time<assesstime);		
    No.impute <- FALSE
    if (sum(yo)==0)	{
        No.impute <- TRUE
        ym <- yo

    #stop("fraction design takes effect when at least one DLT has been observed")
    }
    if (sum(yo)!=0){			
        
        otime = yo*dlt.times+(1-yo)*((current.time-enter.times)*(current.time<assesstime)+tau*(assesstime<=current.time))			
        kmfit = survfit(Surv(otime,yo)~1)	
        #plot(kmfit,xlab="Time (days)",ylab="Survival probability",
        #     main="Kaplan-Meier estimate",cex.lab=1.3,cex.main=1.3,cex.sub=1.3)	
        ym = yo
        
        for (i in 1:length(yo)){
            if (current.time<assesstime[i] & yo[i]==0){
                ym[i]=(kmfit$surv[tail(which(kmfit$time<=(current.time-assesstime[i]+tau+0.00001)),n=1)]- kmfit$surv[tail(which(kmfit$time<=tau),n=1)])/
                    kmfit$surv[tail(which(kmfit$time<=(current.time-assesstime[i]+tau+0.00001)),n=1)]

            }
        }
        
    }
    obsIdxs <- current.time >= assesstime
    obsIdxs[yo==1] <- TRUE
  
    res <- list(y.impute=ym, y.raw=yo, obsIdxs=obsIdxs, No.impute=No.impute)
    res
}

#TITE impuate
TITEImpute.one <- function(followup.times, tau, y, n, prior.paras){
    #args:
    #   followup.times: The follow-up times of the pending patients at the dose level
    #   tau: Observing window size
    #   y: Num of Observed DLT at the dose level 
    #   n: Num of patients with observed results at the dose level
    #   prior.paras: a vector of 2, prior when estimating ptilde
    
    #return: 
    #  ym: imputed y 
    
    p.tilde <- (y+prior.paras[1])/(n+sum(prior.paras))
    ym <- p.tilde * (1-followup.times/tau) /(1-p.tilde)
    #ym[ym >1] <- 1 # trunc the value
    ym
}

TITEImpute <- function(enter.times, dlt.times, current.time, tau, dose.levels, ndose, prior.paras){
    #args:
    # enter.times: The enter times of the patients, a vector 
    # dlt.times: The DLT times of the patients, if no DLT before tau, 0, a vector
    # current.time: Current time point: a value
    # tau: Observing window size
    # dose.levels: dose level for each subject
    # ndose: num of dose levels
    # prior.paras: a vector of 2, prior when estimating ptilde
    
    #return:
    # ym: Imputed y for the patient with no DLT and follow-up time < tau
    
    assesstime = enter.times + tau;	
    followup.times <- current.time - enter.times
    dlt.times[dlt.times==0]= tau+1;
    yo <- (dlt.times<=tau)*(assesstime<=current.time)+(dlt.times<=(current.time-enter.times))*(current.time<assesstime);		
    obsIdxs <- current.time >= assesstime
    obsIdxs[yo==1] <- TRUE
    ym <- yo
    for (i in 1:ndose){
        doseIdxs <- dose.levels == i
        if (sum(1-obsIdxs[doseIdxs]!=0)){
            y <- sum(yo[doseIdxs])
            n <- sum(doseIdxs[obsIdxs==1])
            kpidxs <- doseIdxs & (obsIdxs!=1)
            ym.part <- TITEImpute.one(followup.times[kpidxs], tau, y, n, prior.paras)
            ym[kpidxs] <- ym.part
        }
    }
    res <- list(y.impute=ym, y.raw=yo, obsIdxs=obsIdxs)
    res
    
}


TITEImpute2.one <- function(followup.times, tau, y, n, prior.paras){
    #args:
    #   followup.times: The follow-up times of the pending patients at the dose level
    #   tau: Observing window size
    #   y: Num of Observed DLT at the dose level 
    #   n: Num of patients with observed results at the dose level
    #   prior.paras: a vector of 2, prior when estimating ptilde
    
    #return: 
    #  ym: imputed y 
    
    p.tilde <- (y+prior.paras[1])/(n+sum(prior.paras))
    #ym <- p.tilde * (1-followup.times/tau)
    ym <- p.tilde * (1-followup.times/tau) /((1-p.tilde)+p.tilde * (1-followup.times/tau))
#    ym <- p.tilde * (1-followup.times/tau) /(1-p.tilde)
#    ym[ym >1] <- 1 # trunc the value
    ym
}

TITEImpute2 <- function(enter.times, dlt.times, current.time, tau, dose.levels, ndose, prior.paras){
    #args:
    # enter.times: The enter times of the patients, a vector 
    # dlt.times: The DLT times of the patients, if no DLT before tau, 0, a vector
    # current.time: Current time point: a value
    # tau: Observing window size
    # dose.levels: dose level for each subject
    # ndose: num of dose levels
    # prior.paras: a vector of 2, prior when estimating ptilde
    
    #return:
    # ym: Imputed y for the patient with no DLT and follow-up time < tau
    
    assesstime = enter.times + tau;	
    followup.times <- current.time - enter.times
    dlt.times[dlt.times==0]= tau+1;
    yo <- (dlt.times<=tau)*(assesstime<=current.time)+(dlt.times<=(current.time-enter.times))*(current.time<assesstime);		
    obsIdxs <- current.time >= assesstime
    obsIdxs[yo==1] <- TRUE
    ym <- yo
    for (i in 1:ndose){
        doseIdxs <- dose.levels == i
        if (sum(1-obsIdxs[doseIdxs]!=0)){
            y <- sum(yo[doseIdxs])
            n <- sum(doseIdxs[obsIdxs==1])
            kpidxs <- doseIdxs & (obsIdxs!=1)
            ym.part <- TITEImpute2.one(followup.times[kpidxs], tau, y, n, prior.paras)
            ym[kpidxs] <- ym.part
        }
    }
    res <- list(y.impute=ym, y.raw=yo, obsIdxs=obsIdxs)
    res
    
}
#------------------------------------------------------------------------------------------

# The function is to obtain the DLT results for each subject
#------------------------------------------------------------------------------------------
gen.tite<-function(dist=1, n, pi, tau=1, alpha=0.5){
    #args:
    #   dist: TITE distribution, 1-uniform, 2-weibull, 3-log-log
    #   n: Num of subjects to generate
    #   pi: Target DLT rate, pi=Pr(T<=tau)
    #   tau: Maximal window size
    #   alpha: Parameter for generete time
    #Return:
    #   if no DLT, tox.t=0
    ############ subroutines ############
    weib<-function(n, pi, pihalft)
    {
        ## solve parameters for Weibull given pi=1-S(T) and phalft=1-S(T/2)
        alpha = log(log(1-pi)/log(1-pihalft))/log(2);
        lambda = -log(1-pi)/(tau^alpha);
        t = (-log(runif(n))/lambda)^(1/alpha);
        return(t);
    }
    
    llogit<-function(n, pi, pihalft)
    {
        ## solve parameters for log-logistic given pi=1-S(T) and phalft=1-S(T/2)
        alpha = log((1/(1-pi)-1)/(1/(1-pihalft)-1))/log(2);
        lambda = (1/(1-pi)-1)/(tau^alpha);
        t = ((1/runif(n)-1)/lambda)^(1/alpha);
        return(t);
    }
    ############ end of subroutines ############
    
    
    tox = rep(0, n);
    t.tox = rep(0, n);
    
    #### uniform
    if(dist==1) {  # 50% event in (0, 1/2T)
        tox = rbinom(n, 1, pi);
        ntox.st = sum(tox);
        t.tox[tox==1]=runif(ntox.st, 0, tau);
        t.tox[tox==0]=0;
    }
    #### Weibull
    if(dist==2)
    {
        pihalft = alpha*pi;  # alpha*100% event in (0, 1/2T)
        t.tox = weib(n, pi, pihalft);
        tox[t.tox<=tau]=1;
        ntox.st = sum(tox);
        t.tox[tox==0]=0;
    }
    #### log-logistic
    if(dist==3)
    {
        pihalft = alpha*pi;  # alpha*100% event in (0, 1/2T)
        t.tox = llogit(n, pi, pihalft);
        tox[t.tox<=tau]=1;
        ntox.st = sum(tox);
        t.tox[tox==0]=0;
    }
    return(list(tox=tox, t.tox=t.tox, ntox.st=ntox.st));
}


#------------------------------------------------------------------------------------------

# Below functions are for CFO design
#------------------------------------------------------------------------------------------
# # posterior probability of pj >= phi given data
# post.prob.fn <- function(phi, y, n, alp.prior=0.1, bet.prior=0.1){
#     alp <- alp.prior + y 
#     bet <- bet.prior + n - y
#     1 - pbeta(phi, alp, bet)
# }



# Odds.samples <- function(y1, n1, y2, n2, alp.prior, bet.prior){
#     alp1 <- alp.prior + y1
#     alp2 <- alp.prior + y2
#     bet1 <- alp.prior + n1 - y1
#     bet2 <- alp.prior + n2 - y2
#     sps1 <- c()
#     sps2 <- c()
#     while (length(sps1)<10000){
#         sp1 <- rbeta(1, alp1, bet1)
#         sp2 <- rbeta(1, alp2, bet2)
#         if (sp1 <= sp2){
#             sps1 <- c(sp1, sps1)
#             sps2 <- c(sp2, sps2)
#         }
#     }
#     
#     list(sps1=sps1, sps2=sps2)
# }
# 
# prob.int <- function(phi, y1, n1, y2, n2, alp.prior, bet.prior){
#     alp1 <- alp.prior + y1
#     alp2 <- alp.prior + y2
#     bet1 <- alp.prior + n1 - y1
#     bet2 <- alp.prior + n2 - y2
#     fn.min <- function(x){
#         dbeta(x, alp1, bet1)*(1-pbeta(x, alp2, bet2))
#     }
#     fn.max <- function(x){
#         pbeta(x, alp1, bet1)*dbeta(x, alp2, bet2)
#     }
#     const.min <- integrate(fn.min, lower=0, upper=1)$value
#     const.max <- integrate(fn.max, lower=0, upper=1)$value
#     p1 <- integrate(fn.min, lower=0, upper=phi)$value/const.min
#     p2 <- integrate(fn.max, lower=0, upper=phi)$value/const.max
#     
#     list(p1=p1, p2=p2)
# }
# 
# 
# OR.values <- function(phi, y1, n1, y2, n2, alp.prior, bet.prior, type){
#     ps <- prob.int(phi, y1, n1, y2, n2, alp.prior, bet.prior)
#     if (type=="L"){
#         pC <- 1 - ps$p2
#         pL <- 1 - ps$p1
#         oddsC <- pC/(1-pC)
#         oddsL <- pL/(1-pL)
#         OR <- oddsC*oddsL
#         
#     }else if (type=="R"){
#         pC <- 1 - ps$p1
#         pR <- 1 - ps$p2
#         oddsC <- pC/(1-pC)
#         oddsR <- pR/(1-pR)
#         OR <- (1/oddsC)/oddsR
#     }
#     return(OR)
# }
# 
# All.OR.table <- function(phi, n1, n2, type, alp.prior, bet.prior){
#    ret.mat <- matrix(rep(0, (n1+1)*(n2+1)), nrow=n1+1)
#    for (y1cur in 0:n1){
#        for (y2cur in 0:n2){
#            ret.mat[y1cur+1, y2cur+1] <- OR.values(phi, y1cur, n1, y2cur, n2, alp.prior, bet.prior, type)
#        }
#    }
#    ret.mat
# }
# 
# # compute the marginal prob when lower < phiL/phiC/phiR < upper
# # i.e., Pr(Y=y|lower<phi<upper)
# margin.phi <- function(y, n, lower, upper){
#     C <- 1/(upper-lower)
#     fn <- function(phi) {
#        dbinom(y, n, phi)*C
#     }
#     integrate(fn, lower=lower, upper=upper)$value
# }
# 
# # Obtain the table of marginal distribution of (y1, y2) 
# # after intergrate out (phi1, phi2)
# # under H0 and H1
# # H0: phi1=phi, phi < phi2 < 2phi
# # H1: phi2=phi, 0   < phi1 < phi
# margin.ys.table <- function(n1, n2, phi, hyperthesis){
#     if (hyperthesis=="H0"){
#         p.y1s <- dbinom(0:n1, n1, phi)
#         p.y2s <- sapply(0:n2, margin.phi, n=n2, lower=phi, upper=2*phi)
#     }else if (hyperthesis=="H1"){
#         p.y1s <- sapply(0:n1, margin.phi, n=n1, lower=0, upper=phi)
#         p.y2s <- dbinom(0:n2, n2, phi)
#     }
#     p.y1s.mat <- matrix(rep(p.y1s, n2+1), nrow=n1+1)
#     p.y2s.mat <- matrix(rep(p.y2s, n1+1), nrow=n1+1, byrow=TRUE)
#     margin.ys <- p.y1s.mat * p.y2s.mat
#     margin.ys
# }
# 
# # Obtain the optimal gamma for the hypothesis test
# optim.gamma.fn <- function(n1, n2, phi, type, alp.prior, bet.prior){
#     OR.table <- All.OR.table(phi, n1, n2, type, alp.prior, bet.prior) 
#     ys.table.H0 <- margin.ys.table(n1, n2, phi, "H0")
#     ys.table.H1 <- margin.ys.table(n1, n2, phi, "H1")
#     
#     argidx <- order(OR.table)
#     sort.OR.table <- OR.table[argidx]
#     sort.ys.table.H0 <- ys.table.H0[argidx]
#     sort.ys.table.H1 <- ys.table.H1[argidx]
#     n.tol <- length(sort.OR.table)
#     
#     if (type=="L"){
#         errs <- rep(0, n.tol-1)
#         for (i in 1:(n.tol-1)){
#             err1 <- sum(sort.ys.table.H0[1:i])
#             err2 <- sum(sort.ys.table.H1[(i+1):n.tol])
#             err <- err1 + err2
#             errs[i] <- err
#         }
#         min.err <- min(errs)
#         if (min.err > 1){
#             gam <- 0
#             min.err <- 1
#         }else {
#             minidx <- which.min(errs)
#             gam <- sort.OR.table[minidx]
#         }
#     }else if (type=='R'){
#         errs <- rep(0, n.tol-1)
#         for (i in 1:(n.tol-1)){
#             err1 <- sum(sort.ys.table.H1[1:i])
#             err2 <- sum(sort.ys.table.H0[(i+1):n.tol])
#             err <- err1 + err2
#             errs[i] <- err
#         }
#         min.err <- min(errs)
#         if (min.err > 1){
#             gam <- 0
#             min.err <- 1
#         }else {
#             minidx <- which.min(errs)
#             gam <- sort.OR.table[minidx]
#         }
#         
#     }
#     list(gamma=gam, min.err=min.err)
# }
# 
# make.decision.CFO.fn <- function(phi, cys, cns, alp.prior, bet.prior, cover.doses, diag=FALSE){
#     if (cover.doses[2] == 1){
#         return(1)
#     }else{
#         if (is.na(cys[1]) & (cover.doses[3]==1)){
#             return(2)
#         }else  if (is.na(cys[1]) & (!(cover.doses[3]==1))){
#            gam2 <- optim.gamma.fn(cns[2], cns[3], phi, "R", alp.prior, bet.prior)$gamma 
#            OR.v2 <- OR.values(phi, cys[2], cns[2], cys[3], cns[3], alp.prior, bet.prior, type="R")
#            if (OR.v2>gam2){
#                return(3)
#            }else{
#                return(2)
#            }
#         }else  if (is.na(cys[3]) | (cover.doses[3]==1)){
#            gam1 <- optim.gamma.fn(cns[1], cns[2], phi, "L", alp.prior, bet.prior)$gamma 
#            OR.v1 <- OR.values(phi, cys[1], cns[1], cys[2], cns[2], alp.prior, bet.prior, type="L")
#            if (OR.v1>gam1){
#                return(1)
#            }else{
#                return(2)
#            }
#             
#         }else  if (!(is.na(cys[1]) | is.na(cys[3]) | cover.doses[3]==1)){
#            gam1 <- optim.gamma.fn(cns[1], cns[2], phi, "L", alp.prior, bet.prior)$gamma 
#            gam2 <- optim.gamma.fn(cns[2], cns[3], phi, "R", alp.prior, bet.prior)$gamma 
#            OR.v1 <- OR.values(phi, cys[1], cns[1], cys[2], cns[2], alp.prior, bet.prior, type="L")
#            OR.v2 <- OR.values(phi, cys[2], cns[2], cys[3], cns[3], alp.prior, bet.prior, type="R")
#            v1 <- OR.v1 > gam1
#            v2 <- OR.v2 > gam2
#            if (v1 & !v2){
#                return(1)
#            }else if (!v1 & v2){
#                return(3)
#            }else{
#                return(2)
#            }
#         }
#     }
# }
# 
# overdose.fn <- function(phi, add.args=list()){
#     CV <- add.args$CV
#     if (is.null(CV)){
#         CV <- 0.95
#     }
#     y <- add.args$y
#     n <- add.args$n
#     alp.prior <- add.args$alp.prior
#     bet.prior <- add.args$bet.prior
#     pp <- post.prob.fn(phi, y, n, alp.prior, bet.prior)
#     if ((pp >= CV) & (add.args$n>=3)){
#         return(TRUE)
#     }else{
#         return(FALSE)
#     }
# }

# Propose the next dose level with CFO design
CFOlateonset.next.dose <- function(curDose, phi, tau, impute.method, 
                          enter.times, dlt.times, current.t, doses, tover.doses, simu=FALSE, add.args=list()){
    #args:
    #   curDose: the current dose level
    #   phi: Target DLT rate
    #   tau: maximal observing win size
    #   impute.method: impute method: 1: fractional design, 2: TITE-BOIN way
    #   enter.times: enter.times of each subject
    #   dlt.times: dlt.times of each subject, let 0 if no DLT
    #   current.t: The current time 
    #   doses: Dose level for each subject
    #   tover.doses: over dose index for each dose, 1 over-toxic, 0 safe
    #   simu: Whether simulation or not, if TRUE, return tover.dose also
    #return:
    #   dose: Recommend doso, 0 early stopping
  
    ndose <- length(tover.doses)
    ## Obtain effective results
    if (impute.method == 1){
        impute.res <- fracImpute(enter.times, dlt.times, current.t, tau)
        y.raw <- impute.res$y.raw
        y.impute <- impute.res$y.impute
        if (impute.res$No.impute){
            cn <- sum(doses==curDose)
            # should be cn <- sum(doses==curDose & impute.res$obsIdxs)
            cy <- sum(y.raw[doses==curDose])
        }
        else{
            cn <- sum(doses==curDose)
            cy <- sum(y.impute[doses==curDose])
        }
    }else if(impute.method==2){
        impute.res <-  TITEImpute2(enter.times, dlt.times, current.t, tau, doses, ndose, c(phi/2, 1-phi/2))
        y.raw <- impute.res$y.raw
        y.impute <- impute.res$y.impute
        
        cy <- sum(y.impute[doses==curDose])
        cn <- sum(doses==curDose)
    }
    
    ### early stop
    add.args <- c(list(y=cy, n=cn), add.args)
    if (overdose.fn(phi, add.args)){
        tover.doses[curDose:ndose] <- 1
    }
    
    if (tover.doses[1] == 1){
        dose <- 0
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
        
        dose.chg <- make.decision.CFO.fn(phi, cys, cns, add.args$alp.prior, add.args$bet.prior, cover.doses) - 2
        dose <- dose.chg + curDose
      
    }
    
    if (simu){
        res <- list(dose=dose, tover.doses=tover.doses)
        return(res)
    }else{
        res <- list(dose=dose)
        return(res)
    }
}

#------------------------------------------------------------------------------------------

# Functions below to calculate the delta time to suspend the trial
#------------------------------------------------------------------------------------------
unObsRate.Fn <- function(res.times, enter.times, current.t, curDose, doses){
    followup.ts <- current.t - enter.times
    unObsIdx <- (followup.ts < tau) & (res.times > followup.ts)
    return(mean(unObsIdx[doses==curDose]))
}

suspDelta.Fn <- function(dlt.times, enter.times, current.t, curDose, doses){
    if (sum(doses==curDose)==0){
      # to avoid the case where there is no doses==curDose
      return(0)
    }
    followup.ts <- current.t - enter.times
    res.times <- dlt.times
    res.times[dlt.times==0] <- tau
    unObsIdx <- (followup.ts < tau) & (res.times > followup.ts)
    delta.ts.sorted <- c(0, sort(unique(res.times[unObsIdx] - followup.ts[unObsIdx])))
    for (delta.t in delta.ts.sorted){
        # add 0.0001 for rounding error in R
        unObsRate <- unObsRate.Fn(res.times, enter.times, current.t+delta.t+0.0001, curDose, doses)
        if (unObsRate <= 0.5){
            return(delta.t)
        }
    }
}

#------------------------------------------------------------------------------------------


# main function to RUN simulation
#------------------------------------------------------------------------------------------
Simu.Fn <- function(phi, ps, tau, cohortsize, ncohort,
                          accrual, tite.dist, accrual.dist, design=1,
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
    current.t<- 0
    curDose <- init.dose  #current dose level
    
    if (design==1){
        tover.doses <- rep(0, ndose)
        if (is.null(add.args$impute.method)){
               add.args$impute.method <- 2
        }
    }else if (design==2){
        if (is.null(add.args$crmCI.CV)){
               add.args$crmCI.CV <- 0.95
        }
        if (is.null(add.args$ps.prior)){
               add.args$ps.prior <- getprior(0.05, phi, ceiling(ndose/2), ndose)
        }
      
    }else if (design==3){
        tover.doses <- rep(0, ndose)
        if (is.null(add.args$boinEps)){
               add.args$boinEps <- 0.4
        }
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
        
        
        if (design == 1){
            # CFO
            res <- CFOlateonset.next.dose(curDose, phi, tau, add.args$impute.method, enter.times, 
                             dlt.times, current.t, doses, tover.doses, TRUE, add.args)
            tover.doses <- res$tover.doses
        }else if(design==2){
            #TITE-CRM
            dlt.times1 <- dlt.times
            dlt.times1[dlt.times==0] <- tau + 1
            
            assess.ts <- enter.times + tau;	
            followup.ts <- current.t - enter.times
            y.raw <- (dlt.times1<=tau)*(assess.ts<=current.t)+(dlt.times1<=followup.ts)*(current.t<assess.ts);		
            crm.res <- titecrm(add.args$ps.prior, phi, y.raw, doses, followup=followup.ts, obswin=tau, conf.level=add.args$crmCI.CV)
            
            res <- list()
            if (crm.res$ptoxL[1]>phi){
                res$dose <- 0
            }else{
                if (crm.res$mtd > curDose){
                    res$dose <- curDose + 1
                }else if (crm.res$mtd < curDose){
                    res$dose <- curDose - 1
                }else{
                    res$dose <- curDose 
                }
            }
        }else if (design==3){
           # TITE-BOIN
           if (add.args$suspend){
             # Note that here we assume if suspending, the patient accrual is also suspended.
             # In fact, for a more logic way, we should have two enter.times, 
             # One is the patient-come time, one is the time of patients taking the drug
             # When suspending, patients still come at a given rate.
             # 
              suspend.t <- suspDelta.Fn(dlt.times, enter.times, current.t, curDose, doses)
              current.t <-  suspend.t + current.t
           }
           impute.res <-  TITEImpute(enter.times, dlt.times, current.t, tau, doses, ndose, c(phi/2, 1-phi/2))
    
           y.raw <- impute.res$y.raw
           y.impute <- impute.res$y.impute
           
           cy1 <- sum(y.raw[doses==curDose])
           cn1 <- sum((doses==curDose) & impute.res$obsIdxs)
           ### early stop
           add.args <- c(list(y=cy1, n=cn1), add.args)
           if (overdose.fn(phi, add.args)){
               tover.doses[curDose:ndose] <- 1
           }
    
           cy <- sum(y.impute[doses==curDose])
           cn <- sum(doses==curDose)

           phat <- cy/cn
           phi1 <- phi * (1-add.args$boinEps)
           phi2 <- phi * (1+add.args$boinEps)
           
           lowb <- log((1-phi1)/(1-phi))/log(phi*(1-phi1)/phi1/(1-phi))
           upb <- log((1-phi)/(1-phi2))/log(phi2*(1-phi)/phi/(1-phi2))
           
           res <- list()
           
           doseLoc.addsets <- list()
           for (j in 1:ndose){
              if (j==1)
                cSet <- c(0, 1)
              else if (j==ndose)
                cSet <- c(0, -1)
              else
                cSet <- c(-1, 0, 1)
              doseLoc.addsets[[j]] <- cSet
           }
           
           if (tover.doses[1] == 1){
                 res$dose <- 0
           }else{
             if (tover.doses[curDose]==1){
                 over.addset <- -1
             }else if ((curDose==ndose)|(tover.doses[curDose+1]==1)){
                 over.addset <- c(-1, 0)
             }else{
                 over.addset <- c(-1, 0, 1)
             }
              
             doseLoc.addset <- doseLoc.addsets[[curDose]]
             add.set <- intersect(over.addset, doseLoc.addset)
             if (phat<=lowb){
                 chgidx <- 1
             }else if (phat >= upb){
                 chgidx <- -1
             }else{
                 chgidx <- 0
             }
             
             chgidx <- add.set[which.min(abs(add.set-chgidx))]
             res$dose <- curDose + chgidx
           }
        }else{
            res <- list()
            res$dose <- curDose
        }

        
        
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
        if (design==1 | design==3 | design==4){
            MTD <- select.mtd(phi, tns, tys)$MTD
        }else if (design==2){
            MTD <- res$dose
        }
    }else{
        MTD <- 99
    }
    res <- list(MTD=MTD, dose.ns=tns, DLT.ns=tys, ps=ps, target=phi, total.time=current.t)
    res
}

