# this files contains some utils functions for CFO + (lin and yun, 2020, biostat, LY)
# LY is a impute method 
# It is for revision in PharStat
# LY method is like this
# Suppose at time kappa
# for dose level j, the follow-up time is ui, and tau is the observing win time. 
# dlti=1, if DLT outcome is observed otherwise dlti= 0
# xi=1, if DLT, xi=0 if not DLT in (0, tau)
# tyj = sum dlti xi: num of DLT at time kp
# mj = sum dlti (1-xi): num of non-DLT and ui > tau, at time kp
# tmj = mj + sum (1-dlti) (ui/tau): 
# effective data they use: tyj (DLT) and tmj (no DLT), tnj = tyj + tmj (total)
# LY+CFO, use (tyj, tmj) for CFO

###!!!, I can not make it (Mar 13, 2023)
### In fact, if I use LY's way to do impute, the tyj+tmy can be decimal.
### I can not calculate CV with a decimal value for total sample size.

# RUN simulation for CFO + (lin and yuan, 2020, biostat)
#------------------------------------------------------------------------------------------
LYCFO.Fn <- function(phi, ps, tau, cohortsize, ncohort,
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
  
  tover.doses <- rep(0, ndose)
  
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
    
    
   # CFO next step
   res <- CFOlateonset.next.dose(curDose, phi, tau, add.args$impute.method, enter.times, 
                                    dlt.times, current.t, doses, tover.doses, TRUE, add.args)
   tover.doses <- res$tover.doses
    
    
    
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
    MTD <- select.mtd(phi, tns, tys)$MTD
  }else{
    MTD <- 99
  }
  res <- list(MTD=MTD, dose.ns=tns, DLT.ns=tys, ps=ps, target=phi, total.time=current.t)
  res
}

# Propose the next dose level with CFO design + (Lin and Yuan, 2020, biostat)
LYCFOlateonset.next.dose <- function(curDose, phi, tau, 
                                   enter.times, dlt.times, 
                                   current.t, doses, 
                                   tover.doses, 
                                   simu=FALSE, add.args=list()){
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
  # i.e., impute the pending data
  
  
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