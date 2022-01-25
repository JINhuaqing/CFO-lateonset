source("CFO_tox_utils.R")

make.move.fn <- function(ps, m=10){
    # ps: Output from move.dose.probs.fn 
    # m: number of samples to draw for the majority vote
    # output: 
    #    action: D--1, S--2, E--3
    nlevel <- length(ps)
    cps <- cumsum(ps)
    rvs <- runif(m)
    locs <- sapply(rvs, function(rv)rv<= cps)
    actions <- apply(locs, 2, which.max)
    if (m==1){
        final.action <- actions
    }else{
        res <- sapply(1:nlevel, function(i)actions==i)
        res <- colSums(res)
        final.action <- which.max(res)
    }
    final.action
}


under.eff.fn <- function(phiE, add.args=list()){
    x <- add.args$x
    n <- add.args$n
    alp.prior <- add.args$alp.prior.eff
    bet.prior <- add.args$bet.prior.eff
    ppE <- 1 - post.prob.fn(phiE, x, n, alp.prior, bet.prior)
    if ((ppE >= 0.9) & (n>=3)){
        return(TRUE)
    }else{
        return(FALSE)
    }
}


# The fn to return the probs of the movement based on the efficacy in the admissible set
# In fact, we don't need so complicated form. 
# We can simply select the next dose via comparing x/n for different doses
CFO.Eff.move.probs <- function(ad.xs, ad.ns, alp.prior, bet.prior){
    alps <- ad.xs + alp.prior
    bets <- ad.ns - ad.xs + bet.prior
    nd <- length(ad.xs)
    
    Nsps <- 10000
    sps.list <- list() 
    for (i in 1:nd){
        sps.list[[i]] <- rbeta(Nsps, alps[i], bets[i])
    }
    
    spss <- do.call(rbind, sps.list)
    argMaxs <- apply(spss, 2, which.max)
    probs <- as.vector(table(argMaxs))/Nsps
    
    probs
}

# Simulation function for CFO-phase I/II, allowing skipping dose when dose deescalation
CFO.Eff.simu.fn <- function(phi, phiE, p.true, pE.true, ncohort=10, init.level=1,  cohortsize=3, add.args=list(),
                            ph1=0){
    # phi: Target DIL rate
    # phiE: The minimal efficacy rate, only used for early stop
    # p.true: True DIL rates under the different dose levels
    # pE.true: True efficacy probs under the different dose levels
    # ncohort: The number of cohorts
    # cohortsize: The sample size in each cohort
    # alp.prior, bet.prior: prior parameters
    
    tadd.args <- add.args
    earlystop <- 0
    ndose <- length(p.true)
    cidx <- init.level
    
    tys <- rep(0, ndose) # number of DLT responses for different doses.
    txs <- rep(0, ndose) # number of efficacy responses for different doses.
    tns <- rep(0, ndose) # number of subject for different doses.
    tover.doses <- rep(0, ndose) # Whether each dose is too toxic or not, 1 yes.
    tunder.effs <- rep(0, ndose) # Whether the dose is not efficacious or not, 1 yes
    # if a dose is not efficacious enough or it is too toxic, it is would be eliminated from the admissible set.
    
    for (i in 1:ncohort){
        pc <- p.true[cidx] 
        pEc <- pE.true[cidx] 
        
        # sample from current dose
        cres <- rbinom(cohortsize, 1, pc)
        cEres <- rbinom(cohortsize, 1, pEc)
        
        # update results
        tys[cidx] <- tys[cidx] + sum(cres)
        txs[cidx] <- txs[cidx] + sum(cEres)
        tns[cidx] <- tns[cidx] + cohortsize
        
        #if (i == ncohort){
        #   break()
        #}
        
        cy <- tys[cidx]
        cx <- txs[cidx]
        cn <- tns[cidx]
        
        add.args <- c(list(y=cy, n=cn, x=cx, tys=tys, txs=txs, tns=tns, cidx=cidx), tadd.args)
        
        if (overdose.fn(phi, add.args)){
            tover.doses[cidx:ndose] <- 1
        }
        
        if (under.eff.fn(phiE, add.args)){
            tunder.effs[cidx] <- 1
        }else{
            tunder.effs[cidx] <- 0
        }
        
        
        if ((tover.doses[1] == 1) | (sum(tunder.effs[tover.doses==0])==sum(tover.doses==0))){
            earlystop <- 1
            break()
        }
        
        
        # the results for current 3 dose levels
        if (cidx!=1){
            cys <- tys[(cidx-1):(cidx+1)]
            cns <- tns[(cidx-1):(cidx+1)]
            cover.doses <- tover.doses[(cidx-1):(cidx+1)]
        }else{
            cys <- c(NA, tys[1:(cidx+1)])
            cns <- c(NA, tns[1:(cidx+1)])
            cover.doses <- c(NA, tover.doses[1:(cidx+1)])
        }
        
        # The up.idx of the admissible set
        up.idx <- make.decision.CFO.fn(phi, cys, cns, add.args$alp.prior, add.args$bet.prior, cover.doses) - 2 + cidx
        if (i<=ph1){
            cidx <- up.idx
        }else{
            
            if (up.idx == 1){
                cidx <- 1
            }else{
                low.idx <- 1
                ad.xs <- txs[low.idx:up.idx]
                ad.ns <- tns[low.idx:up.idx]
                if (length(ad.xs)==1){
                    cidx <- low.idx 
                }else{
                    probs <- CFO.Eff.move.probs(ad.xs, ad.ns, add.args$alp.prior.eff, add.args$bet.prior.eff)
                    cidx <- which.max(probs)
                }
                
            }
        }
        
    }
    
    
    if (earlystop==0){
        MTD <- select.mtd(phi, tns, tys)$MTD
        if ( (MTD==99) | (sum(tunder.effs[1:MTD]) == MTD)){
            OBD <- 99
        }else{
            OBD.probs <- CFO.Eff.move.probs(txs[1:MTD], tns[1:MTD], add.args$alp.prior.eff, add.args$bet.prior.eff)
            OBD <- which.max(OBD.probs)
        }
        
    }else{
        OBD <- 99
    }

    list(OBD=OBD, dose.ns=tns, eff.ns=txs, DLT.ns=tys, pE.true=pE.true, min.eff=phiE, 
         p.true=p.true, target=phi, over.doses=tover.doses, under.eff=tunder.effs)
    
}

peestimate<-function(yE, n){
    ndose<-length(yE)
    lik<-rep(0,ndose)
    pe<-(yE+0.05)/(n+0.1)
    p.e<-matrix(NA,ndose,ndose)
    for (i in 1:ndose){
        if (i==1) {
            x<-seq(ndose,1,by=-1)
        } else {
            x<-c(1:(i-1),seq(ndose,i))
        }
        #x<-x
        p.e[i,]<-ufit(pe,lmode=i,x=x,w=n+0.5)[[2]]
        lik[i]<-prod(dbinom(yE,n,p.e[i,]))		
    }
    lik<-lik/sum(lik)
    pe<-t(p.e)%*%lik +0.01*seq(1,ndose)
    return(pe)
}

util.fn <- function(phi, xs, ys, ns, Telimi, Eelimi){
    ndose <- length(xs)
    pT<-(ys+0.05)/(ns+0.1)
    pE<-(xs+0.05)/(ns+0.1)
    pT<-pava(pT,ns+0.1)  +0.001*seq(1,ndose)
    pE<-peestimate(xs,ns)
    
    pT[ns==0]<-20
    pT[Telimi==1]<-20
    pE[ns==0]<-0
    #pE[Eelimi==1]<-0
    #pE[Telimi==1]<-0
    u<-pE-0.33*pT-1.09*(pT>phi) 	
    d_opt<-which.max(u)	
    d_opt
}

