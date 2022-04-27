#rm(list=ls())

library(reshape2)
library(ggplot2)
setwd("C:/Users/JINHU/Documents/ProjectCode/CFO-lateonset")
source("CFO_tox_utils.R")
source("./utilities.R")


# 1. test the value of gamma and OR
y1 <- 1
n1 <- 6
y2 <- 2
n2 <- 6
y3 <- 3
n3 <- 6
phi <- 0.25
alp.prior <- phi
bet.prior <- 1-phi
ORL <- OR.values(phi, y1, n1, y2, n2, alp.prior, bet.prior, type="L");ORL
optim.gamma.fn(n1, n2, phi, type="L", alp.prior, bet.prior)
ORR <- OR.values(phi, y2, n2, y3, n3, alp.prior, bet.prior, type="R");ORR
optim.gamma.fn(n2, n3, phi, type="R", alp.prior, bet.prior)
ORL <- OR.values(phi, y1, n1, y2+1, n2+1, alp.prior, bet.prior, type="L");ORL
optim.gamma.fn(n1, n2+1, phi, type="L", alp.prior, bet.prior)


# 2 draw the plots for gammaL
n1s <- 1:30
n2s <- 1:30
gams.raw.L <- matrix(NA, nrow=30, ncol=30)
for (n1 in n1s){
  for (n2 in n2s){
    print(c(n1, n2))
    gams.raw.L[n1, n2] <- optim.gamma.fn(n1, n2, phi, type="L", alp.prior, bet.prior)$gam
  }
}



gams <- as.data.frame(gams.raw.L)
rownames(gams) <- n1s
colnames(gams) <- n2s
gams$nL <- rownames(gams)

data <- melt(gams, id.vars=c("nL"))
data$nL <- as.numeric(data$nL)
data$nL <- as.factor(data$nL)

names(data) <- c("nL", "nC", "gammaL")
gammaL.new <- data$gammaL
gammaL.new[gammaL.new>1] <- NA
data$gammaL.new <- gammaL.new

p <- ggplot(data=data) + 
  geom_tile(mapping=aes(x=nL, y=nC, fill=gammaL.new, colour="")) + 
  scale_fill_gradient2(expression(gamma[L]), low="blue", mid="#e3c5c5", high="red",  
                       midpoint=0.5, limits=c(0, 1), na.value = "#940000", breaks=c(0, 0.5, 1.0, 1.3, 2.2)) +
  scale_colour_manual(values=NA) +              
  guides(colour=guide_legend(expression(gamma[L]>"1.0"), override.aes=list(colour="#940000", fill="#940000")))  +
  theme(legend.position = "right") + 
  ylab("Num of patients at the current dose") +  xlab("Num of patients at the left neighboring dose") +
  scale_x_discrete(labels=seq(5, 30, 5), breaks=seq(5, 30, 5)) + 
  scale_y_discrete(labels=seq(5, 30, 5), breaks=seq(5, 30, 5)); p

ggsave("./plots/gammaL25.jpg", height=5, width = 6, unit="in")



# 3 draw the plots for gammaR
n1s <- 1:30
n2s <- 1:30
gams.raw.R <- matrix(NA, nrow=30, ncol=30)
for (n1 in n1s){
  for (n2 in n2s){
    print(c(n1, n2))
    gams.raw.R[n1, n2] <- optim.gamma.fn(n1, n2, phi, type="R", alp.prior, bet.prior)$gam
  }
}


gams <- as.data.frame(gams.raw.R)
rownames(gams) <- n1s
colnames(gams) <- n2s
gams$nC <- rownames(gams)

data <- melt(gams, id.vars=c("nC"))
data$nC <- as.numeric(data$nC)
data$nC <- as.factor(data$nC)

names(data) <- c("nC", "nR", "gammaR")
gammaR.new <- data$gammaR
RCV <- 3.0 #0.25
#RCV <- 4.0 #0.2
#RCV <- 2.2 #phi<- 0.3
gammaR.new[gammaR.new>RCV] <- NA
data$gammaR.new <- gammaR.new

p <- ggplot(data=data) + 
  geom_tile(mapping=aes(x=nR, y=nC, fill=gammaR.new, colour="")) + 
  scale_fill_gradient2(expression(gamma[R]), low="blue", mid="#e3c5c5", high="red",  
                       midpoint=1.1, limits=c(0, RCV), na.value = "#940000", 
                       breaks=c(0, 1.0, 2.0, RCV) # 0.2
                       #breaks=c(0, 1.0, 2.0, 3.0, RCV) # 0.2
                       #breaks=c(0, 0.5, 1.0, 1.5, RCV) #0.3
                       ) +
  scale_colour_manual(values=NA) +              
  guides(colour=guide_legend(expression(gamma[R]>"3.0"), override.aes=list(fill="#940000")))  +
  theme(legend.position = "right") + 
  ylab(expression(m[C])) +  xlab(expression(m[R])) +
  ylab("Num of patients at the current dose") +  xlab("Num of patients at the right neighboring dose") +
  scale_x_discrete(labels=seq(5, 30, 5), breaks=seq(5, 30, 5)) + 
  scale_y_discrete(labels=seq(5, 30, 5), breaks=seq(5, 30, 5)) ;p

ggsave("./plots/gammaR25.jpg", height=5, width = 6, unit="in")



