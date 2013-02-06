rm(list=ls())

model <- function(rm1, rm2, rm3, 
                  N1.start, N2.start, N3.start, 
                  K1, K2, K3,
                  evar1, evar2, evar3,
                  beat12, beta13, beta21, beta23, beta31, beta32,
                  time){
  r1 = numeric(time)
  r2 = numeric(time)
  r3 = numeric(time)
  N1 = numeric(time)
  N2 = numeric(time)
  N3 = numeric(time)
  Ntot = numeric(time)
  
  N1[1] <- N1.start
  N2[1] <- N2.start
  N3[1] <- N3.start
  Ntot[1] <- N1[1] + N2[1] +N3[1]
  
  for(t in 2:time){
    #   whitevar1 <- rnorm(1, mean=0)
    whitevar2 <- rnorm(1, mean=0)
    whitevar1 <- rnorm(1,mean=whitevar2)
    whitevar3 <- -1*rnorm(1, mean=whitevar1)
    dnoise1 <- rnorm(1, mean=0)
    dnoise2 <- rnorm(1, mean=0)
    dnoise3 <- rnorm(1, mean=0)
    
    sum.non1 = ((beta12 * N2[t-1])/K2) + ((beta13 * N3[t-1])/K3)
    sum.non2 = ((beta21 * N1[t-1])/K1) + ((beta23 * N3[t-1])/K3)
    sum.non3 = ((beta31 * N1[t-1])/K1) + ((beta32 * N2[t-1])/K2)
    
    r1[t] = rm1 * (1-(N1[t-1]/K1) - sum.non1) + evar1 * whitevar1 + (dvar1 * dnoise1/sqrt(N1[t-1]))
    
    r2[t] = rm2 * (1-(N2[t-1]/K2) - sum.non2) + evar2 * whitevar2 + (dvar2 * dnoise2/sqrt(N2[t-1]))
    
    r3[t] = rm3 * (1-(N3[t-1]/K2) - sum.non3) + evar3 * whitevar3 + (dvar3 * dnoise3/sqrt(N3[t-1]))
    
    
    N1[t] = N1[t-1] + N1[t-1]*r1[t]
    N2[t] = N2[t-1] + N2[t-1]*r2[t]
    N3[t] = N3[t-1] + N3[t-1]*r3[t]
    Ntot[t] = N1[t] + N2[t] + N3[t]
  }
  
  pops <- matrix(ncol=4, nrow=time)
  pops[,1] <- N1
  pops[,2] <- N2
  pops[,3] <- N3
  pops[,4] <- Ntot
  return(pops)
}



rm1 = 0.5
rm2 = 0.8
rm3 = 0.6
K1 = 1000
K2 = 1500
K3 = 1000
evar1 = 0.02
evar2 = 0.1
evar3 = 0.02
dvar1 = 1
dvar2 = 1
dvar3 = 1
beta12 = 0.8
beta13 = 0.08
beta21 = 0.1
beta23 = 0.08
beta31 = 0.1
beta32 = 0.7

time = 600


N1.start <- 200
N2.start <- 1300
N3.start <- 500


model.null <- model(rm1, rm2, rm3, 
                    N1.start, N2.start, N3.start, 
                    K1, K2, K3,
                    evar1, evar2, evar3,
                    beat12, beta13, beta21, beta23, beta31, beta32,
                    time)

N1 <- model.null[,1]
N2 <- model.null[,2]
N3 <- model.null[,3]
Ntot <- model.null[,4]

cv.tot = sd(Ntot)/mean(Ntot)
cv.tot

time.plot = seq(1,time,1)
plot(time.plot, Ntot, type="l", col="grey45", ylim=c(200,2500), ylab="Biomass", xlab="Time")
lines(time.plot, N1, col="lightblue")
lines(time.plot, N2, col="orange")
lines(time.plot, N3, col="purple")
abline(h=mean(Ntot), col="black", lwd=2)
abline(h=(sd(Ntot)+mean(Ntot)), col="black", lty="dashed")
abline(h=(mean(Ntot)-sd(Ntot)), col="black", lty="dashed")

abline(h=mean(N1), col="cadetblue", lwd=2)
abline(h=(sd(N1)+mean(N1)), col="cadetblue", lty="dashed")
abline(h=(mean(N1)-sd(N1)), col="cadetblue", lty="dashed")

abline(h=mean(N2), col="darkorange", lwd=2)
abline(h=(sd(N2)+mean(N2)), col="darkorange", lty="dashed")
abline(h=(mean(N2)-sd(N2)), col="darkorange", lty="dashed")

abline(h=mean(N3), col="purple4", lwd=2)
abline(h=(sd(N3)+mean(N3)), col="purple4", lty="dashed")
abline(h=(mean(N3)-sd(N3)), col="purple4", lty="dashed")
