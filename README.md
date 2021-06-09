# transportation

An R package for statistical modelling and inference in transportation research.
 

## Install

DynamicLatticeBasis is most easily installed from GitHub using devtools:

```
install.packages("devtools")
devtools::install_github("MartinLHazelton/transportat")
```


## Code for simulation results

The code below implements the simulation results in 

Hazelton, M.L., (2021). The dynamics of Stochastic User Equilibrium, submitted for publication.

```
devtools::install_github("MartinLHazelton/transportation")
library(transportation)

# Watling's two route example - logit

N.sim <- 1
N.days <- 2e4
theta <- .5
psi <- .05

ODpair <- c(1,1)  # Index vector - route j services OD pair i indicated by jth element being i
ODdemand <- 20    # Overall demand - later to be scaled by zeta
Alpha <- c(0,10) + 1/0.15      
Beta <- c(10,1)     
pow <- c(4,0)
A <- diag(2)        # Link-path incidence matrix

zeta <- 1     # Demand multiplier
ODdemand <- ODdemand*zeta
Beta <- Beta*zeta     # My functions take in actual demand. To compensate, adjust coefficients b in link cost functions

eg1.SUE <- SUE(ODdemand,ODpair,A=A,Alpha=Alpha,Beta=Beta,pow=pow,theta=theta,verbose=F)
eg1.GSUE <- GSUE(ODdemand,ODpair,A=A,Alpha=Alpha,Beta=Beta,pow=pow,theta=theta,verbose=F) 
GSUE(ODdemand,ODpair,A=A,Alpha=Alpha,Beta=Beta,pow=pow,theta=theta,prob.model="Poisson",verbose=F) 

set.seed(2021)

eg1.sim.TIDDS1 <- rDDSProcess(ODdemand,ODpair,A,N.days=N.days,N.sim=N.sim,x.start=c(20,0),Alpha=Alpha,Beta=Beta,pow=pow,RUM="logit",theta=theta,process="TIDDS1",variable.demand=F)
rowMeans(eg1.sim.TIDDS1$x[1,,])
eg1.sim.TIDDS2 <- rDDSProcess(ODdemand,ODpair,A,N.days=N.days,N.sim=N.sim,x.start=c(20,0),Alpha=Alpha,Beta=Beta,pow=pow,RUM="logit",theta=theta,process="TIDDS2",variable.demand=F)
rowMeans(eg1.sim.TIDDS2$x[1,,])

plot(1:N.days,cumsum(eg1.sim.TIDDS1$x[1,1,])/1:N.days,type="l")

pdf(file="../tex/graphics/fig1a.pdf",height=3,width=6)
par(mar=c(4,4,0,0.5)+0.1)
plot(1:N.days,eg1.sim.TIDDS1$u[1,1,],type="l",ylim=c(13.5,14.2),ylab=expression(paste("Disutility, ",u[1]^t)),xlab=" ")
lines(1:N.days,eg1.sim.TIDDS2$u[1,1,],col=grey(0.7))
abline(h=eg1.SUE$u,lty=2,col=1,lwd=2)
abline(h=eg1.GSUE$u,lty=2,col=grey(0.7),lwd=2)
legend("topright",lty=c(1,1,2,2),col=c("black",grey(0.7)),lwd=c(1,1,2,2),legend=c("TIDDS1","TIDDS2",expression(u[SUE]),expression(u[GSUE])))
dev.off()

pdf(file="../tex/graphics/fig1b.pdf",height=3,width=6)
par(mar=c(4,4,0,0.5)+0.1)
thin <- seq(1,N.days,by=50)
plot(thin,eg1.sim.TIDDS1$x[1,1,thin],type="l",ylab=expression(paste("Flow, ",x[1]^t)),xlab="Time, t",ylim=c(10,23))
lines(thin,eg1.sim.TIDDS2$x[1,1,thin],col=grey(0.7))
abline(h=eg1.SUE$x,lty=2,col=1,lwd=2)
abline(h=eg1.GSUE$x,lty=2,col=grey(0.7),lwd=2)
legend("topright",lty=c(2,2),col=c("black",grey(0.7)),lwd=c(2,2),legend=c(expression(mu[SUE]),expression(mu[GSUE])))
dev.off()

# Figure of 8 network from Hazelton & Watling (Trans Sci, 2004)

A <- matrix(c(0,1,0,0,1,0,0,0,0,1,1,0,0,0,1,0,0,0,0,1,1,0,0,0,0,0,0,1),ncol=4,byrow=T)

# Plotting network

nodes.x <- c(2,7,12,2,7,12)
nodes.y <- c(7,7,7,2,2,2)
link.node.matrix <- rbind(c(1,2),c(1,4),c(2,5),c(3,2),c(3,6),c(4,5),c(6,5))
plotnet(nodes.x,nodes.y,link.node.matrix,ARROW.LENGTH=0.3,link.sep=0,radius=0.7,RightAbove=T)
symbols(nodes.x[c(1,3,5)],nodes.y[c(1,3,5)],circles=rep(0.6,3),add=T,inches=FALSE,lwd=2)

pdf(file="../tex/graphics/fig2.pdf",height=3.5,width=6.5)
par(mar=c(1,1,1,1)*0.1)
nodes.x <- c(2,7,12,2,7,12)
nodes.y <- c(7,7,7,2,2,2)
link.node.matrix <- rbind(c(1,2),c(1,4),c(2,5),c(3,2),c(3,6),c(4,5),c(6,5))
plotnet(nodes.x,nodes.y,link.node.matrix,ARROW.LENGTH=0.3,link.sep=0,radius=0.7,RightAbove=T)
symbols(nodes.x[c(1,3,5)],nodes.y[c(1,3,5)],circles=rep(0.6,3),add=T,inches=FALSE,lwd=2)
dev.off()

Alpha <- rep(10,7)
Beta0 <- rep(1,7)
pow <- rep(4,7)
ODpair <- c(1,1,2,2)
ODdemand0 <- c(1,1)
ODdemand0.rep <- rep(ODdemand0,unname(table(ODpair)))

N.sim <- 1
N.days <- 100000
theta <- 0.7
burnin <- 1000

psi.grid <- seq(0.01,0.99,by=0.02)
zeta.grid <- c(5,50,500)
Eight.RMSE.SUE <- Eight.RMSE.GSUE <- matrix(NA,nrow=length(psi.grid),ncol=length(zeta.grid))
Eight.GSUE.table <- Eight.SUE.table <- matrix(0,nrow=3,ncol=4)


set.seed(2021)

for (i in 1:length(zeta.grid)){
	zeta <- zeta.grid[i]
	ODdemand <- ODdemand0*zeta
	ODdemand.rep <- rep(ODdemand,unname(table(ODpair)))
	Beta <- Beta0*zeta
	Eight.GSUE <- GSUE(ODdemand,ODpair,A=A,Alpha=Alpha,Beta=Beta,pow=pow,theta=theta,verbose=F)
	Eight.SUE <- SUE(ODdemand,ODpair,A=A,Alpha=Alpha,Beta=Beta,pow=pow,theta=theta,verbose=F)
	Eight.GSUE.table[i,] <- Eight.GSUE$x/zeta
 	Eight.SUE.table[i,] <- Eight.SUE$x/zeta
	for (j in 1:length(psi.grid)){
		psi <- psi.grid[j]
		Eight.DDS <- rDDSProcess(ODdemand,ODpair,A,N.days=N.days,N.sim=N.sim,x.start=NA,Alpha=Alpha,Beta=Beta,pow=pow,RUM="logit",theta=theta,psi=psi,process="Standard",variable.demand=F)
		Eight.mean <- rowMeans(Eight.DDS$x[1,,-(1:burnin)])	
		Eight.RMSE.SUE[j,i] <- sqrt(mean((Eight.SUE$x-Eight.mean)^2/ODdemand.rep^2))
		Eight.RMSE.GSUE[j,i] <- sqrt(mean((Eight.GSUE$x-Eight.mean)^2/ODdemand.rep^2))
	}
}

Eight.GSUE.table 
Eight.SUE.table

plot(psi.grid,Eight.RMSE.GSUE[,1],type="l",xlab=expression(psi),ylab="RMSE",lty=1,ylim=c(0,0.06))
lines(psi.grid,Eight.RMSE.SUE[,1],type="l",lty=2)
lines(psi.grid,Eight.RMSE.GSUE[,2],type="l",col=2)
lines(psi.grid,Eight.RMSE.SUE[,2],type="l",col=2,lty=2)
lines(psi.grid,Eight.RMSE.GSUE[,3],type="l",col=3)
lines(psi.grid,Eight.RMSE.SUE[,3],type="l",col=3,lty=2)

pdf(file="../tex/graphics/fig3a.pdf",height=5,width=5)
par(mar=c(4,4,4,0.1))
plot(psi.grid,Eight.RMSE.GSUE[,1],type="l",xlab=expression(paste("Recency parameter ",psi)),ylab="RMSE",lty=1,ylim=c(0,0.06))
lines(psi.grid,Eight.RMSE.SUE[,1],type="l",lty=2)
legend("topleft",lty=1:2,legend=c("GSUE","SUE"))
title(expression(paste("Demand multiplier ",zeta==5)))
dev.off()

pdf(file="../tex/graphics/fig3b.pdf",height=5,width=5)
par(mar=c(4,4,4,0.1))
plot(psi.grid,Eight.RMSE.GSUE[,2],type="l",xlab=expression(paste("Recency parameter ",psi)),ylab="RMSE",lty=1,ylim=c(0,0.06))
lines(psi.grid,Eight.RMSE.SUE[,2],type="l",lty=2)
legend("topleft",lty=1:2,legend=c("GSUE","SUE"))
title(expression(paste("Demand multiplier ",zeta==50)))
dev.off()

pdf(file="../tex/graphics/fig3c.pdf",height=5,width=5)
par(mar=c(4,4,4,0.1))
plot(psi.grid,Eight.RMSE.GSUE[,3],type="l",xlab=expression(paste("Recency parameter ",psi)),ylab="RMSE",lty=1,ylim=c(0,0.06))
lines(psi.grid,Eight.RMSE.SUE[,3],type="l",lty=2)
legend("topleft",lty=1:2,legend=c("GSUE","SUE"))
title(expression(paste("Demand multiplier ",zeta==500)))
dev.off()

```

