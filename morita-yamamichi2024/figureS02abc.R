###=== parameters ===###
rMax1 <- 100.; rMax2 <- 100. # growth rate of species 1 & 2
alpha <- c(1,1) # coefficient of intraspecific competition of species 1 & 2
betaMax1 <- 1.1; betaMax2 <- 1.1  # coefficient of interspecific competition of species 1 & 2
here <- 0.01 # additive genetic variance of species 1 & 2
sigma <- 1
dt <- 0.001

t <- 10
interv <- seq(0,t,by=dt)
interval <- length( interv ) - 1

###=== initial values ===###
x2 = -0.1
x <- c( 0.5, x2 )
xChange01 <- matrix(0,(interval+1),2)
xChange01[1,] <- x
xChange02 <- matrix(0,(interval+1),2)
xChange02[1,] <- x
opt_x <- 0
r1 <- rMax1*exp(-(x[1]-opt_x)^2)
r2 <- rMax2*exp(-(x[2]-opt_x)^2)
r <- c(r1,r2)
K1 <- r1/alpha[1]
K2 <- r2/alpha[2]
beta1 <- betaMax1*exp(-(x[1]-x[2])^2)
beta2 <- betaMax2*exp(-(x[1]-x[2])^2)
beta <- c(beta1,beta2)

nChange01 <- matrix(0,(interval+1),2)
nChange01[1,] <- c(K1, 0) 
nChange02 <- matrix(0,(interval+1),2)
nChange02[1,] <- c(K1, 0)

#########################################
###=== population & trait dynamics ===###
#########################################
source( "~/morita-yamamichi2024/calc_model02.R" )

###--- simulation ---###
sp11 <- 0.2

for (timei in 1:interval) {
  
  if( timei*dt == sp11 ){ nChange01[timei,2] <- nChange01[timei,1]/9 }
  
  ###--- changing intrinsic growth rate & interspecific competition ---###
  r[1] <- rMax1*exp( -(xChange01[timei+1,1]-opt_x)^2 )
  r[2] <- rMax2*exp( -(xChange01[timei+1,2]-opt_x)^2 )
  beta[1] <- betaMax1*exp( -(xChange01[timei+1,1]-xChange01[timei+1,2])^2 )
  beta[2] <- betaMax2*exp( -(xChange01[timei+1,1]-xChange01[timei+1,2])^2 )
  
  ###--- population dynamics & character displacement ---###
  result <- ecd_ad( nChange01[timei,], xChange01[timei,], opt_x, r, rMax1, rMax2, alpha, beta, betaMax1, betaMax2, here, sigma, dt )
  
  nChange01[timei+1,] <- result[c(1,2)]
  xChange01[timei+1,] <- result[c(3,4)]
}

###--- making csv file ---###
###--- figure S02a ---###
fname <- "~/morita-yamamichi2024/result/figureS02a"
pdf( paste( fname, "pdf", sep = "." ), width = 5, height = 3 )
par(mar=c(5, 7, 1, 3)) # margin

yMax <- 100

t1 <- 8
interv01 <- seq(0,t1,by=dt)
interval01 <- length( interv01 ) 

plot( interv01[seq(1,interval01,by=10)], rep(0, length=length(interv01[seq(1,interval01,by=10)])), 
      xlim = c(0,t1), ylim = c(0,yMax), xlab = "", ylab = "",
      lwd = 3, type = "l", lty="dashed", col = "black"
)
par(new=T)
nChange01[1:(sp11/dt),2] <- NA
plot( interv01, nChange01[1:interval01,2], 
      xlim = c(0,t1), ylim = c(0,yMax), 
      xlab = expression(paste("Time, ",t)), 
      ylab = expression(paste("Species density, ", N[i])), 
      lwd = 3, type = "l", col = "blue"
) 
par(new=T)
plot( interv01, nChange01[1:interval01,1], 
      xlim = c(0,t1), ylim = c(0,yMax),
      lwd = 3, type = "l", xlab = "", ylab = "", col = "red"
)

dev.off()
print("Figure S2a completed!")

sp12 <- 0.3

for (timei in 1:interval) {
  
  if( timei*dt == sp12 ){ nChange02[timei,2] <- nChange02[timei,1]/9 }
  
  ###--- changing intrinsic growth rate & interspecific competition ---###
  r[1] <- rMax1*exp( -(xChange02[timei+1,1]-opt_x)^2 )
  r[2] <- rMax2*exp( -(xChange02[timei+1,2]-opt_x)^2 )
  beta[1] <- betaMax1*exp( -(xChange02[timei+1,1]-xChange02[timei+1,2])^2 )
  beta[2] <- betaMax2*exp( -(xChange02[timei+1,1]-xChange02[timei+1,2])^2 )
  
  ###--- population dynamics & character displacement ---###
  result <- ecd_ad( nChange02[timei,], xChange02[timei,], opt_x, r, rMax1, rMax2, alpha, beta, betaMax1, betaMax2, here, sigma, dt )
  
  nChange02[timei+1,] <- result[c(1,2)]
  xChange02[timei+1,] <- result[c(3,4)]
}

###--- figure S02b ---###
fname <- "~/morita-yamamichi2024/result/figureS02b"
pdf( paste( fname, "pdf", sep = "." ), width = 5, height = 3 )
par(mar=c(5, 7, 1, 3)) # margin

yMax <- 100

t2 <- 3 
interv02 <- seq(0,t2,by=dt)
interval02 <- length( interv02 ) 

plot( interv02[seq(1,interval02,by=10)], rep(0, length=length(interv02[seq(1,interval02,by=10)])),
      xlim = c(0,t2), ylim = c(0,yMax), xlab = "", ylab = "",
      lwd = 3, type = "l", lty="dashed", col = "black"
)
par(new=T)
nChange02[1:(sp12/dt),2] <- NA
plot( interv02, nChange02[1:interval02,2], 
      xlim = c(0,t2), ylim = c(0,yMax), 
      xlab = expression(paste("Time, ",t)), 
      ylab = expression(paste("Species density, ", N[i])), 
      lwd = 3, type = "l", col = "blue"
) 
par(new=T)
plot( interv02, nChange02[1:interval02,1], 
      xlim = c(0,t2), ylim = c(0,yMax),
      lwd = 3, type = "l", xlab = "", ylab = "", col = "red"
)

dev.off()
print("Figure S2b completed!")

###--- figure S02c ---###
t <- 8
interv <- seq(0,t,by=dt)
interval <- length( interv ) - 1

fname <- "~/morita-yamamichi2024/result/figureS02c"
pdf( paste( fname, "pdf", sep = "." ), width = 5, height = 3 )
par(mar=c(5, 7, 1, 3)) # margin

yMax <- 0.8

plot( interv[seq(1,interval,by=10)], rep(0, length=length(interv[seq(1,interval,by=10)])), 
      xlim = c(0,t), ylim = c(x2,yMax), xlab = "", ylab = "",
      lwd = 3, type = "l", lty="dashed", col = "green4"
)
par(new=T)
xChange02[1:(sp12/dt),2] <- NA
plot( interv, xChange02[1:(interval+1),2], 
      xlim = c(0,t), ylim = c(x2,yMax), 
      xlab = expression(paste("Time, ",t)), 
      ylab = expression(paste("Competitive traits, ", x[i])), 
      lwd = 3, type = "l", col = "blue"
) 
par(new=T)
xChange01[1:(sp11/dt),2] <- NA
plot( interv, xChange01[1:(interval+1),2], 
      xlim = c(0,t), ylim = c(x2,yMax), 
      xlab = "", ylab = "", 
      lwd = 3, type = "l", lty="dotted", col = "blue"
) 
par(new=T)
plot( interv, xChange02[1:(interval+1),1], 
      xlim = c(0,t), ylim = c(x2,yMax), 
      xlab = "", ylab = "", 
      lwd = 3, type = "l", col = "red"
) 
par(new=T)
plot( interv[seq( 1, (interval+1), by=10 )], xChange01[seq( 1, (interval+1), by=10 ),1], 
      xlim = c(0,t), ylim = c(x2,yMax),
      lwd = 3, type = "l", lty="dashed", xlab = "", ylab = "", col = "red"
)

dev.off()
print("Figure S2c completed!")