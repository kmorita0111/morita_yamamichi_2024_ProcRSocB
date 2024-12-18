###--- parameters ---###
b <- c( 0., 0. ) # maximum strength of reproductive interference
bmax1 <- b[1]
bmax2 <- b[2]
y <- c( 0., 0. ) # initial difference of reproductive traits
alpha <- c( 1., 1. ) # intraspecific resource competition
alphaMax1 <- 1.
alphaMax2 <- 1.
beta <- c( 1.1, 1.1 ) # maximum strength of interspecific resource competition
betaMax1 <- beta[1]
betaMax2 <- beta[2]
x20 <- -0.1
x <- c( 0.5, x20 ) # initial difference of ecological traits
x1_h <- 0 # species 1 optimum ecological trait
x2_h <- 0 # species 2 optimum ecological trait
lambda <- c( 100, 100 ) # fecundity
lambdaMax1 <- lambda[1]
lambdaMax2 <- lambda[2]
N <- c( 10, 99 )
afterN <- c( 0, 0 )
rel_ii <- c( 1, 0. ) # genetic covariance
rel_ij <- c( 0., 0. ) # genetic covariance
G <- 0.01 # additice genetic variance
Sx1 <- 1 # width of distribution of strength of interspecific resource competition
Sx2 <- 1 # width of distribution of strength of intraspecific resource competition
Sy1 <- 1 # width of distribution of strength of reproductive interference
Sy2 <- 1 # width of distribution of strength of reproductive interference
Sl1 <- 1 # width of distribution of strength of reproductive interference
Sl2 <- 1 # width of distribution of strength of reproductive interference
K1 <-(lambda[1]*exp(-(x[1]-x1_h)^2/Sl1) - 1) / alpha[1]
K2 <-(lambda[2]*exp(-(x[2]-x2_h)^2/Sl2) - 1) / alpha[2]

###--- recorder ---###
t1 <- 450 # time interval
Density01 <- matrix( 0, t1+1, 2 )
lambdaChange01 <- matrix( 0, t1, 2 )
alphaChange01 <- matrix( 0, t1, 2 )
betaChange01 <- matrix( 0, t1, 2 )
bChange01 <- matrix( 0, t1, 2 )
xChange01 <- matrix( 0, t1+1, 2 )
xChange01[1,] <- x
yChange01 <- matrix( 0, t1+1, 2 )
yChange01[1,] <- y

###--- simulation of trait & population dynamics ---###
source( "~/morita-yamamichi2024/research_dynamics_optEtrait.R" )

###--- early immigration ---###
xyDensChange01 <- matrix( NA, t1+1, 3 )
DensTraitChange01 <- matrix( NA, t1+1, 4 )
x <- c( 0.5, -0.1 )
N <- c( K1, 0 )

DensTraitChange01 <- research_dynamics_optEtrait( b, bmax1, bmax2, y,
                                                  alpha, alphaMax1, alphaMax2,
                                                  beta, betaMax1, betaMax2, x, x1_h, x2_h,
                                                  lambda, lambdaMax1, lambdaMax2,
                                                  N,
                                                  Gx1=0.01, Gx2=0, Gy1=0, rel_ii, rel_ij,
                                                  Sl1, Sl2, Sx1, Sx2, Sy1, Sy2,
                                                  Density01, 
                                                  lambdaChange01, alphaChange01, betaChange01, bChange01, 
                                                  xChange01, yChange01, 
                                                  xyDensChange01,
                                                  t1, sp=25, mod=3
                                                )
DensTraitChange01[1,] <- c(K1, 0, 0.5, -0.1)

###--- late-immigration ---###
t2 <- 450
Density01 <- matrix( 0, t2+1, 2 )
Density01[1,] <- N
lambdaChange01 <- matrix( 0, t2, 2 )
alphaChange01 <- matrix( 0, t2, 2 )
betaChange01 <- matrix( 0, t2, 2 )
bChange01 <- matrix( 0, t2, 2 )
xChange01 <- matrix( 0, t2+1, 2 )
yChange01 <- matrix( 0, t2+1, 2 )
xyDensChange01 <- matrix( 0, t2+1, 3 )
x <- c( 0.5, -0.1 )
xChange01[1,] <- x
yChange01[1,] <- c(0,0)
N <- c( K1, 0 )
sp2 <- 50

DensTraitChange02 <- matrix( NA, t2+1, 4 )
DensTraitChange02 <- research_dynamics_optEtrait( b, bmax1, bmax2, y,
                                                  alpha, alphaMax1, alphaMax2,
                                                  beta, betaMax1, betaMax2, x, x1_h, x2_h,
                                                  lambda, lambdaMax1, lambdaMax2,
                                                  N,
                                                  Gx1=0.01, Gx2=0, Gy1=0, rel_ii, rel_ij,
                                                  Sl1, Sl2, Sx1, Sx2, Sy1, Sy2,
                                                  Density01, 
                                                  lambdaChange01, alphaChange01, betaChange01, bChange01, 
                                                  xChange01, yChange01, 
                                                  xyDensChange01,
                                                  t2, sp=sp2, mod=3 
                                                )
DensTraitChange02[1,] <- c(K1, 0, 0.5, -0.1)

print( "Simulation completed!" )
###=== graphics for population & trait dynamics ===###
###--- figure (a) ---###
fname <- "~/morita-yamamichi2024/result/figure02a"
pdf( paste( fname, "pdf", sep = "." ), width = 5, height = 3 )
par(mar=c(5, 7, 1, 3)) # margin

tMax = 450
plot( 0:tMax, rep(0, length=(tMax+1)), 
      xlim = c(0,tMax), ylim = c(0,100), 
      type = "l", lty="dashed", lwd=3, col = "black", 
      xlab = "", ylab = ""  
    )
par(new=T)
DensTraitChange01[1:25,2] <- NA
plot( 0:tMax, DensTraitChange01[1:(tMax+1),2], 
      xlim = c(0,tMax), ylim = c(0,100), 
      type = "l", lwd=3, col = "blue", 
      xlab = expression(paste( "Time, ", t)), ylab = expression(paste("Species density, ", N[i])) 
    )
par(new=T)
plot( 0:tMax, DensTraitChange01[1:(tMax+1),1], 
      xlim = c(0,tMax), ylim = c(0,100), 
      type = "l", lwd=3, xlab = "", ylab = "", col = "red" 
    )

dev.off()

###--- figure (b) ---###
fname <- "~/morita-yamamichi2024/result/figure02b"
pdf( paste( fname, "pdf", sep = "." ), width = 5, height = 3 )
par(mar=c(5, 7, 1, 3)) # margin

#tMax = 450
plot( 0:tMax, rep(0, length=(tMax+1)), 
      xlim = c(0,tMax), ylim = c(0,100), 
      type = "l", lty="dashed", lwd=3, col = "black", 
      xlab = "", ylab = ""  
    )
par(new=T)
DensTraitChange02[1:sp2,2] <- NA
plot( 0:tMax, DensTraitChange02[1:(tMax+1),2], 
      xlim = c(0,tMax), ylim = c(0,100), 
      type = "l", lwd=3, col = "blue", 
      xlab = expression(paste( "Time, ", t)), ylab = expression(paste("Species density, ", N[i]))  
    )
par(new=T)
plot( 0:tMax, DensTraitChange02[1:(tMax+1),1], 
      xlim = c(0,tMax), ylim = c(0,100), 
      type = "l", lwd=3, xlab = "", ylab = "", col = "red" 
    )

dev.off()

###--- figure (c) ---###
fname <- "~/morita-yamamichi2024/result/figure02c"
pdf( paste( fname, "pdf", sep = "." ), width = 5, height = 3 )
par(mar=c(5, 7, 1, 3)) # margin

tMax = 450
yMax = 0.8
plot( 0:tMax, rep(0, length=(tMax+1)), 
      xlim = c(0,tMax), ylim = c(-0.1,yMax), 
      type = "l", lty="dashed", lwd=3, col = "green4", 
      xlab = "", ylab = ""  
    )
par(new=T)
DensTraitChange01[1:25,4] <- NA
plot( 0:tMax, DensTraitChange01[1:(tMax+1),4], 
      xlim = c(0,tMax), ylim = c(-0.1,yMax), 
      type = "l", lty="dotted", lwd=3, col = "blue", 
      xlab = expression(paste( "Time, ", t)), ylab = expression(paste("Competitive traits, ", x[i])) 
    )
par(new=T)
DensTraitChange02[1:sp2,4] <- NA
plot( 0:tMax, DensTraitChange02[1:(tMax+1),4], 
      xlim = c(0,tMax), ylim = c(-0.1,yMax), 
      type = "l", lwd=3, col = "blue", 
      xlab = "", ylab = "" 
    )
par(new=T)
plot( 0:tMax, DensTraitChange01[1:(tMax+1),3], 
      xlim = c(0,tMax), ylim = c(-0.1,yMax), 
      type = "l", lty="dashed", lwd=3, xlab = "", ylab = "", col = "red" 
    )
par(new=T)
plot( 0:tMax, DensTraitChange02[1:(tMax+1),3], 
      xlim = c(0,tMax), ylim = c(-0.1,yMax), 
      type = "l", lwd=3, xlab = "", ylab = "", col = "red" 
    )

dev.off()

print( "Figure 2a-c completed!" )
