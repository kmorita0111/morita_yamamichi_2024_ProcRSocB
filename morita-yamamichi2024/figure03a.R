###--- parameters ---###
b <- c( 0., 0. ) # maximum strength of reproductive interference
bmax1 <- b[1]; bmax2 <- b[2]
y <- c( 0., 0. ) # initial difference of reproductive traits
alpha <- c( 1., 1. ) # intraspecific resource competition
alphaMax1 <- 1.; alphaMax2 <- 1.
beta <- c( 1.1, 1.1 ) # maximum strength of interspecific resource competition
betaMax1 <- beta[1]; betaMax2 <- beta[2]
x20 <- -0.1
x <- c( 0.5, x20 ) # initial difference of ecological traits
x1_h <- 0 # species 1 optimum ecological trait
x2_h <- 0 # species 2 optimum ecological trait
lambda <- c( 100, 100 ) # fecundity
lambdaMax1 <- lambda[1]; lambdaMax2 <- lambda[2]
N <- c( 10, 99 ); afterN <- c( 0, 0 )
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

reps01 <- 61; reps02 <- 101 # mesh
xMin <- 0.; xMax <- 60.
yMin <- 0.95; yMax <- 1.05
zMin <- 1.; zMax <- 1.20

list.inv <- seq(xMin,xMax,length=reps01)
list.lambda <- seq(yMin,yMax,length=reps02)
list.alpha <- seq(zMin,zMax,length=reps02)

mat01 <- matrix( NA, reps02, reps01 )
mat02 <- matrix( NA, reps02, reps01 )
mat03 <- matrix( NA, reps02, reps01 )
mat04 <- matrix( NA, reps02, reps01 )
mat05 <- matrix( NA, reps02, reps01 )
mat06 <- matrix( NA, reps02, reps01 )

################################################################################
#
# Figure 3-1 (the change of lambda vs immigration timings)
# (a) betaMax1 = 1.03 (b) betaMax1 = 1.07 (c) betaMax1 = 1.14
#
################################################################################
###--- (a) simulation of trait & population dynamics ---###
source( "~/morita-yamamichi2024/research_dynamics_optEtrait.R" )

t1 <- 10000 # time interval

Density01 <- matrix( 0, t1+1, 2 )
lambdaChange01 <- matrix( 0, t1, 2 )
alphaChange01 <- matrix( 0, t1, 2 )
betaChange01 <- matrix( 0, t1, 2 )
bChange01 <- matrix( 0, t1, 2 )
xChange01 <- matrix( 0, t1+1, 2 )
xChange01[1,] <- x
yChange01 <- matrix( 0, t1+1, 2 )
yChange01[1,] <- y

betaMax1 <- 1.03; betaMax2 <- 1.03
for ( i in 1:reps02 ) {
  for( j in 1:reps01 ){

    xyDensChange01 <- matrix( NA, t1+1, 3 )
    DensTraitChange01 <- matrix( NA, t1+1, 4 )
    
    lambdaMax1 <- list.lambda[i]*lambdaMax2
    x <- c( 0.5, -0.1 )
    K1 <-(lambdaMax1*exp(-(x[1]-x1_h)^2/Sl1) - 1) / alpha[1]
    if( list.inv[j] > 0 ){ N <- c( K1, 0 ) 
    } else if( list.inv[j] == 0 ){ N <- c( K1, K1/9 ) }
    
    mat01[i,j] <- research_dynamics_optEtrait( b, bmax1, bmax2, y,
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
                                               t1, sp=list.inv[j], mod=2
                                              )
    
  }
}
print("Simulation finished - S3a")
###--- (b) simulation of trait & population dynamics ---###
betaMax1 <- 1.07; betaMax2 <- 1.07

Density01 <- matrix( 0, t1+1, 2 )
lambdaChange01 <- matrix( 0, t1, 2 )
alphaChange01 <- matrix( 0, t1, 2 )
betaChange01 <- matrix( 0, t1, 2 )
bChange01 <- matrix( 0, t1, 2 )
xChange01 <- matrix( 0, t1+1, 2 )
xChange01[1,] <- x
yChange01 <- matrix( 0, t1+1, 2 )
yChange01[1,] <- y

for ( i in 1:reps02 ) {
  for( j in 1:reps01 ){
    
    xyDensChange01 <- matrix( NA, t1+1, 3 )
    DensTraitChange01 <- matrix( NA, t1+1, 4 )
    
    lambdaMax1 <- list.lambda[i]*lambdaMax2
    x <- c( 0.5, -0.1 )
    K1 <-(lambdaMax1*exp(-(x[1]-x1_h)^2/Sl1) - 1) / alpha[1]
    if( list.inv[j] > 0 ){ N <- c( K1, 0 ) 
    } else if( list.inv[j] == 0 ){ N <- c( K1, K1/9 ) }
    
    mat02[i,j] <- research_dynamics_optEtrait( b, bmax1, bmax2, y,
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
                                               t1, sp=list.inv[j], mod=2
                                              )
    
  }
}
print("Simulation finished - 3a")
###--- (c) simulation of trait & population dynamics ---###
betaMax1 <- 1.14; betaMax2 <- 1.14

Density01 <- matrix( 0, t1+1, 2 )
lambdaChange01 <- matrix( 0, t1, 2 )
alphaChange01 <- matrix( 0, t1, 2 )
betaChange01 <- matrix( 0, t1, 2 )
bChange01 <- matrix( 0, t1, 2 )
xChange01 <- matrix( 0, t1+1, 2 )
xChange01[1,] <- x
yChange01 <- matrix( 0, t1+1, 2 )
yChange01[1,] <- y

for ( i in 1:reps02 ) {
  for( j in 1:reps01 ){
    
    xyDensChange01 <- matrix( NA, t1+1, 3 )
    DensTraitChange01 <- matrix( NA, t1+1, 4 )
    
    lambdaMax1 <- list.lambda[i]*lambdaMax2
    x <- c( 0.5, -0.1 )
    K1 <-(lambdaMax1*exp(-(x[1]-x1_h)^2/Sl1) - 1) / alpha[1]
    if( list.inv[j] > 0 ){ N <- c( K1, 0 ) 
    } else if( list.inv[j] == 0 ){ N <- c( K1, K1/9 ) }
    
    mat03[i,j] <- research_dynamics_optEtrait( b, bmax1, bmax2, y,
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
                                               t1, sp=list.inv[j], mod=2
    )
    
  }
}
print("Simulation finished - S3c")
################################################################################
#
# Figure 3-2 (the change of alphaij vs immigration timings)
# (a) lambdaMax1 = 97 (b) lambdaMax1 = 100 (c) lambdaMax1 = 104
#
################################################################################
Density01 <- matrix( 0, t1+1, 2 )
lambdaChange01 <- matrix( 0, t1, 2 )
alphaChange01 <- matrix( 0, t1, 2 )
betaChange01 <- matrix( 0, t1, 2 )
bChange01 <- matrix( 0, t1, 2 )
xChange01 <- matrix( 0, t1+1, 2 )
xChange01[1,] <- x
yChange01 <- matrix( 0, t1+1, 2 )
yChange01[1,] <- y

###--- (a) simulation of trait & population dynamics ---###
lambdaMax1 <- 97; lambdaMax2 <- 100

for ( i in 1:reps02 ) {
  for( j in 1:reps01 ){
    
    xyDensChange01 <- matrix( NA, t1+1, 3 )
    DensTraitChange01 <- matrix( NA, t1+1, 4 )
    
    betaMax1 <- list.alpha[i]; betaMax1 <- list.alpha[i]
    x <- c( 0.5, -0.1 )
    K1 <-(lambdaMax1*exp(-(x[1]-x1_h)^2/Sl1) - 1) / alpha[1]
    if( list.inv[j] > 0 ){ N <- c( K1, 0 ) 
    } else if( list.inv[j] == 0 ){ N <- c( K1, K1/9 ) }
    
    mat04[i,j] <- research_dynamics_optEtrait( b, bmax1, bmax2, y,
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
                                               t1, sp=list.inv[j], mod=2
    )
    
  }
}
print("Simulation finished - S3d")
###--- (b) simulation of trait & population dynamics ---###
lambdaMax1 <- 100; lambdaMax2 <- 100

for ( i in 1:reps02 ) {
  for( j in 1:reps01 ){
    
    xyDensChange01 <- matrix( NA, t1+1, 3 )
    DensTraitChange01 <- matrix( NA, t1+1, 4 )
    
    betaMax1 <- list.alpha[i]; betaMax1 <- list.alpha[i]
    x <- c( 0.5, -0.1 )
    K1 <-(lambdaMax1*exp(-(x[1]-x1_h)^2/Sl1) - 1) / alpha[1]
    if( list.inv[j] > 0 ){ N <- c( K1, 0 ) 
    } else if( list.inv[j] == 0 ){ N <- c( K1, K1/9 ) }
    
    mat05[i,j] <- research_dynamics_optEtrait( b, bmax1, bmax2, y,
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
                                               t1, sp=list.inv[j], mod=2
    )
    
  }
}
print("Simulation finished - S3e")
###--- (c) simulation of trait & population dynamics ---###
lambdaMax1 <- 104; lambdaMax2 <- 100

for ( i in 1:reps02 ) {
  for( j in 1:reps01 ){
    
    xyDensChange01 <- matrix( NA, t1+1, 3 )
    DensTraitChange01 <- matrix( NA, t1+1, 4 )

    betaMax1 <- list.alpha[i]; betaMax1 <- list.alpha[i]
    x <- c( 0.5, -0.1 )
    K1 <-(lambdaMax1*exp(-(x[1]-x1_h)^2/Sl1) - 1) / alpha[1]
    if( list.inv[j] > 0 ){ N <- c( K1, 0 ) 
    } else if( list.inv[j] == 0 ){ N <- c( K1, K1/9 ) }
    
    mat06[i,j] <- research_dynamics_optEtrait( b, bmax1, bmax2, y,
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
                                               t1, sp=list.inv[j], mod=2
    )
    
  }
}
print("Simulation finished - S3f")
################################################################################
# figure
################################################################################
###--- loading library ---###
library(maps)
library(fields)
source('~/morita-yamamichi2024/make_color_map.R')
###--- figure ---###
fname <- "~/morita-yamamichi2024/result/figureS03a"
pdf( paste( fname, "pdf", sep = "." ), width = 8, height = 7 )

xlabname <- expression(paste("Species 2's immigration timing, ", t ))
ylabname <- expression(paste("Ratio of maximum growth rate, ", lambda[1][max], '/', lambda[2][max]))

make_color_map( t(mat01), fname, xlabname, ylabname, x[1], Sx1, Sx2, Sl1, Sl2,
                xMax, xMin, yMax, yMin, betaMax1, betaMax2, axis_leg=c(12,10)
              )
par(new=T)
plot( c(25,25), c(yMin,yMax), xlab="", ylab="",  
      xlim=c(xMin,xMax), ylim = c(yMin,yMax), 
      xaxs="i", yaxs="i", xaxt="n", yaxt="n",
      lwd=3, type="l", lty="dashed", col="gray85" )
par(new=T)
plot( c(50,50), c(yMin,yMax), xlab="", ylab="",  
      xlim=c(xMin,xMax), ylim = c(yMin,yMax), 
      xaxs="i", yaxs="i", xaxt="n", yaxt="n",
      lwd=3, type="l", lty="dashed", col="gray85" )

dev.off()
print( "Figure S3-a completed!" )

###--- figure ---###
fname <- "~/morita-yamamichi2024/result/figure03a"
pdf( paste( fname, "pdf", sep = "." ), width = 8, height = 7 )

xlabname <- expression(paste("Species 2's immigration timing, ", t ))
ylabname <- expression(paste("Ratio of maximum growth rate, ", lambda[1][max], '/', lambda[2][max]))

make_color_map( t(mat02), fname, xlabname, ylabname, x[1], Sx1, Sx2, Sl1, Sl2,
                xMax, xMin, yMax, yMin, betaMax1, betaMax2, axis_leg=c(12,10)
              )
par(new=T)
plot( c(25,25), c(yMin,yMax), xlab="", ylab="",  
      xlim=c(xMin,xMax), ylim = c(yMin,yMax), 
      xaxs="i", yaxs="i", xaxt="n", yaxt="n",
      lwd=3, type="l", lty="dashed", col="gray85" )
par(new=T)
plot( c(50,50), c(yMin,yMax), xlab="", ylab="",  
      xlim=c(xMin,xMax), ylim = c(yMin,yMax), 
      xaxs="i", yaxs="i", xaxt="n", yaxt="n",
      lwd=3, type="l", lty="dashed", col="gray85" )

dev.off()
print( "Figure 3a, S3-b completed!" )

###--- figure ---###
fname <- "~/morita-yamamichi2024/result/figureS03c"
pdf( paste( fname, "pdf", sep = "." ), width = 8, height = 7 )

xlabname <- expression(paste("Species 2's immigration timing, ", t ))
ylabname <- expression(paste("Ratio of maximum growth rate, ", lambda[1][max], '/', lambda[2][max]))

make_color_map( t(mat03), fname, xlabname, ylabname, x[1], Sx1, Sx2, Sl1, Sl2,
                xMax, xMin, yMax, yMin, betaMax1, betaMax2, axis_leg=c(12,10)
              )
par(new=T)
plot( c(25,25), c(yMin,yMax), xlab="", ylab="",  
      xlim=c(xMin,xMax), ylim = c(yMin,yMax), 
      xaxs="i", yaxs="i", xaxt="n", yaxt="n",
      lwd=3, type="l", lty="dashed", col="gray85" )
par(new=T)
plot( c(50,50), c(yMin,yMax), xlab="", ylab="",  
      xlim=c(xMin,xMax), ylim = c(yMin,yMax), 
      xaxs="i", yaxs="i", xaxt="n", yaxt="n",
      lwd=3, type="l", lty="dashed", col="gray85" )

dev.off()
print( "Figure S3-c completed!" )

fname <- "~/morita-yamamichi2024/result/figureS03d"
pdf( paste( fname, "pdf", sep = "." ), width = 8, height = 7 )

xlabname <- expression(paste("Species 2's immigration timing, ", t ))
ylabname <- expression(paste("Maximum interspecific resource competition, ", alpha[ij][max]))

make_color_map( t(mat04), fname, xlabname, ylabname, x[1], Sx1, Sx2, Sl1, Sl2,
                xMax, xMin, zMax, zMin, betaMax1, betaMax2, axis_leg=c(12,10) 
              )
par(new=T)
plot( c(25,25), c(zMin,zMax), xlab="", ylab="",  
      xlim=c(xMin,xMax), ylim = c(zMin,zMax), 
      xaxs="i", yaxs="i", xaxt="n", yaxt="n",
      lwd=3, type="l", lty="dashed", col="gray85" )
par(new=T)
plot( c(50,50), c(zMin,zMax), xlab="", ylab="",  
      xlim=c(xMin,xMax), ylim = c(zMin,zMax), 
      xaxs="i", yaxs="i", xaxt="n", yaxt="n",
      lwd=3, type="l", lty="dashed", col="gray85" )

dev.off()
print( "Figure S3-d completed!" )

###--- figure ---###
fname <- "~/morita-yamamichi2024/result/figureS03e"
pdf( paste( fname, "pdf", sep = "." ), width = 8, height = 7 )

xlabname <- expression(paste("Species 2's immigration timing, ", t ))
ylabname <- expression(paste("Maximum interspecific resource competition, ", alpha[ij][max]))

make_color_map( t(mat05), fname, xlabname, ylabname, x[1], Sx1, Sx2, Sl1, Sl2,
                xMax, xMin, zMax, zMin, betaMax1, betaMax2, axis_leg=c(12,10)
              )
par(new=T)
plot( c(25,25), c(zMin,zMax), xlab="", ylab="",  
      xlim=c(xMin,xMax), ylim = c(zMin,zMax), 
      xaxs="i", yaxs="i", xaxt="n", yaxt="n",
      lwd=3, type="l", lty="dashed", col="gray85" )
par(new=T)
plot( c(50,50), c(zMin,zMax), xlab="", ylab="",  
      xlim=c(xMin,xMax), ylim = c(zMin,zMax), 
      xaxs="i", yaxs="i", xaxt="n", yaxt="n",
      lwd=3, type="l", lty="dashed", col="gray85" )

dev.off()
print( "Figure S3-e completed!" )

fname <- "~/morita-yamamichi2024/result/figureS03f"
pdf( paste( fname, "pdf", sep = "." ), width = 8, height = 7 )

xlabname <- expression(paste("Species 2's immigration timing, ", t ))
ylabname <- expression(paste("Maximum interspecific resource competition, ", alpha[ij][max]))

make_color_map( t(mat06), fname, xlabname, ylabname, x[1], Sx1, Sx2, Sl1, Sl2,
                xMax, xMin, zMax, zMin, betaMax1, betaMax2, axis_leg=c(12,10)
              )
par(new=T)
plot( c(25,25), c(zMin,zMax), xlab="", ylab="",  
      xlim=c(xMin,xMax), ylim = c(zMin,zMax), 
      xaxs="i", yaxs="i", xaxt="n", yaxt="n",
      lwd=3, type="l", lty="dashed", col="gray85" )
par(new=T)
plot( c(50,50), c(zMin,zMax), xlab="", ylab="",  
      xlim=c(xMin,xMax), ylim = c(zMin,zMax), 
      xaxs="i", yaxs="i", xaxt="n", yaxt="n",
      lwd=3, type="l", lty="dashed", col="gray85" )

dev.off()
print( "Figure S3-f completed!" )