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
x <- c( 0.1, x20 ) # initial difference of ecological traits
x1_h <- 0 # species 1 optimum ecological trait
x2_h <- 0 # species 2 optimum ecological trait
lambda <- c( 100, 100 ) # fecundity
lambdaMax1 <- lambda[1]
lambdaMax2 <- lambda[2]
#N10 <- 50 # initial density of species 1
N <- c( 10, 99 )
afterN <- c( 0, 0 )
rel_ii <- c( 1, 0. ) # genetic covariance
rel_ij <- c( 0., 0. ) # genetic covariance
agv <- 0.01 # additice genetic variance
Sx1 <- 1 # width of distribution of strength of interspecific resource competition
Sx2 <- 1 # width of distribution of strength of intraspecific resource competition
Sy1 <- 1 # width of distribution of strength of reproductive interference
Sy2 <- 1 # width of distribution of strength of reproductive interference
Sl1 <- 1 # width of distribution of strength of reproductive interference
Sl2 <- 1 # width of distribution of strength of reproductive interference
K1 <-(lambda[1]*exp(-(x[1]-x1_h)^2/Sl1) - 1) / alpha[1]
K2 <-(lambda[2]*exp(-(x[2]-x2_h)^2/Sl2) - 1) / alpha[2]

#############################################
###=== simulation in early immigration ===###
#############################################
###--- trait & population dynamics ---###
t1 <- 700  # time interval
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
x <- c( 0.5, -0.1 )
N <- c( K1, 0 )
xyDensChange01 <- research_dynamics_optEtrait( b, bmax1, bmax2, y,
                                               alpha, alphaMax1, alphaMax2,
                                               beta, betaMax1, betaMax2, x, x1_h, x2_h,
                                               lambda, lambdaMax1, lambdaMax2,
                                               N,
                                               Gx1=agv, Gx2=0, Gy1=0, rel_ii, rel_ij,
                                               Sl1, Sl2, Sx1, Sx2, Sy1, Sy2,
                                               Density01, 
                                               lambdaChange01, alphaChange01, betaChange01, bChange01, 
                                               xChange01, yChange01, 
                                               xyDensChange01,
                                               t1, sp=25, mod=1
)
xyDensChange01[1,] <- c(0.5, 0, 1)

###--- late-immigration ---###
t2 <- 220
Density01 <- matrix( 0, t2+1, 2 )
Density01[1,] <- N
lambdaChange01 <- matrix( 0, t2, 2 )
alphaChange01 <- matrix( 0, t2, 2 )
betaChange01 <- matrix( 0, t2, 2 )
bChange01 <- matrix( 0, t2, 2 )
xChange01 <- matrix( 0, t2+1, 2 )
yChange01 <- matrix( 0, t2+1, 2 )
xyDensChange02 <- matrix( 0, t2+1, 3 )
x <- c( 0.5, -0.1 )
xChange01[1,] <- x
yChange01[1,] <- c(0,0)
N <- c( K1, 0 )
xyDensChange02 <- research_dynamics_optEtrait( b, bmax1, bmax2, y,
                                               alpha, alphaMax1, alphaMax2,
                                               beta, betaMax1, betaMax2, x, x1_h, x2_h,
                                               lambda, lambdaMax1, lambdaMax2,
                                               N,
                                               Gx1=agv, Gx2=0, Gy1=0, rel_ii, rel_ij,
                                               Sl1, Sl2, Sx1, Sx2, Sy1, Sy2,
                                               Density01, 
                                               lambdaChange01, alphaChange01, betaChange01, bChange01, 
                                               xChange01, yChange01, 
                                               xyDensChange01,
                                               t2, sp=50, mod=1 
)
xyDensChange02[1,] <- c( 0.5, 0, 1 )

###--- simultaneous immigration ---###
reps <- 101
xMin <- -0.3; xMax <- 1.3
yMin <- 0.; yMax <- 1.

color_mat00 <- matrix(NA,reps,reps)
t3 <- 10000
xyDensChange03 <- matrix( 0, t3+1, 3 )
Density03 <- matrix( 0, t3+1, 2 )
lambdaChange03 <- matrix( 0, t3, 2 )
alphaChange03 <- matrix( 0, t3, 2 )
betaChange03 <- matrix( 0, t3, 2 )
bChange03 <- matrix( 0, t3, 2 )
xChange03 <- matrix( 0, t3+1, 2 )
yChange03 <- matrix( 0, t3+1, 2 )
yChange03[1,] <- y
xyDensChange03 <- matrix( 0, t3+1, 3 )

for ( i in 1:reps ) {
  print( paste( "i = ", as.character(i) ) )
  
  x[1] <- xMin + (xMax - xMin)/(reps-1)*(i-1); x[2] <- -0.1
  
  K1 <-(lambda[1]*exp(-(x[1]-x1_h)^2/Sl1) - 1) / alpha[1]
  K2 <-(lambda[2]*exp(-(x[2]-x2_h)^2/Sl2) - 1) / alpha[2]
  
  for ( j in 1:reps ) {
    print( paste( "j = ", as.character(j) ) )
    
    mmm <- yMin + (yMax - yMin)/(reps-1)*(j-1)
    if( mmm < 1 ){ 
      N[2] <- K2; N[1] <- mmm/(1-mmm)*N[2]
    } else { 
      N[1] <- K1; N[2] <- 0 
    }
    
    betaMax1 <- 1.1 
    betaMax2 <- 1.1 
    
    lambdaMax2 <- 100 
    lambdaMax1 <- 100 
    
    Density03[1,] <- N
    xChange03[1,] <- x
    xyDensChange03[1,] <- c(x[1],y[1],N[1]/(N[1]+N[2]))
    
    color_mat00[i,j] <- research_dynamics_optEtrait( b, bmax1, bmax2, y,
                                                     alpha, alphaMax1, alphaMax2,
                                                     beta, betaMax1, betaMax2, x, x1_h, x2_h,
                                                     lambda, lambdaMax1, lambdaMax2,
                                                     N,
                                                     Gx1=agv, Gx2=0, Gy1=0, rel_ii, rel_ij,
                                                     Sl1, Sl2, Sx1, Sx2, Sy1, Sy2,
                                                     Density03, 
                                                     lambdaChange03, alphaChange03, betaChange03, bChange03, 
                                                     xChange03, yChange03, 
                                                     xyDensChange03,
                                                     t3, sp=0, mod=2 
    )
  }
}

################################
###=== draw the figure 2d ===###
################################
###--- null-clnie ---###
reps <- 100
minX <- -0.3
maxX <- 1.3
minY <- 0
maxY <- 1
mx <- (maxX-minX)*reps+1
my <- (maxY-minY)*reps+1
x <- seq( minX, maxX, length = mx )
y <- seq( minY, maxY, length = my )
z1 <- matrix( 0, mx, my )
z2 <- matrix( 0, mx, my )
# population density null-cline
for(i in 1:mx){
  for(j in 1:my){
    z1[i,j] <- lambdaMax1*exp(-(x[i]-x1_h)^2)*(1+(y[j]/(1-y[j]))*betaMax2*exp(-(x[i]-x20)^2)) - lambdaMax2*exp(-(x20-x2_h)^2)*((y[j]/(1-y[j]))+betaMax1*exp(-(x[i]-x20)^2))
  }
}

Sx1=1; Sl1=1; alpha1=1; G=0.01
# trait null-cline
for(i in 1:mx){
  for(j in 1:my){
    z2[i,j] <- -2*G*(y[j]*(x[i]-x1_h)+betaMax1*exp(-(x[i]-x20)^2)*(1-y[j])*(x20-x1_h))/(y[j]+betaMax1*exp(-(x[i]-x20)^2)*(1-y[j]))
  }
}

###--- loading library ---###
library(maps)
library(fields)
###--- axis names ---###
fname <- "~/morita-yamamichi2024/result/figure02d"
pdf( paste( fname, "pdf", sep = "." ), width = 8, height = 7 )

xlabname <- expression(paste("Species 1's competitive trait, ", x[1] ))
ylabname <- expression(paste("Sp. 1's frequency in 2 sp., ", N[1]/(N[1]+N[2]) ))
source('~/morita-yamamichi2024/make_color_map.R')

make_color_map( color_mat00, fname, xlabname, ylabname, x[1], Sx1, Sx2, Sl1, Sl2,
                xMax, xMin, yMax, yMin, betaMax1, betaMax2, axis_leg=c(16,10) 
              )
###--- null-cline ---###
par(ps=15)
par(cex=1.5)
par(new=T)
contour( x, y, z1, 
         xlab = "", ylab = "", 
         drawlabels = F, levels=0, 
         col = "black",
         lwd = 3, cex=1, xaxs="i", yaxs="i", xaxt="n", yaxt="n"
)
par(new=T)
contour( x, y, z2, 
         xlab = "", ylab = "", 
         drawlabels = F, levels=0, 
         col = "grey80", lwd = 3,
         xaxs="i", yaxs="i", xaxt="n", yaxt="n"
)
lines(c(x1_h,x1_h),c(minY,maxY),
      col="green4",lwd=3, lty="dashed"
)
lines(c(minX,maxX),c((minY+0.01),(minY+0.01)),
      lwd=3, col="black"
)
lines(c(minX,maxX),c((maxY-0.01),(maxY-0.01)),
      lwd=3, col="black"
)
###--- population & trait dynamics ---###
par(new=T)
plot( xyDensChange01[1:(t1+1),1], xyDensChange01[1:(t1+1),3], 
      xlab = "", ylab = "", 
      xlim = c(xMin,xMax), ylim = c(yMin,yMax), 
      type = "p", pch=16, lwd=3, col = grey(seq(0.8,0.2,length=(t1+1))), 
      xaxs="i", yaxs="i", xaxt="n", yaxt="n"
)
par(new=T)
plot( xyDensChange02[1:(t2+1),1], xyDensChange02[1:(t2+1),3], 
      xlab = "", ylab = "", 
      xlim = c(xMin,xMax), ylim = c(yMin,yMax), 
      type = "p", pch=16, lwd=3, col = grey(seq(0.8,0.2,length=(t2+1))), 
      xaxs="i", yaxs="i", xaxt="n", yaxt="n" 
)

dev.off()

print("Figure2d completed!")