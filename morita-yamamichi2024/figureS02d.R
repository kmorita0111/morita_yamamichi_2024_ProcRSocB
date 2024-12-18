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
xyDensChange01 <- matrix(0,(interval+1),2)
xyDensChange01[1,] <- c(x[1],1)
nChange02 <- matrix(0,(interval+1),2)
nChange02[1,] <- c(K1, 0)
xyDensChange02 <- matrix(0,(interval+1),2)
xyDensChange02[1,] <- c(x[1],1)

#########################################
###=== population & trait dynamics ===###
#########################################
source( "~/morita-yamamichi2024/calc_model02.R" )

###--- figure S02a ---###
sp11 <- 0.2

cnt <- 0
for (timei in 1:interval) {
  
  if( timei*dt >= sp11 && cnt < 1 ){ nChange01[timei,2] <- nChange01[timei,1]/9; cnt <- cnt+1 }
  
  ###--- changing intrinsic growth rate & interspecific competition ---###
  r[1] <- rMax1*exp( -(xChange01[timei+1,1]-opt_x)^2 )
  r[2] <- rMax2*exp( -(xChange01[timei+1,2]-opt_x)^2 )
  beta[1] <- betaMax1*exp( -(xChange01[timei+1,1]-xChange01[timei+1,2])^2 )
  beta[2] <- betaMax2*exp( -(xChange01[timei+1,1]-xChange01[timei+1,2])^2 )
  
  ###--- population dynamics & character displacement ---###
  result <- ecd_ad( nChange01[timei,], xChange01[timei,], opt_x, r, rMax1, rMax2, alpha, beta, betaMax1, betaMax2, here, sigma, dt )
  
  nChange01[timei+1,] <- result[c(1,2)]
  xChange01[timei+1,] <- result[c(3,4)]
  xyDensChange01[timei+1,] <- c( result[3], result[1]/sum(result[c(1,2)]) )
}

###--- late immigration ---###
sp12 <- 0.3

cnt <- 0
for (timei in 1:interval) {
  
  if( timei*dt == sp12 ){ nChange02[timei,2] <- nChange02[timei,1]/9; cnt <- cnt+1 }
  
  ###--- changing intrinsic growth rate & interspecific competition ---###
  r[1] <- rMax1*exp( -(xChange02[timei+1,1]-opt_x)^2 )
  r[2] <- rMax2*exp( -(xChange02[timei+1,2]-opt_x)^2 )
  beta[1] <- betaMax1*exp( -(xChange02[timei+1,1]-xChange02[timei+1,2])^2 )
  beta[2] <- betaMax2*exp( -(xChange02[timei+1,1]-xChange02[timei+1,2])^2 )
  
  ###--- population dynamics & character displacement ---###
  result <- ecd_ad( nChange02[timei,], xChange02[timei,], opt_x, r, rMax1, rMax2, alpha, beta, betaMax1, betaMax2, here, sigma, dt )
  
  nChange02[timei+1,] <- result[c(1,2)]
  xChange02[timei+1,] <- result[c(3,4)]
  xyDensChange02[timei+1,] <- c( result[3], result[1]/sum(result[c(1,2)]) )
}

###--- simultaneous immigration ---###
reps <- 101
xMin <- -0.3; xMax <- 1.3
yMin <- 0.; yMax <- 1.

N <- c(K1, 0)
here <- 0.001
t3 <- 100
interv <- seq(0,t3,by=dt)
interval <- length( interv )

nChange03 <- matrix( 0, interval+1, 2 )
xChange03 <- matrix( 0, interval+1, 2 )
color_mat00 <- matrix(NA,reps,reps)

for ( i in 1:reps ) {
  print( paste( "i = ", as.character(i) ) )
  
  x[1] <- xMin + (xMax - xMin)/(reps-1)*(i-1)
  x[2] <- -0.1
  
  for ( j in 1:reps ) {
    print( paste( "j = ", as.character(j) ) )
    
    mmm <- yMin + (yMax - yMin)/(reps-1)*(j-1)
    if( mmm < 1 ){ 
      N[2] <- (1-mmm)*K2; N[1] <- mmm*K2
    } else { 
      N[1] <- K1; N[2] <- 0 
    }
    
    nChange03[1,] <- N
    xChange03[1,] <- x

    for (timei in 1:interval) {
      
      ###--- changing intrinsic growth rate & interspecific competition ---###
      r[1] <- rMax1*exp( -(xChange03[timei,1]-opt_x)^2 )
      r[2] <- rMax2*exp( -(xChange03[timei,2]-opt_x)^2 )
      beta[1] <- betaMax1*exp( -(xChange03[timei,1]-xChange03[timei,2])^2 )
      beta[2] <- betaMax2*exp( -(xChange03[timei,1]-xChange03[timei,2])^2 )
      
      ###--- population dynamics & character displacement ---###
      result <- ecd_ad( nChange03[timei,], xChange03[timei,], opt_x, r, rMax1, rMax2, alpha, beta, betaMax1, betaMax2, here, sigma, dt )
      
      nChange03[timei+1,] <- result[c(1,2)]
      xChange03[timei+1,] <- result[c(3,4)]
      
    }
    
    sp1freq <- nChange03[(interval+1),1]/sum(nChange03[(interval+1),])
    if( is.nan(sp1freq) ){ print( "NaN!: " ); print( nChange03[1,] ); print( nChange03[(interval+1),] ) }
    color_mat00[i,j] <- sp1freq
  }
}

################################
###=== draw the figure S02d ===###
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
Ntotal <- rMax1
for(i in 1:mx){
  for(j in 1:my){
    
    r[1] <- rMax1*exp( -(x[i]-opt_x)^2 ); r[2] <- rMax2*exp( -(x2-opt_x)^2 )
    beta[1] <- betaMax1*exp( -(x[i]-x2)^2 ); beta[2] <- betaMax2*exp( -(x[i]-x2)^2 )

    z1[i,j] <- Ntotal*y[j]*(1-y[j])*( (r[1]-r[2])/Ntotal - (alpha[1]-beta[2])*y[j] - (beta[1]-alpha[2])*(1-y[j]) )
  }
}

Sx1=1
Sl1=1
# trait null-cline
for(i in 1:mx){
  for(j in 1:my){
    z2[i,j] <- -2*(x[i]-opt_x)*rMax1*exp( -(x[i]-opt_x)^2 ) + 2*(x[i]-x2)*betaMax1*exp( -(x[i]-x2)^2 )*Ntotal*(1-y[j])
  }
}

###--- loading library ---###
library(maps)
library(fields)
###--- axis names ---###
fname <- "~/morita-yamamichi2024/result/figureS02d"
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
lines(c(opt_x,opt_x),c(minY,maxY),
      col="green4",lwd=3, lty="dashed"
)
lines(c(minX,maxX),c((minY+0.01),(minY+0.01)),
      lwd=3, col="black"#, lty="dashed"
)
lines(c(minX,maxX),c((maxY-0.01),(maxY-0.01)),
      lwd=3, col="black"#, lty="dashed"
)
###--- population & trait dynamics ---###
t1 <- 7; t2 <- 2

par(new=T)
plot( xyDensChange01[1:(t1/dt),1], xyDensChange01[1:(t1/dt),2], 
      xlab = "", ylab = "", 
      xlim = c(xMin,xMax), ylim = c(yMin,yMax), 
      type = "p", pch=16, lwd=3, col = grey(seq(0.8,0.2,length=(t1/dt))), 
      xaxs="i", yaxs="i", xaxt="n", yaxt="n"
)
par(new=T)
plot( xyDensChange02[1:(t2/dt),1], xyDensChange02[1:(t2/dt),2], 
      xlab = "", ylab = "", 
      xlim = c(xMin,xMax), ylim = c(yMin,yMax), 
      type = "p", pch=16, lwd=3, col = grey(seq(0.8,0.2,length=(t2/dt))), 
      xaxs="i", yaxs="i", xaxt="n", yaxt="n" 
)

dev.off()

print("Figure S02d completed!")
