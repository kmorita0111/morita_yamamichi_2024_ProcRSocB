##-- parameters --##
t1 <- 10000  # time interval
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
G <- 0.01 # additice genetic variance
Sx1 <- 1 # width of distribution of strength of interspecific resource competition of sp.1
Sx2 <- 1 # width of distribution of strength of intraspecific resource competition of sp.2
Sy1 <- 1 # width of distribution of strength of reproductive interference of sp.1
Sy2 <- 1 # width of distribution of strength of reproductive interference of sp.2
Sl1 <- 1 # width of distribution of strength of reproductive interference of sp.1
Sl2 <- 1 # width of distribution of strength of reproductive interference of sp.2
K1 <-(lambda[1]*exp(-(x[1]-x1_h)^2/Sl1) - 1) / alpha[1]
K2 <-(lambda[2]*exp(-(x[2]-x2_h)^2/Sl2) - 1) / alpha[2]

Density01 <- matrix( 0, t1+1, 2 )
lambdaChange01 <- matrix( 0, t1, 2 )
alphaChange01 <- matrix( 0, t1, 2 )
betaChange01 <- matrix( 0, t1, 2 )
bChange01 <- matrix( 0, t1, 2 )
xChange01 <- matrix( 0, t1+1, 2 )
xChange01[1,] <- x
yChange01 <- matrix( 0, t1+1, 2 )
yChange01[1,] <- y
xyDensChange01 <- matrix( 0, t1+1, 3 )

###--- recorder ---###
reps <- 101
xMin <- 1.
xMax <- 1.1
yMin <- 0.95
yMax <- 1.05

color_mat01 <- 0
color_mat02 <- 0
color_mat03 <- 0
color_mat04 <- 0
color_mat002 <- matrix(NA,reps,reps)
color_mat003 <- matrix(NA,reps,reps)

##################################################
####=== simulation in supplementary figures ===###
##################################################
source( "~/morita-yamamichi2024/research_dynamics_optEtrait.R" )

Gx1 <- 0.01
Gx2 <- 0.

for ( i in 1:reps ) {
  print( paste( "i = ", as.character(i) ) )
  
  betaMax1 <- xMin + (xMax - xMin)/(reps-1)*(i-1)
  betaMax2 <- xMin + (xMax - xMin)/(reps-1)*(i-1)
  
  for ( j in 1:reps ) {
    print( paste( "j = ", as.character(j) ) )
    
    lambdaMax1 <- lambdaMax2*(yMin + (yMax - yMin)/reps*(j-1)) 
    lambdaMax2 <- 100 
    
    ###=== research dynamcis ===###
    x <- c(0,-0.1)
    N1 <- lambdaMax1*exp(-(x[1]-x1_h)^2)/alpha[1]; N2<-N1/99
    N <- c( N1, N2 )
    Density01[1,] <- N
    xChange01[1,] <- x
    xyDensChange01[1,] <- c(x[1],y[1],N[1]/(N[1]+N[2]))
    
    color_mat01 <- research_dynamics_optEtrait( b, bmax1, bmax2, y,
                                                alpha, alphaMax1, alphaMax2,
                                                beta, betaMax1, betaMax2, x, x1_h, x2_h,
                                                lambda, lambdaMax1, lambdaMax2,
                                                N,
                                                Gx1, Gx2, Gy1=0, rel_ii, rel_ij,
                                                Sl1, Sl2, Sx1, Sx2, Sy1, Sy2,
                                                Density01, 
                                                lambdaChange01, alphaChange01, betaChange01, bChange01, 
                                                xChange01, yChange01, 
                                                xyDensChange01,
                                                t1, sp=1, mod=2 
    )
    
    if( color_mat01 < 0.95 ){ color_mat003[i,j] <- 0 }
    
    
    x <- c(1,-0.1)
    N2<-(lambdaMax2*exp(-(x[2]-x2_h)^2/Sl2)-1)/alpha[2]; N1<-N2/4
    N <- c( N1, N2 )
    Density01[1,] <- N
    xChange01[1,] <- x
    xyDensChange01[1,] <- c(x[1],y[1],N[1]/(N[1]+N[2]))
    
    color_mat02 <- research_dynamics_optEtrait( b, bmax1, bmax2, y,
                                                alpha, alphaMax1, alphaMax2,
                                                beta, betaMax1, betaMax2, x, x1_h, x2_h,
                                                lambda, lambdaMax1, lambdaMax2,
                                                N,
                                                Gx1, Gx2, Gy1=0, rel_ii, rel_ij,
                                                Sl1, Sl2, Sx1, Sx2, Sy1, Sy2,
                                                Density01, 
                                                lambdaChange01, alphaChange01, betaChange01, bChange01, 
                                                xChange01, yChange01, 
                                                xyDensChange01,
                                                t1, sp=1, mod=2 
    )
    
    if( color_mat02 > 0.95 ){ color_mat003[i,j] <- 0.9 }
    
  }
}

###########################################
###=== Drawing color map in figure 3 ===###
###########################################
###--- load libraries ---###
#library(maps)
#library(fields)
source('~/morita-yamamichi2024/make_color_map_fig02.R')

###--- figure S4 ---###
fname <- "~/morita-yamamichi2024/result/figureS04"
xlabname <- expression(paste("Maximum interspecific resource competition, ", alpha[ij][max]))
ylabname <- expression(paste("Ratio of maximum growth rate, ", lambda[1][max], '/', lambda[2][max]))

pdf( paste( fname, "pdf", sep = "." ), width = 8, height = 7 )
par(mar=c(5, 7, 2, 3)) # margin

###=== boundary of extinction of either species ===### 
source( "~/morita-yamamichi2024/extinctBoundary.R" )

reps01 = (xMax-xMin)*100*100+1
reps02 = (yMax-yMin)*100*100+1

row_extBound <- (xMin + (xMax-xMin)*(which( extBound > 0, arr.ind = TRUE )[,1]-1)/(reps01-1))
low_extBound <- (yMin + (yMax-yMin)*(which( extBound > 0, arr.ind = TRUE )[,2]-1)/(reps02-1))
simu_x <- c( 1, row_extBound[-14] )
simu_y <- c( 1, low_extBound[-14] )

###--- Simulation result ---###
make_color_map( color_mat003, fname, xlabname, ylabname, x[1], Sx1, Sx2, Sl1, Sl2,
                xMax, xMin, yMax, yMin, betaMax1, betaMax2, axis_leg=c(10,10)
)
###--- always coexistence ---###
par(new=T)
alphaijmax <- seq(xMin,1/yMin,length=1000) 
plot( alphaijmax, 1/alphaijmax, 
      xlim = c(xMin,xMax), ylim = c(yMin,yMax), 
      xlab = "", ylab = "", type = "l", lwd = 2, 
      xaxs="i", yaxs="i", xaxt="n", yaxt="n"
)
###--- extinction of sp.2 ---###
par(new=T)
plot( simu_x, simu_y, xlim = c(xMin,xMax), ylim = c(yMin,yMax),
      xlab="", ylab="", type="l", lty="solid", lwd=1, 
      xaxs="i", yaxs="i", xaxt="n", yaxt="n" 
)

dev.off()

print( "Figure S5 completed!" )