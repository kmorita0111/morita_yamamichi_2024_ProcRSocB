sim_traj00 = function( xeq1, xeq2, betaMax1, coex, res_ext, minX, maxX, minY, maxY ){
  
  x10 <- 0.1
  x2 <- -0.1
  opt_x <- 0
  
  lambdaMax <- 100
  
  trait <- seq( x10, xeq1, length=1001 )
  plot( 1-betaMax1*exp( -(x2-trait)^2 ), rep( ((lambdaMax-1)/(lambdaMax-1)), length=1001 ), 
        xlab="", ylab="", xlim = c(minX,maxX), ylim = c(minY,maxY),
        type = "l", lty = "dashed", col = "black", 
        lwd = 2, xaxs="i", yaxs="i"
      )
  points(  1-betaMax1*exp( -(x2-x10)^2 ), ((lambdaMax-1)/(lambdaMax-1)), 
           xlim = c(minX,maxX), ylim = c(minY,maxY), col = "blue", pch = 19 
      ) 
  if( coex ){
    
    points( 1-betaMax1*exp( -(x2-xeq1)^2 ), ((lambdaMax-1)/(lambdaMax-1)), 
            xlim = c(minX,maxX), ylim = c(minY,maxY), col = "black", pch = 19 
          ) 
  } 

}
### discrete-time 2 competitors model 
sim_traj01 = function( xeq1, xeq2, betaMax1, coex, res_ext, minX, maxX, minY, maxY ){
  
  x10 <- 0.1
  x2 <- -0.1
  opt_x <- 0
  
  lambdaMax <- 100
  
  trait <- seq( opt_x, xeq1, length=1001 )
  plot( 1-betaMax1*exp( -(x2-trait)^2 ), ((lambdaMax*exp( -(trait-opt_x)^2 )-1)/(lambdaMax*exp( -(x2-opt_x)^2 )-1)), 
        xlab="", ylab="", xlim = c(minX,maxX), ylim = c(minY,maxY),
        type = "l", lty = "dashed", col = "black", 
        lwd = 2, xaxs="i", yaxs="i"
      )
  points(  1-betaMax1*exp( -(x2-x10)^2 ), ((lambdaMax*exp( -(x2-opt_x)^2 )-1)/(lambdaMax*exp( -(x2-opt_x)^2 )-1)), 
           xlim = c(minX,maxX), ylim = c(minY,maxY), col = "blue", pch = 19 
        ) 
  if( coex ){
    
    points( 1-betaMax1*exp( -(x2-xeq1)^2 ), ((lambdaMax*exp( -(xeq1-opt_x)^2 )-1)/(lambdaMax*exp( -(x2-opt_x)^2 )-1)), 
            xlim = c(minX,maxX), ylim = c(minY,maxY), col = "black", pch = 19 
          ) 
  } 
  if( res_ext ){
    
    points( 1-betaMax1*exp( -(x2-xeq2)^2 ), ((lambdaMax*exp( -(xeq2-opt_x)^2)-1)/(lambdaMax*exp( -(x2-opt_x)^2 )-1)), 
            xlim = c(minX,maxX), ylim = c(minY,maxY), col = "black", pch = 19 
          ) 
  }
  
}
### continuous-time 2 competitors model 
sim_traj02 = function( xeq1, xeq2, betaMax1, coex, res_ext, minX, maxX, minY, maxY ){
  
  x10 <- 0.1
  x2 <- -0.1
  opt_x <- 0
  
  rMax <- 100
  
  trait <- seq( opt_x, xeq1, length=1001 )
  plot( 1-betaMax1*exp( -(x2-trait)^2 ), (rMax*exp( -(trait-opt_x)^2 )/(rMax*exp( -(x2-opt_x)^2 ))), 
        xlab="", ylab="", xlim = c(minX,maxX), ylim = c(minY,maxY),
        type = "l", lty = "dashed", col = "black", 
        lwd = 2, xaxs="i", yaxs="i"
      )
  points(  1-betaMax1*exp( -(x2-x10)^2 ), (rMax*exp( -(x2-opt_x)^2 )/(rMax*exp( -(x2-opt_x)^2 ))), 
           xlim = c(minX,maxX), ylim = c(minY,maxY), col = "blue", pch = 19 
        ) 
  if( coex ){
    
    points( 1-betaMax1*exp( -(x2-xeq1)^2 ), (rMax*exp( -(xeq1-opt_x)^2 )/(rMax*exp( -(x2-opt_x)^2 ))), 
            xlim = c(minX,maxX), ylim = c(minY,maxY), col = "black", pch = 19 
          ) 
  } 
  if( res_ext ){
    
    points( 1-betaMax1*exp( -(x2-xeq2)^2 ), (rMax*exp( -(xeq2-opt_x)^2 )/(rMax*exp( -(x2-opt_x)^2 ))), 
            xlim = c(minX,maxX), ylim = c(minY,maxY), col = "black", pch = 19 
          ) 
  }
  
}
### continuous-time 2-cosumer-2-resource model 
sim_traj03 = function( xeq1, xeq2, betaMax1, coex, res_ext, minX, maxX, minY, maxY ){
  
  csum <- 0.25; c11 <- seq( xeq2, xeq1, length=1000 ); c12 <- csum - c11
  c2 <- c(0.1,0.15)
  e1 <- c(0.06,0.1); e2 <- c(0.1,0.06)
  m1 <- 6*(c11 - 0.1)^2 + 0.47; m2 <- 0.5
  r <- c(10,10); K <- c(100,100)
  
  B1 <- e1[1]*c11*K[1]  
  B2 <- e1[2]*c12*K[2]  
  B3 <- e2[1]*c2[1]*K[1]  
  B4 <- e2[2]*c2[2]*K[2]  
  
  alpha11 <- (B1*c11/r[1] + B2*c12/r[2])/(B1+B2-m1)
  alpha12 <- (B1*c2[1]/r[1] + B2*c2[2]/r[2])/(B1+B2-m1)
  alpha21 <- (B3*c11/r[1] + B4*c12/r[2])/(B3+B4-m2)
  alpha22 <- (B3*c2[1]/r[1] + B4*c2[2]/r[2])/(B3+B4-m2)
  
  nd <- 1 - sqrt( (alpha12*alpha21)/(alpha11*alpha22) )
  fd <- sqrt( (alpha21*alpha22)/(alpha11*alpha12) )
  
  plot( nd, fd, 
        xlab="", ylab="", xlim = c(minX,maxX), ylim = c(minY,maxY),
        type = "l", lty = "dashed", col = "black", 
        lwd = 2, xaxs="i", yaxs="i"
      )
  
  ### initial state
  c11 <- 0.15; c12 <- csum - c11
  m1 <- 6*(c11 - 0.1)^2 + 0.47
  
  B1 <- e1[1]*c11*K[1]  
  B2 <- e1[2]*c12*K[2]  
  B3 <- e2[1]*c2[1]*K[1]  
  B4 <- e2[2]*c2[2]*K[2]  
  
  alpha11 <- (B1*c11/r[1] + B2*c12/r[2])/(B1+B2-m1)
  alpha12 <- (B1*c2[1]/r[1] + B2*c2[2]/r[2])/(B1+B2-m1)
  alpha21 <- (B3*c11/r[1] + B4*c12/r[2])/(B3+B4-m2)
  alpha22 <- (B3*c2[1]/r[1] + B4*c2[2]/r[2])/(B3+B4-m2)
  
  nd <- 1 - sqrt( (alpha12*alpha21)/(alpha11*alpha22) )
  fd <- sqrt( (alpha21*alpha22)/(alpha11*alpha12) )
  
  points( nd, fd, xlim = c(minX,maxX), ylim = c(minY,maxY), col = "blue", pch = 19 ) 
  
  if( coex ){
    
    c11 <- xeq1; c12 <- csum - c11
    m1 <- 6*(c11 - 0.1)^2 + 0.47
    
    B1 <- e1[1]*c11*K[1]  
    B2 <- e1[2]*c12*K[2]  
    B3 <- e2[1]*c2[1]*K[1]  
    B4 <- e2[2]*c2[2]*K[2]
    
    alpha11 <- (B1*c11/r[1] + B2*c12/r[2])/(B1+B2-m1)
    alpha12 <- (B1*c2[1]/r[1] + B2*c2[2]/r[2])/(B1+B2-m1)
    alpha21 <- (B3*c11/r[1] + B4*c12/r[2])/(B3+B4-m2)
    alpha22 <- (B3*c2[1]/r[1] + B4*c2[2]/r[2])/(B3+B4-m2)
    
    nd <- 1 - sqrt( (alpha12*alpha21)/(alpha11*alpha22) )
    fd <- sqrt( (alpha21*alpha22)/(alpha11*alpha12) )
    
    points( nd, fd, 
            xlim = c(minX,maxX), ylim = c(minY,maxY), col = "black", pch = 19 
          ) 
  } 
  if( res_ext ){
    
    c11 <- xeq2; c12 <- csum - c11
    m1 <- 6*(c11 - 0.1)^2 + 0.47
    
    B1 <- e1[1]*c11*K[1]  
    B2 <- e1[2]*c12*K[2]  
    B3 <- e2[1]*c2[1]*K[1]  
    B4 <- e2[2]*c2[2]*K[2]
    
    alpha11 <- (B1*c11/r[1] + B2*c12/r[2])/(B1+B2-m1)
    alpha12 <- (B1*c2[1]/r[1] + B2*c2[2]/r[2])/(B1+B2-m1)
    alpha21 <- (B3*c11/r[1] + B4*c12/r[2])/(B3+B4-m2)
    alpha22 <- (B3*c2[1]/r[1] + B4*c2[2]/r[2])/(B3+B4-m2)
    
    nd <- 1 - sqrt( (alpha12*alpha21)/(alpha11*alpha22) )
    fd <- sqrt( (alpha21*alpha22)/(alpha11*alpha12) )
    
    points( nd, fd, 
            xlim = c(minX,maxX), ylim = c(minY,maxY), col = "black", pch = 19 
          ) 
    
  }
  
}

niche_fitnessDiff = function(){

  reps <- 1000
  minX <- -1
  maxX <- 1.
  minY <- 0.01
  maxY <- 1.2
  
  mx <- (maxX-minX)*reps+1
  my <- (maxY-minY)*reps+1
  x <- seq( minX, maxX, length = mx )
  y <- seq( minY, maxY, length = my )
  z1 <- matrix( 0, length(x), length(y) )
  z2 <- matrix( 0, length(x), length(y) )
  # population density null-cline
  for(i in 1:mx){
    for(j in 1:my){
      z1[i,j] <- 1 - x[i] - y[j]
    }
  }
  par(ps=15)
  par(cex=1.5)
  contour( x, y, z1, 
           xlab = expression(paste("Niche difference, ", 1-rho)), 
           ylab = expression(paste("Competitive ability difference,", f[1]/f[2])), 
           drawlabels = F, levels=0, 
           col = "black", lwd = 3, cex=1, xaxs="i", yaxs="i"
          )
  
  for(i in 1:mx){
    for(j in 1:my){
      z2[i,j] <- y[j] - 1/(1-x[i])
    }
  }
  par(new=T)
  contour( x, y, z2, 
           xlab = "", ylab = "", 
           drawlabels = F, levels=0, 
           col = "black", lwd = 3, xaxs="i", yaxs="i"
          )
  
  
  ###=== coexistence or extinction region ===###
  int <- 0

  poly1 <- seq( minX, int, length=1000 )
  poly2 <- seq( int, maxX, length=1000 )
  poly3 <- seq( 1-1/minY, int, length=1000 )
  poly4 <- seq( int, 1-1/maxY, length=1000 )
  poly5 <- seq( 1-1/maxY, maxX, length=1000 )
  
  px <- c( poly1, poly4, rev(c(poly1,poly4)) )
  py <- c( 1-poly1, 1/(1-poly4), rep( maxY, length=length(rev(c(poly1,poly4))) ) )
  polygon(px,py,col="pink")
  
  px <- c( poly2, rev(poly5), rev(poly4) )
  py <- c( 1-poly2, rep(maxY,length=length(poly5)), 1/(1-rev(poly4)) )
  polygon(px,py,col="palegreen")
  
  px <- c( poly3, poly2, rev(c(poly3,poly2)) )
  py <- c( 1/(1-poly3), 1-poly2, rep(minY,length=length(rev(c(poly3,poly2)))) )
  polygon(px,py,col="cyan")
  
  ###=== simulation trajectory ===###
  par(new=T)
  sim_traj02( 0.81, 0, 1.1, TRUE, TRUE, minX, maxX, minY, maxY )
  par(new=T)
  sim_traj02( 2.2, 0, 1.5, FALSE, TRUE, minX, maxX, minY, maxY )
  par(new=T)
  sim_traj02( 0.15, 0, 0.9, TRUE, FALSE, minX, maxX, minY, maxY )
  

}

###=== graphic ===###
fname <- "~/morita-yamamichi2024/result/figureS14"
pdf( paste( fname, "pdf", sep = "." ), width = 8, height = 7 )
par(mar=c(5, 7, 1, 3)) # margin
par(oma = c(1, 0.1, 1, 0.5)) # outer margin

niche_fitnessDiff()

dev.off()
print( "supplementary figure S14 completed!" )