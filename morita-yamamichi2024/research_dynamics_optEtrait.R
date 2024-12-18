################################################################################
# population dynamics
################################################################################
research_dynamics_optEtrait = function( b, bmax1, bmax2, y,
                                        alpha, alphaMax1, alphaMax2,
                                        beta, betaMax1, betaMax2, x, x1_h, x2_h,
                                        lambda, lambdaMax1, lambdaMax2,
                                        N,
                                        Gx1, Gx2, Gy1, rel_ii, rel_ij,
                                        Sl1, Sl2, Sx1, Sx2, Sy1, Sy2, 
                                        Density01, 
                                        lambdaChange01, alphaChange01, 
                                        betaChange01, bChange01, 
                                        xChange01, yChange01,
                                        xyDensChange01,
                                        t, sp, mod  
                                      ){
  
  # counting loop
  count <- 0

  #population dynamics
  calcN1 <- function( N, lambda, b, alpha, beta ){
    
    if( N[1] <= 0 ) return( 0 )
    else if( N[1] > 0 ){
      k <- N[1] / ( N[1] + b[1] * N[2] )
      l <- N[1] / ( alpha[1] * N[1] + beta[1] * N[2] ) # approximation 
      return( lambda[1] * ( k * l ) )
    }
  }
  calcN2 <- function( N, lambda, b, alpha, beta ){
    
    if( N[2] <= 0 ) return( 0 )
    else if( N[2] > 0 ){
      k <- N[2] / ( N[2] + b[2] * N[1] )
      l <- N[2] / ( alpha[2] * N[2] + beta[2] * N[1] ) # approximation 
      return( lambda[2] * ( k * l ) )
    }
  }
  # character displacement
  calcx1 = function( G, N, x, x1_h, lambda, alpha, beta, betaMax1, Sx1, Sl1 ){ # ecological CD
    
    k <- (1+alpha[1]*N[1])*Sx1*(x[1]-x1_h)
    l <- beta[1] * N[2] * (-Sl1*x[1] + Sl1*x[2] + Sx1*x[1] - Sx1*x1_h)
    m <- ( (alpha[1]*N[1]) + beta[1] * N[2] ) * Sl1 * Sx1 # approximation 
    return( -2*G*(k+l)/m )
  }
  calcx2 = function( G, N, x, x2_h, lambda, alpha, beta, betaMax2, Sx2, Sl2 ){ # ecological CD
    
    k <- -2 * (1+alpha[2]*N[2])*Sx2*(x[2]-x2_h)
    l <- 2 * beta[2] * N[1] * (Sl2*x[2] - Sl2*x[1] - Sx2*x[2] + Sx2*x2_h)
    m <- ( (1+alpha[2]*N[2]) + beta[2] * N[1] ) * Sl2 * Sx2 # modified 
    
    return( G * (k+l)/m )
  }
  calcy1 = function( G, N, y, b, S ){ # reproductive CD
    
    k <- 2.0 * N[2] * ( y[1] - y[2] ) / S
    e <- b[1]
    l <- (N[1] + b[1]*N[2])
    
    return( G * k * e / l )
  }
  calcy2 <- function( G, N, y2, alpha, beta, betaMax1, S ){ # reproductive CD
    
    k <- 2.0 * N[1] * ( y[2] - y[1] ) / S
    e <- b[2]
    l <- (N[2] + b[2]*N[1])
    
    return( G * k * e / l )
  }
  
  invN2 <- 0
  if( mod>3 ){
    invN2 <- N[2]; N[2] <- 0
  } 
    
  for( i in 1:t ){
    
    if( mod>3 ){
      if( sp >= 1 ){
        if( i == sp ){
          N[2] <- invN2
        }
      }
    } else {
      if( sp >= 1 ){
        if( i == sp ){
          N[2] <- N[1]/9 
        }
      }
    }

    # cahnging interspecific interactions
    alpha[1] <- alphaMax1 ## species 1 intraspecific competition with no character displacement
    alpha[2] <- alphaMax2 ## species 2 intraspecific competition with no character displacement
    beta[1] <- betaMax1 * exp( -1.0 * (x[1]-x[2])^2 /Sx1 ) ## species 1 interspecific competition 
    beta[2] <- betaMax2 * exp( -1.0 * (x[1]-x[2])^2 /Sx2 ) ## species 2 interspecific competition 
    b[1] <- bmax1 * exp( -1.0 * (y[1]-y[2])^2 /Sy1 ) ## species 1 reproductive interference
    b[2] <- bmax2 * exp( -1.0 * (y[1]-y[2])^2 /Sy2 ) ## species 2 reproductive interference
    lambda[1] <- lambdaMax1*exp( -1.0 * (x[1]-x1_h)^2 /Sl1 ) ## species 1 fecundity reduction
    lambda[2] <- lambdaMax2*exp( -1.0 * (x[2]-x2_h)^2 /Sl2 ) ## species 2 fecundity reduction
    
    alphaChange01[i,] <- alpha
    betaChange01[i,] <- beta
    bChange01[i,] <- b
    lambdaChange01[i,] <- lambda
    
    N1 <- calcN1( N, lambda, b, alpha, beta )
    N2 <- calcN2( N, lambda, b, alpha, beta )
    
    ###=== threshold ===###
    if( mod == 2 ){
      if( i > 100 ){
        if( N1 <= 0.0 ){ return( 0 ) }
        else if( N2 <= 0.0 ){ return( 1 ) }
      }
    }

    N[1] <- N1
    N[2] <- N2
    Density01[i+1,] <- N

    Wx1 <- calcx1(Gx1, N, x, x1_h, lambda, alpha, beta, betaMax1, Sx1, Sl1 )
    x[1] <- x[1] + rel_ii[1] * Wx1
    x[2] <- x[2]
    xChange01[i+1,] <- x
    yChange01[i+1,] <- y
    xyDensChange01[i+1,] <- c(x[1],y[1],N1/(N1+N2))
    
  }
  
 if( mod == 1 ){
    return( xyDensChange01 )
  } else if( mod == 2 ){
    sp1freq <- Density01[i+1,1]/(Density01[i+1,1]+Density01[i+1,2])
    return( sp1freq )
  } else if( mod == 3 ){
    return( cbind( Density01, xChange01 ) )
  } else if( mod == 4 ){
    sp1freq <- Density01[i+1,1]/(Density01[i+1,1]+Density01[i+1,2])
    return( sp1freq )
  } else if( mod == 5 ){
    return( cbind( Density01, xChange01 ) )
  }

}
################################################################################
# frequency dynamics
################################################################################
frequency_dynamics_optEtrait = function( b, bmax1, bmax2, y,
                                         alpha, alphaMax1, alphaMax2,
                                         beta, betaMax1, betaMax2, x, x1_h, x2_h,
                                         lambda, lambdaMax1, lambdaMax2,
                                         M,
                                         Gx1, Gx2, Gy1, rel_ii, rel_ij,
                                         Sl1, Sl2, Sx1, Sx2, Sy1, Sy2, 
                                         frequency01, 
                                         lambdaChange01, alphaChange01, 
                                         betaChange01, bChange01, 
                                         xChange01, yChange01,
                                         xyDensChange01,
                                         t, sp, mod  
){
  
  # counting loop
  count <- 0
  
  #frequency dynamics
  calcM1 <- function( M, lambda, b, alpha, beta ){
    
    if( M[1] <= 0 ) return( 0 )
    else if( N[1] > 0 ){
      k <- lambda[1]*M[1] / ( M[1] + beta[2]*M[2] )
      l <- lambda[1]*M[1] / ( M[1] + beta[2]*M[2] ) + lambda[2]*M[2] / ( beta[1]*M[1] + M[2] ) # approximation 
      return( k/l )
    }
  }
  # character displacement
  calcx1 = function( G, M, x, x1_h, lambda, alpha, beta, betaMax1, Sx1, Sl1 ){ # ecological CD
    
    k <- (M[1] + beta[2]*M[2])/lambda[1]
    l <- beta[1] * M[2] * lambda[1] * (x[2] - x[1]) / ( M[1] + beta[2]*M[2] )^2
    m <- lambda[1] * (x[1] - x1_h) / ( M[1] + beta[2]*M[2] )
    return( -2*G*k*l*m )
  }
  for( i in 1:t ){
    
    if( sp >= 1 ){
      if( i == sp ){
        M[2] <- M[1]/9 
      }
    }
    
    # cahnging interspecific interactions
    alpha[1] <- alphaMax1 ## species 1 intraspecific competition with no character displacement
    alpha[2] <- alphaMax2 ## species 2 intraspecific competition with no character displacement
    beta[1] <- betaMax1 * exp( -1.0 * (x[1]-x[2])^2 /Sx1 ) ## species 1 interspecific competition 
    beta[2] <- betaMax2 * exp( -1.0 * (x[1]-x[2])^2 /Sx2 ) ## species 2 interspecific competition 
    b[1] <- bmax1 * exp( -1.0 * (y[1]-y[2])^2 /Sy1 ) ## species 1 reproductive interference
    b[2] <- bmax2 * exp( -1.0 * (y[1]-y[2])^2 /Sy2 ) ## species 2 reproductive interference
    lambda[1] <- lambdaMax1*exp( -1.0 * (x[1]-x1_h)^2 /Sl1 ) ## species 1 fecundity reduction
    lambda[2] <- lambdaMax2*exp( -1.0 * (x[2]-x2_h)^2 /Sl2 ) ## species 2 fecundity reduction
    
    alphaChange01[i,] <- alpha
    betaChange01[i,] <- beta
    bChange01[i,] <- b
    lambdaChange01[i,] <- lambda
    
    M1 <- calcM1( M, lambda, b, alpha, beta )
    print(M1)

    M[1] <- M1; M[2] <- 1-M1
    frequency01[i+1,] <- M
    
    Wx1 <- calcx1(Gx1, M, x, x1_h, lambda, alpha, beta, betaMax1, Sx1, Sl1 )
    x[1] <- x[1] + rel_ii[1] * Wx1
    x[2] <- x[2]
    xChange01[i+1,] <- x
    yChange01[i+1,] <- y

  }
  
  cbind( frequency01, xChange01 )  
}