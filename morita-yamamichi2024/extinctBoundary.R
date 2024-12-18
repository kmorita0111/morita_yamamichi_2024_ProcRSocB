extinctBoundary = function( iMin, iMax, jMin, jMax, reps01, reps02 ){
  
  f = function( x ){
    
    ###--- population dynamics ---###
    frac011 = 2*exp(x^2 + x2^2 - 4*x2*xopt + 2*xopt^2 - 4*x*(x2 + xopt))
    frac012 = betaMax1*lambdaMax2*exp(2*x2*(x + xopt))
    frac013 = alpha1*exp(2*(x^2 + x2*xopt))*lambdaMax2*(x - x2) + betaMax2*exp(2*x*(x2 + xopt))*lambdaMax1*(x - xopt)
    frac014 = alpha2*lambdaMax1*exp(2*(x2^2 + x*xopt))
    frac015 = -betaMax2*exp(2*x*(x2 + xopt))*lambdaMax1*(x - x2) - alpha1*exp(2*(x^2 + x2*xopt))*lambdaMax2*(x - xopt)
    
    frac01 = frac011*( frac012*frac013 + frac014*frac015 )
    
    
    frac021 = alpha2*lambdaMax1*exp((x-x2)^2 + (x2-xopt)^2)
    frac022 = betaMax2*lambdaMax1*exp((x2 - xopt)^2)
    frac023 = lambdaMax2*exp((x-xopt)^2)*(-betaMax1+alpha1*exp((x-x2)^2))
    
    frac02 = ( frac021 - frac022 + frac023 )^2
    
    
    ###--- trait dynamics ---###
    frac03 = alpha1*betaMax1*exp((x - x2)^2)*(x2 - xopt)*( 1 + 2*x^2 + 2*x2*xopt - 2*x*(x2 + xopt) )
    frac04 = ( alpha1*exp((x - x2)^2)*(x - xopt) + betaMax1*(-x2 + xopt) )^2
    
    
    return( frac01/frac02 - frac03/frac04 )
  }
  
  g = function( x ){
    
    ###--- population dynamics ---###
    frac01 = alpha2*lambdaMax1*exp((x-x2)^2 + (x2-xopt)^2) - betaMax1*lambdaMax2*exp((x - xopt)^2)
    
    frac021 = alpha2*lambdaMax1*exp((x-x2)^2 + (x2-xopt)^2) - betaMax2*lambdaMax1*exp((x2 - xopt)^2)
    frac022 = exp((x-xopt)^2)*( -betaMax1 + alpha1*exp((x - x2)^2) )*lambdaMax2
    frac02 = frac021 + frac022
    
    ###--- trait dynamics ---###
    frac03 = betaMax1*(x2 - xopt)
    frac04 = -alpha1*exp((x-x2)^2)*(x - xopt) + betaMax1*(x2 - xopt)
    
    return( frac01/frac02 - frac03/frac04 )
  }
  
  ###--- loading library ---###
  library(rootSolve)
  
  ###=== main ===###
  extBound <- matrix(NA,reps01,reps02)
  
  x2 <- -0.1; xopt <- 0 
  alpha1 <- 1; alpha2 <- 1
  lambdaMax2 <- 100
  
  for ( j in 1:reps02 ) { # 1:reps
    
    # print( paste( "j =", as.character(j) ) )
    
    lambdaMax1 <- lambdaMax2*( jMin + (jMax-jMin)*(j-1)/(reps02-1) )
    
    count <- 1
    
    for ( i in 1:reps01 ) { # 1:reps
      
      betaMax1 <- iMin + (iMax-iMin)*(i-1)/(reps01-1)
      betaMax2 <- iMin + (iMax-iMin)*(i-1)/(reps01-1)
      
      if( (lambdaMax1/lambdaMax2 - 0.5454*betaMax1) < 0.4445 ){ break } # prevent to pick up unintended point 
      
      res01 <- uniroot.all( f, lower = 0, upper = 2 ) 
      res02 <- uniroot.all( g, lower = 0, upper = 2 ) 
      
      len_res01 <- length( res01 )
      len_res02 <- length( res02 )
      
      if( len_res01 > 0 ){ # 0 2
        if( len_res02 > 0 ){ # 0 2
          
          for ( num_res01 in 1:len_res01 ) { # 1:len_res01, c(2,3)
            for ( num_res02 in 1:len_res02 ) { # 1:len_res02, c(2,3)
              
              if( 0 < res01[num_res01] && res01[num_res01] < 1 ){ 
                if( 0 < res02[num_res02] && res02[num_res02] < 1 ){ 
                  
                  if( abs(res01[num_res01] - res02[num_res02]) < 0.01 ){ 
                    extBound[i,j] <- 1
                    
                    print( "===" )
                    print( paste("bataMax1", as.character(betaMax1)) )
                    print( paste("lambdaMax1", as.character(lambdaMax1)) )
                    print( paste("-> x =", as.character(res01[num_res01])) )
                    
                  }
                  
                }
              }
              
            }
          }
          
        }
      }
      
    }
  }
  
  return( extBound )  
}