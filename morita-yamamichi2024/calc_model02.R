populDynam = function( N, r, alpha, beta, h ){
  
  k1 <- c(NA, NA, NA, NA)
  k2 <- c(NA, NA, NA, NA)
  
  k1[1] <- N[1]*(r[1] - alpha[1]*N[1] - beta[1]*N[2])*h
  k2[1] <- N[2]*(r[2] - alpha[2]*N[2] - beta[2]*N[1])*h
  
  k1[2] <- (N[1]+k1[1]/2)*(r[1] - alpha[1]*(N[1]+k1[1]/2) - beta[1]*(N[2]+k2[1]/2))*h
  k2[2] <- (N[2]+k2[1]/2)*(r[2] - alpha[2]*(N[2]+k2[1]/2) - beta[2]*(N[1]+k1[1]/2))*h
  
  k1[3] <- (N[1]+k1[2]/2)*(r[1] - alpha[1]*(N[1]+k1[2]/2) - beta[1]*(N[2]+k2[2]/2))*h
  k2[3] <- (N[2]+k2[2]/2)*(r[2] - alpha[2]*(N[2]+k2[2]/2) - beta[2]*(N[1]+k1[2]/2))*h
  
  k1[4] <- (N[1]+k1[3])*(r[1] - alpha[1]*(N[1]+k1[3]) - beta[1]*(N[2]+k2[3]))*h
  k2[4] <- (N[2]+k2[3])*(r[2] - alpha[2]*(N[2]+k2[3]) - beta[2]*(N[1]+k1[3]))*h
  
  dN1 <- (k1[1]+2*k1[2]+2*k1[3]+k1[4])/6
  dN2 <- (k2[1]+2*k2[2]+2*k2[3]+k2[4])/6
  
  return( c((N[1]+dN1), (N[2]+dN2)) )
}
traitDiver = function( N, x, r, alpha, beta, betaMax1, betaMax2, here, sigma, h ){
  
  k1 <- c(NA, NA, NA, NA)
  k2 <- c(0, 0, 0, 0)
  
  k1[1] <- -2*betaMax1*exp(-(-x[1]+x[2])^2)*(-x[1]+x[2])*N[2]*h
  #k2[1] <- -2*betaMax2*exp(-(x[1]-x[2])^2)*(x[1]-x[2])*N[1]*h
  
  k1[2] <- -2*betaMax1*exp(-(-(x[1]+k1[1]/2)+(x[2]+k2[1]/2))^2)*(-(x[1]+k1[1]/2)+(x[2]+k2[1]/2))*N[2]*h
  #k2[2] <- -2*betaMax2*exp(-((x[1]+k1[1]/2)-(x[2]+k2[1]/2))^2)*((x[1]+k1[1]/2)-(x[2]+k2[1]/2))*N[1]*h
  
  k1[3] <- -2*betaMax1*exp(-(-(x[1]+k1[2]/2)+(x[2]+k2[2]/2))^2)*(-(x[1]+k1[2]/2)+(x[2]+k2[2]/2))*N[2]*h
  #k2[3] <- -2*betaMax2*exp(-((x[1]+k1[2]/2)-(x[2]+k2[2]/2))^2)*((x[1]+k1[2]/2)-(x[2]+k2[2]/2))*N[1]*h
  
  k1[4] <- -2*betaMax1*exp(-(-(x[1]+k1[3])+(x[2]+k2[3]))^2)*(-(x[1]+k1[3])+(x[2]+k2[3]))*N[2]*h
  #k2[4] <- -2*betaMax2*exp(-((x[1]+k1[3])-(x[2]+k2[3]))^2)*((x[1]+k1[3])-(x[2]+k2[3]))*N[1]*h
  
  dx1 <- (k1[1]+2*k1[2]+2*k1[3]+k1[4])/6
  dx2 <- (k2[1]+2*k2[2]+2*k2[3]+k2[4])/6
  
  return( c((x[1]+here*sigma*dx1), (x[2]+here*sigma*dx2)) )
}

ecd_ad = function( N, x, opt_x, r, rMax1, rMax2, alpha, beta, betaMax1, betaMax2, here, sigma, dt ){
  
  kn1 <- c(0, 0, 0, 0)
  kn2 <- c(0, 0, 0, 0)
  kx1 <- c(0, 0, 0, 0)
  kx2 <- c(0, 0, 0, 0)
  
  ### k1 ###
  kn1[1] <- N[1]*(r[1] - alpha[1]*N[1] - beta[1]*N[2])*dt
  kn2[1] <- N[2]*(r[2] - alpha[2]*N[2] - beta[2]*N[1])*dt
  kx1[1] <- (2*exp(-(opt_x-x[1])^2)*rMax1*(opt_x-x[1]) + 2*betaMax1*exp(-(x[1]-x[2])^2)*(x[1]-x[2])*N[2])*dt
  #kx2[1] <- (2*exp(-(opt_x-x[2])^2)*rMax2*(opt_x-x[2]) + 2*betaMax2*exp(-(x[1]-x[2])^2)*(x[2]-x[1])*N[1])*dt
  
  ### k2 ###
  r[1] <- rMax1*exp( -((x[1]+kx1[1]/2)-opt_x)^2 )
  r[2] <- rMax2*exp( -(x[2]-opt_x)^2 )
  beta[1] <- betaMax1*exp( -((x[1]+kx1[1]/2)-x[2])^2 )
  beta[2] <- betaMax2*exp( -((x[1]+kx1[1]/2)-x[2])^2 )
  
  kn1[2] <- (N[1]+kn1[1]/2)*(r[1] - alpha[1]*(N[1]+kn1[1]/2) - beta[1]*(N[2]+kn2[1]/2))*dt
  kn2[2] <- (N[2]+kn2[1]/2)*(r[2] - alpha[2]*(N[2]+kn2[1]/2) - beta[2]*(N[1]+kn1[1]/2))*dt
  kx1[2] <- (2*exp(-(opt_x-(x[1]+kx1[1]/2))^2)*rMax1*(opt_x-(x[1]+kx1[1]/2)) + 2*betaMax1*exp(-((x[1]+kx1[1]/2)-(x[2]+kx2[1]/2))^2)*((x[1]+kx1[1]/2)-(x[2]+kx2[1]/2))*(N[2]+kn2[1]/2))*dt
  #kx2[2] <- (2*exp(-(opt_x-(x[2]+kx2[1]/2))^2)*rMax2*(opt_x-(x[2]+kx2[1]/2)) + 2*betaMax2*exp(-((x[1]+kx1[1]/2)-(x[2]+kx2[1]/2))^2)*((x[2]+kx2[1]/2)-(x[1]+kx1[1]/2))*N[1])*dt

  ### k3 ###
  r[1] <- rMax1*exp( -((x[1]+kx1[2]/2)-opt_x)^2 )
  r[2] <- rMax2*exp( -(x[2]-opt_x)^2 )
  beta[1] <- betaMax1*exp( -((x[1]+kx1[2]/2)-x[2])^2 )
  beta[2] <- betaMax2*exp( -((x[1]+kx1[2]/2)-x[2])^2 )
  
  kn1[3] <- (N[1]+kn1[2]/2)*(r[1] - alpha[1]*(N[1]+kn1[2]/2) - beta[1]*(N[2]+kn2[2]/2))*dt
  kn2[3] <- (N[2]+kn2[2]/2)*(r[2] - alpha[2]*(N[2]+kn2[2]/2) - beta[2]*(N[1]+kn1[2]/2))*dt
  kx1[3] <- (2*exp(-(opt_x-(x[1]+kx1[2]/2))^2)*rMax1*(opt_x-(x[1]+kx1[2]/2)) + 2*betaMax1*exp(-((x[1]+kx1[2]/2)-(x[2]+kx2[2]/2))^2)*((x[1]+kx1[2]/2)-(x[2]+kx2[2]/2))*(N[2]+kn2[2]/2))*dt
  #kx2[3] <- (2*exp(-(opt_x-(x[2]+kx2[2]/2))^2)*rMax2*(opt_x-(x[2]+kx2[2]/2)) + 2*betaMax2*exp(-((x[1]+kx1[2]/2)-(x[2]+kx2[2]/2))^2)*((x[2]+kx2[2]/2)-(x[1]+kx1[2]/2))*N[1])*dt

  ### k4 ###
  r[1] <- rMax1*exp( -((x[1]+kx1[3])-opt_x)^2 )
  r[2] <- rMax2*exp( -(x[2]-opt_x)^2 )
  beta[1] <- betaMax1*exp( -((x[1]+kx1[3])-x[2])^2 )
  beta[2] <- betaMax2*exp( -((x[1]+kx1[3])-x[2])^2 )
  
  kn1[4] <- (N[1]+kn1[3])*(r[1] - alpha[1]*(N[1]+kn1[3]) - beta[1]*(N[2]+kn2[3]))*dt
  kn2[4] <- (N[2]+kn2[3])*(r[2] - alpha[2]*(N[2]+kn2[3]) - beta[2]*(N[1]+kn1[3]))*dt
  kx1[4] <- (2*exp(-(opt_x-(x[1]+kx1[3]))^2)*rMax1*(opt_x-(x[1]+kx1[3])) + 2*betaMax1*exp(-((x[1]+kx1[3])-(x[2]+kx2[3]))^2)*((x[1]+kx1[3])-(x[2]+kx2[3]))*(N[2]+kn2[3]))*dt
  #kx2[4] <- (2*exp(-(opt_x-(x[2]+kx2[3]))^2)*rMax2*(opt_x-(x[2]+kx2[3])) + 2*betaMax2*exp(-((x[1]+kx1[3])-(x[2]+kx2[3]))^2)*((x[2]+kx2[3])-(x[1]+kx1[3]))*N[1])*dt
  
  ### weighted average ###
  dN1 <- (kn1[1]+2*kn1[2]+2*kn1[3]+kn1[4])/6
  dN2 <- (kn2[1]+2*kn2[2]+2*kn2[3]+kn2[4])/6
  dx1 <- (kx1[1]+2*kx1[2]+2*kx1[3]+kx1[4])/6
  dx2 <- (kx2[1]+2*kx2[2]+2*kx2[3]+kx2[4])/6
  
  
  return( c( (N[1]+dN1), (N[2]+dN2), (x[1]+here*sigma*dx1), (x[2]+here*sigma*dx2) ) )
}
