###=== parameters setting ===###
cii0 <- 0.20
wii0 <- 0.06
m10 <- 6*(cii0-0.1)^2 + 0.47
m20 <- 0.5

b <- c(1,1)
c1sum <- 0.25; c2sum <- 0.25
c1 <- c(cii0, c1sum - cii0)
c2 <- c(0.10, 0.15)
w1 <- c(wii0,0.1) 
w2 <- c(0.1,wii0)
m <- c(m10,m20)
r <- c(10,10)
K <- c(100,100)

here = 0.01 # additive genetic variance

t <- 2000; h <- 0.001 # simulation duration and time interval
t.inv01 <- 13; t.inv02 <- 15  # immigration timing

interv <- seq(0,t,by=h)

path = "~/morita-yamamichi2024"

###=== initial values ===###
C10 <- 59.5745; C20 <- 0
R10 <- 40.4255; R20 <- 10.6383

nChange01 <- matrix(NA,length(interv),4)
nChange01[1,] <- c(C10,0,R10,R20)
nChange02 <- matrix(NA,length(interv),4)
nChange02[1,] <- c(C10,0,R10,R20)

c1Change01 <- matrix(NA,length(interv),2)
c1Change01[1,] <- c1
c2Change01 <- matrix(NA,length(interv),2)
c2Change01[1,] <- c2
c1Change02 <- matrix(NA,length(interv),2)
c1Change02[1,] <- c1
c2Change02 <- matrix(NA,length(interv),2)
c2Change02[1,] <- c2

m1Change01 <- matrix(NA,length(interv),2)
m1Change01[1,] <- m
m1Change02 <- matrix(NA,length(interv),2)
m1Change02[1,] <- m

###=== simulation ===###
source( paste(path, "calc_S12.R", sep = "/") )

###--- additive genetic variance ---###
print(nChange01[1,])
print(nChange02[1,])
print(c1Change01[1,])
print(c2Change01[1,])
print(c1Change02[1,])
print(c2Change02[1,])
print(m1Change01[1,])
print(m1Change02[1,])
print(b)

for ( timei in 1:(length(interv)-1) ) {
  
  ###=== invasion timing ===###
  if( timei == t.inv01/h ){ nChange01[timei,2] <- 0.1*nChange01[timei,1] }
  if( timei == t.inv02/h ){ nChange02[timei,2] <- 0.1*nChange02[timei,1] }
  
  ###=== 2 consumer & 2 resource population dynamics ===###
  nChange01[timei+1,] <- con2_res2( nChange01[timei,c(1,2)], nChange01[timei,c(3,4)], b, c1Change01[timei,], c2Change01[timei,], w1, w2, m1Change01[timei,], r, K, h )
  nChange02[timei+1,] <- con2_res2( nChange02[timei,c(1,2)], nChange02[timei,c(3,4)], b, c1Change02[timei,], c2Change02[timei,], w1, w2, m1Change02[timei,], r, K, h )
  
  ###=== consumption rate evolutionary dynamics ===###
  #- species 1 -#
  c1Change01[timei+1,] <- traitDynam03( nChange01[timei,c(1,2)], nChange01[timei,c(3,4)], b, c1Change01[timei,], c1sum, c2Change01[timei,], c2sum, w1, w2, m1Change01[timei,], n=1, r, K, here, sigma=1, h )
  c1[1] <- c1Change01[timei+1,1]
  c1[2] <- (c1sum - c1[1])
  c1Change01[timei+1,2] <- c1[2]
  c2Change01[timei+1,] <- c2Change01[timei,]
  #- species 2 -#
  c1Change02[timei+1,] <- traitDynam03( nChange02[timei,c(1,2)], nChange02[timei,c(3,4)], b, c1Change02[timei,], c1sum, c2Change02[timei,], c2sum, w1, w2, m1Change02[timei,], n=1, r, K, here, sigma=1, h )
  c1[1] <- c1Change02[timei+1,1]
  c1[2] <- (c1sum - c1[1])
  c1Change02[timei+1,2] <- c1[2]
  c2Change02[timei+1,] <- c2Change02[timei,]
  
  ###=== trade-off between c11 & m1 ===###
  #- species 1 -#
  m[1] <- 6*(c1Change01[timei+1,1]-0.1)^2 + 0.47; m[2] <- 0.5
  m1Change01[timei+1,] <- m
  
  #- species 2 -#
  m[1] <- 6*(c1Change02[timei+1,1]-0.1)^2 + 0.47; m[2] <- 0.5
  m1Change02[timei+1,] <- m
  
}
print(nChange01[timei+1,])
print(nChange02[timei+1,])

print( "Simulation completed!" )
###=== graphics ===###
source( paste(path, "graphic_S12.R", sep = "/") )

#--- population (figureS10a) ---#
interv <- seq(0,200,h); timei <- 200/h + 1

fname=paste(path, "result/figureS10a", sep = "/")
pdf(paste( fname, "pdf", sep = "." ), width=8, height=6 )
par(ps=15)
par(mar=c(5.5, 6.0, 4, 2)) # margin

graphicPopulDynam( nChange01[1:length(interv),], t.inv01, interv, interv[timei], leg=c(200,120), addL=TRUE )

dev.off()

#--- population (figureS10c) ---#
#--- update interval ---###
interv <- seq(0,2000,by=h); timei <- (length(interv)-1)

fname=paste(path, "result/figureS10c", sep = "/")
pdf(paste( fname, "pdf", sep = "." ), width=8, height=6 )
par(ps=15)
par(mar=c(5.5, 6.0, 4, 2)) # margin

graphicPopulDynam( nChange01[1:length(interv),], t.inv01, interv, interv[timei], leg=c(200,120), addL=FALSE )

dev.off()

#--- trait dynamics (figureS10e)  ---#
fname=paste(path, "result/figureS10e", sep = "/")
pdf(paste( fname, "pdf", sep = "." ), width=8, height=6 )
par(ps=15)
par(mar=c(5.5, 6.0, 4, 2)) # margin

graphicTraitDynam( c1Change01[1:length(interv),], c2Change01[1:length(interv),], t.inv01, interv, interv[timei], leg=c(0,0), addL=FALSE )

dev.off()

#--- population (figureS10b) ---#
interv <- seq(0,200,h); timei <- 200/h + 1

fname=paste(path, "result/figureS10b", sep = "/")
pdf(paste( fname, "pdf", sep = "." ), width=8, height=6 )
par(ps=15)
par(mar=c(5.5, 6.0, 4, 2)) # margin

graphicPopulDynam( nChange02[1:length(interv),], t.inv02, interv, interv[timei], leg=c(200,120), addL=FALSE )

dev.off()

#--- population (figureS10d) ---#
interv <- seq(0,500,h); timei <- 500/h + 1

fname=paste(path, "result/figureS10d", sep = "/")
pdf(paste( fname, "pdf", sep = "." ), width=8, height=6 )
par(ps=15)
par(mar=c(5.5, 6.0, 4, 2)) # margin

graphicPopulDynam( nChange02[1:length(interv),], t.inv02, interv, interv[timei], leg=c(200,120), addL=FALSE )

dev.off()

#--- trait dynamics (figureS10f)  ---#
fname=paste(path, "result/figureS10f", sep = "/")
pdf(paste( fname, "pdf", sep = "." ), width=8, height=6 )
par(ps=15)
par(mar=c(5.5, 6.0, 4, 2)) # margin

graphicTraitDynam( c1Change02[1:length(interv),], c2Change02[1:length(interv),], t.inv02, interv, interv[timei], leg=c(0,0), addL=FALSE )

dev.off()

print( "Supplementary figure S12 completed!" )