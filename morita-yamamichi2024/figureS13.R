###--- parameters ---###
csum <- 0.25; c2 <- c(0.1, 0.15)
wii <- 0.06; wij <- 0.1
K <- 100; m2 <- 0.5

xMin <- 0.1; xMax <- 0.23
yMin <- 0.4; yMax <- 0.6
###-- intersects ---###
intsct <- c( c2[1]*csum/(c2[1] + c2[2]), (c2[1]/c2[2])*csum/((wii/wij)^2+(c2[1]/c2[2])) )

ccc1 <- -c2[2]*csum*K*wii^2+2*csum*m2*wij-3*c2[1]*csum*K*wij^2-c2[2]*wii*yMin + c2[1]*wij*yMin
ccc2 <- -4*csum*(m2*(wii + wij) -  2*K*(c2[2]*wii^2 + c2[1]*wij^2))*(csum*wij*(m2 - c2[1]*K*wij) - c2[2]*wii*yMin)
ccc3 <- c2[2]*wii*(csum*K*wii + yMin) + wij*(-2*csum*m2 + 3*c2[1]*csum*K*wij - c2[1]*yMin)
ccc4 <- 2*m2*(wii + wij) - 4*K*(c2[2]*wii^2 + c2[1]*wij^2)
m1.yMin01 <- (ccc1 + sqrt( ccc2 + ccc3^2 ))/ccc4

ddd1 <- c2[2]*csum*K*wii^2 - 2*csum*m2*wij + 3*c2[1]*csum*K*wij^2 + c2[2]*wii*yMin - c2[1]*wij*yMin
ddd2 <- -4*csum*(m2*(wii + wij) - 2*K*(c2[2]*wii^2 + c2[1]*wij^2))*(csum*wij*(m2 - c2[1]*K*wij) - c2[2]*wii*yMin)
ddd3 <- c2[2]*wii*(csum*K*wii + yMin) + wij*(-2*csum*m2 + 3*c2[1]*csum*K*wij - c2[1]*yMin)
ddd4 <- 2*m2*(wii + wij) - 4*K*(c2[2]*wii^2 + c2[1]*wij^2)
m1.yMin02 <- -(ddd1 + sqrt( ddd2 + ddd3^2 ))/ddd4

eee1 <- c2[2]*csum*wij*(m2 - c2[1]*K*wij) + c2[1]^2*wij*(csum*K*wij - yMax) - c2[2]^2*wii*yMax
eee2 <- c2[1]^2*K*wij^2 + c2[2]*(-c2[2]*K*wii^2 + m2*wij) - c2[1]*(m2*wii + c2[2]*K*(-wii^2 + wij^2))
m2.yMax <- eee1/eee2

fff1 <- c2[2]*csum*wij*(m2 - c2[1]*K*wij) + c2[1]^2*wij*(csum*K*wij - yMin) - c2[2]^2*wii*yMin
fff2 <- c2[1]^2*K*wij^2 + c2[2]*(-c2[2]*K*wii^2 + m2*wij) - c2[1]*(m2*wii + c2[2]*K*(-wii^2 + wij^2))
m2.yMin <- fff1/fff2

###
###--- variables ---###
xLen <- 10000

###--- conditions 1 ---###
c11 <- seq(m1.yMin01, intsct[2],length=xLen)

aaa1 <- c11*K*wii
aaa2 <- (-c11+csum)*K*wij
aaa3 <- m2 - K*(c2[2]*wii+c2[1]*wij)
aaa4 <- -2*(csum/c11) + (csum/c11)^2 + (wii/wij+1)
aaa5 <- -1 + (csum/c11) + (c2[1]/c2[2])*(wij/wii)
aaa6 <- c11/c2[2]
aaa7 <- wij/wii

m11 <- aaa1 + aaa2 + aaa3*aaa4/aaa5*aaa6*aaa7

###--- conditions 2 ---###
c11 <- seq(intsct[2], m1.yMin02,length=xLen)

aaa1 <- c11*K*wii
aaa2 <- (-c11+csum)*K*wij
aaa3 <- m2 - K*(c2[2]*wii+c2[1]*wij)
aaa4 <- -2*(csum/c11) + (csum/c11)^2 + (wii/wij+1)
aaa5 <- -1 + (csum/c11) + (c2[1]/c2[2])*(wij/wii)
aaa6 <- c11/c2[2]
aaa7 <- wij/wii

m12 <- aaa1 + aaa2 + aaa3*aaa4/aaa5*aaa6*aaa7

###--- conditions 3 ---###
c11 <- seq(m2.yMin, intsct[2],length=xLen)

bbb1 <- c11*(K*wii- K*wij)
bbb2 <- csum*K*wij
bbb3 <- m2 - c2[2]*K*wii - c2[1]*K*wij
bbb4 <- (c2[1]/c2[2])*wii + (-1+(csum/c11))*wij
bbb5 <- wii + (c2[1]/c2[2])^2*wij
bbb6 <- c11/c2[2]

m21 <- bbb1 + bbb2 + bbb3*bbb4/bbb5*bbb6

###--- conditions 4 ---###
c11 <- seq(intsct[2], m2.yMax, length=xLen)

bbb1 <- c11*(K*wii- K*wij)
bbb2 <- csum*K*wij
bbb3 <- m2 - c2[2]*K*wii - c2[1]*K*wij
bbb4 <- (c2[1]/c2[2])*wii + (-1+(csum/c11))*wij
bbb5 <- wii + (c2[1]/c2[2])^2*wij
bbb6 <- c11/c2[2]

m22 <- bbb1 + bbb2 + bbb3*bbb4/bbb5*bbb6

###--- graphic ---###
###--- figure 1 ---###
fname="~/morita-yamamichi2024/result/figureS13"
pdf(paste( fname, "pdf", sep = "." ), width=7, height=7 )
par(ps=15)
par(mar=c(5.5, 6.0, 4, 2)) # margin

yMin <- 0.45; yMax <- 0.55
plot( c(0,0),c(0,0),
      xlab = expression(paste("Consumption rate, ", c[11] == c[sum] - c[12])), 
      ylab = expression(paste("Mortality rate, ", m[1])),
      cex.lab=1.3, 
      xlim = c(xMin,xMax), ylim = c(yMin,yMax), 
      type="n", xaxt="n",yaxt="n", xaxs="i", yaxs="i", 
    )

x1 <- seq(m1.yMin01, intsct[2], length=xLen)
x2 <- seq(intsct[2], m2.yMax, length=xLen)

xxx1 <- c( c(xMin,m1.yMin01), x1, x2, c(m2.yMax,intsct[2]), c(xMin,xMin) )
yyy1 <- c( c(yMin,yMin), m11, m22, c(yMax,yMax), c(yMax,yMin) )
polygon(xxx1, yyy1, col="deepskyblue")

x3 <- seq(intsct[2], m1.yMin02,length=xLen)

xxx2 <- c( x3, c(m1.yMin02,xMax), c(xMax,m2.yMax), rev(x2) )
yyy2 <- c( m12, c(yMin,yMin), c(yMax,yMax), rev(m22) )
polygon(xxx2, yyy2, col="green2")

x4 <- seq(m2.yMin, intsct[2], length=xLen)

xxx3 <- c( x4, x3, c(m1.yMin02,m2.yMin) )
yyy3 <- c( m21, m12, c(yMin,yMin)  )
polygon(xxx3, yyy3, col="red1")

par(new=T)
plot( seq(xMin,xMax, length=xLen), 0.47 + 6*(seq(xMin,xMax, length=xLen)-0.1)^2,
      xlim = c(xMin,xMax), ylim = c(yMin,yMax), 
      xlab = "", ylab = "",
      type = "l", lwd=3, xaxt="n",yaxt="n", xaxs="i", yaxs="i"
    )

axis(side=1, at=seq(xMin,xMax,by=0.01), labels=seq(xMin,xMax,by=0.01), cex.axis=1.3 ) # x-axis
axis(side=2, at=seq(yMin,yMax,by=0.01), labels=seq(yMin,yMax,by=0.01), cex.axis=1.3) # y-axis

x.fp <- 0.112719
points( x.fp, 0.47 + 6*(x.fp-0.1)^2, pch=21, cex=3, col="white", bg="black", lwd=1.5 )

x.fp <- 0.194749
points( x.fp, 0.47 + 6*(x.fp-0.1)^2, pch=21, cex=3, col="white", bg="black", lwd=1.5 )

x.fp <- 0.129926
points( x.fp, 0.47 + 6*(x.fp-0.1)^2, pch=21, cex=3, col="black", bg="white", lwd=1.5 )

dev.off()

print("Supplementary figure S13 completed!")