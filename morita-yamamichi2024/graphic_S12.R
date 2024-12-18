graphicPopulDynam <- function( Density, t.inv, interv, t, leg, addL ){
  
  x <- Density[,1]
  y <- Density[,2]
  r1 <- Density[,3]
  r2 <- Density[,4]
  
  m <- max(max(x),max(y),max(r1),max(r2))
  
  plot( interv, rep(0,length(interv)), 
        xlim = c(0,t), ylim = c(0,m),
        lwd = 3, type = "l", lty = "dashed", xlab = "", ylab = "", col = "black"
      )
  par(new=T)
  plot( interv, r1, 
        xlim = c(0,t), ylim = c(0,m),
        lwd = 3, type = "l", xlab = "", ylab = "", col = "green3"
      )
  par(new=T)
  plot( interv, r2, 
        xlim = c(0,t), ylim = c(0,m),
        lwd = 3, type = "l", xlab = "", ylab = "", col = "darkorange"
      )
  par(new=T)
  plot( interv, x, xlim = c(0,t), ylim = c(0,m), 
        xlab = "Time", ylab = "Species density", 
        lwd = 3, type = "l", col = "red"
      ) 
  par(new=T)
  plot( interv, c(rep(NA,length=t.inv/h), y[(t.inv/h+1):length(interv)]), 
        xlim = c(0,t), ylim = c(0,m),
        lwd = 3, type = "l", xlab = "", ylab = "", col = "blue"
      )

  if( addL == TRUE ){
      legend( "topright",
              legend = c("Consumer 1","Consumer 2", "Resource 1","Resource 2" ),
              col = c("red","blue","green3","darkorange"),
              lwd = c(2,2,2,2), lty = c(1,1,1,1)
      )
  }
}

graphicTraitDynam <- function( trait1, trait2, t.inv, interv, t, leg, addL ){
  
  x1 <- trait1[,1]
  x2 <- trait1[,2]
  y1 <- trait2[,1]
  y2 <- trait2[,2]
  
  m <- max(max(max(x1),max(y1)),max(max(x2),max(y2)))
  yMin <- min(min(min(x1),min(y1)), min(min(x2),min(y2)))
  
  plot( interv, c(rep(NA,length=t.inv/h), y1[(t.inv/h+1):length(interv)]),
        xlim = c(0,t), ylim = c(yMin,m),
        lwd = 3, lty="dashed", type = "l", xlab = "", ylab = "", col = "blue"
  )
  par(new=T)
  plot( interv, x2, xlim = c(0,t), ylim = c(yMin,m), 
        xlab = "Time", ylab = expression(paste("Consumption rate, ", c[ij])), 
        lwd = 3, lty="dashed", type = "l", col = "red"
  )
  par(new=T)
  plot( interv, c(rep(NA,length=t.inv/h), y2[(t.inv/h+1):length(interv)]),
        xlim = c(0,t), ylim = c(yMin,m),
        lwd = 3, type = "l", xlab = "", ylab = "", col = "blue"
  )
  par(new=T)
  plot( interv, x1, xlim = c(0,t), ylim = c(yMin,m), 
        xlab = "", ylab = "", 
        lwd = 3, type = "l", col = "red"
  )
  if( addL == TRUE ){
      legend( "topright",
              legend = c("c11","c22","c12","c21"),
              col = c("red","blue","red","blue"), lwd = c(2,2,2,2), lty = c(1,1,2,2)
            )
  }
}