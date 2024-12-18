drawImage <- function( xlabname, ylabname, mainname, color_mat, xMax, xMin, yMax, yMin, axis_leg ){ #, ext.sp2
  
  par(ps=20)
  
  xAxis<-axis_leg[1]
  yAxis<-axis_leg[2]
  colnames <- paste0("",c(seq(xMin, xMax, by=((xMax-xMin)/xAxis) )) )
  rownames <- paste0("",c(seq(yMin, yMax, by=((yMax-yMin)/yAxis) )) )

  original <- colorRampPalette(c("blue", "black", "red"))
  image( color_mat, zlim=c(0,1), xaxt="n", yaxt="n", 
         xlab=xlabname, ylab=ylabname, main=mainname, 
         col=original(10)
       )

  axis(1,c(seq( 0, 1, by=1/xAxis )),labels=colnames)
  axis(2,c(seq( 0, 1, by=1/yAxis )),labels=rownames)
}
###--- main ---###
make_color_map <- function( color_mat, fname, xlabname, ylabname, x1, Sx1, Sx2, Sl1, Sl2, #N_var, 
                            xMax, xMin, yMax, yMin, betaMax1, betaMax2, axis_leg
                          ){
 
  mainname <- ""

  par(mar=c(5, 7, 2, 3)) # margin
  print( paste( "Max of Sp.1 frequency", as.character(max(color_mat)), sep = " = " ) )
  drawImage(xlabname, ylabname, mainname, color_mat, xMax, xMin, yMax, yMin, axis_leg)

}