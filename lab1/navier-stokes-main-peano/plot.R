plotea_snap <- function(dd,main=main,zlim=zlim)
{
   ndim= sqrt(length(dd))
   dd2=matrix(dd,nrow=ndim)
   x = 1:ndim
   y = 1:ndim
   image(x,y,dd2, col=topo.colors(60),main=main,zlim=zlim)
}
plota <-function()
{
     for(i in seq(0,2044,by=1))
     {
       ff=sprintf("results/results-headless-peano-%04d-sorted.dat", i)
       d1 <- read.table(ff)
       ff=sprintf("frames/headless-peano-%04d.png", i)
       png(ff, width = 800, height = 800)
       plotea_snap(d1$V3,main="Peano-Hilbert order",zlim=c(0,3))
       dev.off()
     
       ff=sprintf("results/results-headless-%04d.dat", i)
       d2 <- read.table(ff)
       ff=sprintf("frames/headless-original-%04d.png", i)
       png(ff, width = 800, height = 800)
       plotea_snap(d2$V3,main="Matrix order",zlim=c(0,3))
       dev.off();
     
     
       #readline(prompt="Press [enter] to continue")
       
     }
}

