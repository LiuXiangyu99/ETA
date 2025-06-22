ETAmeas<-function(x,y,t=100)
{
  n <- length(y)
  xcon<-x[order(y,decreasing=T)]   ### concomitant of decreasing Y-o.s.
  xrk<-n+1-rank(xcon)              ### reverse-rank of xcon: e.g., max has rank 1
  
  s<-min(t,n-1)                    ### number of ETA computed
  TailETA<-array(s)               ### Initialize ETA 1:s
  for (k in 1:s)
  {
    pxy=0
    qxy=numeric()
    for (i in 1:k)
    {
      pxy=pxy+sum(pmax(k+1-pmax(xrk[i],xrk[1:k]),0))  ### compute sum(k-max(R_i, R_(1:k))
    }
    TailETA[k]<-3*pxy/k^3
  }  
  return(TailETA)
}



PlotETA<-function(x,y,t=200, cap="Tail ETA(X|Y) and ETA(Y|X)",ybot=0,ytop=0)
{
  ts<-min(100,t)   # starting index
  ex<-1          # ylimit addition on graph
  xy<-ETAmeas(x,y,t)
  yx<-ETAmeas(y,x,t)
  
  valmax<-max(max(xy[ts:t],yx[ts:t]),ytop)
  valmin<-max(min(xy[ts:t],yx[ts:t]),ybot)
  
  plot(ts:t,xy[ts:t],main=cap,col="red",ylab = "value",xlim=c(ts,t),ylim=c(max(valmin-ex,0),min(valmax+ex,1)),cex.main=1.2,xlab="order statistics",type="l",lwd=1.5)
  lines(ts:t,yx[ts:t], col="blue",lty=1,lwd=1.5)
  lines(ts:t,rep(0,length(ts:t)), col="darkgreen",lty=2,lwd=1.5)
}


PlotDelta<-function(x,y,t=200, cap="Tail Delta(X,Y)")
{
  ts<-min(100,t)   # starting index
  xy<-ETAmeas(x,y,t)
  yx<-ETAmeas(y,x,t)
  
  plot(ts:t,xy[ts:t]-yx[ts:t],main=cap,col="red",ylab = "value",xlim=c(ts,t),ylim=c(-0.2,0.2),cex.main=1.2,xlab="order statistics",type="l",lwd=1.5)
  lines(ts:t,rep(0,length(ts:t)), col="blue",lty=2,lwd=1.5)
}


ETAboot1meas<-function(x,y,m,t=100)
{
  n <- length(y)
    ## need to rewrite fully
  xcon<-x[order(y,decreasing=T)]   ### concomitant of decreasing Y-o.s.
  xrk<-n+1-rank(xcon)              ### reverse-rank of xcon: e.g., max has rank 1
  
  s<-min(t,n-1)                    ### number of ETA computed
  TailETA<-array(s)                 ### Initialize ETA 1:s
  for (k in 1:s)
  {
    pxy=0
    qxy=numeric()
    for (i in 1:k)
    {
      pxy=pxy+sum(pmax(k+1-pmax(xrk[i],xrk[1:k]),0))  ### compute sum(k-max(R_i, R_(1:k))
    }
    TailETA[k]<-3*pxy/k^3
  }  
  ma<-m/m[n]
  rx<-n+1-rank(x)
  Rmx<-(cumsum(ma[order(x,decreasing=T)])-ma[order(x,decreasing=T)])[rx]
  ry<-n+1-rank(y)
  Rmy<-(cumsum(ma[order(y,decreasing=T)])-ma[order(y,decreasing=T)])[ry]
  
  TailETAboot<-array(s) 
  for (k in 1:s)
  {
    pxy=0
    qxy=numeric()
    for (i in 1:k)
    {
      pxy=pxy+sum(pmax(k+1-pmax(xrk[i],xrk[1:k]),0))  ### compute sum(k-max(R_i, R_(1:k))
    }
    TailETA[k]<-3*pxy/k^3
  }  
  
  return(TailETA)
}

