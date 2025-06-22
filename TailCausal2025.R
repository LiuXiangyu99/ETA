q99<-function(x)
{
  return(quantile(x,0.99))
}

q95<-function(x)
{
  return(quantile(x,0.95))
}

q05<-function(x)
{
  return(quantile(x,0.05))
}

q01<-function(x)
{
  return(quantile(x,0.01))
}





ETAmeas<-function(x,y,t=100)
{
  n <- length(y)
  xcon<-x[order(y,decreasing=T)]   ### concomitant of decreasing Y-o.s.
  xrk<-n+1-rank(xcon)              ### reverse-rank of xcon: e.g., max has rank 1
  
  s<-min(t,n-1)                    ### number of ETA computed
  TailETA<-array(s)                ### Initialize ETA 1:s
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



PlotETA<-function(x,y,tstart=100,tend=200, cap="Tail ETA(X|Y) and ETA(Y|X)",ybot=0,ytop=1)
{
  ts<-min(tstart,tend)   # starting index
  ex<-1          # ylimit addition on graph
  xy<-ETAmeas(x,y,tend)
  yx<-ETAmeas(y,x,tend)
  
  valmax<-max(max(xy[ts:tend],yx[ts:tend]),ytop)
  valmin<-max(min(xy[ts:tend],yx[ts:tend]),ybot)
  
  plot(ts:tend,xy[ts:tend],main=cap,col="red",ylab = expression(eta), xlim=c(ts,tend),ylim=c(max(valmin-ex,0),min(valmax+ex,1)),cex.main=1.2,xlab="order statistics",type="l",lwd=1.5)
  lines(ts:tend,yx[ts:tend], col="blue",lty=1,lwd=1.5)
  lines(ts:tend,rep(0,length(ts:tend)), col="darkgreen",lty=2,lwd=1.5)
}


PlotDelta<-function(x,y, tstart = 100, tend = 200, cap="Tail Delta(X,Y)",lt=0.5) {
  ts <- min(tstart, tend)  # starting index
  xy <- ETAmeas(x, y, tend)
  yx <- ETAmeas(y, x, tend)
  a <- abs(max(xy - yx))   # max absolute difference (not currently used)

# Plot the difference between ETA(X|Y) and ETA(Y|X)
  plot(ts:tend, xy[ts:tend] - yx[ts:tend],
     main = cap,
     col = "red",
     ylab = expression(Delta),
     xlim = c(ts, tend),
     ylim = c(-lt, lt),
     cex.main = 1.2,
     xlab = "order statistics",
     type = "l",
     lwd = 1.5)

# Add horizontal line at zero
  lines(ts:tend, rep(0, length(ts:tend)), col = "blue", lty = 2, lwd = 1.5)
}








ETAmeas_boot<-function(x,y,t=100,B=100,gap=10)
{
  n<- length(x)                    ### length of x and y
  sindex<-seq(gap,min(t,n-1),gap) ### indices of ETA computed with as given gap
  s<-length(sindex)               ### number of ETA computed with as given gap
  
  xcon<-x[order(y,decreasing=T)]   ### x-concomitant of decreasing Y-o.s.
  xconmat<-matrix(rep(xcon,n),n,n,byrow=T) ### repeat x-concomitant n times
  xcomp<-matrix(as.integer(xconmat>t(xconmat)),n,n) ### aiding computation of rank of x-concomitant
  
  ETAboot<-matrix(0,B,s)
  
  for (b in 1:B)
  {
  w<-rexp(n)                        ### generate n exponential random variables
  a<-w/(mean(w))                    ### compute exp/mean(exp)
  acon<-a[order(y,decreasing=T)]    ### a-concomitant of decreasing Y-o.s.
  axrk<-apply(acon*xcomp,1,sum)    ### rank of bootstrapped version of x-concomitant multiplied with a_i
  aconsumtot<-cumsum(acon)             ### needed to find the amount upto which to sum (upto t)
  temp<-floor(aconsumtot)
  tau<-array(0,s)
  for (i in sindex)
  {tau[i]<-sum(temp<i)}
  
  l<-1
  for (k in sindex)
  {
    pxy=0
    for (i in 1:tau[k])
    {
      pxy=pxy+sum(acon[i]*acon[1:tau[k]]*pmax(k-pmax(axrk[i],axrk[1:tau[k]],0)))  ### compute suma[i]a[j](k-max(aR_i, aR_(1:k))
    }
    ETAboot[b,l]<-3*pxy/k^3
    l<-l+1
  }
  }
  return(ETAboot)
  
}





DELTAmeas_boot<-function(x,y,tstart=100,tend=400,B=100,gap=10)
{
  n<- length(x)                    ### length of x and y
  sindex<-seq(tstart,min(tend,n-1),gap) ### indices of ETA computed with as given gap
  s<-length(sindex)               ### number of ETA computed with as given gap
  
  xcon<-x[order(y,decreasing=T)]   ### x-concomitant of decreasing Y-o.s.
  xconmat<-matrix(rep(xcon,n),n,n,byrow=T) ### repeat x-concomitant n times
  xcomp<-matrix(as.integer(xconmat>t(xconmat)),n,n) ### aiding computation of rank of x-concomitant
  
  ycon<-y[order(x,decreasing=T)]   ### y-concomitant of decreasing X-o.s.
  yconmat<-matrix(rep(ycon,n),n,n,byrow=T) ### repeat y-concomitant n times
  ycomp<-matrix(as.integer(yconmat>t(yconmat)),n,n) ### aiding computation of rank of y-concomitant
  
  DELTAboot<-matrix(0,B,s)
  ETAbooty<-matrix(0,B,s)
  ETAbootx<-matrix(0,B,s)
  for (b in 1:B)
  {
    w<-rexp(n)                        ### generate n exponential random variables
    a<-w/(mean(w))                    ### compute exp/mean(exp)
    
    acon_y<-a[order(y,decreasing=T)]    ### a-concomitant of decreasing Y-o.s.
    axrk_y<-apply(acon_y*xcomp,1,sum)    ### rank of bootstrapped version of x-concomitant multiplied with a_i
    acon_y_sumtot<-cumsum(acon_y)             ### needed to find the amount up to which to sum (upto t)
    temp_y<-floor(acon_y_sumtot)
    tau_y<-array(0,s)
    for (i in sindex)
    {tau_y[i]<-sum(temp_y<i)}
    
    acon_x<-a[order(x,decreasing=T)]    ### a-concomitant of decreasing X-o.s.
    ayrk_x<-apply(acon_x*ycomp,1,sum)    ### rank of bootstrapped version of y-concomitant multiplied with a_i
    acon_x_sumtot<-cumsum(acon_x)             ### needed to find the amount up to which to sum (upto t)
    temp_x<-floor(acon_x_sumtot)
    tau_x<-array(0,s)
    for (i in sindex)
    {tau_x[i]<-sum(temp_x<i)}
    
    l<-1
    for (k in sindex)
    { 
      {
      pxy=0
      for(i in 1:tau_y[k])
      {
        pxy=pxy+sum(acon_y[i]*acon_y[1:tau_y[k]]*pmax(k-pmax(axrk_y[i],axrk_y[1:tau_y[k]]),0))  ### compute suma[i]a[j](k-max(aR_i, aR_(1:k))
      }
      ETAbooty[b,l]<-3*pxy/k^3
      }
      {
      pyx=0
      for (j in 1:tau_x[k])
      {
        pyx=pyx+sum(acon_x[i]*acon_x[1:tau_x[k]]*pmax(k-pmax(ayrk_x[i],ayrk_x[1:tau_x[k]]),0))  ### compute suma[i]a[j](k-max(aR_i, aR_(1:k))
      }
      ETAbootx[b,l]<-3*pyx/k^3
      }
      l<-l+1
    }
  }
  DELTAboot<- ETAbooty-ETAbootx
  
  return(list(ETAbooty,ETAbootx,DELTAboot))
  
}