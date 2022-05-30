#'NPMLE_BSP' function returns information concerning the given random samples' density function 'p(x)', where smooth spline works as the weak learner.
# x: the observed univariate random samples
# df: degree of freedom, which controls the model complexity of single weak learners
# L: the size of grid 
# M: the number of boosting iterations
# mu: learning rate for single weak learner
# density.return: if its value is 'True', the probability density function is returned (rl$density)
# rl$Gs: is the value of ln p(x)
# rl$gs: is the first-order derivatives of ln p(x)
# rl$gps: is the second-order derivatives of ln p(x)
# cite:  YunPeng Li, ZhaoHui Ye, "Boosting in Univariate Nonparametric Maximum Likelihood Estimation", in IEEE Signal Processing Letters, vol. 28, pp. 623-627, 2021, doi: 10.1109/LSP.2021.3065881.


NPMLE_BSP<-function (x, df = 3, L = 500, order = 1, widen = 1.2,M=5, mu=1,density.return = TRUE, 
                     ...) 
{ 
  ylim.scale <- function (ylim, scale = 0) 
  {
    scale2 <- diff(ylim)
    if (scale2 < scale) 
      rep(mean(ylim), 2) + ((ylim - mean(ylim)) * scale)/scale2
    else ylim
  }
  
  mu0<-mean(x)
  sd0<-sd(x)
  x <- drop(scale(x))
  n <- length(x)
  mu=1.0
  
  rangex <- range(x)
  if (order == 1) 
    rx <- rangex
  else {
    rx <- sort(x)[c(order, n - order + 1)]
  }
  rx <- ylim.scale(rx, diff(rx) * widen)
  xg <- seq(from = rx[1], to = rx[2], length = L)
  gaps <- diff(rx)/(L - 1)
  xcuts <- c(min(rangex[1], rx[1]) - gaps/2, xg[-L] + gaps/2, 
             max(rangex[2], rx[2]) + gaps/2)
  
  yg <- as.vector(table(cut(x, xcuts)))/n
  
  
  
  
  nx <-length(xg)
  eps =1e-14
  tol =1e-9
  
  
  zg<-array(gaps,dim=c(nx))
  
  
  
  
  wg<-array(0,dim=c(nx))
  fg<-array(0,dim=c(nx))
  rg<-array(0,dim=c(nx))
  tg<-array(0,dim=c(nx))
  
  
  
  
  
  
  wg<-array(1.0/nx,dim=c(nx))
  
  rg<-(yg-wg)/(wg+eps)
  
  
  
  mdi.fit <- smooth.spline(xg,rg, 0.5*wg, df)
  tg <-mu*predict(mdi.fit,xg,deriv=0)$y
  fg<-fg+tg
  
  Gs <- mu*predict(mdi.fit, x, deriv = 0)$y
  gs <- mu*predict(mdi.fit, x, deriv = 1)$y
  gps <- mu*predict(mdi.fit, x, deriv = 2)$y
  
  
  for(i in 2:M){
    if(i>M)
      break

    
    wg<-wg*exp(tg)
    rg<-(yg-wg)/(wg+eps)

    
    mdi.fit <- smooth.spline(xg, rg, 0.5*wg, df)
    tg <-mu*predict(mdi.fit,xg,deriv=0)$y
    fg<-fg+tg
    
    Gs <- Gs+mu*predict(mdi.fit, x, deriv = 0)$y
    gs <- gs+mu*predict(mdi.fit, x, deriv = 1)$y
    gps <- gps+mu*predict(mdi.fit, x, deriv = 2)$y
  }
  
  
  partition <- gaps*sum(exp(fg))*sd0
  if (density.return) {
    density = list(x = xg*sd0+mu0, y = exp(fg))
    density$y <- density$y/partition
  }
  
  
  rl = list(Gs = Gs/partition, gs = gs/partition, gps = gps/partition)
  
  if (density.return) 
    rl$density = density
  rl
  
}