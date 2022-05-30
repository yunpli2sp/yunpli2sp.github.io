#'NPMLE_BGK' function returns information concerning the given random samples' density function 'p(x)', where Gaussian kernel works as the weak learner.
# x: the observed univariate random samples
# df: degree of freedom, which controls the model complexity of single weak learners
# L: the size of grid 
# M: the number of boosting iterations
# mu: learning rate for single weak learner
# density.return: if its value is 'True', the probability density function is returned (rl$density)
# bw: bandwidth for Gaussian kernel 
# rl$Gs: is the value of ln p(x)
# rl$gs: is the first-order derivatives of ln p(x)
# rl$gps: is the second-order derivatives of ln p(x)
# cite:  YunPeng Li, ZhaoHui Ye, "Boosting in Univariate Nonparametric Maximum Likelihood Estimation", in IEEE Signal Processing Letters, vol. 28, pp. 623-627, 2021, doi: 10.1109/LSP.2021.3065881.


NPMLE_BGK<-function (x, df=3,bw=0.5,L = 500, order = 1, widen = 1.2,M=5, mu=1,density.return = TRUE,
                   ...){
  
  
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
  
  # calculate the Lagrange multiplier \lambda with given degree of freedom df
  caldf <- function(X,W,df=3,lambda = 0,maxit=1000){
    Z <- t(X)%*%W%*%X
    I <- diag(ncol(X))
    it <- 0
    d <- svd(Z)$d
    
    tol <- 1e-14
    grad <- 99
    hessian <- 99
    while(it<maxit){
      it <- it+1
      grad <- sum(d/(d+lambda))-df
      #if(it==1) print(paste("lambda =", lambda,"df =", grad+df));
      hessian <- -1.0*sum(d/(d+lambda)^2)
      if(grad<1e-9) break
      lambda <- lambda - grad/hessian
    }
    lambda
  }
  

  
  
  #Gaussian kernel function
  gausskernel <- function(x1,x2,h){
    m <- length(x1)
    n <- length(x2)
    xnew0 <- array(0,dim=c(m,n))
    xnew1 <- array(0,dim = c(m,n))
    xnew2 <- array(0,dim = c(m,n))
    for(i in 1:m){
      for(j in 1:n){
        xnew0[i,j] <- exp(-0.5*((x1[i]-x2[j])/h)^2)
        xnew1[i,j] <- xnew0[i,j]*(-1*((x1[i]-x2[j])/(h^2)))
        xnew2[i,j] <- xnew0[i,j]*(((x1[i]-x2[j])^2/h^4)-1.0/h^2)
      }
    }
    
    rl <- list(x=xnew0,der1=xnew1,der2=xnew2)
    rl
  }
  
  
  
  
  
  n <- length(x)
  

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
  nx <-length(xg)
  yg <- array(as.vector(table(cut(x, xcuts)))/n,dim = c(nx))
  
  h = bw;
 
  
  xgk <- gausskernel(x1=xg,x2=xg,h=h)$x
  xrl <- gausskernel(x1=x,x2=xg,h=h)
  x0 <- xrl$x
  derx1 <- xrl$der1
  derx2 <- xrl$der2
  
  
  

  
  eps =1e-14
  tol =1e-9
  
  
  zg<-array(gaps,dim=c(nx))
  
  
  
  
  
  wg<-array(0,dim=c(nx))
  fg<-array(0,dim=c(nx))
  rg<-array(0,dim=c(nx))  
  tg<-array(0,dim=c(nx))  
  xg <- as.array(xg)
  
  wg<-array(1.0/nx,dim=c(nx))
  

  
  rg<-(yg-wg)/(wg+eps)
  
  #lambda <- caldf(X=xgk,W=diag(wg),df=df)
  lambda <- 10000 # to avoid the calculation of specific df, we fix lambda as extremely large value
  
  glm.fit <- glmnet(xgk, rg, family="gaussian", alpha=0,weights = wg,lambda = lambda)
  tg <-array(predict(glm.fit, xgk,type ="link"),dim = c(nx))
  
 
  fg<-fg+tg
  betas <- coef(glm.fit)
  intercept <- betas[1]
  betas <- betas[2:(nx+1)]

  
  Gs <-  mu*(intercept +x0%*%betas)
  gs <-  mu*(derx1%*%betas)
  gps <- mu*(derx2%*%betas)
  

  
  
  
  
  for(i in 2:M){
    if(i>M)
      break
 
    
    wg<-wg*exp(tg)
    rg<-(yg-wg)/(wg+eps)
 
    glm.fit <- glmnet(xgk, rg, family="gaussian",weights = wg,alpha=0,lambda = lambda)
    tg <-array(predict(glm.fit, xgk,type ="link"),dim = c(nx))

    fg<-fg+tg
    
    betas <- coef(glm.fit)
    intercept <- betas[1]
    betas <- betas[2:(nx+1)]
    Gs <- Gs + mu*(intercept +x0%*%betas)
    gs <- gs + mu*(derx1%*%betas)
    gps <- gps + mu*(derx2%*%betas)
    
    
  }
  
  
  
  partition <- gaps*sum(exp(fg))*sd0
 
  if (density.return) {
    density = list(x = xg*sd0+mu0, y = exp(fg))
    density$y <-density$y/partition
  }
  
  
  
  
  
  rl = list(Gs = Gs/partition, gs = gs/partition, gps = gps/partition)
  if (density.return) 
    rl$density = density
  rl
  
}