#'Gmdi' function returns information concerning the given random samples' minimum discrimination information (with respect to the standard Gaussian distribution).
# x: the observed univariate random samples
# L: the size of grid 
# rl$Gs: is the value of f(x) 
# rl$gs: is the first-order derivatives of f(x)
# rl$gps: is the second-order derivatives of f(x)
# cite: YunPeng Li, "Second-order Approximation of Minimum Discrimination Information in Independent Component Analysis", in IEEE Signal Processing Letters, vol. 29, pp. 334-338, 2022, doi: 10.1109/LSP.2021.3135193.

Gmdi<-function (x, L = 500, order = 1, widen = 1.2, 
                ...) 
{
  mu0<-mean(x)
  sd0<-sd(x)
  x <- drop(scale(x))
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
  ys <- as.vector(table(cut(x, xcuts)))/n
  
  eps =1e-9
  gxg <- dnorm(xg)
  weights<-0.5*gaps*gxg
  zs<-((ys-gaps*gxg)/(gaps*gxg+eps))
  
  kernel_fun <- function(x1){
    m <- length(x1)
    n <- 2
    xnew0 <- array(0,dim=c(m,2))
    xnew1 <- array(0,dim = c(m,2))
    xnew2 <- array(0,dim = c(m,2))
    
    A <- exp(-0.5*x1^2)
    B <- x1*exp(-0.5*x1^2)
    C <- (x1^2)*exp(-0.5*x1^2)
    D <- (x1^3)*exp(-0.5*x1^2)
    
    
    
    xnew0[,1] <- B 
    xnew0[,2] <- A
    
    xnew1[,1] <- A - C
    xnew1[,2] <- -1.0*B
      
    xnew2[,1] <- -3.0*B + D
    xnew2[,2] <- C - A
    
    rl <- list(x=xnew0,der1=xnew1,der2=xnew2)
    rl
  }
  
  K1 <- kernel_fun(xg)
  
  W <- diag(weights)
  H <- K1$x
  
  beta = solve(t(H)%*%H,t(H)%*%t(W)%*%scale(zs,scale = FALSE))  
 
  
  K2 <- kernel_fun(x)
  
  Gs <- mean(zs) + K2$x%*%beta
  gs <- K2$der1%*%beta
  gps <- K2$der2%*%beta
   
  rl = list(Gs = Gs, gs = gs, gps = gps)
  
  rl
}
