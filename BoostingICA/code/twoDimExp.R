#'twoDimExp' function implements the two-dimensional ICA experiments in past researches[1][2].
# N: sample size
# p: the dimension of the observed mixtures
# m: number of repetitions
# maxit: the maximum iterations for ICA methods
# rseed: random seed
# reference: [1]F. Bach and M. Jordan, “Kernel independent component analysis,” Journal of Machine Learning Research, vol. 3, pp. 1–48, 03 2003.
#            [2]T. Hastie and R. Tibshirani, “Independent components analysis throughproduct density estimation,” in Advances in Neural Information Pro-
#             cessing Systems, S. Becker, S. Thrun, and K. Obermayer, Eds., vol. 15.MIT Press, 2003, pp. 665–672.
# cite: YunPeng Li, ZhaoHui Ye, "Boosting Independent Component Analysis", in IEEE Signal Processing Letters
#       YunPeng Li, "Second-order Approximation of Minimum Discrimination Information in Independent Component Analysis", in IEEE Signal Processing Letters, vol. 29, pp. 334-338, 2022, doi: 10.1109/LSP.2021.3135193.

twoDimExp<-function (N=1024,p=2,m=20,maxit=20,rseed=12345,...){
  set.seed(rseed)
  k<-18
  n<-5
  dists<-c("a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r")
  amaris<-array(0,dim=c(n,k,m))     # Amari metrics
  SIR<-array(0,dim=c(n,k,m))        # Single -to-Interference-Ratio (SIR)
  timecost<-array(0,dim=c(n,k,m))   # CPU elapsed time
  
  row_names <- c("F-G0","F-G1","PICA","B-SP","MICA")
  row.names(amaris) <- row_names
  row.names(SIR) <- row_names
  row.names(timecost) <- row_names
  
  rls <- list(
    rl1= function(x,W0) {
      time_start=proc.time()
      W = ProDenICA(x=x, W0=W0,Gfunc=G0,trace=FALSE,maxit=maxit)$W
      time_end =proc.time()
      rl =list(W=W,time=(time_end-time_start)[3][1])
      rl
    },
    rl2= function(x,W0) {
      time_start=proc.time()
      W = ProDenICA(x=x, W0=W0,Gfunc=G1,trace=FALSE,maxit=maxit)$W
      time_end =proc.time()
      rl =list(W=W,time=(time_end-time_start)[3][1])
      rl
    },
    rl3= function(x,W0) {
      time_start=proc.time()
      W = ProDenICA(x=x, W0=W0,Gfunc=GPois,trace=FALSE,maxit=maxit)$W
      time_end =proc.time()
      rl =list(W=W,time=(time_end-time_start)[3][1])
      rl
    },
    rl4= function(x,W0) {
      time_start=proc.time()
      W = ProDenICA(x=x, W0=W0,Gfunc=NPMLE_BSP,trace=FALSE,maxit=maxit)$W
      time_end =proc.time()
      rl =list(W=W,time=(time_end-time_start)[3][1])
      rl
    },
   
    rl5= function(x,W0) {
      time_start=proc.time()
      W = ProDenICA(x=x, W0=W0,Gfunc=Gmdi,trace=FALSE,maxit=maxit)$W
      time_end =proc.time()
      rl =list(W=W,time=(time_end-time_start)[3][1])
      rl
    }
  )
  
  
  
  for (i in 1:k){
    
    for (j in 1:m){
      print(paste("Distribution",i,"m",j))
      dist =dists[i]
      A0<-mixmat(p)
      s<-scale(cbind(rjordan(dist,N),rjordan(dist,N)))
      
      x <- s %*% A0
 
      #Whiten the data
      x <- scale(x, TRUE, FALSE)
      #print("##")
      sx <- svd(x)	### orthogonalization function
      x <- sqrt(N) * sx$u
      target <- solve(A0)
      target <- diag(sx$d) %*% t(sx$v) %*% target/sqrt(N)
      W0 <- matrix(rnorm(p*p), p, p)
      W0 <- ICAorthW(W0)
      
      time_start = Sys.time();
      cur_rls <- lapply(rls, function(f) f(x=x,W0=W0))
      time_end =Sys.time();
      print(paste("total time consumed: ",(time_end-time_start)))
      
      amaris[1,i,j]<-amari(cur_rls$rl1$W,target)
      amaris[2,i,j]<-amari(cur_rls$rl2$W,target)
      amaris[3,i,j]<-amari(cur_rls$rl3$W,target)
      amaris[4,i,j]<-amari(cur_rls$rl4$W,target)
      amaris[5,i,j]<-amari(cur_rls$rl5$W,target)
     
      
      timecost[1,i,j]<-cur_rls$rl1$time
      timecost[2,i,j]<-cur_rls$rl2$time
      timecost[3,i,j]<-cur_rls$rl3$time
      timecost[4,i,j]<-cur_rls$rl4$time
      timecost[5,i,j]<-cur_rls$rl5$time
    
      
      print(paste("add time consumed: ",(sum(timecost[,i,j]))))
      
      SIR[1,i,j]<-SIR(s,x%*%cur_rls$rl1$W)
      SIR[2,i,j]<-SIR(s,x%*%cur_rls$rl2$W)
      SIR[3,i,j]<-SIR(s,x%*%cur_rls$rl3$W)
      SIR[4,i,j]<-SIR(s,x%*%cur_rls$rl4$W)
      SIR[5,i,j]<-SIR(s,x%*%cur_rls$rl5$W)
      
      
    }
  }
  k <-18
  amaris <- 100.*amaris
  amaris.mean =apply(amaris,c(1,2),mean)
  amaris.sd   =apply(amaris,c(1,2),sd)
  
  SIR.mean =apply(SIR,c(1,2),mean)
  SIR.sd   =apply(SIR,c(1,2),sd)
  
  timecost.mean =apply(timecost,c(1,2),mean)
  timecost.sd   =apply(timecost,c(1,2),sd)
  
  
  
 
  
  
  showtext_auto(enable = TRUE)
  font_add('SimSun', 'simsun.ttc')
  
  y_min = 0.9*min(amaris.mean)
  y_max = 1.1*max(amaris.mean)
  
  plot(1:k,amaris.mean[1,],type="b",lty=1,pch=1,col="red",xaxt="n",xlab="Distribution",ylab="Amari Distance from True W",ylim=c(y_min,y_max))
  lines(1:k,amaris.mean[2,],type="b",lty=1,,pch=2,col="blue")
  lines(1:k,amaris.mean[3,],type="b",lty=1,,pch=3,col="green")
  lines(1:k,amaris.mean[4,],type="b",lty=1,,pch=4,col="pink")
  lines(1:k,amaris.mean[5,],type="b",lty=1,,pch=5,col="black")
  
  axis(1, at = 1:k, labels = letters[1:18])
  legend("topleft", inset=.05,lty=c(1,2,3,4,5),lwd=c(2,2,2,2,2),legend=c("F-G0","F-G1","PICA","B-SP","MICA"),pch=c(1,1,1,1,1,1),
         col=c("red", "blue","green","pink","black"))
  
  
  y_min = 0.9*min(SIR.mean)
  y_max = 1.1*max(SIR.mean)
  
  plot(1:k,SIR.mean[1,],type="b",lty=1,pch=1,col="red",xaxt="n",xlab="Distribution",ylab="SIR",ylim=c(y_min,y_max))
  lines(1:k,SIR.mean[2,],type="b",lty=1,,pch=2,col="blue")
  lines(1:k,SIR.mean[3,],type="b",lty=1,,pch=3,col="green")
  lines(1:k,SIR.mean[4,],type="b",lty=1,,pch=4,col="pink")
  lines(1:k,SIR.mean[5,],type="b",lty=1,,pch=5,col="black")
  
  axis(1, at = 1:k, labels = letters[1:18])
  legend("topleft", inset=.05,lty=c(1,2,3,4,5),lwd=c(2,2,2,2,2),legend=c("F-G0","F-G1","PICA","B-SP","MICA"),pch=c(1,1,1,1,1,1),
         col=c("red", "blue","green","pink","black"))
  
  
  y_min = 0.9*min(timecost.mean)
  y_max = 1.1*max(timecost.mean)
  
  plot(1:k,timecost.mean[1,],type="b",lty=1,pch=1,col="red",xaxt="n",xlab="Distribution",ylab="CPU elapsed time",ylim=c(y_min,y_max))
  lines(1:k,timecost.mean[2,],type="b",lty=1,,pch=2,col="blue")
  lines(1:k,timecost.mean[3,],type="b",lty=1,,pch=3,col="green")
  lines(1:k,timecost.mean[4,],type="b",lty=1,,pch=4,col="pink")
  lines(1:k,timecost.mean[5,],type="b",lty=1,,pch=5,col="black")
  
  axis(1, at = 1:k, labels = letters[1:18])
  legend("topleft", inset=.05,lty=c(1,2,3,4,5),lwd=c(2,2,2,2,2),legend=c("F-G0","F-G1","PICA","B-SP","MICA"),pch=c(1,1,1,1,1,1),
         col=c("red", "blue","green","pink","black"))
  
  rl<-list(amaris=amaris,time=timecost,SIR=SIR)
  rl
  #print(apply(amaris.mean, c(1,2),mean))
} 