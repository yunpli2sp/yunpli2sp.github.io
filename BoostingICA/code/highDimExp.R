#'twoDimExp' function implements the high dimensional ICA experiments in past researches[1][2].
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


highDimExp<-function (N=1024,p=2,m=100,maxit=20, rseed=12345,...){
  time_start2 = proc.time()
  set.seed(rseed)
  
  n<-5
  dists<-c("a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r")
  
  
  amaris<-array(0,dim=c(n,m))
  SIR<-array(0,dim=c(n,m))
  timecost<-array(0,dim=c(n,m))
  
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
  
  for (i in 1:m){
    time_start3 = proc.time()
    print(paste("iteration ",i))
    dist =sample(dists,p)
    A0<-mixmat(p)
    s<-rjordan(dist[1],N)
    for(d in 2:p){
      s<-cbind(s,rjordan(dist[d],N))
    }
    s <- scale(s)
    x <- s %*% A0
    #print("dim(x)")
    #print(dim(x))
    ###Whiten the data
    x <- scale(x, TRUE, FALSE)
    #print("##")
    sx <- svd(x)	### orthogonalization function
    x <- sqrt(N) * sx$u
    target <- solve(A0)
    target <- diag(sx$d) %*% t(sx$v) %*% target/sqrt(N)
    W0 <- matrix(rnorm(p*p), p, p)
    W0 <- ICAorthW(W0)
    time_end3 =proc.time()
    print(paste("prepare time consumed: ",(time_end3-time_start3)[3][1]))
    
    time_start = proc.time()
    cur_rls <- lapply(rls, function(f) f(x=x,W0=W0))
    time_end =proc.time()
    print(paste("time consumed: ",(time_end-time_start)[3][1]))
    
    
    amaris[1,i]<-amari(cur_rls$rl1$W,target)
    amaris[2,i]<-amari(cur_rls$rl2$W,target)
    amaris[3,i]<-amari(cur_rls$rl3$W,target)
    amaris[4,i]<-amari(cur_rls$rl4$W,target)
    amaris[5,i]<-amari(cur_rls$rl5$W,target)

    
    timecost[1,i]<-cur_rls$rl1$time
    timecost[2,i]<-cur_rls$rl2$time
    timecost[3,i]<-cur_rls$rl3$time
    timecost[4,i]<-cur_rls$rl4$time
    timecost[5,i]<-cur_rls$rl5$time

    

    
    SIR[1,i]<-SIR(s,x%*%cur_rls$rl1$W)
    SIR[2,i]<-SIR(s,x%*%cur_rls$rl2$W)
    SIR[3,i]<-SIR(s,x%*%cur_rls$rl3$W)
    SIR[4,i]<-SIR(s,x%*%cur_rls$rl4$W)
    SIR[5,i]<-SIR(s,x%*%cur_rls$rl5$W)

  }
  
  amaris <- amaris*100
  
  cl<-c(rep(row_names[1],m),rep(row_names[2],m),rep(row_names[3],m),rep(row_names[4],m),rep(row_names[5],m))
  am<-c(amaris[1,],amaris[2,],amaris[3,],amaris[4,],amaris[5,])
  cmp<-data.frame(id=cl,val=am)
  means <- with(cmp, reorder(id,val,mean))
  boxplot(val ~ means, data=cmp,xlab=NULL,ylab="Amari Distance from True W",main=paste("dimension: ",p))
  
  cl<-c(rep(row_names[1],m),rep(row_names[2],m),rep(row_names[3],m),rep(row_names[4],m),rep(row_names[5],m))
  am<-c(SIR[1,],SIR[2,],SIR[3,],SIR[4,],SIR[5,])
  cmp<-data.frame(id=cl,val=am)
  means <- with(cmp, reorder(id,val,mean))
  boxplot(val ~ means, data=cmp,xlab=NULL,ylab="SIR",main=paste("dimension: ",p))
  
  cl<-c(rep(row_names[1],m),rep(row_names[2],m),rep(row_names[3],m),rep(row_names[4],m),rep(row_names[5],m))
  am<-c(timecost[1,],timecost[2,],timecost[3,],timecost[4,],timecost[5,])
  cmp<-data.frame(id=cl,val=am)
  means <- with(cmp, reorder(id,val,mean))
  boxplot(val ~ means, data=cmp,xlab=NULL,ylab="CPU elapsed time",main=paste("dimension: ",p))


  rl<-list(amaris=amaris,time=timecost,SIR=SIR)
  rl
  
  #print(apply(amaris.mean, c(1,2),mean))
} 