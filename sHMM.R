# Implementation of the model proposed in the paper: 
#"A Scalable and Interpretable Supervised Hidden Markov Model 
# for Tractor Operational State Prediction"


sHMM=function(train,test,s.train,s.test,modello,p)
{
  k=length(table(s.train))
  
  # PHASE 1: MODEL ESTIMATION ON THE TRAINING SET
  # Response variable (state to be predicted):
  s.train= train$OptState  
  s.test=test$OptState
  
  # Total number of rows in the training set:
  ntrain=nrow(train)
  
  # Number of data loggers in the test set:
  ntest=nrow(test)
  
  
  # Estimation of parameters according to variable type and the state:
  
  parameter.s=list()
  parameters=list()
  
  for (h in 1:k) {sub=train[s.train==h,]
  
  for(j in 1:p)  {if (modello[j]=="binary") parameter.s[[j]]=sum(sub[,j])/sum(sub$w)
  if (modello[j]=="gaussian"){
    mu=sum(sub[,j])/sum(sub$w)
    sigma2=mean((sub[,j]-sub$w*mu)^2/sub$w)
    sigma2=sqrt(sigma2)
    parameter.s[[j]]=c(mu,sigma2)}
  }
  
  
  parameters[[h]]=parameter.s              
  }
  
  
  
  # PHASE 2: ESTIMATION OF INITIAL AND TRANSITION PROBABILITIES:
  
  #INITIAL PROBABILITIES:
  
  pi.0h =prop.table(table(train[train$Tempo==1,]$OptState))
  
  #TRANSITION PROBABILITIES:
  
  quali.dl=unique(train$Data.logger)
  n.dl=length(quali.dl)
  num=array(0,c(k,k,n.dl))
  
  for (i in 1:n.dl) {
    s.sub.train=s.train[train$Data.logger==quali.dl[i]]
    quanti=length(s.sub.train)
    if (quanti>1) for(h in 1:k) for(r in 1:k) {
      for(tt in 2:quanti) if (s.sub.train[tt-1]==h & s.sub.train[tt]==r) num[h,r,i]=num[h,r,i]+1
    }
  }
  
  
  num=apply(num,c(1,2),sum)
  den=matrix(rowSums(num),k,k)
  pi.hr=num/den
  
  
  # PHASE 3: PREDICTION ON THE TEST SET
  
  ### theta is a vector containing means and standard deviations
  
  dens.gaussian= function(y,theta,w){
    mu<-theta[1]*w
    sigma<-theta[2]*sqrt(w)
    f.y <- dnorm(y, mean=mu, sd=sigma, log=TRUE)
    return(f.y)
  }
  
  
  dens.binary=function(y,theta,w){
    f.y <- dbinom(y,w,prob=theta,log=TRUE)
    return(f.y)
  }
  
  
  # Variable-wise density computed for each possible state:
  
  dens=array(0,c(ntest,p,k))
  
  for (h in 1:k) for(j in 1:p){
    if (modello[j]=="binary") dens[,j,h]=dens.binary(test[,j],parameters[[h]][[j]],test$w)
    if (modello[j]=="gaussian") dens[,j,h]=dens.gaussian(test[,j],parameters[[h]][[j]],test$w)
  }
  
  px.y=apply(dens,c(1,3),sum) 
  px.y=exp(px.y)
  px.y=ifelse(px.y==0,exp(-740),px.y)
  px.y=as.matrix(px.y)
  
  
  ## alpha ##
  
  quali.dl=unique(test$Data.logger)
  n.dl=length(quali.dl)
  py.x=matrix(0,ntest,k) 
  
  for (i in 1:n.dl) {
    index=test$Data.logger==quali.dl[i]
    alpha.bw=matrix(0,sum(index),k)
    beta.bw=matrix(0,sum(index),k)  
    px.y.sub=px.y[index,]
    py.x.sub=matrix(0,sum(index),k)
    if(is.vector(px.y.sub)) px.y.sub=t(as.matrix(px.y.sub))
    
    foo <-(px.y.sub[1,]*pi.0h) 
    foo=ifelse(foo==0,10^(-300),foo)
    sumfoo <- sum(foo) 
    lscale <- log(sumfoo)
    foo <- foo/sumfoo
    alpha.bw[1,] <- log(foo)+lscale
    
    for (h in 2:sum(index))  {
      if(nrow(px.y.sub)!=1) {
        foo = (foo%*%(pi.hr))*px.y.sub[h,] 
        foo=ifelse(foo==0,10^(-300),foo)
        sumfoo <- sum(foo)
        lscale <- lscale+log(sumfoo)
        foo <- foo/sumfoo
        alpha.bw[h,] <- log(foo)+lscale
      }
    }
    
    beta.bw[sum(index),] <- rep(0,k)
    foo <- rep (1/k,k)
    lscale <- log(k)
    for (h in (sum(index)-1):1) {
      if(nrow(px.y.sub)!=1){ 
        foo <- (pi.hr)%*%( px.y.sub[h+1,]*foo) 
        foo=ifelse(foo==0,10^(-300),foo)
        beta.bw[h,] <- log(foo)+lscale
        sumfoo <- sum(foo)
        foo <- foo/sumfoo
        lscale <- lscale+log(sumfoo)}
    }
    
    
    for (h in 1:sum(index)) {
      cc <- max(alpha.bw[h,]+beta.bw[h,])
      den <- cc+log(sum(exp(alpha.bw[h,]+beta.bw[h,]-cc)))
      py.x.sub[h,]<-alpha.bw[h,]+beta.bw[h,] -den 
      py.x.sub[h,]<-exp(py.x.sub[h,])}
    
    py.x[index,]=py.x.sub
  }
  
  prev.test=apply(py.x,1,which.max)
  
  
  misc=function(a1,a2)
  {
    library(mclust)
    return(classError(a1,a2)$error)
  }
  
  mr=misc(prev.test,s.test)
  
  return(mr)
}

