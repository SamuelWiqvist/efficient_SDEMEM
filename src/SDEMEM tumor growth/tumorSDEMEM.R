#Forward simulation (SDE)
analyticsim=function(T=20,x0=c(150,75,75),theta=c(5.8/20,1.1/sqrt(20),1.8/20,1.5/sqrt(20)))
{
 d=length(x0)
 x=matrix(0,nrow=T+1,ncol=d)
 x[1,]=x0 
 for (i in 2:(T+1)) {
  x[i,2]=x[i-1,2]*exp(rnorm(1,theta[1],theta[2]))
  x[i,3]=x[i-1,3]*exp(rnorm(1,-theta[3],theta[4]))
  x[i,1]=x[i,2]+x[i,3]
  
 }
 return(x)
}

#Forward simulation (ODE)
analyticsimODE=function(T=20,x0=c(150,75,75),theta=c(5.8/20,1.8/20))
{
 d=length(x0)
 x=matrix(0,nrow=T+1,ncol=d)
 times=seq(0,T,1)
 x[,2]=x0[2]*exp(theta[1]*times)
 x[,3]=x0[3]*exp(-theta[2]*times)
 x[,1]=x[,2]+x[,3] 
 return(x)
}

#Helper functions for inference

rmvn=function(m,V) #Multivariate normal simulation
{
  p=length(m)
  z=rnorm(p)
  return(m+t(chol(V))%*%z)
}

logprior=function(param,a=0,b=10)
{
 #assumes params follow a N(a,b^2) a priori
 return(sum(dnorm(param,a,b,log=T)))
}

#Systematic resampling
sysresamp2=function(wts,N,uni)
{
 vec=rep(0,N)
 wsum=sum(wts)
 k=1
 u=uni/N
 wsumtarg=u
 wsumcurr=wts[k]/wsum
 delta=1/N
 for(i in 1:N)
 {
  while (wsumcurr<wsumtarg)
  {
   k=k+1
   wsumcurr=wsumcurr+wts[k]/wsum
  }   
  vec[i]=k 
  wsumtarg=wsumtarg+delta
 }
 return(vec)
}

euclid=function(vec)
{
return(sqrt(sum(vec^2)))
}

#returns indices of Euclidean sorted (smallest to largest) vectors 
esort=function(xmat)
{
 N=dim(xmat)[1]
 indices=rep(0,N)
 iset=1:N
 indices[1]=which.min(xmat[,1])
 for(j in 2:N)
 {
  xstar=xmat[indices[j-1],]
  iset=iset[iset!=indices[(j-1)]]
  dist=rep(0,length(iset))
  for(i in 1:length(iset))
  {
   dist[i]=euclid(xstar-xmat[iset[i],])
  }
  ind=which.min(dist)
  #print(xmat)
  indices[j]=iset[ind]
 }
 return(indices)
}

#Forward simulation (SDE) - supply noise increments (for PMMH)
analyticsimInc=function(x0=c(75,75),theta=c(5.8/20,1.1/sqrt(20),1.8/20,1.5/sqrt(20)),bmvec)
{
 #bmvec is a length d vector of N(0,1) quantities
 d=length(x0)
 x=rep(0,d)
 x[1]=x0[1]*exp(theta[1]+theta[2]*bmvec[1])
 x[2]=x0[2]*exp(-theta[3]+theta[4]*bmvec[2])
 return(x)
}

#Weighted resampling
wrInc=function(N=100,x0=matrix(75,ncol=2,nrow=100),yT=5.02,theta=c(5.8/20,1.1/sqrt(20),1.8/20,1.5/sqrt(20)),sige=0.2,bmmat,uni)
{
  #bmmat is N*d matrix of N(0,1) quantities
  d=dim(x0)[2] #cols = no. components 
  mat=matrix(0,nrow=N,ncol=d) #store values at end here
  wts=rep(0,N)
  flag=0
  #propose and calculate weight
  if(N==1){
    end=analyticsimInc(x0,theta,bmmat)
    wts=exp(dnorm(yT,log(sum(end)),sige,log=T))
    mat[1,]=end
    indices=1
  }else{
   for(i in 1:N)
   {
     end=analyticsimInc(x0[i,],theta,bmmat[i,])
     if(is.nan(end[1])==TRUE){
      end[1]=0; flag=1
      }
      if(is.nan(end[2])==TRUE){
       end[2]=0; flag=1
      }
     wts[i]=exp(dnorm(yT,log(sum(end)),sige,log=T))
     mat[i,]=end
   }
   if(flag==0){
    sorted=esort(mat) #Euclidean sorting
   }else{
    sorted=seq(1,N,1);
   }
   mat=mat[sorted,]; wts=wts[sorted]
   #systematic resampling
   if(sum(wts)==0){
    indices=seq(1,N,1); flag=1
    }else{
   indices=sysresamp2(wts,N,uni)
   }
 } 
 return(list(mat[indices,],mean(wts),mat,flag))
}

#Estimate of log-likelihood under SDE (for CPMMH)
llike=function(N,xi,x,param,sige,bmarray,u)
{
 #returns estimated loglike using sequential weighted resampling
 n=length(x) #no. obs
 loglikecan=0
 for(j in 1:(n-1))
   {
    wrstuff=wrInc(N,xi,x[j+1],param,sige,bmarray[,,j],u[j])
    if(wrstuff[[4]]==1)
    {
     loglikecan= -1e10
     j=n
    }else{
     loglikecan=loglikecan+log(wrstuff[[2]])
     xi=wrstuff[[1]]
    } 
  }
 loglikecan
}

#Log-likelihood under ODE
llikeODE=function(x0,x,param,sige)
{
 #returns loglike 
 n=length(x) #no. obs
 sum(dnorm(log(analyticsimODE((n-1),x0,param)[,1]),x,sige,log=TRUE)) 
}

#######
#CPMMH#
#######

#Main Gibbs function
corPMMHrand=function(iters=10, bvec=rep(-2,4), dvec=rep(10,4), gvec=rep(2,4), hvec=rep(0.2,4), sigma,N=200,x0=matrix(75,ncol=2,nrow=200),x=datasim,init=c(5.8/20,1.1/sqrt(20),1.8/20,1.5/sqrt(20)),sige=0.2,tuneSige=0.01,rho=0.0)
{
  n=dim(x)[1] #no. obs
  M=dim(x)[2] #no. mice
  d=dim(x0)[2] #cols = no. components #known and same initial condition for each mouse
  bmarray=array(rnorm(N*d*(n-1)*M,0,1),dim=c(N,d,n-1,M)) #4d array of innovations
  umat=matrix(runif((n-1)*M),ncol=M,nrow=(n-1))
  p=length(init)  #no. params per mouse
  mat=array(0,dim=c(iters,(p+1),M)) #Store samples here
  curr=log(init) #current param value (log-scale)
  loglikecurr=rep(-1000000,M) #log-likelihood of current value for each mouse (to accept first candidate)
  loglikecan=rep(0,M); loglikecurrTMP=rep(0,M)
  for(j in 1:M){
    mat[1,,j]=c(curr,0) #params on log-scale and candidate log-like 
  }
  muMat=matrix(0,nrow=iters,ncol=p) #mean of log\theta_1
  tauMat=matrix(0,nrow=iters,ncol=p) #precision
  
  sigevec=rep(0,iters)
  
  muMat[1,]=c(log(5.8/20),log(1.1/sqrt(20)),log(1.8/20),log(1.5/sqrt(20))); tauMat[1,]=rep(10,4)
   
  sigevec[1]=log(sige)
  count=rep(0,M); count2=0

  for (i in 2:iters) {
    
    for(j in 1:M)
    {
      #update mouse specific params
      sige=exp(sigevec[i-1]) # current error sd value
      curr=as.vector(mat[(i-1),1:p,j])
      can = rmvn(curr,sigma[[j]]) #assuming different tuning matrix for each mouse
      #can=curr
      bmarrayprop=rho*bmarray[,,,j]+sqrt(1-rho^2)*rnorm(N*d*(n-1),0,1) #update bm increments
      uprop=pnorm(rho*qnorm(umat[,j])+sqrt(1-rho^2)*rnorm(n-1)) #update uniforms for resampling step
      loglikecan[j]=llike(N,x0,x[,j],exp(can),sige,bmarrayprop,uprop)
      laprob = loglikecan[j]-loglikecurr[j]+logprior(can,muMat[i-1,],1/sqrt(tauMat[i-1,]))-logprior(curr,muMat[i-1,],1/sqrt(tauMat[i-1,]))
      u = runif(1)
      if (log(u) < laprob)
      { 
        curr = can
        loglikecurr[j]=loglikecan[j] 
        bmarray[,,,j]=bmarrayprop
        umat[,j]=uprop
        count[j]=count[j]+1  #track no. acceptances
      }   
      mat[i,,j] = c(curr,loglikecan[j])
    }
    #update sige
    canSige=rnorm(1,log(sige),sqrt(tuneSige))
    #canSige=log(sige) 
    for(j in 1:M)
    {
     curr=as.vector(mat[i,1:p,j])
     loglikecan[j]=llike(N,x0,x[,j],exp(curr),exp(canSige),bmarray[,,,j],umat[,j])
     #uncomment for updating u - only use for rho=0 case (don't need to track)
     #bmarrayprop=rho*bmarray[,,,j]+sqrt(1-rho^2)*rnorm(N*d*(n-1),0,1) #update bm increments
     #uprop=pnorm(rho*qnorm(umat[,j])+sqrt(1-rho^2)*rnorm(n-1)) #update uniforms for resampling step
     #loglikecan[j]=llike(N,x0,x[,j],exp(curr),exp(canSige),bmarrayprop,uprop)
    }
    laprob = sum(loglikecan)-sum(loglikecurr)+logprior(canSige,0,1)-logprior(log(sige),0,1)
 
    u = runif(1)
    if (log(u) < laprob)
     { 
      sige=exp(canSige)
      loglikecurr=loglikecan
      count2=count2+1
     }
     sigevec[i]=log(sige)
  
    #update top level params
    for(k in 1:p)
    {
     muMat[i,k]=rnorm(1,((dvec[k]*bvec[k])+(tauMat[i-1,k]*sum((mat[i,k,]))))/(dvec[k]+M*tauMat[i-1,k]),sqrt(1/(dvec[k]+M*tauMat[i-1,k])))
     tauMat[i,k]=rgamma(1,gvec[k]+M/2,hvec[k]+0.5*(sum(((mat[i,k,])-muMat[i,k])^2)))  
    }
 
  }
  print(count/(iters-1)); print(count2/iters)
  return(list(mat, muMat, tauMat, sigevec))
}

#Read in data
data=scan("vol.dat")
datasim=matrix(data,ncol=10,byrow=TRUE)


#For pilot run
sigma=diag(c(0.01,0.1,0.01,0.1))
vars = list(sigma, sigma, sigma, sigma, sigma, sigma, sigma, sigma, sigma, sigma)
set.seed(1)
system.time(test<-corPMMHrand(iters=500,N=7,x=datasim,sigma=vars,x0=matrix(75,ncol=2,nrow=7),sige=0.2,tuneSige=0.02,rho=0.0))
pmat=test[[1]]
var(pmat[,5,1])
rhoL=cor(pmat[1:998,5,1],pmat[2:999,5,1])
2.16^2/(1-rhoL^2)


sT1=var(pmat[200:9999,1:4,1])
vt1 =((2.56^2)/4)*sT1
sT2=var(pmat[200:9999,1:4,2])
vt2 =((2.56^2)/4)*sT2
sT3=var(pmat[200:9999,1:4,3])
vt3 =((2.56^2)/4)*sT3
sT4=var(pmat[200:9999,1:4,4])
vt4 =((2.56^2)/4)*sT4
sT5=var(pmat[200:9999,1:4,5])
vt5 =((2.56^2)/4)*sT5

sT6=var(pmat[200:9999,1:4,6])
vt6 =((2.56^2)/4)*sT6
sT7=var(pmat[200:9999,1:4,7])
vt7 =((2.56^2)/4)*sT7
sT8=var(pmat[200:9999,1:4,8])
vt8 =((2.56^2)/4)*sT8
sT9=var(pmat[200:9999,1:4,9])
vt9 =((2.56^2)/4)*sT9
sT10=var(pmat[200:9999,1:4,10])
vt10 =((2.56^2)/4)*sT10

vars2 = list(vt1, vt2, vt3, vt4, vt5, vt6, vt7, vt8, vt9, vt10)
set.seed(1)
system.time(test<-corPMMHrand(iters=10000,N=8,x=datasim,sigma=vars2,x0=matrix(75,ncol=2,nrow=8),sige=0.2,tuneSige=0.017,rho=0.999))

################
#Euler-Maruyama#
################

#Forward simulation (SDE) - supply noise increments (for PMMH) - Euler-Maruyama
analyticsimIncEM=function(x0=c(75,75),theta=c(5.8/20,1.1/sqrt(20),1.8/20,1.5/sqrt(20)),bmmat,deltat)
{
 #bmmat is a d*n vector of N(0,deltat) quantities
 n = 1/deltat
 x=x0
 tol=1e-4
 for(i in 1:n)
 {
 x[1]=x[1]+(theta[1]+theta[2]^2/2)*x[1]*deltat+theta[2]*x[1]*bmmat[1,i]
 x[2]=x[2]+(-theta[3]+theta[4]^2/2)*x[2]*deltat+theta[4]*x[2]*bmmat[2,i]
 if(x[1]<tol){x[1]=tol}
 if(x[2]<tol){x[2]=tol} 	
 }
 return(x)
}


#Weighted resampling
wrIncEM=function(N=100,x0=matrix(75,ncol=2,nrow=100),yT=5.02,theta=c(5.8/20,1.1/sqrt(20),1.8/20,1.5/sqrt(20)),sige=0.2,bmarray,uni,deltat)
{
  #bmarray is N*d*n array of N(0,deltat) quantities
  d=dim(x0)[2] #cols = no. components 
  mat=matrix(0,nrow=N,ncol=d) #store values at end here
  wts=rep(0,N)
  flag=0
  #propose and calculate weight
  if(N==1){
    end=analyticsimIncEM(x0,theta,bmarray,deltat)
    wts=exp(dnorm(yT,log(sum(end)),sige,log=T))
    mat[1,]=end
    indices=1
  }else{
   for(i in 1:N)
   {
     end=analyticsimIncEM(x0[i,],theta,bmarray[i,,],deltat)
     if(is.nan(end[1])==TRUE){
      end[1]=0; flag=1
      }
      if(is.nan(end[2])==TRUE){
       end[2]=0; flag=1
      }
     wts[i]=exp(dnorm(yT,log(sum(end)),sige,log=T))
     mat[i,]=end
   }
   if(flag==0){
    sorted=esort(mat) #Euclidean sorting
   }else{
    sorted=seq(1,N,1);
   }
   mat=mat[sorted,]; wts=wts[sorted]
   #systematic resampling
   if(sum(wts)==0){
    indices=seq(1,N,1); flag=1
    }else{
   indices=sysresamp2(wts,N,uni)
   }
 } 
 return(list(mat[indices,],mean(wts),mat,flag))
}

#Estimate of log-likelihood under SDE (for CPMMH) -- EM
llikeEM=function(N,xi,x,param,sige,bmarray,u,deltat)
{
 #returns estimated loglike using sequential weighted resampling
 n=length(x) #no. obs
 loglikecan=0
 for(j in 1:(n-1))
   {
    wrstuff=wrIncEM(N,xi,x[j+1],param,sige,bmarray[,,,j],u[j],deltat)
    if(wrstuff[[4]]==1)
    {
     loglikecan= -1e10
     j=n
    }else{
     loglikecan=loglikecan+log(wrstuff[[2]])
     xi=wrstuff[[1]]
    } 
  }
 loglikecan
}

#Main Gibbs function
corPMMHrandEM=function(iters=10, bvec=rep(-2,4), dvec=rep(10,4), gvec=rep(2,4), hvec=rep(0.2,4), sigma,N=200,x0=matrix(75,ncol=2,nrow=200),x=datasim,init=c(5.8/20,1.1/sqrt(20),1.8/20,1.5/sqrt(20)),sige=0.2,tuneSige=0.01,rho=0.0,deltat=0.2)
{
  n2=1/deltat #no. increments per time interval, per particle, per dimension, per mouse
  n=dim(x)[1] #no. obs
  M=dim(x)[2] #no. mice
  d=dim(x0)[2] #cols = no. components #known and same initial condition for each mouse
  bmarray=array(rnorm(N*d*n2*(n-1)*M,0,sqrt(deltat)),dim=c(N,d,n2,n-1,M)) #5d array of innovations
  umat=matrix(runif((n-1)*M),ncol=M,nrow=(n-1))
  p=length(init)  #no. params per mouse
  mat=array(0,dim=c(iters,(p+1),M)) #Store samples here
  curr=log(init) #current param value (log-scale)
  loglikecurr=rep(-1000000,M) #log-likelihood of current value for each mouse (to accept first candidate)
  loglikecan=rep(0,M); loglikecurrTMP=rep(0,M)
  for(j in 1:M){
    mat[1,,j]=c(curr,0) #params on log-scale and candidate log-like 
  }
  muMat=matrix(0,nrow=iters,ncol=p) #mean of log\theta_1
  tauMat=matrix(0,nrow=iters,ncol=p) #precision
  
  sigevec=rep(0,iters)
  
  muMat[1,]=c(log(5.8/20),log(1.1/sqrt(20)),log(1.8/20),log(1.5/sqrt(20))); tauMat[1,]=rep(10,4)
   
  sigevec[1]=log(sige)
  count=rep(0,M); count2=0

  for (i in 2:iters) {
    
    for(j in 1:M)
    {
      #update mouse specific params
      sige=exp(sigevec[i-1]) # current error sd value
      curr=as.vector(mat[(i-1),1:p,j])
      can = rmvn(curr,sigma[[j]]) #assuming different tuning matrix for each mouse
      bmarrayprop=rho*bmarray[,,,,j]+sqrt(1-rho^2)*rnorm(N*d*n2*(n-1),0,sqrt(deltat)) #update bm increments
      uprop=pnorm(rho*qnorm(umat[,j])+sqrt(1-rho^2)*rnorm(n-1)) #update uniforms for resampling step
      loglikecan[j]=llikeEM(N,x0,x[,j],exp(can),sige,bmarrayprop,uprop,deltat)
      laprob = loglikecan[j]-loglikecurr[j]+logprior(can,muMat[i-1,],1/sqrt(tauMat[i-1,]))-logprior(curr,muMat[i-1,],1/sqrt(tauMat[i-1,]))
      u = runif(1)
      if (log(u) < laprob)
      { 
        curr = can
        loglikecurr[j]=loglikecan[j] 
        bmarray[,,,,j]=bmarrayprop
        umat[,j]=uprop
        count[j]=count[j]+1  #track no. acceptances
      }   
      mat[i,,j] = c(curr,loglikecan[j])
    }
    #update sige
    canSige=rnorm(1,log(sige),sqrt(tuneSige))
    for(j in 1:M)
    {
     curr=as.vector(mat[i,1:p,j])
     loglikecan[j]=llikeEM(N,x0,x[,j],exp(curr),exp(canSige),bmarray[,,,,j],umat[,j],deltat)
     #uncomment for updating u - only use for rho=0 case (don't need to track)
     #bmarrayprop=rho*bmarray[,,,,j]+sqrt(1-rho^2)*rnorm(N*d*n2*(n-1),0,sqrt(deltat)) #update bm increments
     #uprop=pnorm(rho*qnorm(umat[,j])+sqrt(1-rho^2)*rnorm(n-1)) #update uniforms for resampling step
     #loglikecan[j]=llikeEM(N,x0,x[,j],exp(curr),exp(canSige),bmarrayprop,uprop,deltat)
    }
    laprob = sum(loglikecan)-sum(loglikecurr)+logprior(canSige,0,1)-logprior(log(sige),0,1)
 
    u = runif(1)
    if (log(u) < laprob)
     { 
      sige=exp(canSige)
      loglikecurr=loglikecan
      count2=count2+1
     }
     sigevec[i]=log(sige)
  
    #update top level params
    for(k in 1:p)
    {
     muMat[i,k]=rnorm(1,((dvec[k]*bvec[k])+(tauMat[i-1,k]*sum((mat[i,k,]))))/(dvec[k]+M*tauMat[i-1,k]),sqrt(1/(dvec[k]+M*tauMat[i-1,k])))
     tauMat[i,k]=rgamma(1,gvec[k]+M/2,hvec[k]+0.5*(sum(((mat[i,k,])-muMat[i,k])^2)))  
    }
 
  }
  print(count/(iters-1)); print(count2/iters)
  return(list(mat, muMat, tauMat, sigevec))
}

set.seed(1)
system.time(test<-corPMMHrand(iters=100,N=10,x=datasim,sigma=vars2,x0=matrix(75,ncol=2,nrow=10),sige=0.2,tuneSige=0.017,rho=0.999))
set.seed(1)
system.time(test<-corPMMHrandEM(iters=100,N=30,x=datasim,sigma=vars2,x0=matrix(75,ncol=2,nrow=30),sige=0.2,tuneSige=0.017,rho=0.0,deltat=0.2))



#####
#ODE#
#####

#Main Gibbs function
corMHrandODE=function(iters=10, bvec=rep(-2,2), dvec=rep(10,2), gvec=rep(2,2), hvec=rep(0.2,2), sigma,x0=c(150,75,75),x=datasim,init=c(5.8/20,1.8/20),sige=0.2,tuneSige=0.01)
{
  n=dim(x)[1] #no. obs
  M=dim(x)[2] #no. mice
  d=dim(x0)[2] #cols = no. components #known and same initial condition for each mouse
  p=length(init)  #no. params per mouse
  mat=array(0,dim=c(iters,(p+1),M)) #Store samples here
  curr=log(init) #current param value (log-scale)
  loglikecurr=rep(-1000000,M) #log-likelihood of current value for each mouse (to accept first candidate)
  loglikecan=rep(0,M); loglikecurrTMP=rep(0,M)
  for(j in 1:M){
    mat[1,,j]=c(curr,0) #params on log-scale and candidate log-like 
  }
  muMat=matrix(0,nrow=iters,ncol=p) #mean of log\theta_1
  tauMat=matrix(0,nrow=iters,ncol=p) #precision
  
  sigevec=rep(0,iters)
  
  muMat[1,]=c(log(5.8/20),log(1.8/20)); tauMat[1,]=rep(10,2)
   
  sigevec[1]=log(sige)
  count=rep(0,M); count2=0

  for (i in 2:iters) {
    
    for(j in 1:M)
    {
      #update mouse specific params
      sige=exp(sigevec[i-1]) # current error sd value
      curr=as.vector(mat[(i-1),1:p,j])
      can = rmvn(curr,sigma[[j]]) #assuming different tuning matrix for each mouse
      loglikecan[j]=llikeODE(x0,x[,j],exp(can),sige)
      laprob = loglikecan[j]-loglikecurr[j]+logprior(can,muMat[i-1,],1/sqrt(tauMat[i-1,]))-logprior(curr,muMat[i-1,],1/sqrt(tauMat[i-1,]))
      u = runif(1)
      if (log(u) < laprob)
      { 
        curr = can
        loglikecurr[j]=loglikecan[j] 
        count[j]=count[j]+1  #track no. acceptances
      }   
      mat[i,,j] = c(curr,loglikecurr[j])
    }
    #update sige
    canSige=rnorm(1,log(sige),sqrt(tuneSige))
    for(j in 1:M)
    {
     curr=as.vector(mat[i,1:p,j])
     loglikecan[j]=llikeODE(x0,x[,j],exp(curr),exp(canSige))
    }
    laprob = sum(loglikecan)-sum(loglikecurr)+logprior(canSige,0,1)-logprior(log(sige),0,1)
 
    u = runif(1)
    if (log(u) < laprob)
     { 
      sige=exp(canSige)
      loglikecurr=loglikecan
      count2=count2+1
     }
     sigevec[i]=log(sige)
  
    #update top level params
    for(k in 1:p)
    {
     muMat[i,k]=rnorm(1,((dvec[k]*bvec[k])+(tauMat[i-1,k]*sum((mat[i,k,]))))/(dvec[k]+M*tauMat[i-1,k]),sqrt(1/(dvec[k]+M*tauMat[i-1,k])))
     tauMat[i,k]=rgamma(1,gvec[k]+M/2,hvec[k]+0.5*(sum(((mat[i,k,])-muMat[i,k])^2)))  
    }
 
  }
  print(count/(iters-1)); print(count2/iters)
  return(list(mat, muMat, tauMat, sigevec))
}


#####
#LNA#
#####

#Helper functions

sqrtmat=function(V)
{
 spec=svd(V)
 return(spec$u%*%diag(sqrt(spec$d))%*%t(spec$u))
}

rmvn2=function(m,V) #Use SVD rather than Cholesky
{
  p=length(m)
  z=rnorm(p)
  return(m+sqrtmat(V)%*%z)
}

alphaL=function(x,theta) #Drift of log process
{
 c((theta[1]+0.5*theta[2]^2)*exp(x[2]-x[1])+(-theta[3]+0.5*theta[4]^2)*exp(x[3]-x[1])-0.5*((theta[2]^2)*exp(2*(x[2]-x[1])) + (theta[4]^2)*exp(2*(x[3]-x[1]))),theta[1],-theta[3])
}

betaL=function(x,theta) #Diffusion matrix of log process
{
 mat=matrix(0,ncol=3,nrow=3,byrow=T)
 mat[1,1]=(theta[2]^2)*exp(2*(x[2]-x[1]))+(theta[4]^2)*exp(2*(x[3]-x[1]))
 mat[2,2]=theta[2]^2
 mat[3,3]=theta[4]^2
 mat[1,2]=(theta[2]^2)*exp(2*(x[2]-x[1])); mat[2,1]=mat[1,2]
 mat[1,3]=(theta[4]^2)*exp(2*(x[3]-x[1])); mat[3,1]=mat[1,3]
 return(mat)
}

Fmat=function(x=c(log(150),log(75),log(75)),theta=c(5.8/20,1.1/sqrt(20),1.8/20,1.5/sqrt(20))) #Jacobian matrix
{
 F=matrix(0,ncol=3,nrow=3)
 F[1,]=c(-(theta[1]+0.5*theta[2]^2)*exp(x[2]-x[1])-(-theta[3]+0.5*theta[4]^2)*exp(x[3]-x[1])+(theta[2]^2)*exp(2*(x[2]-x[1]))+(theta[4]^2)*exp(2*(x[3]-x[1])),(theta[1]+0.5*theta[2]^2)*exp(x[2]-x[1])-(theta[2]^2)*exp(2*(x[2]-x[1])),(-theta[3]+0.5*theta[4]^2)*exp(x[3]-x[1])-(theta[4]^2)*exp(2*(x[3]-x[1])))
F
}

tstep3=function(T=1,dt=0.1,initz=c(log(150),log(75),log(75)),initV=matrix(0,ncol=3,nrow=3),theta=c(5.8/20,1.1/sqrt(20),1.8/20,1.5/sqrt(20)),afun=alphaL,bfun=betaL)
{
 n=T/dt
 d=length(initz)
 z=matrix(0,nrow=d,ncol=n+1) #store z solution here
 V=array(0,dim=c(d,d,n+1)) #store V solution here
 P=array(0,dim=c(d,d,n+1)) #store P solution here
 #Initialise
 z[,1]=initz; V[,,1]=initV; P[,,1]=diag(rep(1,d))
 for(i in 2:(n+1))
 { 
  Fm=Fmat(z[,i-1],theta)
  bmat=bfun(z[,i-1],theta)
  z[,i]=z[,i-1]+afun(z[,i-1],theta)*dt
  P[,,i]=P[,,i-1]+dt*Fm%*%P[,,i-1]
  V[1,1,i]=V[1,1,i-1]+dt*(2*(V[1,1,i-1]*Fm[1,1]+V[1,2,i-1]*Fm[1,2]+V[1,3,i-1]*Fm[1,3])+bmat[1,1])
  V[2,2,i]=V[2,2,i-1]+dt*bmat[2,2]; V[3,3,i]=V[3,3,i-1]+dt*bmat[3,3] 
  V[1,2,i]=V[1,2,i-1]+dt*(Fm[1,1]*V[1,2,i-1]+Fm[1,2]*V[2,2,i-1]+bmat[1,2])
  V[1,3,i]=V[1,3,i-1]+dt*(Fm[1,1]*V[1,3,i-1]+Fm[1,3]*V[3,3,i-1]+bmat[1,3])
  V[2,1,i]=V[1,2,i]; V[3,1,i]=V[1,3,i]
 }
list(z,V,P)
}

filter=function(ODEdt=0.2,param=c(5.8/20,1.1/sqrt(20),1.8/20,1.5/sqrt(20)),sige=0.2,a=log(c(150,75,75)),C=diag(c(0,0,0)),F=c(1,0,0),x=data,afun=alphaL,bfun=betaL,flag=0)
{
 ODEend=1/ODEdt +1
 n=length(x) #no. obs
 xsamp=matrix(0,ncol=n,nrow=3)
 mt=rep(0,3); Vt=matrix(0, ncol=3, nrow=3); Pt=diag(c(1,1,1))
 mtA=array(0,dim=c(3,1,n)); VtA=array(0,dim=c(3,3,n)); #Store m and V here
 atA=array(0,dim=c(3,1,n)); CtA=array(0,dim=c(3,3,n)); #Store a and C here
 PtA=array(0,dim=c(3,3,n));
 ll=0
 #update marginal likelihood
 ll=ll+dnorm(x[1],t(F)%*%a,sqrt(t(F)%*%C%*%F+sige^2),log=T)
 #update posterior
 a=a+C%*%F%*%solve(t(F)%*%C%*%F+sige^2)%*%(x[1]-t(F)%*%a)
 C=C-C%*%F%*%solve(t(F)%*%C%*%F+sige^2)%*%t(F)%*%C 
 atA[,,1]=a; CtA[,,1]=C
 #loop
 for(i in 1:(n-1))
 {
  #update m and V
  lnastep=tstep3(1,ODEdt,a,C,param,afun,bfun)     
  mt=lnastep[[1]][,ODEend]
  Vt=lnastep[[2]][,,ODEend]
  Pt=lnastep[[3]][,,ODEend]  
  mtA[,,(i+1)]=mt; VtA[,,(i+1)]=Vt; PtA[,,(i+1)]=Pt	
  #update marginal likelihood
  ll=ll+dnorm(x[i+1],t(F)%*%mt,sqrt(t(F)%*%Vt%*%F+sige^2),log=T)
  #update posterior
  a=mt+Vt%*%F%*%solve(t(F)%*%Vt%*%F+sige^2)%*%(x[i+1]-t(F)%*%mt)
  C=Vt-Vt%*%F%*%solve(t(F)%*%Vt%*%F+sige^2)%*%t(F)%*%Vt  
  atA[,,(i+1)]=a; CtA[,,(i+1)]=C
 }
 if(flag>0) #Backward sweep
 {
  #sample end
  xsamp[,n]=rmvn2(a,C)
  #Backward sampler
  for(i in (n-1):1)
  {
   #calculate ahat and Chat
   ahat = atA[,,i] + CtA[,,(i)]%*%t(PtA[,,(i+1)])%*%solve(VtA[,,(i+1)])%*%(xsamp[,i+1]-mtA[,,(i+1)])
   Chat = CtA[,,i] - CtA[,,i]%*%t(PtA[,,(i+1)])%*%solve(VtA[,,(i+1)])%*%PtA[,,(i+1)]%*%CtA[,,(i)]
   xsamp[,i]=rmvn2(ahat,Chat)
  }   
 }
 return(list(ll,xsamp))
}

#Main Gibbs function
LNArand=function(iters=10, bvec=rep(-2,4), dvec=rep(10,4), gvec=rep(2,4), hvec=rep(0.2,4),sigma,x=datasim,init=c(5.8/20,1.1/sqrt(20),1.8/20,1.5/sqrt(20)),sige=0.2,tuneSige=0.01,a=log(c(150,75,75)),C=diag(c(0,0,0)),F=c(1,0,0),afun=alphaL,bfun=betaL)
{
  n=dim(x)[1] #no. obs
  M=dim(x)[2] #no. mice
  p=length(init)  #no. params per mouse
  mat=array(0,dim=c(iters,(p+1),M)) #Store samples here
  curr=log(init) #current param value (log-scale)
  loglikecurr=rep(-1000000,M) #log-likelihood of current value for each mouse (to accept first candidate)
  loglikecan=rep(0,M); loglikecurrTMP=rep(0,M)
  for(j in 1:M){
    mat[1,,j]=c(curr,0) #params on log-scale and candidate log-like 
  }
  muMat=matrix(0,nrow=iters,ncol=p) #mean of log\theta_1
  tauMat=matrix(0,nrow=iters,ncol=p) #precision
  sigevec=rep(0,iters)
  muMat[1,]=c(log(5.8/20),log(1.1/sqrt(20)),log(1.8/20),log(1.5/sqrt(20))); tauMat[1,]=rep(10,4)
  sigevec[1]=log(sige)
  count=rep(0,M); count2=0

  for (i in 2:iters) {
    
    for(j in 1:M)
    {
      #update mouse specific params
      sige=exp(sigevec[i-1]) # current error sd value
      curr=as.vector(mat[(i-1),1:p,j])
      can = rmvn(curr,sigma[[j]]) #assuming different tuning matrix for each mouse
      loglikecan[j]=filter(0.2,exp(can),sige,a,C,F,x[,j],afun,bfun,0)[[1]]  
      if(is.na(loglikecan[j]))
      {
        loglikecan[j]=-1e10 #reject
      }
      laprob = loglikecan[j]-loglikecurr[j]+logprior(can,muMat[i-1,],1/sqrt(tauMat[i-1,]))-logprior(curr,muMat[i-1,],1/sqrt(tauMat[i-1,]))
      u = runif(1)
      if (log(u) < laprob)
      { 
        curr = can
        loglikecurr[j]=loglikecan[j] 
        count[j]=count[j]+1  #track no. acceptances
      }   
      mat[i,,j] = c(curr,loglikecurr[j])
    }
    #update sige
    canSige=rnorm(1,log(sige),sqrt(tuneSige))
    for(j in 1:M)
    {
     curr=as.vector(mat[i,1:p,j])
     loglikecan[j]=filter(0.2,exp(curr),exp(canSige),a,C,F,x[,j],afun,bfun,0)[[1]]
     if(is.na(loglikecan[j]))
      {
        loglikecan[j]=-1e10 #reject
      }
     }
    laprob = sum(loglikecan)-sum(loglikecurr)+logprior(canSige,0,1)-logprior(log(sige),0,1)
 
    u = runif(1)
    if (log(u) < laprob)
     { 
      sige=exp(canSige)
      loglikecurr=loglikecan
      count2=count2+1
     }
     sigevec[i]=log(sige)
  
    #update top level params
    for(k in 1:p)
    {
     muMat[i,k]=rnorm(1,((dvec[k]*bvec[k])+(tauMat[i-1,k]*sum((mat[i,k,]))))/(dvec[k]+M*tauMat[i-1,k]),sqrt(1/(dvec[k]+M*tauMat[i-1,k])))
     tauMat[i,k]=rgamma(1,gvec[k]+M/2,hvec[k]+0.5*(sum(((mat[i,k,])-muMat[i,k])^2)))  
    }
 
  }
  print(count/(iters-1)); print(count2/iters)
  return(list(mat, muMat, tauMat, sigevec))
}
