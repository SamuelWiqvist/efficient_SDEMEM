
rmvn=function(m,V)
{
  p=length(m)
  z=rnorm(p)
  return(m+t(chol(V))%*%z)
}

#params (beta,gamma,delta,tau)

alpha=function(x,theta)
{
  c((theta[1]+0.5*theta[2]^2)*x[1],(-theta[3]+0.5*theta[4]^2)*x[2])
}

beta=function(x,theta)
{
  mat=matrix(0,ncol=2,nrow=2,byrow=T)
  mat[1,1]=(theta[2]*x[1])^2
  mat[2,2]=(theta[4]*x[2])^2
  return(mat)
}


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
  x
}

#vola=analyticsim()
#plot(ts(log(vola),start=0,deltat=1))

#generate multiple trajectories of log(Vt) and look at distribution against time

analyticsimMany=function(N=1000,T=20,x0=c(150,75,75),theta=c(5.8/20,1.1/sqrt(20),1.8/20,1.5/sqrt(20)))
{
  mat=matrix(0,nrow=T+1,ncol=N)
  for (i in 1:N) 
  {
    mat[,i]=log(analyticsim(T,x0,theta)[,1])
  }
  mat
}

#amat=analyticsimMany(N=10000)
#ame=apply(amat,1,mean)
#alq=apply(amat,1,quantile,0.025)
#auq=apply(amat,1,quantile,0.975)
#plot(ts(ame,start=0,deltat=1),ylim=c(2.5,12))
#lines(ts(alq,start=0,deltat=1),col=2)
#lines(ts(auq,start=0,deltat=1),col=2)

###Inference###

#Simulate some data

set.seed(1)
#vola=analyticsim()
#data=log(vola[,1])+rnorm(21,0,0.2)
#plot(ts(log(vola[,1]),start=0,deltat=1))
#lines(ts(data,start=0,deltat=1),col=2)


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
  vec
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
  indices
}

analyticsimInc=function(x0=c(75,75),theta=c(5.8/20,1.1/sqrt(20),1.8/20,1.5/sqrt(20)),bmvec)
{
  #bmvec is a length d vector of N(0,1) quantities
  d=length(x0)
  x=rep(0,d)
  x[1]=x0[1]*exp(theta[1]+theta[2]*bmvec[1])
  x[2]=x0[2]*exp(-theta[3]+theta[4]*bmvec[2])
  
  x
}


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
    #if(sum(is.nan(mat))>0)
    #{
    #flag=1; sorted=1:N
    #}else{
    if(flag==0){
      sorted=esort(mat) #Euclidean sorting
    }else{
      sorted=seq(1,N,1);
    }
    #}
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

#test

#analyticsimInc(bmvec=rnorm(2))

#test=wrInc(bmmat=matrix(rnorm(100*2),ncol=2,nrow=100),uni=runif(1))


logprior=function(param,a=0,b=10)
{
  #assumes params follow a N(a,b^2) a priori
  return(sum(dnorm(param,a,b,log=T)))
}

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
      #print(x0)
      #print(log(wrstuff[[2]]))
    } 
  }
  loglikecan
}

#system.time(test<-corPMMH(iters=10000,N=6,x=data,sigma=sigtune,rho=0.999))

#system.time(test<-corPMMH(iters=10000,N=6,x=data,sigma=sigtune,rho=0.999))

############
#Try a single site update on the latent states. Use LNA to initialise.
############
####LNA#####
############

sqrtmat=function(V)
{
  spec=svd(V)
  spec$u%*%diag(sqrt(spec$d))%*%t(spec$u)
}

rmvn2=function(m,V)
{
  p=length(m)
  z=rnorm(p)
  return(m+sqrtmat(V)%*%z)
}

alphaL=function(x,theta)
{
  c((theta[1]+0.5*theta[2]^2)*exp(x[2]-x[1])+(-theta[3]+0.5*theta[4]^2)*exp(x[3]-x[1])-0.5*((theta[2]^2)*exp(2*(x[2]-x[1])) + (theta[4]^2)*exp(2*(x[3]-x[1]))),theta[1],-theta[3])
}

betaL=function(x,theta)
{
  mat=matrix(0,ncol=3,nrow=3,byrow=T)
  mat[1,1]=(theta[2]^2)*exp(2*(x[2]-x[1]))+(theta[4]^2)*exp(2*(x[3]-x[1]))
  mat[2,2]=theta[2]^2
  mat[3,3]=theta[4]^2
  mat[1,2]=(theta[2]^2)*exp(2*(x[2]-x[1])); mat[2,1]=mat[1,2]
  mat[1,3]=(theta[4]^2)*exp(2*(x[3]-x[1])); mat[3,1]=mat[1,3]
  return(mat)
}

Fmat=function(x=c(log(150),log(75),log(75)),theta=c(5.8/20,1.1/sqrt(20),1.8/20,1.5/sqrt(20)))
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
  z[,1]=initz
  V[,,1]=initV
  P[,,1]=diag(rep(1,d))
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
    
    #print(V[,,i-1]%*%t(Fm)+bfun(z[,i-1],theta)+Fm%*%V[,,i-1])
    #print(z[,i])
    #print(V[,,i])
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

#test=filter(ODEdt=0.2,param=c(5.8/20,1.1/sqrt(20),1.8/20,1.5/sqrt(20)),sige=0.2)

loglike=function(param=c(5.8/20,1.1/sqrt(20),1.8/20,1.5/sqrt(20)),x=(test[[2]][2:3,]))
{
  #x on log scale
  n=dim(x)[2]
  sum(dlnorm(exp(x[1,2:n])/exp(x[1,1:(n-1)]),param[1],param[2],log=TRUE))+sum(dlnorm(exp(x[2,2:n])/exp(x[2,1:(n-1)]),-param[3],param[4],log=TRUE))
}

#loglike(param=c(5.8/20,1.1/sqrt(20),1.8/20,1.5/sqrt(20)),x=(test[[2]][2:3,]))

#Generate data

#Simulating 10 units each with 21 observations

set.seed(3)

theta1vec=exp(rnorm(10, log(5.8/20), sqrt(0.1)))
theta2vec=exp(rnorm(10, log(1.1/sqrt(20)), sqrt(0.1)))
theta3vec=exp(rnorm(10, log(1.8/20), sqrt(0.1)))
theta4vec=exp(rnorm(10, log(1.5/sqrt(20)), sqrt(0.1)))


M=10
datasim=matrix(0,nrow=21,ncol=M)
datasim2=matrix(0,nrow=21,ncol=M)

for(i in 1:M){
  datasim2[,i]=log(analyticsim(theta=c(theta1vec[i], theta2vec[i], theta3vec[i], theta4vec[i]))[,1])
  rnorm(1,0,0.2)
  #datasim[,i]=datasim2[,i]+rnorm(21,0,0.2) #sd=0.2
} 
for(i in 1:M){
  datasim[,i]=datasim2[,i]+rnorm(21,0,0.2) #sd=0.2
} 

plot(datasim[,1], ylim=c(5, 15), type = "l")
lines(datasim2[,1], type = "l",lty=2) #without noise
for(i in 2:M){
  lines(datasim[,i], col = i)
  lines(datasim2[,i], col = i, lty=2) #without noise
}


########################################
#####Random effects on all 4 parameters#
########################################

#######
#CPMMH#
#######


corPMMHrand=function(iters=10, bvec=rep(-2,4), dvec=rep(10,4), gvec=rep(2,4), hvec=rep(0.2,4), sigma=diag(c(0.01,0.1,0.01,0.1)),N=200,x0=matrix(75,ncol=2,nrow=200),x=datasim,init=c(5.8/20,1.1/sqrt(20),1.8/20,1.5/sqrt(20)),sige=0.2,tuneSige=0.01,rho=0.0)
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
      mat[i,,j] = c(curr,loglikecurr[j])
    }
    #update sige
    canSige=rnorm(1,log(sige),sqrt(tuneSige))
    for(j in 1:M)
    {
      curr=as.vector(mat[i,1:p,j])
      loglikecan[j]=llike(N,x0,x[,j],exp(curr),exp(canSige),bmarray[,,,j],umat[,j])
      #uncomment for updating u - only use for rho=0 case
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

sigma=diag(c(0.01,0.1,0.01,0.1))
vars = list(sigma, sigma, sigma, sigma, sigma, sigma, sigma, sigma, sigma, sigma)

set.seed(1)
#set N and rho accordingly
system.time(test<-corPMMHrand(iters=10000,N=30,x=datasim,sigma=vars,x0=matrix(75,ncol=2,nrow=30),sige=0.2,tuneSige=0.017,rho=0.0))

pmat=test[[1]]
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
#set N and rho accordingly 
system.time(test1<-corPMMHrand(iters=500000,N=30,x=datasim,sigma=vars2,x0=matrix(75,ncol=2,nrow=30),sige=0.2,tuneSige=0.017,rho=0.0))
pmat=test1[[1]]
mus = test1[[2]]
taus=test1[[3]]
sige = test1[[4]]

#run plot sheet
for(i in 1:10){
  print(effectiveSize(pmat[,,i]))
}
mus = pmat2[[2]]
effectiveSize(mus)
effectiveSize(taus)
effectiveSize(sige)

#####
#LNA#
#####

LNArand=function(iters=10, bvec=rep(-2,4), dvec=rep(10,4), gvec=rep(2,4), hvec=rep(0.2,4), sigma=diag(c(0.01,0.1,0.01,0.1)),x=datasim,init=c(5.8/20,1.1/sqrt(20),1.8/20,1.5/sqrt(20)),sige=0.2,tuneSige=0.01,a=log(c(150,75,75)),C=diag(c(0,0,0)),F=c(1,0,0),afun=alphaL,bfun=betaL)
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

set.seed(1)
system.time(test3<-LNArand(iters=1000,x=datasim,sigma=vars2,sige=0.2,tuneSige=0.017))
pmat3=test3[[1]]

###################
#Data Augmentation#
###################
#Target joint posterior of params and latent x given data y

gibbsRand=function(iters=1000,bvec=rep(-2,4), dvec=rep(10,4), gvec=rep(2,4), hvec=rep(0.2,4), sigma=diag(c(0.01,0.1,0.01,0.1)), x0=c(log(75),log(75)),x=datasim,init=c(5.8/20,1.1/sqrt(20),1.8/20,1.5/sqrt(20)),sige=0.2,tuneSige=0.01)
{ 
  n=dim(x)[1] #no. obs
  M=dim(x)[2] #no. mice
  p=length(init)  #no. params per mouse
  matThet=array(0,dim=c(iters,p,M))  #store parameter values here
  matX=array(0,dim=c(2,n,M,iters)) #store process values here 
  curr=log(init) #current param value (log-scale)
  matThet[1,,]=rep(curr,M) #params on log-scale 
  #Initialise matX[,,k,1] with draws from LNA
  for(k in 1:M){
    test=filter(0.2,init,sige,log(c(150,75,75)),diag(c(0,0,0)),c(1,0,0),x[,k],alphaL,betaL,1)
    matX[,,k,1]=test[[2]][2:3,]; matX[,1,k,1]=x0
  }
  
  muMat=matrix(0,nrow=iters,ncol=p) #mean of log\theta_1
  tauMat=matrix(0,nrow=iters,ncol=p) #precision
  sigevec=rep(0,iters)
  muMat[1,]=c(log(5.8/20),log(1.1/sqrt(20)),log(1.8/20),log(1.5/sqrt(20))); tauMat[1,]=rep(10,4)
  sigevec[1]=log(sige)
  
  for (i in 2:iters) 
  {
    for(k in 1:M)
    { 
      xcurr = matX[,,k,(i-1)]
      curr=matThet[(i-1),,k]
      ###Update parameters 
      #propose
      can = rmvn(curr,sigma[[k]])
      
      #calculate likelihood at proposed and current values 
      laprob = loglike(exp(can),xcurr)-loglike(exp(curr),xcurr)+logprior(can,muMat[i-1,],1/sqrt(tauMat[i-1,]))-logprior(curr,muMat[i-1,],1/sqrt(tauMat[i-1,]))
      #store in matThet
      if(log(runif(1))<laprob)
      {
        curr = can  #with correct probability, chain moves
      }
      matThet[i,,k]=curr #store
      
      ###Update latent process between observations - use a loop
      
      for(j in 2:(n-1))
      {
        #propose half way between neighbouring log values - independence proposal
        xprop = rnorm(2,0.5*(xcurr[,j-1]+xcurr[,j+1]),c(exp(curr[2]),exp(curr[4])))
        propcan = sum(dnorm(xprop,0.5*(xcurr[,j-1]+xcurr[,j+1]),c(exp(curr[2]),exp(curr[4])),log=TRUE))
        propcurr = sum(dnorm(xcurr[,j],0.5*(xcurr[,j-1]+xcurr[,j+1]),c(exp(curr[2]),exp(curr[4])),log=TRUE))
        targetcan = sum(dnorm(xprop,xcurr[,j-1]+c(exp(curr[1]),-exp(curr[3])),c(exp(curr[2]),exp(curr[4])),log=TRUE))+sum(dnorm(xcurr[,j+1],xprop+c(exp(curr[1]),-exp(curr[3])),c(exp(curr[2]),exp(curr[4])),log=TRUE))+dnorm(x[j,k],log(sum(exp(xprop))),sige,log=TRUE)
        targetcurr = sum(dnorm(xcurr[,j],xcurr[,j-1]+c(exp(curr[1]),-exp(curr[3])),c(exp(curr[2]),exp(curr[4])),log=TRUE))+sum(dnorm(xcurr[,j+1],xcurr[,j]+c(exp(curr[1]),-exp(curr[3])),c(exp(curr[2]),exp(curr[4])),log=TRUE))+dnorm(x[j,k],log(sum(exp(xcurr[,j]))),sige,log=TRUE)
        laprob = targetcan-targetcurr+propcurr-propcan
        
        if(log(runif(1))<laprob)
        {
          xcurr[,j]=xprop 
        }
      }
      #Update end - forward simulate
      xprop=rnorm(2,xcurr[,n-1]+c(exp(curr[1]),-exp(curr[3])),c(exp(curr[2]),exp(curr[4])))
      laprob = dnorm(x[n,k],log(sum(exp(xprop))),sige,log=TRUE) - dnorm(x[n,k],log(sum(exp(xcurr[,n]))),sige,log=TRUE) 
      if(log(runif(1))<laprob)
      {
        xcurr[,n]=xprop 
      }
      matX[,,k,i]=xcurr #store  
    }
    #update sige
    canSige=rnorm(1,log(sige),sqrt(tuneSige))
    loglikecurr=0; loglikecan=0
    for(k in 1:M)
    {
      xcurr = matX[,,k,i]
      for( j in 1:n)
      {
        loglikecan=loglikecan+dnorm(x[j,k],log(sum(exp(xcurr[,j]))),exp(canSige),log=TRUE)
        loglikecurr=loglikecurr+dnorm(x[j,k],log(sum(exp(xcurr[,j]))),sige,log=TRUE)
      }
    }
    laprob = loglikecan-loglikecurr+logprior(canSige,0,1)-logprior(log(sige),0,1)
    u = runif(1)
    if (log(u) < laprob)
    { 
      sige=exp(canSige)
    }
    sigevec[i]=log(sige)
    
    #update top level params
    for(k in 1:p)
    {
      muMat[i,k]=rnorm(1,((dvec[k]*bvec[k])+(tauMat[i-1,k]*sum((matThet[i,k,]))))/(dvec[k]+M*tauMat[i-1,k]),sqrt(1/(dvec[k]+M*tauMat[i-1,k])))
      tauMat[i,k]=rgamma(1,gvec[k]+M/2,hvec[k]+0.5*(sum(((matThet[i,k,])-muMat[i,k])^2)))  
    }
  }
  
  return(list(matThet,matX,muMat,tauMat,sigevec))
}

set.seed(1)
system.time(test4<-gibbsRand(iters=100,x=datasim,sigma=vars2,sige=0.2,tuneSige=0.017))



