LR_int=function(y1,len1,l1){
  if(y1>0 & y1<=l1){
    a=c(.Machine$double.eps,l1)
  }else{
    k=as.integer((y1-l1)/len1)+1
    a=c(l1+((k-1)*len1),l1+(k*len1))
  }
  return(a)
}
data_0_BC=function(n,b0,b1,b2,g1,g2,lam,cenrate){ 
  # lam is baseline exponential parameter 
  L=rep(NA,n)
  R=rep(NA,n)
  d=rep(NA,n)
  x1=rbinom(n=n,size=1,prob=0.5) # binary covariate
  x2=runif(n,min = 0.1,max = 20) # continuous covariate
  phi=exp(b0+(b1*x1)+(b2*x2))
  U=runif(n,min=0,max=1)
  C=rexp(n,rate=cenrate)
  p0 = exp(-phi)
  count.obs=0
  count.cure=0
  for(i in 1:n){
    if(U[i]<=p0[i]){
      L[i]=C[i]
      R[i]=Inf
      d[i]=0
      count.cure=count.cure+1
    }else{
      U1 = runif(1,min=0,max=1)
      y=-exp(-(g1*x1[i])-(g2*x2[i]))
      *(1/lam)*log(1+(exp(-b0-(b1*x1[i])-(b2*x2[i]))
                      *log(p0[i]+((1-p0[i])*U1))))
      t=min(y,C[i])
      if(t==C[i]){    
        L[i]=C[i]
        R[i]=Inf
        d[i]=0
      }else{
        len=runif(1,0.2,0.7)
        l=runif(1,0,1)
        ans=LR_int(t,len,l)
        L[i]=ans[1]
        R[i]=ans[2]
        d[i]=1
        count.obs = count.obs + 1
      }# end of inner else
    }# end of outer else
  }# end of for
  print("count.cure")
  print(count.cure)
  print("count.obs")
  print(count.obs)
  return(data.frame(L,R,d,x1,x2))
}
# calling the function
# g1 should be negative and g2 should be positive
# b1 is negative and b2 is positive
data = data_0_BC
(n=200,b0=0.6,b1=-1.5,b2=0.1,g1=-1.2,g2=0.1,lam=0.5,cenrate=0.2) 
\end{lstlisting}}


\subsubsection*{A4.1.2. \underline{Case $\alpha \in (0,1]$}}

\begin{lstlisting}[language=R]
LR_int=function(y1,len1,l1){
  if(y1>0 & y1<=l1){
    a=c(.Machine$double.eps,l1)
  }else{
    k=as.integer((y1-l1)/len1)+1
    a=c(l1+((k-1)*len1),l1+(k*len1))
  }
  return(a)
}
data_gen_BC=function(n,alpha,b0,b1,b2,g1,g2,lam,cenrate){ 
  # alpha is the BC index parameter, lam is baseline exponential parameter 
  L=rep(NA,n)
  R=rep(NA,n)
  d=rep(NA,n)
  x1=rbinom(n=n,size=1,prob=0.5) # binary covariate
  x2=runif(n,min = 0.1,max = 20) # continuous covariate
  phi=exp(b0+(b1*x1)+(b2*x2))/(1+(alpha*exp(b0+(b1*x1)+(b2*x2))))
  U=runif(n,min=0,max=1)
  C=rexp(n,rate=cenrate)
  p0 = (1-(alpha*phi))^(1/alpha)
  count.obs=0
  count.cure=0
  for(i in 1:n){
    if(U[i]<=p0[i]){
      L[i]=C[i]
      R[i]=Inf
      d[i]=0
      count.cure=count.cure+1
    }else{
      U1 = runif(1,min=0,max=1)
      y=-exp(-(g1*x1[i])-(g2*x2[i]))
      *(1/lam)
      *log(((alpha*phi[i])+((p0[i]+((1-p0[i])*U1))^alpha)-1)/(alpha*phi[i]))
      t=min(y,C[i])
      if(t==C[i]){    
        L[i]=C[i]
        R[i]=Inf
        d[i]=0
      }else{
        len=runif(1,0.2,0.7)
        l=runif(1,0,1)
        ans=LR_int(t,len,l)
        L[i]=ans[1]
        R[i]=ans[2]
        d[i]=1
        count.obs = count.obs + 1
      }# end of inner else
    }# end of outer else
  }# end of for
  print("count.cure")
  print(count.cure)
  print("count.obs")
  print(count.obs)
  return(data.frame(L,R,d,x1,x2))
}
# calling the function


data = data_gen_BC
(n=200,alpha=0.75,b0=0.6,b1=-1.5,b2=0.1,g1=-1.2,g2=0.1,lam=0.5,cenrate=0.2)