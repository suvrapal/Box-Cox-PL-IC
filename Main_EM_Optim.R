###################################################
#######   Calling necessary libraries
###################################################
library(optimx)
library(survival)
library(nloptr)
library(numDeriv)
library(MASS)
library(pracma)
library(matrixcalc)

#############################################
#########  Data pre-processing 
#############################################
data2=data ### Reallocation data so as to keep original generated data intact

### Separating the censored data from the uncensored data
inc=which(data2$d==1)
ic=which(data2$d==0)

### Sorting time points excluding the 
### "Inf" time points since this will be required for the choice of cut points
timepts=sort(unique(c(data2$L,data2$R)))
timn=length(timepts)
timepts=timepts[-timn]
plot(timepts,type="l")
hist(timepts)
data3=data2 ### Reallocation data so as to keep data2 intact

### Some more pre-processing
data3$Timept1=data3$L
data3$Timept2=data3$R
data3=data3[,c(4,5,6,7)]

### Defining the X and Z vectors
Xvec=as.matrix(data3[,c(1,2)])
Zvec=cbind(1,Xvec)


######################################
####  Some intermediate functions
######################################

cush=function(x){ 
  for(i in 1:length(x)){
    if(x[i]<=10^(-15)){x[i]=10^(-15)
    }
  }
  return(x)
}### Replaces values less than 10^(-15) with 10^(-15) 


cush3=function(x,imp){
  lenx=length(x);val=rep(NA,lenx)
  for(i in 1:lenx){
    if(x[i]>=0){val[i]=runif(1,(1-imp)*x[i],(1+imp)*x[i])}
    else{val[i]=runif(1,(1+imp)*x[i],(1-imp)*x[i])}
  }
  return(val)
}### Used for generating initial parameter values within a given interval


###########################################################
### Choice of cut points and hazard values at those points
###########################################################
ctppsi_func=function(num_lines){
  ctp_choice=c(0,quantile(timepts,c(0,0.25,0.33,0.5,0.66,0.75,0.95)))
  psi_choice=c(0.01,0.01,0.01,0.02,0.02,0.1,0.3,1.0)
  if(num_lines==1){ctp_cons=ctp_choice[c(1,8)];
  psi_cons=psi_choice[c(1,8)]}
  if(num_lines==2){ctp_cons=ctp_choice[c(1,5,8)];
  psi_cons=psi_choice[c(1,6,8)]}
  if(num_lines==3){ctp_cons=ctp_choice[c(1,2,5,8)];
  psi_cons=psi_choice[c(1,3,6,8)]}
  if(num_lines==4){ctp_cons=ctp_choice[c(1,3,5,6,8)];
  psi_cons=psi_choice[c(1,3,5,6,8)]}
  #if(num_lines==5){ctp_cons=ctp_choice[c(1,2,3,4,5,8)];
  psi_cons=psi_choice[c(1,2,3,4,5,8)]}
#if(num_lines==6){ctp_cons=ctp_choice[c(1,2,3,4,5,6,8)];
psi_cons=psi_choice[c(1,2,3,4,5,6,8)]}
return(c(ctp_cons,psi_cons))
}

###################################
#### Initial value allocation  
###################################

### initial values for regression parameter related to Z
betinit=round(cush3(bt,imp),3) 

### initial values for regression parameter related to X
gaminit=round(cush3(gt,imp),3) 


######### initial values related to  shape and scale 
######### for Weibull baseline of susceptible dist
laminit2=round(cush3(c(1,lamt),imp),3)


################################################################
###########   All necessary user-defined functions 
################################################################

################################################################################
##########  Defining optimization functions and calculation of std errors
################################################################################
pla_optim=function(num_lines){
  
  ### Cut points and initial estimates of hazard
  chopts=ctppsi_func(num_lines)[1:(num_lines+1)]
  
  ### First one hazard and 2 onward are incremental hazards
  psinit=ctppsi_func(num_lines)[(num_lines+2):(2*num_lines+2)]  
  
  thetain=round(c(psinit,gaminit,betinit),3)
  
  ### The main PLA function
  pla=function(tt,psi,ctp){
    NN=length(psi)
    hj=rep(0,NN)
    cj=rep(0,NN)
    fj=rep(0,NN)
    sj=rep(0,NN)
    ej=rep(0,NN)
    II=rep(0,NN)
    
    for(ii in 2:NN){
      sj[ii]=(psi[ii]-psi[ii-1])/(ctp[ii]-ctp[ii-1])
      cj[ii]=psi[ii]-ctp[ii]*sj[ii]
      hj[ii]=cj[ii]+sj[ii]*tt
      ej[ii]=cj[ii]*(min(tt,ctp[ii])-ctp[ii-1])
      fj[ii]=ej[ii]+(sj[ii]/2)*((min(tt,ctp[ii]))^2-(ctp[ii-1])^2)
    }
    
    for(ii in 2:NN){
      if(tt>=ctp[ii-1]&tt<ctp[ii]){II[ii]=1}
    }
    
    iii=which(II==1)
    
    if(is.integer(iii) && length(iii) == 0){
      hval=psi[NN]
      Hval=sum(fj[1:NN])
      if(Hval<0){Hval=0}
      Sval=exp(-Hval)
      Fval=1-Sval
      fval=hval*Sval
      return(c(hval,Hval,Sval,Fval,fval))
    }
    else{    
      hval=hj[iii]
      Hval=sum(fj[1:iii])
      if(Hval<0){Hval=0}
      Sval=exp(-Hval)
      Fval=1-Sval
      fval=hval*Sval
      return(c(hval,Hval,Sval,Fval,fval))
    }
  }
  
  
  ### Complete data log-likelihood function with 
  ### baseline hazard being estimated by PLA
  
  comloglik=function(theta){
    psidf1=theta[1:(length(psinit))];
    psidf=rep(NA,length(psidf1));
    psidf[1]=(psidf1[1])^2
    for(uu in 2:length(psidf1)){
      psidf[uu]=(psidf1[uu])^2+psidf[uu-1]
    }
    gamth=theta[(length(psinit)+1):(length(psinit)+length(gaminit))];
    beth=theta[(length(psinit)+length(gaminit)+1):(length(theta))]
    valmatl=matrix(NA,nrow(data3),5);valmatr=matrix(NA,nrow(data3),5)
    for(ii in 1:nrow(data3)){
      valmatl[ii,]=pla(data3$Timept1[ii],psidf,chopts)
      valmatr[ii,]=pla(data3$Timept2[ii],psidf,chopts)
    }
    #### alpha=0
    if(alp==0){
      thfunvec=rep(NA,nrow(data3))
      Spvecl=rep(NA,nrow(data3))
      Spvecr=rep(NA,nrow(data3))
      p0vec=rep(NA,nrow(data3))
      Suvecl=rep(NA,nrow(data3))  
      for(ii in 1:nrow(data3)){
        thfunvec[ii]=exp((Zvec[ii,])%*%(beth))
        Spvecl[ii]=exp(-thfunvec[ii]*valmatl[ii,4])
        Spvecr[ii]=exp(-thfunvec[ii]*valmatr[ii,4])
        p0vec[ii]=exp(-thfunvec[ii])
        Suvecl[ii]=(Spvecl[ii]-p0vec[ii])/(1-p0vec[ii])
      }
    }
    #### alpha in (0,1]
    if(alp>0&&alp<=1){
      thfunvec=rep(NA,nrow(data3))
      intml=rep(NA,nrow(data3))
      intmr=rep(NA,nrow(data3))
      Spvecl=rep(NA,nrow(data3))
      Spvecr=rep(NA,nrow(data3))
      p0vec=rep(NA,nrow(data3))
      Suvecl=rep(NA,nrow(data3))
      for(ii in 1:nrow(data3)){
        thfunvec[ii]=exp((Zvec[ii,])%*%(beth))/(1+(alp*exp(Zvec[ii,]%*%(beth))))
        intml[ii]=1-((valmatl[ii,3])^(exp(Xvec[ii,]%*%gamth)))
        intmr[ii]=1-((valmatr[ii,3])^(exp(Xvec[ii,]%*%gamth)))
        Spvecl[ii]=(1-(alp*(thfunvec[ii])*intml[ii]))^(1/alp)
        Spvecr[ii]=(1-(alp*(thfunvec[ii])*intmr[ii]))^(1/alp)
        p0vec[ii]=(1-alp*(thfunvec[ii]))^(1/alp)
        Suvecl[ii]=(Spvecl[ii]-p0vec[ii])/(1-p0vec[ii])
      }
    }
    finval=sum(log(cush(Spvecl[inc]-Spvecr[inc])))
    +sum((1-w0)*log(cush(p0vec[ic])))
    +sum(w0*log(cush((1-p0vec[ic])*Suvecl[ic])))
    return(-finval)
  }
  
  
  
  ### Observed data log-likelihood function with 
  ### baseline hazard being estimated by PLA
  loglikf=function(theta){
    psidf1=theta[1:(length(psinit))];
    psidf=rep(NA,length(psidf1));
    psidf[1]=(psidf1[1])^2
    for(uu in 2:length(psidf1)){
      psidf[uu]=(psidf1[uu])^2+psidf[uu-1]
    }
    gamth=theta[(length(psinit)+1):(length(psinit)+length(gaminit))];
    beth=theta[(length(psinit)+length(gaminit)+1):(length(theta))]
    valmatl=matrix(NA,nrow(data3),5);valmatr=matrix(NA,nrow(data3),5) 
    for(ii in 1:nrow(data3)){
      valmatl[ii,]=pla(data3$Timept1[ii],psidf,chopts)
      valmatr[ii,]=pla(data3$Timept2[ii],psidf,chopts)
    }
    #### alpha=0
    if(alp==0){
      thfunvec=rep(NA,nrow(data3))
      Spvecl=rep(NA,nrow(data3))
      Spvecr=rep(NA,nrow(data3))
      p0vec=rep(NA,nrow(data3))
      Suvecl=rep(NA,nrow(data3))  
      for(ii in 1:nrow(data3)){
        thfunvec[ii]=exp((Zvec[ii,])%*%(beth))
        Spvecl[ii]=exp(-thfunvec[ii]*valmatl[ii,4])
        Spvecr[ii]=exp(-thfunvec[ii]*valmatr[ii,4])
        p0vec[ii]=exp(-thfunvec[ii])
        Suvecl[ii]=(Spvecl[ii]-p0vec[ii])/(1-p0vec[ii])
      }
    }
    #### alpha in (0,1]
    if(alp>0&&alp<=1){
      thfunvec=rep(NA,nrow(data3))
      intml=rep(NA,nrow(data3))
      intmr=rep(NA,nrow(data3))
      Spvecl=rep(NA,nrow(data3))
      Spvecr=rep(NA,nrow(data3))
      p0vec=rep(NA,nrow(data3))
      Suvecl=rep(NA,nrow(data3))
      for(ii in 1:nrow(data3)){
        thfunvec[ii]=exp((Zvec[ii,])%*%(beth))/(1+(alp*exp(Zvec[ii,]%*%(beth))))
        intml[ii]=1-((valmatl[ii,3])^(exp(Xvec[ii,]%*%gamth)))
        intmr[ii]=1-((valmatr[ii,3])^(exp(Xvec[ii,]%*%gamth)))
        Spvecl[ii]=(1-(alp*(thfunvec[ii])*intml[ii]))^(1/alp)
        Spvecr[ii]=(1-(alp*(thfunvec[ii])*intmr[ii]))^(1/alp)
        p0vec[ii]=(1-alp*(thfunvec[ii]))^(1/alp)
        Suvecl[ii]=(Spvecl[ii]-p0vec[ii])/(1-p0vec[ii])
      }
    }   
    finval=sum(log(cush((Spvecl[inc]-Spvecr[inc]))))+sum(log(cush(Spvecl[ic])))
    return(finval)
  }
  
  
  ### Observed data log-likelihood function for Hessian  
  loglikfH=function(theta){   
    psidf=theta[1:(length(psinit))];
    gamth=theta[(length(psinit)+1):(length(psinit)+length(gaminit))];
    beth=theta[(length(psinit)+length(gaminit)+1):(length(theta))]
    valmatl=matrix(NA,nrow(data3),5);valmatr=matrix(NA,nrow(data3),5)
    for(ii in 1:nrow(data3)){
      valmatl[ii,]=pla(data3$Timept1[ii],psidf,chopts)
      valmatr[ii,]=pla(data3$Timept2[ii],psidf,chopts)
    }  
    #### alpha=0
    if(alp==0){
      thfunvec=rep(NA,nrow(data3))
      Spvecl=rep(NA,nrow(data3))
      Spvecr=rep(NA,nrow(data3))
      p0vec=rep(NA,nrow(data3))
      Suvecl=rep(NA,nrow(data3))  
      for(ii in 1:nrow(data3)){
        thfunvec[ii]=exp((Zvec[ii,])%*%(beth))
        Spvecl[ii]=exp(-thfunvec[ii]*valmatl[ii,4])
        Spvecr[ii]=exp(-thfunvec[ii]*valmatr[ii,4])
        p0vec[ii]=exp(-thfunvec[ii])
        Suvecl[ii]=(Spvecl[ii]-p0vec[ii])/(1-p0vec[ii])
      }
    }
    #### alpha in (0,1]
    if(alp>0&&alp<=1){
      thfunvec=rep(NA,nrow(data3))
      intml=rep(NA,nrow(data3))
      intmr=rep(NA,nrow(data3))
      Spvecl=rep(NA,nrow(data3))
      Spvecr=rep(NA,nrow(data3))
      p0vec=rep(NA,nrow(data3))
      Suvecl=rep(NA,nrow(data3))
      for(ii in 1:nrow(data3)){
        thfunvec[ii]=exp((Zvec[ii,])%*%(beth))/(1+(alp*exp(Zvec[ii,]%*%(beth))))
        intml[ii]=1-((valmatl[ii,3])^(exp(Xvec[ii,]%*%gamth)))
        intmr[ii]=1-((valmatr[ii,3])^(exp(Xvec[ii,]%*%gamth)))
        Spvecl[ii]=(1-(alp*(thfunvec[ii])*intml[ii]))^(1/alp)
        Spvecr[ii]=(1-(alp*(thfunvec[ii])*intmr[ii]))^(1/alp)
        p0vec[ii]=(1-alp*(thfunvec[ii]))^(1/alp)
        Suvecl[ii]=(Spvecl[ii]-p0vec[ii])/(1-p0vec[ii])
      }
    }   
    finval=sum(log(cush((Spvecl[inc]-Spvecr[inc]))))+sum(log(cush(Spvecl[ic])))
    return(finval)
  }
  
  ########################################
  ### Main Loop and optimization #########
  ########################################
  nalr=length(alpr)
  ESTV=matrix(NA,nalr,length(thetain))
  LLV=rep(NA,nalr)
  
  for (jj in 1:nalr){
    alp=alpr[jj];
    ctt=0;eps0=0.01;eps=10;
    theta0=thetain
    while(eps>eps0){
      psidf1=theta0[1:(length(psinit))];
      psidf=rep(NA,length(psidf1));
      psidf[1]=(psidf1[1])^2
      for(uu in 2:length(psidf1)){
        psidf[uu]=(psidf1[uu])^2+psidf[uu-1]
      }
      gamth=theta0[(length(psinit)+1):(length(psinit)+length(gaminit))];
      beth=theta0[(length(psinit)+length(gaminit)+1):(length(theta0))]
      valmatl=matrix(NA,nrow(data3),5);valmatr=matrix(NA,nrow(data3),5)
      for(ii in 1:nrow(data3)){
        valmatl[ii,]=pla(data3$Timept1[ii],psidf,chopts)
        valmatr[ii,]=pla(data3$Timept2[ii],psidf,chopts)
      }
      #### alpha=0
      if(alp==0){
        thfunvec=rep(NA,nrow(data3))
        Spvecl=rep(NA,nrow(data3))
        Spvecr=rep(NA,nrow(data3))
        p0vec=rep(NA,nrow(data3))
        Suvecl=rep(NA,nrow(data3))  
        for(ii in 1:nrow(data3)){
          thfunvec[ii]=exp((Zvec[ii,])%*%(beth))
          Spvecl[ii]=exp(-thfunvec[ii]*valmatl[ii,4])
          Spvecr[ii]=exp(-thfunvec[ii]*valmatr[ii,4])
          p0vec[ii]=exp(-thfunvec[ii])
          Suvecl[ii]=(Spvecl[ii]-p0vec[ii])/(1-p0vec[ii])
        }
      }
      #### alpha in (0,1]
      if(alp>0&&alp<=1){
        thfunvec=rep(NA,nrow(data3))
        intml=rep(NA,nrow(data3))
        intmr=rep(NA,nrow(data3))
        Spvecl=rep(NA,nrow(data3))
        Spvecr=rep(NA,nrow(data3))
        p0vec=rep(NA,nrow(data3))
        Suvecl=rep(NA,nrow(data3))
        for(ii in 1:nrow(data3)){
          thfunvec[ii]=exp((Zvec[ii,])%*%(beth))/(1+(alp*exp(Zvec[ii,]%*%(beth))))
          intml[ii]=1-((valmatl[ii,3])^(exp(Xvec[ii,]%*%gamth)))
          intmr[ii]=1-((valmatr[ii,3])^(exp(Xvec[ii,]%*%gamth)))
          Spvecl[ii]=(1-(alp*(thfunvec[ii])*intml[ii]))^(1/alp)
          Spvecr[ii]=(1-(alp*(thfunvec[ii])*intmr[ii]))^(1/alp)
          p0vec[ii]=(1-alp*(thfunvec[ii]))^(1/alp)
          Suvecl[ii]=(Spvecl[ii]-p0vec[ii])/(1-p0vec[ii])
        }
      }
      w0=((1-p0vec[ic])*Suvecl[ic])/Spvecl[ic]
      obj1=optim(theta0,comloglik,method="Nelder-Mead")
      theta1=obj1$par
      eps=sqrt(sum((theta1-theta0)^2))
      #eps=max(abs((theta1-theta0)/theta0))
      theta0=theta1;
      ctt=ctt+1;
    }
    ESTV[jj,]=theta0
    LLV[jj]=loglikf(theta0)  
  }
  
  ########################################
  ######   Var-Cov matrix ################
  ########################################
  idx=which.max(LLV)
  maxl=LLV[idx]
  estm=ESTV[idx,]
  alp=alpr[idx]
  estm2=estm
  estm2[1]=(estm[1])^2
  for(uu in 2:(num_lines+1)){
    estm2[uu]=(estm[uu])^2+estm2[uu-1]
  }
  estm3=round(estm2,4)
  H1=-(hessian(loglikfH,estm3))
  H2=(pinv(H1))
  hdiag=diag(H2)
  seval=round(sqrt(hdiag),4)
  lcl=estm3-1.96*seval
  ucl=estm3+1.96*seval
  retval=list()
  retval[[1]]=estm3
  retval[[2]]=seval
  retval[[3]]=lcl
  retval[[4]]=ucl
  retval[[5]]=c(alp,maxl)
  return(retval)
}

#############################################
###### Calling the main function ############
#############################################

### The function call below returns parameter estimates and associated SEs 
### with the estimating technique discussed in the main article where 
### the baseline hazard function is approximated by 3 lines 

pla_optim(num_lines=3) 