###########Simulation for SRAB,SBB, SRAH,SBH procedures for bernoulli distribution#############

alpha=.05                      ###Fixing the value of alpha
beta=.1                        ###Fixing the value of beta
mu1=.3                         ###mean of control
mu2=.4                         ###mean of experimental treatment 1
mu3=.4                        ###mean of experimental treatment 2
mu4=.6                         ###mean of experimental treatment 3
sigma1=sqrt(mu1*(1-mu1))       ###standard deviation of control
sigma2=sqrt(mu2*(1-mu2))       ###standard deviation of experimental treatment 1
sigma3=sqrt(mu3*(1-mu3))       ###standard deviation of experimental treatment 2
sigma4=sqrt(mu4*(1-mu4))       ###standard deviation of experimental treatment 3


 
delta0=0.1                   ##Fixing the value of delta0
delta1=0.3                   ## Fixing the value of delta1

delta=delta1-delta0         ##Defining delta

n0=5                       ##Fixing pilot sample size 'K'
k=(delta0*qnorm(1-beta/3)+delta1*qnorm(1-alpha/3))/(qnorm(1-alpha/3)+qnorm(1-beta/3))      ##cutoff for bonferroni's procedure

iterations = 1000        ##Number of iterations         
###Parallelization
library(parallel)        
library(doParallel)
library(foreach)
nc=detectCores()
cl=makeCluster(nc-2)
registerDoParallel(cl)
###Simulation code for the SRAB procedure
resultbonopt=foreach (i=1:iterations,.combine=rbind)%dopar%{  
  set.seed(i)
  x1=c(1,0,rbinom(n0,1,mu1))  ######generating pilot observations from control, initial 1,0 added to get sample mean in (0,1) 
  x4=c(1,0,rbinom(n0,1,mu4))  ###generating pilot observations from experimental treatment 3, initial 1,0 added to get sample mean in (0,1) 
  x2=c(1,0,rbinom(n0,1,mu2))  ###generating pilot observations from experimental treatment 1, initial 1,0 added to get sample mean in (0,1) 
  x3=c(1,0,rbinom(n0,1,mu3))  ###generating pilot observations from experimental treatment 2, initial 1,0 added to get sample mean in (0,1) 
  if((mean(x1)*(1-mean(x1))/length(x1))+(mean(x2)*(1-mean(x2))/length(x2))<=(delta^2)/(qnorm(1-alpha/3)+qnorm(1-beta/3))^2 
     && (mean(x1)*(1-mean(x1))/length(x1))+(mean(x3)*(1-mean(x3))/length(x3))<=(delta^2)/(qnorm(1-alpha/3)+qnorm(1-beta/3))^2
     && (mean(x1)*(1-mean(x1))/length(x1))+(mean(x4)*(1-mean(x4))/length(x4))<=(delta^2)/(qnorm(1-alpha/3)+qnorm(1-beta/3))^2){
    stopcrmet=TRUE   ###checking for stopping criteria 
  }else {
    stopcrmet=FALSE
  }
  ###running a while loop until the stopping criteria is reached
  while(!stopcrmet){
    pi1=sqrt(mean(x1)*(1-mean(x1))/(1-mean(x1)))*
      sqrt(mean(x2)*(1-mean(x2))*(1-mean(x2))+
             mean(x3)*(mean(x3))*(1-mean(x3))+
             mean(x4)*(mean(x4))*(1-mean(x4)))
    pi2=mean(x2)*(1-mean(x2))
    pi3=mean(x3)*(1-mean(x3))
    pi4=mean(x4)*(1-mean(x4))
    
    pi11=pi1/(pi1+pi2+pi3+pi4)
    pi22=pi2/(pi1+pi2+pi3+pi4)
    pi33=pi3/(pi1+pi2+pi3+pi4)
    pi44=pi4/(pi1+pi2+pi3+pi4)
    p=c(pi11,pi22,pi33,pi44)
    
    u=sample(c(1,2,3,4),1,prob=p)
    if(u==1){
      x1=c(x1,rbinom(1,1,mu1))
      x2=x2
      x3=x3
      x4=x4
    }else if(u==2){
      x1=x1
      x2=c(x2,rbinom(1,1,mu2))
      x3=x3
      x4=x4
    } else if (u==3) {
      x2=x2
      x1=x1
      x3=c(x3,rbinom(1,1,mu3))
      x4=x4
    } else {
      x2=x2
      x1=x1
      x3=x3
      x4=c(x4,rbinom(1,1,mu4))
    }
    if((mean(x1)*(1-mean(x1))/length(x1))+(mean(x2)*(1-mean(x2))/length(x2))<=(delta^2)/(qnorm(1-alpha/3)+qnorm(1-beta/3))^2 
       && (mean(x1)*(1-mean(x1))/length(x1))+(mean(x3)*(1-mean(x3))/length(x3))<=(delta^2)/(qnorm(1-alpha/3)+qnorm(1-beta/3))^2
       && (mean(x1)*(1-mean(x1))/length(x1))+(mean(x4)*(1-mean(x4))/length(x4))<=(delta^2)/(qnorm(1-alpha/3)+qnorm(1-beta/3))^2){
      stopcrmet=TRUE ###checking for stopping criteria
    }
    
  }
  N=length(x1)+length(x2)+length(x3)+length(x4)           ###calculating the sample size          
  TF=sum(x1==0)+sum(x2==0)+sum(x3==0)+sum(x4==0)          ###calculating the number of treatment failures           
  ###making decisions for three experimental treatments
  if(mean(x2)-mean(x1)>k){
    r2=1
    a2=0
  }else{
    r2=0
    a2=1
  }
  if(mean(x3)-mean(x1)>k){
    r3=1
    a3=0
  }else{
    r3=0
    a3=1
  }
  if(mean(x4)-mean(x1)>k){
    r4=1
    a4=0
  }else{
    r4=0
    a4=1
  }
  c(N,TF,max(r2,r3),max(a4))  ### change it according to the definition of type-I and type-II errors from true means. 
  ####Third component corresponds to type-I error and fourth component corresponds to type-II error. Note that for global null there is  no type-II error.
}
##Results for SRAB procedure
ASN=colMeans(resultbonopt)[1]
ANF=colMeans(resultbonopt)[2]
FWER_I=colMeans(resultbonopt)[3]
FWER_II=colMeans(resultbonopt)[4]
se_samplesize=sd(resultbonopt[,1])/sqrt(iterations)
se_treatmentfailures=sd(resultbonopt[,2])/sqrt(iterations)

#####################Simulation for SBB Procedure#################
stopCluster(cl)
nc=detectCores()
cl=makeCluster(nc-2)
registerDoParallel(cl)

resultbonbal=foreach(i=1:iterations,.combine=rbind)%dopar%{  
  set.seed(i)
  
  x1=c(1,0,rbinom(n0,1,mu1))              ###generating pilot observations from control, initial 1,0 added to get sample mean in (0,1) 
  x4=c(1,0,rbinom(n0,1,mu4))              ###generating pilot observations from experimental treatment 3, initial 1,0 added to get sample mean in (0,1) 
  x2=c(1,0,rbinom(n0,1,mu2))              ###generating pilot observations from experimental treatment 1, initial 1,0 added to get sample mean in (0,1) 
  x3=c(1,0,rbinom(n0,1,mu3))              ###generating pilot observations from experimental treatment 2, initial 1,0 added to get sample mean in (0,1) 
  p=c(1/4,1/4,1/4,1/4)                   ###equal randomization
  
  if((mean(x1)*(1-mean(x1))/length(x1))+(mean(x2)*(1-mean(x2))/length(x2))<=(delta^2)/(qnorm(1-alpha/3)+qnorm(1-beta/3))^2 
     && (mean(x1)*(1-mean(x1))/length(x1))+(mean(x3)*(1-mean(x3))/length(x3))<=(delta^2)/(qnorm(1-alpha/3)+qnorm(1-beta/3))^2
     && (mean(x1)*(1-mean(x1))/length(x1))+(mean(x4)*(1-mean(x4))/length(x4))<=(delta^2)/(qnorm(1-alpha/3)+qnorm(1-beta/3))^2){
    stopcrmet=TRUE ##checking for stopping criteria
  }else {
    stopcrmet=FALSE
  }
  while(!stopcrmet){
    u=sample(c(1,2,3,4),1,prob=p)
    ###generating new observation and updating the information
    if(u==1){
      x1=c(x1,rbinom(1,1,mu1))
      x3=x3
      x2=x2
      x4=x4
    }else if(u==2 ){
      x3=x3
      x2=c(x2,rbinom(1,1,mu2))
      x1=x1
      x4=x4
    } else if (u==3) {
      x2=x2
      x1=x1
      x3=c(x3,rbinom(1,1,mu3))
      x4=x4
    }else  {
      x2=x2
      x1=x1
      x4=c(x4,rbinom(1,1,mu4))
      x3=x3
    }
    
    if((mean(x1)*(1-mean(x1))/length(x1))+(mean(x2)*(1-mean(x2))/length(x2))<=(delta^2)/(qnorm(1-alpha/3)+qnorm(1-beta/3))^2 
       && (mean(x1)*(1-mean(x1))/length(x1))+(mean(x3)*(1-mean(x3))/length(x3))<=(delta^2)/(qnorm(1-alpha/3)+qnorm(1-beta/3))^2
       && (mean(x1)*(1-mean(x1))/length(x1))+(mean(x4)*(1-mean(x4))/length(x4))<=(delta^2)/(qnorm(1-alpha/3)+qnorm(1-beta/3))^2){
      stopcrmet=TRUE    ###checking for stopping criteria
    }
  }
  
  N=length(x1)+length(x2)+length(x3)+length(x4)                            #######calculating total sample size
  TF=sum(x1==0)+sum(x2==0)+sum(x3==0)+sum(x4==0)                           #######calculating the number of treatment failures
  ###making decisions for three experimental treatments
  if(mean(x2)-mean(x1)>k){
    r2=1
    a2=0
  }else{
    r2=0
    a2=1
  }
  if(mean(x3)-mean(x1)>k){
    r3=1
    a3=0
  }else{
    r3=0
    a3=1
  }
  if(mean(x4)-mean(x1)>k){
    r4=1
    a4=0
  }else{
    r4=0
    a4=1
  }
  
  c(N,TF,max(r2,r3),max(a4)) ### change it according to the definition of type-I and type-II errors from true means. 
  ####Third component corresponds to type-I error and fourth component corresponds to type-II error.Note that for global null there is  no type-II error.
}
##Results for SBB procedure
ASN=colMeans(resultbonbal)[1]
ANF=colMeans(resultbonbal)[2]
FWER_I=colMeans(resultbonbal)[3]
FWER_II=colMeans(resultbonbal)[4]
se_samplesize=sd(resultbonbal[,1])/sqrt(iterations)
se_treatmentfailures=sd(resultbonbal[,2])/sqrt(iterations)


#######################Simulation for the SRAH procedure##########

n0=5                              ###Fixing the value of pilot sample size 'K'

iterations = 1000
library(parallel)
library(doParallel)
library(foreach)
nc=detectCores()
cl=makeCluster(nc-2)
registerDoParallel(cl)

resultholmopt=foreach (i=1:iterations,.combine=rbind,.packages = c("NlcOptim","MASS"))%dopar%{  
  
  set.seed(i)
  x1=c(1,0,rbinom(n0,1,mu1))                     ###generating pilot observations from control, initial 1,0 added to get sample mean in (0,1)
  x4=c(1,0,rbinom(n0,1,mu4))                     ###generating pilot observations from experimental treatment 3, initial 1,0 added to get sample mean in (0,1) 
  x2=c(1,0,rbinom(n0,1,mu2))                     ###generating pilot observations from experimental treatment 1, initial 1,0 added to get sample mean in (0,1)
  x3=c(1,0,rbinom(n0,1,mu3))                     ###generating pilot observations from experimental treatment 2, initial 1,0 added to get sample mean in (0,1)
  d=sort(c(mean(x2),mean(x3),mean(x4)),decreasing = TRUE)       ###ordering of sample means
  ##breking ties for equal means
  if(mean(x4)!=mean(x2) && mean(x2)!=mean(x3) && mean(x4)!=mean(x3)){
    d2=which(d==mean(x2))
    d3=which(d==mean(x3))
    d4=which(d==mean(x4))
  } else if (mean(x4)==mean(x2) && mean(x3)>mean(x4)){
    d3=1
    u=1+rank(c(mean(x4),mean(x2)),ties.method= "random")
    d2=u[1]
    d4=u[2]
  }else if (mean(x4)==mean(x2) && mean(x4)>mean(x3)){
    d3=3
    u=rank(c(mean(x2),mean(x4)),ties.method= "random")
    d2=u[1]
    d4=u[2]
  }else if (mean(x3)==mean(x4) && mean(x2)>mean(x3)){
    d2=1
    u=1+rank(c(mean(x3),mean(x4)),ties.method= "random")
    d3=u[1]
    d4=u[2]
  }else if (mean(x3)==mean(x4) && mean(x2)<mean(x3)){
    d2=3
    u=rank(c(mean(x3),mean(x4)),ties.method= "random")
    d3=u[1]
    d4=u[2]
  }else if (mean(x2)==mean(x3) && mean(x2)>mean(x4)){
    d4=3
    u=rank(c(mean(x2),mean(x3)),ties.method= "random")
    d2=u[1]
    d3=u[2]
  }else if (mean(x3)==mean(x2) && mean(x4)>mean(x2)){
    d4=1
    u=1+rank(c(mean(x3),mean(x2)),ties.method= "random")
    d2=u[1]
    d3=u[2]
  }else{
    u=rank(c(mean(x3),mean(x2),mean(x1)),ties.method= "random")
    d4=u[1]
    d2=u[2]
    d3=u[3]
  }
  
  
  
  if((mean(x1)*(1-mean(x1))/length(x1))+(mean(x2)*(1-mean(x2))/length(x2))<=delta^2/(qnorm(1-alpha/(3-d2+1))+
                                                                                     qnorm(1-beta/d2))^2
     
     && (mean(x1)*(1-mean(x1))/length(x1))+(mean(x3)*(1-mean(x3))/length(x3))<=delta^2/(qnorm(1-alpha/(3-d3+1))+
                                                                                        qnorm(1-beta/d3))^2
     
     && (mean(x1)*(1-mean(x1))/length(x1))+(mean(x4)*(1-mean(x4))/length(x4))<=delta^2/(qnorm(1-alpha/(3-d4+1))+
                                                                                        qnorm(1-beta/d4))^2){
    stopcrmet=TRUE ###checking for stopping criteria
  }else{
    stopcrmet=FALSE
  }
  while(!stopcrmet){
    d=sort(c(mean(x2),mean(x3),mean(x4)),decreasing = TRUE)                   ###ordering of sample means
    ###breaking ties for equal means
    if(mean(x4)!=mean(x2) && mean(x2)!=mean(x3) && mean(x4)!=mean(x3)){
      d2=which(d==mean(x2))
      d3=which(d==mean(x3))
      d4=which(d==mean(x4))
    } else if (mean(x4)==mean(x2) && mean(x3)>mean(x4)){
      d3=1
      u=1+rank(c(mean(x4),mean(x2)),ties.method= "random")
      d2=u[1]
      d4=u[2]
    }else if (mean(x4)==mean(x2) && mean(x4)>mean(x3)){
      d3=3
      u=rank(c(mean(x2),mean(x4)),ties.method= "random")
      d2=u[1]
      d4=u[2]
    }else if (mean(x3)==mean(x4) && mean(x2)>mean(x3)){
      d2=1
      u=1+rank(c(mean(x3),mean(x4)),ties.method= "random")
      d3=u[1]
      d4=u[2]
    }else if (mean(x3)==mean(x4) && mean(x2)<mean(x3)){
      d2=3
      u=rank(c(mean(x3),mean(x4)),ties.method= "random")
      d3=u[1]
      d4=u[2]
    }else if (mean(x2)==mean(x3) && mean(x2)>mean(x4)){
      d4=3
      u=rank(c(mean(x2),mean(x3)),ties.method= "random")
      d2=u[1]
      d3=u[2]
    }else if (mean(x3)==mean(x2) && mean(x4)>mean(x2)){
      d4=1
      u=1+rank(c(mean(x3),mean(x2)),ties.method= "random")
      d2=u[1]
      d3=u[2]
    }else{
      u=rank(c(mean(x3),mean(x2),mean(x1)),ties.method= "random")
      d4=u[1]
      d2=u[2]
      d3=u[3]
    }
    
    ####defining objective function
    obj=function(x){
      return(x[1]*(1-mean(x1))+x[2]*(1-mean(x2))
             +x[3]*(1-mean(x3))+x[4]*(1-mean(x4)))
    }
    ##defining constraints
    budgetconstraint<-function(x){
      f=NULL
      f=rbind(f,mean(x2)*(1-mean(x2))/x[2]+mean(x1)*(1-mean(x1))/x[1]-delta^2/(qnorm(1-alpha/(3-d2+1))+
                                                                                 qnorm(1-beta/d2))^2)
      f=rbind(f,mean(x3)*(1-mean(x3))/x[3]+mean(x1)*(1-mean(x1))/x[1]-delta^2/(qnorm(1-alpha/(3-d3+1))+
                                                                                 qnorm(1-beta/d3))^2)
      f= rbind(f,mean(x4)*(1-mean(x4))/x[4]+mean(x1)*(1-mean(x1))/x[1]-delta^2/(qnorm(1-alpha/(3-d4+1))+
                                                                                  qnorm(1-beta/d4))^2)
      f=rbind(f,1-x[1])
      f=rbind(f,1-x[2])
      f=rbind(f,1-x[3])
      f=rbind(f,1-x[4])
      
      return(list(ceq = NULL, c = f))
    }
    
    x0=c(10,10,10,10)            ####initial value for numerical optimization
    a=solnl(x0,obj,budgetconstraint,tolX = 1e-2,
            tolFun = 1e-5, tolCon = 1e-5, maxnFun = 1e+10,maxIter = 8000)$par
    pi=a/sum(a)                ####obtaining response-adaptive allocation rule
    
    u=sample(c(1,2,3,4),1,prob=pi)      ####response-adaptive randomization
    ##generating observation and updating the information
    if(u==1){
      x1=c(x1,rbinom(1,1,mu1))
      x2=x2
      x3=x3
      x4=x4
    }else if(u==2){
      x1=x1
      x2=c(x2,rbinom(1,1,mu2))
      x3=x3
      x4=x4
    } else if (u==3) {
      x2=x2
      x1=x1
      x3=c(x3,rbinom(1,1,mu3))
      x4=x4
    } else {
      x2=x2
      x1=x1
      x3=x3
      x4=c(x4,rbinom(1,1,mu4))
      
    }
    d=sort(c(mean(x2),mean(x3),mean(x4)),decreasing = TRUE)               ### ordering the sample means
    ###breaking ties for equal means
    if(mean(x4)!=mean(x2) && mean(x2)!=mean(x3) && mean(x4)!=mean(x3)){
      d2=which(d==mean(x2))
      d3=which(d==mean(x3))
      d4=which(d==mean(x4))
    } else if (mean(x4)==mean(x2) && mean(x3)>mean(x4)){
      d3=1
      u=1+rank(c(mean(x4),mean(x2)),ties.method= "random")
      d2=u[1]
      d4=u[2]
    }else if (mean(x4)==mean(x2) && mean(x4)>mean(x3)){
      d3=3
      u=rank(c(mean(x2),mean(x4)),ties.method= "random")
      d2=u[1]
      d4=u[2]
    }else if (mean(x3)==mean(x4) && mean(x2)>mean(x3)){
      d2=1
      u=1+rank(c(mean(x3),mean(x4)),ties.method= "random")
      d3=u[1]
      d4=u[2]
    }else if (mean(x3)==mean(x4) && mean(x2)<mean(x3)){
      d2=3
      u=rank(c(mean(x3),mean(x4)),ties.method= "random")
      d3=u[1]
      d4=u[2]
    }else if (mean(x2)==mean(x3) && mean(x2)>mean(x4)){
      d4=3
      u=rank(c(mean(x2),mean(x3)),ties.method= "random")
      d2=u[1]
      d3=u[2]
    }else if (mean(x3)==mean(x2) && mean(x4)>mean(x2)){
      d4=1
      u=1+rank(c(mean(x3),mean(x2)),ties.method= "random")
      d2=u[1]
      d3=u[2]
    }else{
      u=rank(c(mean(x3),mean(x2),mean(x1)),ties.method= "random")
      d4=u[1]
      d2=u[2]
      d3=u[3]
    }
    
    if((mean(x1)*(1-mean(x1))/length(x1))+(mean(x2)*(1-mean(x2))/length(x2))<=delta^2/(qnorm(1-alpha/(3-d2+1))+
                                                                                       qnorm(1-beta/d2))^2
       
       && (mean(x1)*(1-mean(x1))/length(x1))+(mean(x3)*(1-mean(x3))/length(x3))<=delta^2/(qnorm(1-alpha/(3-d3+1))+
                                                                                          qnorm(1-beta/d3))^2
       
       && (mean(x1)*(1-mean(x1))/length(x1))+(mean(x4)*(1-mean(x4))/length(x4))<=delta^2/(qnorm(1-alpha/(3-d4+1))+
                                                                                          qnorm(1-beta/d4))^2){
      stopcrmet=TRUE       ###checking the stopping criteria
    }
    
  }
  N=length(x1)+length(x2)+length(x3)+length(x4)                                 ##calculating total sample size
  TF= sum(x1==0)+sum(x2==0)+sum(x3==0)+sum(x4==0)                               ##calculating number of treatment failures
  d=sort(c(mean(x2),mean(x3),mean(x4)),decreasing = TRUE)                       ##ordering sample means
  ###breaking ties for equal means
  if(mean(x4)!=mean(x2) && mean(x2)!=mean(x3) && mean(x4)!=mean(x3)){
    d2=which(d==mean(x2))
    d3=which(d==mean(x3))
    d4=which(d==mean(x4))
  } else if (mean(x4)==mean(x2) && mean(x3)>mean(x4)){
    d3=1
    u=1+rank(c(mean(x4),mean(x2)),ties.method= "random")
    d2=u[1]
    d4=u[2]
  }else if (mean(x4)==mean(x2) && mean(x4)>mean(x3)){
    d3=3
    u=rank(c(mean(x2),mean(x4)),ties.method= "random")
    d2=u[1]
    d4=u[2]
  }else if (mean(x3)==mean(x4) && mean(x2)>mean(x3)){
    d2=1
    u=1+rank(c(mean(x3),mean(x4)),ties.method= "random")
    d3=u[1]
    d4=u[2]
  }else if (mean(x3)==mean(x4) && mean(x2)<mean(x3)){
    d2=3
    u=rank(c(mean(x3),mean(x4)),ties.method= "random")
    d3=u[1]
    d4=u[2]
  }else if (mean(x2)==mean(x3) && mean(x2)>mean(x4)){
    d4=3
    u=rank(c(mean(x2),mean(x3)),ties.method= "random")
    d2=u[1]
    d3=u[2]
  }else if (mean(x3)==mean(x2) && mean(x4)>mean(x2)){
    d4=1
    u=1+rank(c(mean(x3),mean(x2)),ties.method= "random")
    d2=u[1]
    d3=u[2]
  }else{
    u=rank(c(mean(x3),mean(x2),mean(x1)),ties.method= "random")
    d4=u[1]
    d2=u[2]
    d3=u[3]
  }
  
  ###making decisions for three equal treatments
  
  if(mean(x2)-mean(x1)>(delta0*qnorm(1-beta/d2)+
                        delta1*qnorm(1-alpha/(3-d2+1)))/
     (qnorm(1-alpha/(3-d2+1))+
      qnorm(1-beta/d2))){
    r2=1
    a2=0
  }else{
    r2=0
    a2=1
  }
  if(mean(x3)-mean(x1)>(delta0*qnorm(1-beta/d3)+
                        delta1*qnorm(1-alpha/(3-d3+1)))/
     (qnorm(1-alpha/(3-d3+1))+
      qnorm(1-beta/d3))){
    r3=1
    a3=0
  }else{
    r3=0
    a3=1
  }
  if(mean(x4)-mean(x1)>(delta0*qnorm(1-beta/d4)+
                        delta1*qnorm(1-alpha/(3-d4+1)))/
     (qnorm(1-alpha/(3-d4+1))+
      qnorm(1-beta/d4))){
    r4=1
    a4=0
  }else{
    r4=0
    a4=1
  }
  c(N,TF,max(r2,r3),max(a4)) ### change it according to the definition of type-I and type-II errors from true means. 
  ####Third component corresponds to type-I error and fourth component corresponds to type-II error. Note that for global null there is  no type-II error.
  
}
##Results for SRAH procedure
ASN=colMeans(resultholmopt)[1]
ANF=colMeans(resultholmopt)[2]
FWER_I=colMeans(resultholmopt)[3]
FWER_II=colMeans(resultholmopt)[4]
se_samplesize=sd(resultholmopt[,1])/sqrt(iterations)
se_treatmentfailures=sd(resultholmopt[,2])/sqrt(iterations)
##################Simulation for SBH procedure######################
stopCluster(cl)
nc=detectCores()
cl=makeCluster(nc-2)
registerDoParallel(cl)

resultholmbal=foreach (i=1:iterations,.combine=rbind,.packages = c("NlcOptim","MASS"))%dopar%{  
  set.seed(i)
  x1=c(1,0,rbinom(n0,1,mu1))  ######generating pilot observations from control, initial 1,0 added to get sample mean in (0,1) 
  x4=c(1,0,rbinom(n0,1,mu4))  ###generating pilot observations from experimental treatment 3, initial 1,0 added to get sample mean in (0,1) 
  x2=c(1,0,rbinom(n0,1,mu2))  ###generating pilot observations from experimental treatment 1, initial 1,0 added to get sample mean in (0,1) 
  x3=c(1,0,rbinom(n0,1,mu3))  ###generating pilot observations from experimental treatment 2, initial 1,0 added to get sample mean in (0,1) 
  d=sort(c(mean(x2),mean(x3),mean(x4)),decreasing = TRUE)         ###ordering the sample means
  ### breaking the ties for equal means
  if(mean(x4)!=mean(x2) && mean(x2)!=mean(x3) && mean(x4)!=mean(x3)){
    d2=which(d==mean(x2))
    d3=which(d==mean(x3))
    d4=which(d==mean(x4))
  } else if (mean(x4)==mean(x2) && mean(x3)>mean(x4)){
    d3=1
    u=1+rank(c(mean(x4),mean(x2)),ties.method= "random")
    d2=u[1]
    d4=u[2]
  }else if (mean(x4)==mean(x2) && mean(x4)>mean(x3)){
    d3=3
    u=rank(c(mean(x2),mean(x4)),ties.method= "random")
    d2=u[1]
    d4=u[2]
  }else if (mean(x3)==mean(x4) && mean(x2)>mean(x3)){
    d2=1
    u=1+rank(c(mean(x3),mean(x4)),ties.method= "random")
    d3=u[1]
    d4=u[2]
  }else if (mean(x3)==mean(x4) && mean(x2)<mean(x3)){
    d2=3
    u=rank(c(mean(x3),mean(x4)),ties.method= "random")
    d3=u[1]
    d4=u[2]
  }else if (mean(x2)==mean(x3) && mean(x2)>mean(x4)){
    d4=3
    u=rank(c(mean(x2),mean(x3)),ties.method= "random")
    d2=u[1]
    d3=u[2]
  }else if (mean(x3)==mean(x2) && mean(x4)>mean(x2)){
    d4=1
    u=1+rank(c(mean(x3),mean(x2)),ties.method= "random")
    d2=u[1]
    d3=u[2]
  }else{
    u=rank(c(mean(x3),mean(x2),mean(x1)),ties.method= "random")
    d4=u[1]
    d2=u[2]
    d3=u[3]
  }
  
  
  if((mean(x1)*(1-mean(x1))/length(x1))+(mean(x2)*(1-mean(x2))/length(x2))<=delta^2/(qnorm(1-alpha/(3-d2+1))+
                                                                                     qnorm(1-beta/d2))^2
     
     && (mean(x1)*(1-mean(x1))/length(x1))+(mean(x3)*(1-mean(x3))/length(x3))<=delta^2/(qnorm(1-alpha/(3-d3+1))+
                                                                                        qnorm(1-beta/d3))^2
     
     && (mean(x1)*(1-mean(x1))/length(x1))+(mean(x4)*(1-mean(x4))/length(x4))<=delta^2/(qnorm(1-alpha/(3-d4+1))+
                                                                                        qnorm(1-beta/d4))^2){
    stopcrmet=TRUE           ###checking for stopping criteria
  }else{
    stopcrmet=FALSE
  }
  while(!stopcrmet){
    pi=c(1/4,1/4,1/4,1/4)
    u=sample(c(1,2,3,4),1,prob=pi)        ###equal randomization
    ###generating new observation and updating the information
    if(u==1){
      x1=c(x1,rbinom(1,1,mu1))
      x2=x2
      x3=x3
      x4=x4
    }else if(u==2){
      x1=x1
      x2=c(x2,rbinom(1,1,mu2))
      x3=x3
      x4=x4
    } else if (u==3) {
      x2=x2
      x1=x1
      x3=c(x3,rbinom(1,1,mu3))
      x4=x4
    } else {
      x2=x2
      x1=x1
      x3=x3
      x4=c(x4,rbinom(1,1,mu4))
      
    }
    
    d=sort(c(mean(x2),mean(x3),mean(x4)),decreasing = TRUE)        ###ordering the sample means
    ### breaking the ties for equal means
    if(mean(x4)!=mean(x2) && mean(x2)!=mean(x3) && mean(x4)!=mean(x3)){
      d2=which(d==mean(x2))
      d3=which(d==mean(x3))
      d4=which(d==mean(x4))
    } else if (mean(x4)==mean(x2) && mean(x3)>mean(x4)){
      d3=1
      u=1+rank(c(mean(x4),mean(x2)),ties.method= "random")
      d2=u[1]
      d4=u[2]
    }else if (mean(x4)==mean(x2) && mean(x4)>mean(x3)){
      d3=3
      u=rank(c(mean(x2),mean(x4)),ties.method= "random")
      d2=u[1]
      d4=u[2]
    }else if (mean(x3)==mean(x4) && mean(x2)>mean(x3)){
      d2=1
      u=1+rank(c(mean(x3),mean(x4)),ties.method= "random")
      d3=u[1]
      d4=u[2]
    }else if (mean(x3)==mean(x4) && mean(x2)<mean(x3)){
      d2=3
      u=rank(c(mean(x3),mean(x4)),ties.method= "random")
      d3=u[1]
      d4=u[2]
    }else if (mean(x2)==mean(x3) && mean(x2)>mean(x4)){
      d4=3
      u=rank(c(mean(x2),mean(x3)),ties.method= "random")
      d2=u[1]
      d3=u[2]
    }else if (mean(x3)==mean(x2) && mean(x4)>mean(x2)){
      d4=1
      u=1+rank(c(mean(x3),mean(x2)),ties.method= "random")
      d2=u[1]
      d3=u[2]
    }else{
      u=rank(c(mean(x3),mean(x2),mean(x1)),ties.method= "random")
      d4=u[1]
      d2=u[2]
      d3=u[3]
    }
    
    
    if((mean(x1)*(1-mean(x1))/length(x1))+(mean(x2)*(1-mean(x2))/length(x2))<=delta^2/(qnorm(1-alpha/(3-d2+1))+
                                                                                       qnorm(1-beta/d2))^2      &&
       (mean(x1)*(1-mean(x1))/length(x1))+(mean(x3)*(1-mean(x3))/length(x3))<=delta^2/(qnorm(1-alpha/(3-d3+1))+
                                                                                       qnorm(1-beta/d3))^2   &&
       (mean(x1)*(1-mean(x1))/length(x1))+(mean(x4)*(1-mean(x4))/length(x4))<=delta^2/(qnorm(1-alpha/(3-d4+1))+
                                                                                       qnorm(1-beta/d4))^2){
      stopcrmet=TRUE ###checking for stopping criteria
    }
    
  }
  N=length(x1)+length(x2)+length(x3)+length(x4)                       ###calculating the sample size
  TF= sum(x1==0)+sum(x2==0)+sum(x3==0)+sum(x4==0)                     ###calculating the number of treatment failures
  d=sort(c(mean(x2),mean(x3),mean(x4)),decreasing = TRUE)            ###ordering the sample means
  ### breaking the ties for equal means
  if(mean(x4)!=mean(x2) && mean(x2)!=mean(x3) && mean(x4)!=mean(x3)){
    d2=which(d==mean(x2))
    d3=which(d==mean(x3))
    d4=which(d==mean(x4))
  } else if (mean(x4)==mean(x2) && mean(x3)>mean(x4)){
    d3=1
    u=1+rank(c(mean(x4),mean(x2)),ties.method= "random")
    d2=u[1]
    d4=u[2]
  }else if (mean(x4)==mean(x2) && mean(x4)>mean(x3)){
    d3=3
    u=rank(c(mean(x2),mean(x4)),ties.method= "random")
    d2=u[1]
    d4=u[2]
  }else if (mean(x3)==mean(x4) && mean(x2)>mean(x3)){
    d2=1
    u=1+rank(c(mean(x3),mean(x4)),ties.method= "random")
    d3=u[1]
    d4=u[2]
  }else if (mean(x3)==mean(x4) && mean(x2)<mean(x3)){
    d2=3
    u=rank(c(mean(x3),mean(x4)),ties.method= "random")
    d3=u[1]
    d4=u[2]
  }else if (mean(x2)==mean(x3) && mean(x2)>mean(x4)){
    d4=3
    u=rank(c(mean(x2),mean(x3)),ties.method= "random")
    d2=u[1]
    d3=u[2]
  }else if (mean(x3)==mean(x2) && mean(x4)>mean(x2)){
    d4=1
    u=1+rank(c(mean(x3),mean(x2)),ties.method= "random")
    d2=u[1]
    d3=u[2]
  }else{
    u=rank(c(mean(x3),mean(x2),mean(x1)),ties.method= "random")
    d4=u[1]
    d2=u[2]
    d3=u[3]
  }
  ##making decisions for three experimental treatments
  
  if(mean(x2)-mean(x1)>(delta0*qnorm(1-beta/d2)+
                        delta1*qnorm(1-alpha/(3-d2+1)))/
     (qnorm(1-alpha/(3-d2+1))+
      qnorm(1-beta/d2))){
    r2=1
    a2=0
  }else{
    r2=0
    a2=1
  }
  if(mean(x3)-mean(x1)>(delta0*qnorm(1-beta/d3)+
                        delta1*qnorm(1-alpha/(3-d3+1)))/
     (qnorm(1-alpha/(3-d3+1))+
      qnorm(1-beta/d3))){
    r3=1
    a3=0
  }else{
    r3=0
    a3=1
  }
  if(mean(x4)-mean(x1)>(delta0*qnorm(1-beta/d4)+
                        delta1*qnorm(1-alpha/(3-d4+1)))/
     (qnorm(1-alpha/(3-d4+1))+
      qnorm(1-beta/d4))){
    r4=1
    a4=0
  }else{
    r4=0
    a4=1
  }
  
  c(N,TF,max(r2,r3),max(a4))   ### change it according to the definition of type-I and type-II errors from true means. 
  ####Third component corresponds to type-I error and fourth component corresponds to type-II error. Note that for global null there is  no type-II error.
}
##Results for SBH procedure
ASN=colMeans(resultholmbal)[1]
ANF=colMeans(resultholmbal)[2]
FWER_I=colMeans(resultholmbal)[3]
FWER_II=colMeans(resultholmbal)[4]
se_samplesize=sd(resultholmbal[,1])/sqrt(iterations)
se_treatmentfailures=sd(resultholmbal[,2])/sqrt(iterations)



