####Simullation studies for the SRAB, SBB, SRAH, SBH procedures for t5 distribution####


alpha=.05       ###Fixing the value of alpha
beta=.1         ###Fixing the value of beta
v=5             ##degree of t-distribution fixed as 5
mu1=1           ###mean of the control
mu2=1           ###mean of the experimental treatment 1
mu3=2         ###mean of the experimental treatment 1
mu4=2          ###mean of the experimental treatment 1
sigma1=1.5      ###scale parameter of the control
sigma2=1.5      ###scale parameter of the experimental treatment 1
sigma3=1.5      ###scale parameter of the experimental treatment 2
sigma4=1.5      ###scale parameter of the experimental treatment 3

h=1.5             ### fixing the threshold value which defines success and failure


delta0=0       ###fixing the value of delta0
delta1=1       ###fixing the value of delta1

delta=delta1-delta0     
n0=5                ###fixing the value of pilot sample size 'K'
k=(delta0*qnorm(1-beta/3)+delta1*qnorm(1-alpha/3))/(qnorm(1-alpha/3)+qnorm(1-beta/3))      ###cutoff for bonferroni's procedure

iterations = 1000                         ###number of Monte Carlo runs

###########Simulation for SRAB procedure#######
library(parallel)
library(doParallel)
library(foreach)
nc=detectCores()
cl=makeCluster(nc-1)
registerDoParallel(cl)

resultbonopt=foreach (i=1:iterations,.combine=rbind,.packages = c("smoothmest","NlcOptim","MASS"))%dopar%{  
  
  
  set.seed(i)
  x1=sigma1*rt(n0,v)+mu1      ###generating pilot observations from control 
  x4=sigma4*rt(n0,v)+mu4     ###generating pilot observations from experimental treatment 3
  x2=sigma2*rt(n0,v)+mu2      ###generating pilot observations from experimental treatment 1
  x3=sigma3*rt(n0,v)+mu3     ###generating pilot observations from experimental treatment 2
  if((var(x1)/length(x1))+(var(x2)/length(x2))<=(delta^2)/(qnorm(1-alpha/3)+qnorm(1-beta/3))^2 
     && (var(x1)/length(x1))+(var(x3)/length(x3))<=(delta^2)/(qnorm(1-alpha/3)+qnorm(1-beta/3))^2
     && (var(x1)/length(x1))+(var(x4)/length(x4))<=(delta^2)/(qnorm(1-alpha/3)+qnorm(1-beta/3))^2){
    stopcrmet=TRUE     ###checking for stopping criteria
  }else {
    stopcrmet=FALSE
  }
  ###running a while loop until the stopping criteria is reached 
  while(!stopcrmet){
    ###calculating steps for response-adaptive allocation rules
    mu1hat=mean(x1)    
    mu2hat=mean(x2)
    mu3hat=mean(x3)
    mu4hat=mean(x4)
    var1=var(x1)
    var2=var(x2)
    var3=var(x3)
    var4=var(x4)
    sigma1hat=sqrt(var1*(v-2)/v)
    sigma2hat=sqrt(var2*(v-2)/v)
    sigma3hat=sqrt(var3*(v-2)/v)
    sigma4hat=sqrt(var4*(v-2)/v)
    psi1hat=pt((h-mu1hat)/sigma1hat,df=v)
    psi2hat=pt((h-mu2hat)/sigma2hat,df=v)
    psi3hat=pt((h-mu3hat)/sigma3hat,df=v)
    psi4hat=pt((h-mu4hat)/sigma4hat,df=v)
    
    
    pi1=sqrt(var1/psi1hat)*
      sqrt(var2*psi2hat+
             var3*psi3hat+
             var4*psi4hat)
    pi2=var2
    pi3=var3
    pi4=var4
    
    pi11=pi1/(pi1+pi2+pi3+pi4)       ###response-adaptive allocation probability for control
    pi22=pi2/(pi1+pi2+pi3+pi4)       ###response-adaptive allocation probability for experimental treatment 1
    pi33=pi3/(pi1+pi2+pi3+pi4)        ###response-adaptive allocation probability for experimental treatment 2
    pi44=pi4/(pi1+pi2+pi3+pi4)         ###response-adaptive allocation probability for experimental treatment 3
    p=c(pi11,pi22,pi33,pi44)
    u=sample(c(1,2,3,4),1,prob=p)     ###response-adaptive randomization
    ###generating new observation and updating the information
    if(u==1){
      x1=c(x1,sigma1*rt(1,v)+mu1)
      x2=x2
      x3=x3
      x4=x4
    }else if(u==2){
      x1=x1
      x2=c(x2,sigma2*rt(1,v)+mu2)
      x3=x3
      x4=x4
    } else if (u==3) {
      x2=x2
      x1=x1
      x3=c(x3,sigma3*rt(1,v)+mu3)
      x4=x4
    } else {
      x2=x2
      x1=x1
      x3=x3
      x4=c(x4,sigma4*rt(1,v)+mu4)
    }
    var1=var(x1)
    var2=var(x2)
    var3=var(x3)
    var4=var(x4)
    if(var1/length(x1)+(var2/length(x2))<=(delta^2)/(qnorm(1-alpha/3)+qnorm(1-beta/3))^2 
       && (var1/length(x1))+(var3/length(x3))<=(delta^2)/(qnorm(1-alpha/3)+qnorm(1-beta/3))^2
       && (var1/length(x1))+(var4/length(x4))<=(delta^2)/(qnorm(1-alpha/3)+qnorm(1-beta/3))^2){
      stopcrmet=TRUE       ###checking for stopping criteria 
    }
    
    
  }
  N=length(x1)+length(x2)+length(x3)+length(x4)              ###calculating total sample size
  TF=sum(x1<h)+sum(x2<h)+sum(x3<h)+sum(x4<h)                 ###calculating the number of treatment failures 
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
  
  c(N,TF,max(r2),max(a3,a4))### change it according to the definition of type-I and type-II errors from true means. 
  ####Third component corresponds to type-I error and fourth component corresponds to type-II error. Note that for global null there is  no type-II error.
  
}
##Results for SRAB procedure
ASN=colMeans(resultbonopt)[1]
ANF=colMeans(resultbonopt)[2]
FWER_I=colMeans(resultbonopt)[3]
FWER_II=colMeans(resultbonopt)[4]
se_samplesize=sd(resultbonopt[,1])/sqrt(iterations)
se_treatmentfailures=sd(resultbonopt[,2])/sqrt(iterations)




#############Simulation for SBB##############

stopCluster(cl)
library(parallel)
library(doParallel)
library(foreach)
nc=detectCores()
cl=makeCluster(nc-1)
registerDoParallel(cl)

resultbonbal=foreach (i=1:iterations,.combine=rbind,.packages = c("smoothmest","NlcOptim","MASS"))%dopar%{  
  
  
  set.seed(i)
  x1=sigma1*rt(n0,v)+mu1      ###generating pilot observations from control 
  x4=sigma4*rt(n0,v)+mu4     ###generating pilot observations from experimental treatment 3
  x2=sigma2*rt(n0,v)+mu2      ###generating pilot observations from experimental treatment 1
  x3=sigma3*rt(n0,v)+mu3     ###generating pilot observations from experimental treatment 2
  if((var(x1)/length(x1))+(var(x2)/length(x2))<=(delta^2)/(qnorm(1-alpha/3)+qnorm(1-beta/3))^2 
     && (var(x1)/length(x1))+(var(x3)/length(x3))<=(delta^2)/(qnorm(1-alpha/3)+qnorm(1-beta/3))^2
     && (var(x1)/length(x1))+(var(x4)/length(x4))<=(delta^2)/(qnorm(1-alpha/3)+qnorm(1-beta/3))^2){
    stopcrmet=TRUE      ###checking for stopping criteria    
  }else {
    stopcrmet=FALSE
  }
  ###running a while loop until the stopping criteria is reached 
  while(!stopcrmet){
    
    pi11=1/4
    pi22=1/4
    pi33=1/4
    pi44=1/4
    p=c(pi11,pi22,pi33,pi44)         
    u=sample(c(1,2,3,4),1,prob=p)       ###equal probability randomization
    ###generating new observation and updating the information 
    if(u==1){
      x1=c(x1,sigma1*rt(1,v)+mu1)
      x2=x2
      x3=x3
      x4=x4
    }else if(u==2){
      x1=x1
      x2=c(x2,sigma2*rt(1,v)+mu2)
      x3=x3
      x4=x4
    } else if (u==3) {
      x2=x2
      x1=x1
      x3=c(x3,sigma3*rt(1,v)+mu3)
      x4=x4
    } else {
      x2=x2
      x1=x1
      x3=x3
      x4=c(x4,sigma4*rt(1,v)+mu4)
    }
    var1=var(x1)
    var2=var(x2)
    var3=var(x3)
    var4=var(x4)
    if(var1/length(x1)+(var2/length(x2))<=(delta^2)/(qnorm(1-alpha/3)+qnorm(1-beta/3))^2 &&
       (var1/length(x1))+(var3/length(x3))<=(delta^2)/(qnorm(1-alpha/3)+qnorm(1-beta/3))^2  &&
       (var1/length(x1))+(var4/length(x4))<=(delta^2)/(qnorm(1-alpha/3)+qnorm(1-beta/3))^2){
      stopcrmet=TRUE      ####checking for stopping criteria
    }
    
  }
  N=length(x1)+length(x2)+length(x3)+length(x4)                ###calculating the sample size
  TF=sum(x1<h)+sum(x2<h)+sum(x3<h)+sum(x4<h)                   ###calculating the treatment failures
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
  
  c(N,TF,max(r2),max(a3,a4)) ### change it according to the definition of type-I and type-II errors from true means. 
  ####Third component corresponds to type-I error and fourth component corresponds to type-II error. Note that for global null there is  no type-II error.
  
} 
##Results for SBB procedure
ASN=colMeans(resultbonbal)[1]
ANF=colMeans(resultbonbal)[2]
FWER_I=colMeans(resultbonbal)[3]
FWER_II=colMeans(resultbonbal)[4]
se_samplesize=sd(resultbonbal[,1])/sqrt(iterations)
se_treatmentfailures=sd(resultbonbal[,2])/sqrt(iterations)


########Simulation for SRAH procedure##################
stopCluster(cl)
library(parallel)
library(doParallel)
library(foreach)
nc=detectCores()
cl=makeCluster(nc-1)
registerDoParallel(cl)
resultholmopt=foreach (i=1:iterations,.combine=rbind,.packages = c("NlcOptim","MASS","smoothmest"))%dopar%{  
  
  set.seed(i)
  x1=sigma1*rt(n0,v)+mu1      ###generating pilot observations from control 
  x4=sigma4*rt(n0,v)+mu4     ###generating pilot observations from experimental treatment 3
  x2=sigma2*rt(n0,v)+mu2      ###generating pilot observations from experimental treatment 1
  x3=sigma3*rt(n0,v)+mu3     ###generating pilot observations from experimental treatment 2
  d=sort(c(mean(x2),mean(x3),mean(x4)),decreasing = TRUE)       ###ordering the sample means
  ####breaking ties for equal means
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
  
  if((var(x1)/length(x1))+(var(x2)/length(x2))<=delta^2/(qnorm(1-alpha/(3-d2+1))+
                                                   qnorm(1-beta/d2))^2
     
     && (var(x1)/length(x1))+(var(x3)/length(x3))<=delta^2/(qnorm(1-alpha/(3-d3+1))+
                                                      qnorm(1-beta/d3))^2
     
     && (var(x1)/length(x1))+(var(x4)/length(x4))<=delta^2/(qnorm(1-alpha/(3-d4+1))+
                                                      qnorm(1-beta/d4))^2){
    stopcrmet=TRUE       ###checking for stopping criteria
  }else{
    stopcrmet=FALSE
  }
  ###running a while loop until the stopping criteria is reached 
  while(!stopcrmet){
    mu1hat=mean(x1)  
    mu2hat=mean(x2)
    mu3hat=mean(x3)
    mu4hat=mean(x4)
    var1=var(x1)
    var2=var(x2)
    var3=var(x3)
    var4=var(x4)
    sigma1hat=sqrt(var1*(v-2)/v)  ###estimating sigma1
    sigma2hat=sqrt(var2*(v-2)/v)  ###estimating sigma2
    sigma3hat=sqrt(var3*(v-2)/v)  ###estimating sigma3
    sigma4hat=sqrt(var4*(v-2)/v)  ###estimating sigma4
    d=sort(c(mean(x2),mean(x3),mean(x4)),decreasing = TRUE)           ###ordering the sample means
    ####breaking ties for equal means             
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
    
    psi1hat=pt((h-mu1hat)/sigma1hat,df=v)
    psi2hat=pt((h-mu2hat)/sigma2hat,df=v)
    psi3hat=pt((h-mu3hat)/sigma3hat,df=v)
    psi4hat=pt((h-mu4hat)/sigma4hat,df=v)
    ###defining the objective function 
    obj=function(x){
      return(x[1]*psi1hat+x[2]*psi2hat
             +x[3]*psi3hat+x[4]*psi4hat)
    }
    ###defining the constraints
    budgetconstraint<-function(x){
      f=NULL
      f=rbind(f,var2/x[2]+var1/x[1]-delta^2/(qnorm(1-alpha/(3-d2+1))+
                                               qnorm(1-beta/d2))^2)
      f=rbind(f,var3/x[3]+var1/x[1]-delta^2/(qnorm(1-alpha/(3-d3+1))+
                                               qnorm(1-beta/d3))^2)
      f= rbind(f,var4/x[4]+var1/x[1]-delta^2/(qnorm(1-alpha/(3-d4+1))+
                                                qnorm(1-beta/d4))^2)
      f=rbind(f,1-x[1])
      f=rbind(f,1-x[2])
      f=rbind(f,1-x[3])
      f=rbind(f,1-x[4])
      
      return(list(ceq = NULL, c = f))
    }
    
    x0=c(10,10,10,10)           ####initial value for numerical optimization 
    a=solnl(x0,obj,budgetconstraint,tolX = 1e-2,
            tolFun = 1e-5, tolCon = 1e-5, maxnFun = 1e+10,maxIter = 8000)$par
    pi=a/sum(a)                           ###obtaining response-adaptive allocation rule by numerical optimization
    
    u=sample(c(1,2,3,4),1,prob=pi)        ###response-adaptive randomization 
    ###generating new observation and updating the information 
    if(u==1){
      x1=c(x1,sigma1*rt(1,v)+mu1)
      x2=x2
      x3=x3
      x4=x4
    }else if(u==2){
      x1=x1
      x2=c(x2,sigma2*rt(1,v)+mu2)
      x3=x3
      x4=x4
    } else if (u==3) {
      x2=x2
      x1=x1
      x3=c(x3,sigma3*rt(1,v)+mu3)
      x4=x4
    } else {
      x2=x2
      x1=x1
      x3=x3
      x4=c(x4,sigma4*rt(1,v)+mu4)
      
    }
    
    var1=var(x1)
    var2=var(x2)
    var3=var(x3)
    var4=var(x4)
    d=sort(c(mean(x2),mean(x3),mean(x4)),decreasing = TRUE) ###ordering the sample means
    ####breaking ties for equal means 
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
    
    if((var(x1)/length(x1))+(var(x2)/length(x2))<=delta^2/(qnorm(1-alpha/(3-d2+1))+
                                                     qnorm(1-beta/d2))^2
       
       && (var(x1)/length(x1))+(var(x3)/length(x3))<=delta^2/(qnorm(1-alpha/(3-d3+1))+
                                                        qnorm(1-beta/d3))^2
       
       && (var(x1)/length(x1))+(var(x4)/length(x4))<=delta^2/(qnorm(1-alpha/(3-d4+1))+
                                                        qnorm(1-beta/d4))^2){
      stopcrmet=TRUE          ###checking for stopping criteria 
    }
    
  }
  N=length(x1)+length(x2)+length(x3)+length(x4)                                ###calculating the sample size
  TF= sum(x1<h)+sum(x2<h)+sum(x3<h)+sum(x4<h)                                  ### calculating the number of treatment failures
  d=sort(c(mean(x2),mean(x3),mean(x4)),decreasing = TRUE)             ###ordering the sample means
  ####breaking ties for equal means 
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
  ###making decisions for three experimental treatments
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
  
  c(N,TF,max(r2),max(a3,a4))         ### change it according to the definition of type-I and type-II errors from true means. 
  ####Third component corresponds to type-I error and fourth component corresponds to type-II error. Note that for global null there is  no type-II error.
  
}
##Results for SRAH procedure
ASN=colMeans(resultholmopt)[1]
ANF=colMeans(resultholmopt)[2]
FWER_I=colMeans(resultholmopt)[3]
FWER_II=colMeans(resultholmopt)[4]
se_samplesize=sd(resultholmopt[,1])/sqrt(iterations)
se_treatmentfailures=sd(resultholmopt[,2])/sqrt(iterations)

########################Simulation for SBH Procedure#########################
stopCluster(cl)
iterations = 1000
library(parallel)
library(doParallel)
library(foreach)
nc=detectCores()
cl=makeCluster(nc-2)
registerDoParallel(cl)

resultholmbal=foreach (i=1:iterations,.combine=rbind,.packages = c("NlcOptim","MASS","smoothmest"))%dopar%{  
  
  set.seed(i)
  x1=sigma1*rt(n0,v)+mu1      ###generating pilot observations from control 
  x4=sigma4*rt(n0,v)+mu4     ###generating pilot observations from experimental treatment 3
  x2=sigma2*rt(n0,v)+mu2      ###generating pilot observations from experimental treatment 1
  x3=sigma3*rt(n0,v)+mu3     ###generating pilot observations from experimental treatment 2
  d=sort(c(mean(x2),mean(x3),mean(x4)),decreasing = TRUE)         ###ordering the sample means
  ####breaking ties for equal means 
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
  
  if((var(x1)/length(x1))+(var(x2)/length(x2))<=delta^2/(qnorm(1-alpha/(3-d2+1))+
                                                   qnorm(1-beta/d2))^2
     
     && (var(x1)/length(x1))+(var(x3)/length(x3))<=delta^2/(qnorm(1-alpha/(3-d3+1))+
                                                      qnorm(1-beta/d3))^2
     
     && (var(x1)/length(x1))+(var(x4)/length(x4))<=delta^2/(qnorm(1-alpha/(3-d4+1))+
                                                      qnorm(1-beta/d4))^2){
    stopcrmet=TRUE   ###checking for stopping criteria 
  }else{
    stopcrmet=FALSE
  }
  ###running a while loop until the stopping criteria is reached 
  while(!stopcrmet){
    
    
    pi=c(1,1,1,1)/4
    u=sample(c(1,2,3,4),1,prob=pi)       ###equal probability randomization
    ###generting new observation and updating the information
    if(u==1){
      x1=c(x1,sigma1*rt(1,v)+mu1)
      x2=x2
      x3=x3
      x4=x4
    }else if(u==2){
      x1=x1
      x2=c(x2,sigma2*rt(1,v)+mu2)
      x3=x3
      x4=x4
    } else if (u==3) {
      x2=x2
      x1=x1
      x3=c(x3,sigma3*rt(1,v)+mu3)
      x4=x4
    } else {
      x2=x2
      x1=x1
      x3=x3
      x4=c(x4,sigma4*rt(1,v)+mu4)
      
    }
    
    var1=var(x1)
    var2=var(x2)
    var3=var(x3)
    var4=var(x4)
    d=sort(c(mean(x2),mean(x3),mean(x4)),decreasing = TRUE) ###ordering the sample means
    ####breaking ties for equal means 
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
    
    
    if((var1/length(x1))+(var2/length(x2))<=delta^2/(qnorm(1-alpha/(3-d2+1))+
                                                     qnorm(1-beta/d2))^2
       
       && (var1/length(x1))+(var3/length(x3))<=delta^2/(qnorm(1-alpha/(3-d3+1))+
                                                        qnorm(1-beta/d3))^2
       
       && (var1/length(x1))+(var4/length(x4))<=delta^2/(qnorm(1-alpha/(3-d4+1))+
                                                        qnorm(1-beta/d4))^2){
      stopcrmet=TRUE      ##checking for stopping criteria
    }
    
  }
  N=length(x1)+length(x2)+length(x3)+length(x4)                      ###calculating the sample size
  TF= sum(x1<h)+sum(x2<h)+sum(x3<h)+sum(x4<h)                        ### calculating the number of treatment failures
  d=sort(c(mean(x2),mean(x3),mean(x4)),decreasing = TRUE)             ###ordering the sample means
  ####breaking ties for equal means 
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
  ###making decisions for three experimental treatments 
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
  
  c(N,TF,max(r2),max(a3,a4))  ### change it according to the definition of type-I and type-II errors from true means. 
  ####Third component corresponds to type-I error and fourth component corresponds to type-II error. Note that for global null there is  no type-II error.
  
}
##Results for SBH procedure
ASN=colMeans(resultholmbal)[1]
ANF=colMeans(resultholmbal)[2]
FWER_I=colMeans(resultholmbal)[3]
FWER_II=colMeans(resultholmbal)[4]
se_samplesize=sd(resultholmbal[,1])/sqrt(iterations)
se_treatmentfailures=sd(resultholmbal[,2])/sqrt(iterations)




