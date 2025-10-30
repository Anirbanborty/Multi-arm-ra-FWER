## Simulation Code for SRAB, SBB, SRAH, SBH Procedures for Normal Distribution


###Fix the values of alpha,beta, delta0,delta1, delta=delta1-delta0
alpha=.05
beta=.1
delta0=0
delta1=1
delta=delta1-delta0

mu1=1   ###mean of control treatment
mu2=1   ###mean of experimental treatment 1
mu3=2   ###mean of experimental treatment 2
mu4=2   ###mean of experimental treatment 3
h=1     ###threshold value which defines success and failure
sigma1=2   ###standard deviation of control treatment
sigma2=2   ###standard deviation of experimental treatment 1
sigma3=2   ###standard deviation of experimental treatment 2
sigma4=2   ###standard deviation of experimental treatment 3
n0=5       ### Pilot sample size K for each treatment 
k=(delta0*qnorm(1-beta/3)+delta1*qnorm(1-alpha/3))/(qnorm(1-alpha/3)+qnorm(1-beta/3))   ###Cutoff for Bonferroni's procedure

iterations = 1000          ######Number of MonteCarlo Runs
###Parallelization 
library(parallel)
library(doParallel)
library(foreach)
nc=detectCores()
cl=makeCluster(nc-2)
registerDoParallel(cl)
#######Simulation for the SRAB procedure
resultbonopt=foreach (i=1:iterations,.combine=rbind)%dopar%{  
  set.seed(i)
  x1=rnorm(n0,mu1,sigma1)   ####generating pilot observations from control
  x4=rnorm(n0,mu4,sigma4)   ####generating pilot observations from experimental treatment 3
  x2=rnorm(n0,mu2,sigma2)   ####generating pilot observations from experimental treatment 1
  x3=rnorm(n0,mu3,sigma3)   ####generating pilot observations from experimental treatment 2
  if( (var(x1)/length(x1))+(var(x2)/length(x2))<=(delta^2)/(qnorm(1-alpha/3)+qnorm(1-beta/3))^2 
     && (var(x1)/length(x1))+(var(x3)/length(x3))<=(delta^2)/(qnorm(1-alpha/3)+qnorm(1-beta/3))^2
     && (var(x1)/length(x1))+(var(x4)/length(x4))<=(delta^2)/(qnorm(1-alpha/3)+qnorm(1-beta/3))^2){
    stopcrmet=TRUE        ### checking for stopping criteria 
  }else {
    stopcrmet=FALSE
  }
  
  while(!stopcrmet){
    ### defining allocation rule
    p=c(sqrt(sd(x1)^2/pnorm((h-mean(x1))/sd(x1)))*sqrt(pnorm((h-mean(x2))/sd(x2))*var(x2)+
                                                         pnorm((h-mean(x3))/sd(x3))*var(x3)+pnorm((h-mean(x4))/sd(x4))*var(x4)),
        var(x2),var(x3),var(x4))/
      (sqrt(sd(x1)^2/pnorm((h-mean(x1))/sd(x1)))*sqrt(pnorm((h-mean(x2))/sd(x2))*var(x2)+
                                                        pnorm((h-mean(x3))/sd(x3))*var(x3)+pnorm((h-mean(x4))/sd(x4))*var(x4))+
         var(x2)+var(x3)+var(x4))
    u=sample(c(1,2,3,4),1,prob=p)    ####sampling by response-adaptive randomization
    ###getting new observation and updating accordingly
    if(u==1){
      x1=c(x1,rnorm(1,mu1,sigma1))
      x2=x2
      x3=x3
      x4=x4
    }else if(u==2){
      x1=x1
      x2=c(x2,rnorm(1,mu2,sigma2))
      x3=x3
      x4=x4
    } else if (u==3) {
      x2=x2
      x1=x1
      x3=c(x3,rnorm(1,mu3,sigma3))
      x4=x4
    } else {
      x2=x2
      x1=x1
      x3=x3
      x4=c(x4,rnorm(1,mu4,sigma4))
    }
    if((var(x1)/length(x1))+(var(x2)/length(x2))<=(delta^2)/(qnorm(1-alpha/3)+qnorm(1-beta/3))^2 
       && (var(x1)/length(x1))+(var(x3)/length(x3))<=(delta^2)/(qnorm(1-alpha/3)+qnorm(1-beta/3))^2
       && (var(x1)/length(x1))+(var(x4)/length(x4))<=(delta^2)/(qnorm(1-alpha/3)+qnorm(1-beta/3))^2){
      stopcrmet=TRUE ### checking for stopping criteria 
    }
    
  }
  N=length(x1)+length(x2)+length(x3)+length(x4)       ###Calculating total sample size
  TF=(sum(x1<h)+sum(x2<h)+sum(x3<h)+sum(x4<h))        ######Calculating total number of failures
  ### decisions for three treatments
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
  c(N,TF,max(r2),max(a3,a4))      ### change it according to the definition of type-I and type-II errors from true means. 
  ####Third component corresponds to type-I error and fourth component corresponds to type-II error.Note that for global null there is  no type-II error.
}
##Results for SRAB procedure
ASN=colMeans(resultbonopt)[1]
ANF=colMeans(resultbonopt)[2]
FWER_I=colMeans(resultbonopt)[3]
FWER_II=colMeans(resultbonopt)[4]
se_samplesize=sd(resultbonopt[,1])/sqrt(iterations)
se_treatmentfailures=sd(resultbonopt[,2])/sqrt(iterations)

######Simulation for the SBB procedure
stopCluster(cl)
iterations = 1000
library(parallel)
library(doParallel)
library(foreach)
nc=detectCores()
cl=makeCluster(nc-2)
registerDoParallel(cl)

resultbonbal=foreach (i=1:iterations,.combine=rbind)%dopar%{  
  set.seed(i)
  x1=rnorm(n0,mu1,sigma1)   ####generating pilot observations from control
  x4=rnorm(n0,mu4,sigma4)   ####generating pilot observations from experimental treatment 3
  x2=rnorm(n0,mu2,sigma2)   ####generating pilot observations from experimental treatment 1
  x3=rnorm(n0,mu3,sigma3)   ####generating pilot observations from experimental treatment 2
  if((var(x1)/length(x1))+(var(x2)/length(x2))<=(delta^2)/(qnorm(1-alpha/3)+qnorm(1-beta/3))^2 
     && (var(x1)/length(x1))+(var(x3)/length(x3))<=(delta^2)/(qnorm(1-alpha/3)+qnorm(1-beta/3))^2
     && (var(x1)/length(x1))+(var(x4)/length(x4))<=(delta^2)/(qnorm(1-alpha/3)+qnorm(1-beta/3))^2){
    stopcrmet=TRUE      ###checking for stopping criteria
  }else {
    stopcrmet=FALSE
  }
  while(!stopcrmet){
    p=c(1/4,1/4,1/4,1/4)       ###Randomization with equal probability 
    u=sample(c(1,2,3,4),1,prob=p)
    ###Getting new observation and updating the information
    if(u==1){
      x1=c(x1,rnorm(1,mu1,sigma1))
      x2=x2
      x3=x3
      x4=x4
    }else if(u==2){
      x1=x1
      x2=c(x2,rnorm(1,mu2,sigma2))
      x3=x3
      x4=x4
    } else if (u==3) {
      x2=x2
      x1=x1
      x3=c(x3,rnorm(1,mu3,sigma3))
      x4=x4
    } else {
      x2=x2
      x1=x1
      x3=x3
      x4=c(x4,rnorm(1,mu4,sigma4))
    }
    if((var(x1)/length(x1))+(var(x2)/length(x2))<=(delta^2)/(qnorm(1-alpha/3)+qnorm(1-beta/3))^2 
       && (var(x1)/length(x1))+(var(x3)/length(x3))<=(delta^2)/(qnorm(1-alpha/3)+qnorm(1-beta/3))^2
       && (var(x1)/length(x1))+(var(x4)/length(x4))<=(delta^2)/(qnorm(1-alpha/3)+qnorm(1-beta/3))^2){
      stopcrmet=TRUE ###checking for stopping criteria
    }
    
  }
  N=length(x1)+length(x2)+length(x3)+length(x4)            ###calculating total sample sizes
  TF=(sum(x1<h)+sum(x2<h)+sum(x3<h)+sum(x4<h))             ###calculating number of treatment failures
  #####decision for three experimental treatments
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
  ####Third component corresponds to type-I error and fourth component corresponds to type-II error.Note that for global null there is  no type-II error.
}
##Results for SBB procedure
ASN=colMeans(resultbonbal)[1]
ANF=colMeans(resultbonbal)[2]
FWER_I=colMeans(resultbonbal)[3]
FWER_II=colMeans(resultbonbal)[4]
se_samplesize=sd(resultbonbal[,1])/sqrt(iterations)
se_treatmentfailures=sd(resultbonbal[,2])/sqrt(iterations)

########Simulation Code for SRAH Procedure############
stopCluster(cl)
n0=5                #####Pilot sample size "K'
iterations = 1000
library(parallel)
library(doParallel)
library(foreach)
nc=detectCores()
cl=makeCluster(nc-2)
registerDoParallel(cl)

resultopt=foreach (i=1:iterations,.combine=rbind,.packages = c("NlcOptim","MASS"))%dopar%{  
  
  set.seed(i)
  x1=rnorm(n0,mu1,sigma1)   ####generating pilot observations from control
  x4=rnorm(n0,mu4,sigma4)   ####generating pilot observations from experimental treatment 3
  x2=rnorm(n0,mu2,sigma2)   ####generating pilot observations from experimental treatment 1
  x3=rnorm(n0,mu3,sigma3)   ####generating pilot observations from experimental treatment 2
  
  ###oredring of the treatments based on sample means
  d=sort(c(mean(x2),mean(x3),mean(x4)),decreasing = TRUE)
  ##breaking tie randomly for equal sample means
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
    stopcrmet=TRUE        ###checking for stopping criteria
  }else{
    stopcrmet=FALSE
  }
  while(!stopcrmet){
    d=sort(c(mean(x2),mean(x3),mean(x4)),decreasing = TRUE)
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
    
    #####defining objective function
    obj=function(x){
      return(x[1]*pnorm((h-mean(x1))/sd(x1))+x[2]*pnorm((h-mean(x2))/sd(x2))
             +x[3]*pnorm((h-mean(x3))/sd(x3))+x[4]*pnorm((h-mean(x4))/sd(x4)))
    }
    ###defing constraints
    budgetconstraint<-function(x){
      f=NULL
      
      
      
      f=rbind(f,var(x2)/x[2]+var(x1)/x[1]-delta^2/(qnorm(1-alpha/(3-d2+1))+
                                                     qnorm(1-beta/d2))^2)
      f=rbind(f,var(x3)/x[3]+var(x1)/x[1]-delta^2/(qnorm(1-alpha/(3-d3+1))+
                                                     qnorm(1-beta/d3))^2)
      f= rbind(f,var(x4)/x[4]+var(x1)/x[1]-delta^2/(qnorm(1-alpha/(3-d4+1))+
                                                      qnorm(1-beta/d4))^2)
      f=rbind(f,1-x[1])
      f=rbind(f,1-x[2])
      f=rbind(f,1-x[3])
      f=rbind(f,1-x[4])
      
      return(list(ceq = NULL, c = f))
    }
    ###initial value for optimization 
    x0=c(10,10,10,10)
    a=solnl(x0,obj,budgetconstraint,tolX = 1e-2,
            tolFun = 1e-5, tolCon = 1e-5, maxnFun = 1e+10,maxIter = 8000)$par     ###numerical optimization
    pi=a/sum(a)   ####obtaining response-adaptive allocation rules
    
    u=sample(c(1,2,3,4),1,prob=pi)      ###randomization using response-adaptive allocation probabilities
    #### Generating new observation and updating the information
    if(u==1){
      x1=c(x1,rnorm(1,mu1,sigma1))
      x2=x2
      x3=x3
      x4=x4
    }else if(u==2){
      x1=x1
      x2=c(x2,rnorm(1,mu2,sigma2))
      x3=x3
      x4=x4
    } else if (u==3) {
      x2=x2
      x1=x1
      x3=c(x3,rnorm(1,mu3,sigma3))
      x4=x4
    } else {
      x2=x2
      x1=x1
      x3=x3
      x4=c(x4,rnorm(1,mu4,sigma4))
      
    }
    d=sort(c(mean(x2),mean(x3),mean(x4)),decreasing = TRUE)
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
      stopcrmet=TRUE ####checking for stopping criteria
    }
    
  }
  N=length(x1)+length(x2)+length(x3)+length(x4)      ####calculating total sample size
  TF=(sum(x1<h)+sum(x2<h)+sum(x3<h)+sum(x4<h))       ####calculating number of failures
  ####ordering the sample means of experimental treatments
  d=sort(c(mean(x2),mean(x3),mean(x4)),decreasing = TRUE)
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
  
  c(N,TF,max(r2),max(a4,a3))  ### change it according to the definition of type-I and type-II errors from true means. 
  ####Third component corresponds to type-I error and fourth component corresponds to type-II error. Note that for global null there is  no type-II error.
  
}
##Results for SRAH procedure
ASN=colMeans(resultopt)[1]
ANF=colMeans(resultopt)[2]
FWER_I=colMeans(resultopt)[3]
FWER_II=colMeans(resultopt)[4]
se_samplesize=sd(resultopt[,1])/sqrt(iterations)
se_treatmentfailures=sd(resultopt[,2])/sqrt(iterations)

########Simulation Code for SBH Procedure############

stopCluster(cl)
iterations = 1000
library(parallel)
library(doParallel)
library(foreach)
nc=detectCores()
cl=makeCluster(nc-2)
registerDoParallel(cl)

resultbal=foreach (i=1:iterations,.combine=rbind,.packages = c("NlcOptim","MASS"))%dopar%{  
  
  set.seed(i)
  x1=rnorm(n0,mu1,sigma1)   ####generating pilot observations from control
  x4=rnorm(n0,mu4,sigma4)   ####generating pilot observations from experimental treatment 3
  x2=rnorm(n0,mu2,sigma2)   ####generating pilot observations from experimental treatment 1
  x3=rnorm(n0,mu3,sigma3)   ####generating pilot observations from experimental treatment 2
  
  d=sort(c(mean(x2),mean(x3),mean(x4)),decreasing = TRUE)         ###ordering the sample means of experimental treatment
  ####breaking tie for equal means using randomization
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
    stopcrmet=TRUE    ###checking stopping criteria
  }else{
    stopcrmet=FALSE
  }
  while(!stopcrmet){
    pi=c(1/4,1/4,1/4,1/4)        ###equal randomization
    
    u=sample(c(1,2,3,4),1,prob=pi)
    ###drawing observation by equal probability allocation and updating the information
    if(u==1){
      x1=c(x1,rnorm(1,mu1,sigma1))
      x2=x2
      x3=x3
      x4=x4
    }else if(u==2){
      x1=x1
      x2=c(x2,rnorm(1,mu2,sigma2))
      x3=x3
      x4=x4
    } else if (u==3) {
      x2=x2
      x1=x1
      x3=c(x3,rnorm(1,mu3,sigma3))
      x4=x4
    } else {
      x2=x2
      x1=x1
      x3=x3
      x4=c(x4,rnorm(1,mu4,sigma4))
      
    }
    d=sort(c(mean(x2),mean(x3),mean(x4)),decreasing = TRUE) ###ordering the sample means of experimental treatment
    ####breaking tie for equal means using randomization
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
      stopcrmet=TRUE     ####checking for stopping criteria
    }
    
  }
  N=length(x1)+length(x2)+length(x3)+length(x4)                 ###calculating the sample size
  TF=(sum(x1<h)+sum(x2<h)+sum(x3<h)+sum(x4<h))                  ###calculating the treatment failures
  d=sort(c(mean(x2),mean(x3),mean(x4)),decreasing = TRUE) ###ordering the sample means of experimental treatment
  ####breaking tie for equal means using randomization
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
  
  #### making decisions for the experimental treatments
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
  
  c(N,TF,max(r2),max(a3,a4)) ### change it according to the definition of type-I and type-II errors from true means. 
  ####Third component corresponds to type-I error and fourth component corresponds to type-II error. Note that for global null there is  no type-II error.
  
}
##Results for SBH procedure
ASN=colMeans(resultbal)[1]
ANF=colMeans(resultbal)[2]
FWER_I=colMeans(resultbal)[3]
FWER_II=colMeans(resultbal)[4]
se_samplesize=sd(resultbal[,1])/sqrt(iterations)
se_treatmentfailures=sd(resultbal[,2])/sqrt(iterations)




