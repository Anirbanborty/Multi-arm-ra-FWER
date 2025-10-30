###Simulation for SRAH* procedure###
mu1=0      ###mean of control
mu2=0      ###mean of experimental treatment 1
mu3=0      ###mean of experimental treatment 2
mu4=1      ###mean of experimental treatment 3
alpha=.05  ###Fixing the value of alpha
beta=.1    ###Fixing the value of beta
sigma1=1               ###variance of control
sigma2=1              ###variance of experimental treatment 1
sigma3=1              ###variance of experimental treatment 2
sigma4=1              ###variance of experimental treatment 1
var1=sigma1^2
var2=sigma2^2
var3=sigma3^2
var4=sigma4^2


delta0=0                   ###Fixing the value of delta0
delta1=1                   ###Fixing the value of delta1
delta=delta1-delta0

psi1=1
psi2=1
psi3=1
psi4=1


n0=5                      ###Fixing the pilot sample size 'K'
iterations = 10^5        ###number of Monte Carlo runs
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
  
  
  if((1/length(x1))+(1/length(x2))<=delta^2/(qnorm(1-alpha/(3-d2+1))+
                                             qnorm(1-beta/d2))^2
     
     && (1/length(x1))+(1/length(x3))<=delta^2/(qnorm(1-alpha/(3-d3+1))+
                                                qnorm(1-beta/d3))^2
     
     && (1/length(x1))+(1/length(x4))<=delta^2/(qnorm(1-alpha/(3-d4+1))+
                                                qnorm(1-beta/d4))^2){
    stopcrmet=TRUE        ###checking for stopping criteria
  }else{
    stopcrmet=FALSE
  }
  ###running a while loop until stopping criteria or maximum allowed sample size is reached
  while(!stopcrmet){
    d=sort(c(mean(x2),mean(x3),mean(x4)),decreasing = TRUE)      ###ordering the sample means 
    ### breaking the ties for equal means
    if(mean(x4)!=mean(x2) && mean(x2)!=mean(x3) && mean(x4)!=mean(x3)){
      d2=which(d==mean(x2))
      d3=which(d==mean(x3))
      d4=which(d==mean(x4))
    } else if (mean(x4)==mean(x2) && mean(x3)>mean(x4)){
      d3=1
      u=1+rank(c(mean(x4),mean(x2)),ties.method= "random")
      d4=u[1]
      d2=u[2]
    }else if (mean(x3)==mean(x2) && mean(x4)>mean(x3)){
      d4=1
      u=1+rank(c(mean(x3),mean(x2)),ties.method= "random")
      d2=u[1]
      d3=u[2]
    }else if (mean(x3)==mean(x4) && mean(x2)>mean(x3)){
      d2=1
      u=1+rank(c(mean(x3),mean(x4)),ties.method= "random")
      d4=u[1]
      d3=u[2]
    }else if (mean(x3)==mean(x2) && mean(x4)<mean(x3)){
      d4=3
      u=rank(c(mean(x3),mean(x2)),ties.method= "random")
      d2=u[1]
      d3=u[2]
    }else if (mean(x4)==mean(x2) && mean(x4)>mean(x3)){
      d3=3
      u=rank(c(mean(x4),mean(x2)),ties.method= "random")
      d2=u[1]
      d4=u[2]
    }else if (mean(x3)==mean(x2) && mean(x4)>mean(x2)){
      d2=3
      u=rank(c(mean(x3),mean(x2)),ties.method= "random")
      d4=u[1]
      d3=u[2]
    }else{
      u=rank(c(mean(x3),mean(x2),mean(x1)),ties.method= "random")
      d4=u[1]
      d2=u[2]
      d3=u[3]
    }
    ###defining the objective function 
    obj=function(x){
      return(x[1]+x[2]
             +x[3]+x[4])
    }
    ###defining the constraints 
    budgetconstraint<-function(x){
      f=NULL
      f=rbind(f,sigma2^2/x[2]+sigma1^2/x[1]-delta^2/(qnorm(1-alpha/(3-d2+1))+
                                                       qnorm(1-beta/d2))^2)
      f=rbind(f,sigma3^2/x[3]+sigma1^2/x[1]-delta^2/(qnorm(1-alpha/(3-d3+1))+
                                                       qnorm(1-beta/d3))^2)
      f= rbind(f,sigma4^2/x[4]+sigma1^2/x[1]-delta^2/(qnorm(1-alpha/(3-d4+1))+
                                                        qnorm(1-beta/d4))^2)
      f=rbind(f,1-x[1])
      f=rbind(f,1-x[2])
      f=rbind(f,1-x[3])
      f=rbind(f,1-x[4])
      
      return(list(ceq = NULL, c = f))
    }
    
    x0=c(10,10,10,10)          ###initial value for numerical optimization
    a=solnl(x0,obj,budgetconstraint,tolX = 1e-2,
            tolFun = 1e-5, tolCon = 1e-5, maxnFun = 1e+10,maxIter = 8000)$par
    pi=a/sum(a)     ###obtaining response-adaptive allocation rule  
    
    u=sample(c(1,2,3,4),1,prob=pi) ### response-adaptive randomization 
    ###generating new observation and updating the information 
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
    d=sort(c(mean(x2),mean(x3),mean(x4)),decreasing = TRUE)  ###ordering the sample means 
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
    
    
    if(((1/length(x1))+(1/length(x2))<=delta^2/(qnorm(1-alpha/(3-d2+1))+
                                                qnorm(1-beta/d2))^2
        
        && (1/length(x1))+(1/length(x3))<=delta^2/(qnorm(1-alpha/(3-d3+1))+
                                                   qnorm(1-beta/d3))^2
        
        && (1/length(x1))+(1/length(x4))<=delta^2/(qnorm(1-alpha/(3-d4+1))+
                                                   qnorm(1-beta/d4))^2)
       |(length(x1)+length(x2)+length(x3)+length(x4)==85)){
      stopcrmet=TRUE    ###checking for stopping criteria 
    }
    
  }
  ###making decisions for experimental treatments based on the reason behind stopping (maximum allowable sample size is reached or stopping criteria is satisfied)
  
  if((1/length(x1))+(1/length(x2))<=delta^2/(qnorm(1-alpha/(3-d2+1))+
                                             qnorm(1-beta/d2))^2
     
     && (1/length(x1))+(1/length(x3))<=delta^2/(qnorm(1-alpha/(3-d3+1))+
                                                qnorm(1-beta/d3))^2
     
     && (1/length(x1))+(1/length(x4))<=delta^2/(qnorm(1-alpha/(3-d4+1))+
                                                qnorm(1-beta/d4))^2){
    s=1
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
  } else {
    s=0
    v2=(mean(x2)-mean(x1))/sqrt(1/length(x1)+1/length(x2))
    v3=(mean(x3)-mean(x1))/sqrt(1/length(x1)+1/length(x3))
    v4=(mean(x4)-mean(x1))/sqrt(1/length(x1)+1/length(x4))
    d=sort(c(v2,v3,v4),decreasing = TRUE)
    if(v4!=v2 && v2!=v3 && v4!=v3){
      d2=which(d==v2)
      d3=which(d==v3)
      d4=which(d==v4)
    } else if (v4==v2 && v3>v4){
      d3=1
      u=1+rank(c(v4,v2),ties.method= "random")
      d2=u[1]
      d4=u[2]
    }else if (v4==v2 && v4>v3){
      d3=3
      u=rank(c(v2,v4),ties.method= "random")
      d2=u[1]
      d4=u[2]
    }else if (v3==v4 && v2>v3){
      d2=1
      u=1+rank(c(v3,v4),ties.method= "random")
      d3=u[1]
      d4=u[2]
    }else if (v3==v4 && v2<v3){
      d2=3
      u=rank(c(v3,v4),ties.method= "random")
      d3=u[1]
      d4=u[2]
    }else if (v2==v3 && v2>v4){
      d4=3
      u=rank(c(v2,v3),ties.method= "random")
      d2=u[1]
      d3=u[2]
    }else if (v3==v2 && v4>v2){
      d4=1
      u=1+rank(c(v3,v2),ties.method= "random")
      d2=u[1]
      d3=u[2]
    }else{
      u=rank(c(v3,v2,v1),ties.method= "random")
      d4=u[1]
      d2=u[2]
      d3=u[3]
    }
    
    if(v2>qnorm(1-alpha/(3-d2+1))){
      r2=1
      a2=0
    }else{
      r2=0
      a2=1
    }
    if(v3>qnorm(1-alpha/(3-d3+1))){
      r3=1
      a3=0
    }else{
      r3=0
      a3=1
    }
    if(v4>qnorm(1-alpha/(3-d4+1))){
      r4=1
      a4=0
    }else{
      r4=0
      a4=1
    }
  }
  sample_size=length(x1)+length(x2)+length(x3)+length(x4)             ###calculating sample size
  c(s,max(r2,r3),max(r4),sample_size)          ### change it according to the definition of type-I error and the disjunctive power from true means. 
  ####Second component corresponds to type-I error and third component corresponds to disjunctive power. Note that for global null there is  no power. 
  
}
ASN=colMeans(resultopt)[4]
FWER_I=colMeans(resultopt)[2]
FWER_II=colMeans(resultopt)[3]
se_samplesize=sd(resultopt[,4])/sqrt(iterations)
PES=colMeans(resultopt)[1]
