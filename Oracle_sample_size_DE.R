#####Oracle Sample Size calculation for DE####
alpha=0.05           ###Fix the value of alpha
beta=.1              ###Fix the value of beta
mu1=1                ##mean of the control 
mu2=1                ##mean of experimental treatment 1
mu3=2                ##mean of experimental treatment 2
mu4=2                ##mean of experimental treatment 3
b1=1            ##scale parameter of control 
b2=1              ##scale parameter of experimental treatment 1
b3=1              ##scale parameter of experimental treatment 2
b4=1              ##scale parameter of experimental treatment 3
sigma1=b1*sqrt(2)    
sigma2=b2*sqrt(2)
sigma3=b3*sqrt(2)
sigma4=b4*sqrt(2)
var1=sigma1^2           ##variance of control
var2=sigma2^2           ##variance of experimental treatment 1
var3=sigma3^2           ##variance of experimental treatment 2
var4=sigma4^2           ##variance of experimental treatment 3

h=1.5                     ##fixing the value of h
delta0=0                ##fixing the value of delta0
delta1=1                 ##fixing the value of delta1
delta=delta1-delta0
###calculating allocation proportions 
if(h<=mu1){
  psi1=.5*exp((h-mu1)/b1)
} else {
  psi1=1-.5*exp(-(h-mu1)/b1)
}
if(h<=mu2){
  psi2=.5*exp((h-mu2)/b2)
} else {
  psi2=1-.5*exp(-(h-mu2)/b2)
}
if(h<=mu3){
  psi3=.5*exp((h-mu3)/b3)
} else {
  psi3=1-.5*exp(-(h-mu3)/b3)
}
if(h<=mu4){
  psi4=.5*exp((h-mu4)/b4)
} else {
  psi4=1-.5*exp(-(h-mu4)/b4)
}

pi1=sqrt(var1/psi1)*
  sqrt(var2*psi2+
         var3*psi3+
         var4*psi4)
pi2=var2
pi3=var3
pi4=var4

pi11=pi1/(pi1+pi2+pi3+pi4)       ###allocation proportion of control
pi22=pi2/(pi1+pi2+pi3+pi4)         ###allocation proportion of experimental treatment 1
pi33=pi3/(pi1+pi2+pi3+pi4)         ###allocation proportion of experimental treatment 2
pi44=pi4/(pi1+pi2+pi3+pi4)         ###allocation proportion of experimental treatment 3
p=c(pi11,pi22,pi33,pi44)

pi11=pi1/(pi1+pi2+pi3+pi4)
pi22=pi2/(pi1+pi2+pi3+pi4)
pi33=pi3/(pi1+pi2+pi3+pi4)
pi44=pi4/(pi1+pi2+pi3+pi4)

a11=(qnorm(1-alpha/(3))+qnorm(1-beta/3))^2*((sigma2^2)/pi22+(sigma1^2)/pi11)/delta^2
b11=(qnorm(1-alpha/(3))+qnorm(1-beta/3))^2*((sigma3^2)/pi33+(sigma1^2)/pi11)/delta^2
c11=(qnorm(1-alpha/(3))+qnorm(1-beta/3))^2*((sigma4^2)/pi44+(sigma1^2)/pi11)/delta^2
n_SRAB=a11      ##  oracle sample size for SRAB. We report it as the smallest integer greater than n_SRAB.Here a11,b11,c11 should be same. 

####SBB###


a22=(qnorm(1-alpha/(3))+qnorm(1-beta/3))^2*((sigma2^2)+(sigma1^2))*4/delta^2
b22=(qnorm(1-alpha/(3))+qnorm(1-beta/3))^2*((sigma3^2)+(sigma1^2))*4/delta^2
c22=(qnorm(1-alpha/(3))+qnorm(1-beta/3))^2*((sigma4^2)+(sigma1^2))*4/delta^2
n_SBB=max(a22,b22,c22)  ## oracle sample size for SBB. We report it as the smallest integer greater than n_SBB.




#####True sample size for holm####
#####optimal###
mu1=1
mu2=1
mu3=2
mu4=2
b1=1
b2=1
b3=1
b4=1
sigma1=b1*sqrt(2)
sigma2=b2*sqrt(2)
sigma3=b3*sqrt(2)
sigma4=b4*sqrt(2)
var1=sigma1^2
var2=sigma2^2
var3=sigma3^2
var4=sigma4^2

h=1.5
delta0=0
delta1=1
delta=delta1-delta0
if(h<=mu1){
  psi1=.5*exp((h-mu1)/b1)
} else {
  psi1=1-.5*exp(-(h-mu1)/b1)
}
if(h<=mu2){
  psi2=.5*exp((h-mu2)/b2)
} else {
  psi2=1-.5*exp(-(h-mu2)/b2)
}
if(h<=mu3){
  psi3=.5*exp((h-mu3)/b3)
} else {
  psi3=1-.5*exp(-(h-mu3)/b3)
}
if(h<=mu4){
  psi4=.5*exp((h-mu4)/b4)
} else {
  psi4=1-.5*exp(-(h-mu4)/b4)
}
###defining objective function
obj=function(x){
  return(x[1]*psi1+x[2]*psi2
         +x[3]*psi3+x[4]*psi4)
}
###defining the order of the true means
d4=1
d3=2
d2=3
###defining constraints 
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
library(NlcOptim)
x0=c(10,10,10,10)            ###initial value for numerical optimization 
a=solnl(x0,obj,budgetconstraint,tolX = 1e-2,
        tolFun = 1e-5, tolCon = 1e-5, maxnFun = 1e+10,maxIter = 8000)$par
pi=a/sum(a)                  ### obtaing optimal allocation rule from numerical optimization  

a33=(qnorm(1-alpha/(3-d2+1))+qnorm(1-beta/d2))^2*((sigma2^2)/pi[2]+(sigma1^2)/pi[1])/delta^2
b33=(qnorm(1-alpha/(3-d3+1))+qnorm(1-beta/d3))^2*((sigma3^2)/pi[3]+(sigma1^2)/pi[1])/delta^2
c33=(qnorm(1-alpha/(3-d4+1))+qnorm(1-beta/d4))^2*((sigma4^2)/pi[4]+(sigma1^2)/pi[1])/delta^2
n_SRAH=a33   ##  Oracle sample size for SRAH. We report it as the smallest integer greater than n_SRAH. Here a3,b3,c3 should be same.

### Calculation for SBH 

a44=(qnorm(1-alpha/(3-d2+1))+qnorm(1-beta/d2))^2*(sigma2^2+sigma1^2)*4/delta^2
b44=(qnorm(1-alpha/(3-d3+1))+qnorm(1-beta/d3))^2*(sigma3^2+sigma1^2)*4/delta^2
c44=(qnorm(1-alpha/(3-d4+1))+qnorm(1-beta/d4))^2*(sigma4^2+sigma1^2)*4/delta^2

n_SBH=max(a44,b44,c44)  ## oracle sample size for SBH. We report it as the smallest integer greater than n_SBH.



