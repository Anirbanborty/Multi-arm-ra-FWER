###Calculation of oracle sample size for Bernoulli distribution###
alpha=.05       ###Fixing the value of alpha
beta=.1         ###Fixing the value of beta
mu1=.2          ###mean of control 
mu2=.3         ###mean of experimental treatment 1
mu3=.5          ###mean of experimental treatment 2 
mu4=.5           ###mean of experimental treatment 1
delta0=.1         ##fixing the value of delta0
delta1=.3         ##fixing the value of delta1
delta=delta1-delta0
sigma1=sqrt(mu1*(1-mu1))        ###sd of control 
sigma2=sqrt(mu2*(1-mu2))        ###sd of experimental treatment 1
sigma3=sqrt(mu3*(1-mu3))         ###sd of experimental treatment 2 
sigma4=sqrt(mu4*(1-mu4))         ###sd of experimental treatment 3  
pi1=sqrt(mu1*(1-mu1)/(1-mu1))*
  sqrt(mu2*(1-mu2)*(1-mu2)+
         mu3*(1-mu3)*(1-mu3)+
         mu4*(1-mu4)*(1-mu4))
pi2=mu2*(1-mu2)
pi3=mu3*(1-mu3)
pi4=mu4*(1-mu4)

pi11=pi1/(pi1+pi2+pi3+pi4)       ###allocation proportion of the control 
pi22=pi2/(pi1+pi2+pi3+pi4)       ###allocation proportion of the experimental treatment 1
pi33=pi3/(pi1+pi2+pi3+pi4)       ###allocation proportion of the experimental treatment 2
pi44=pi4/(pi1+pi2+pi3+pi4)       ###allocation proportion of the experimental treatment 3
p=c(pi11,pi22,pi33,pi44)

a11=(qnorm(1-alpha/(3))+qnorm(1-beta/3))^2*((sigma2^2)/pi22+(sigma1^2)/pi11)/delta^2
b11=(qnorm(1-alpha/(3))+qnorm(1-beta/3))^2*((sigma3^2)/pi33+(sigma1^2)/pi11)/delta^2
c11=(qnorm(1-alpha/(3))+qnorm(1-beta/3))^2*((sigma4^2)/pi44+(sigma1^2)/pi11)/delta^2
n_SRAB=a11      ##  oracle sample size for SRAB. We report it as the smallest integer greater than n_SRAB.Here a11,b11,c11 should be same. 
##############SBB####################
a22=(qnorm(1-alpha/(3))+qnorm(1-beta/3))^2*((sigma2^2)+(sigma1^2))*4/delta^2
b22=(qnorm(1-alpha/(3))+qnorm(1-beta/3))^2*((sigma3^2)+(sigma1^2))*4/delta^2
c22=(qnorm(1-alpha/(3))+qnorm(1-beta/3))^2*((sigma4^2)+(sigma1^2))*4/delta^2
n_SBB=max(a22,b22,c22)  ## oracle sample size for SBB. We report it as the smallest integer greater than n_SBB.


#################Holm################
#####True sample size for holm####
#####optimal###

###defining the objective function
obj=function(x){
  return(x[1]*(1-mu1)+x[2]*(1-mu2)
         +x[3]*(1-mu3)+x[4]*(1-mu4))
}
###defining the order of true means 
b4=1
b3=2
b2=3
###defining the constraints 
budgetconstraint<-function(x){
  f=NULL
  
  
  
  f=rbind(f,sigma2^2/x[2]+sigma1^2/x[1]-delta^2/(qnorm(1-alpha/(3-b2+1))+
                                                   qnorm(1-beta/b2))^2)
  f=rbind(f,sigma3^2/x[3]+sigma1^2/x[1]-delta^2/(qnorm(1-alpha/(3-b3+1))+
                                                   qnorm(1-beta/b3))^2)
  f= rbind(f,sigma4^2/x[4]+sigma1^2/x[1]-delta^2/(qnorm(1-alpha/(3-b4+1))+
                                                    qnorm(1-beta/b4))^2)
  f=rbind(f,1-x[1])
  f=rbind(f,1-x[2])
  f=rbind(f,1-x[3])
  f=rbind(f,1-x[4])
  
  return(list(ceq = NULL, c = f))
}
library(NlcOptim)
x0=c(10,10,10,10)          #### initial value for numerical optimization 
a=solnl(x0,obj,budgetconstraint,tolX = 1e-2,
        tolFun = 1e-5, tolCon = 1e-5, maxnFun = 1e+10,maxIter = 8000)$par
pi=a/sum(a)               ###obtaining optimal allocation proportion 

a33=(qnorm(1-alpha/(3-b2+1))+qnorm(1-beta/b2))^2*((sigma2^2)/pi[2]+(sigma1^2)/pi[1])/delta^2
b33=(qnorm(1-alpha/(3-b3+1))+qnorm(1-beta/b3))^2*((sigma3^2)/pi[3]+(sigma1^2)/pi[1])/delta^2
c33=(qnorm(1-alpha/(3-b4+1))+qnorm(1-beta/b4))^2*((sigma4^2)/pi[4]+(sigma1^2)/pi[1])/delta^2

n_SRAH=a33    ##  Oracle sample size for SRAH. We report it as the smallest integer greater than n_SRAH. Here a33,b33,c33 should be same.

### Calculation for SBH 


a43=(qnorm(1-alpha/(3-b2+1))+qnorm(1-beta/b2))^2*(sigma2^2+sigma1^2)*4/delta^2
b44=(qnorm(1-alpha/(3-b3+1))+qnorm(1-beta/b3))^2*(sigma3^2+sigma1^2)*4/delta^2
c44=(qnorm(1-alpha/(3-b4+1))+qnorm(1-beta/b4))^2*(sigma4^2+sigma1^2)*4/delta^2

n_SBH=max(a44,b44,c44)  ## oracle sample size for SBH. We report it as the smallest integer greater than n_SBH.














