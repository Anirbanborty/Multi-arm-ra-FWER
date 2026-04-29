#####Calculation of Oracle sample size for Normal####
alpha=0.05  ###Fix the value of alpha
beta=0.1    ###Fix the value of beta
mu1=1       ###mean of control 
mu2=1       ###mean of experimental treatment 1
mu3=2      ###mean of experimental treatment 2
mu4=2       ###mean of experimental treatment 3

sigma1=2  ###SD of control
sigma2=2  ### SD of experimental treatment 1
sigma3=2   ### SD of experimental treatment 2
sigma4=2    ### SD of experimental treatment 3
var1=sigma1^2  ###### variance of control
var2=sigma2^2  ###### variance of experimental treatment 1
var3=sigma3^2  ###### variance of experimental treatment 2
var4=sigma4^2  ###### variance of experimental treatment 3

h=1           ###fixing the value of h
delta0=0      ###fixing the value of delta0
delta1=1      ###fixing the value of delta1
delta=delta1-delta0

psi1=pnorm((h-mu1)/sigma1)
psi2=pnorm((h-mu2)/sigma2)
psi3=pnorm((h-mu3)/sigma3)
psi4=pnorm((h-mu4)/sigma4)

###calculating allocation proportions 
pi1=sqrt(var1/psi1)*
  sqrt(var2*psi2+
         var3*psi3+
         var4*psi4)
pi2=var2
pi3=var3
pi4=var4

pi11=pi1/(pi1+pi2+pi3+pi4)      ####allocation probability of the control 
pi22=pi2/(pi1+pi2+pi3+pi4)      ########allocation probability of the experimental treatment 1 
pi33=pi3/(pi1+pi2+pi3+pi4)      ########allocation probability of the experimental treatment 2 
pi44=pi4/(pi1+pi2+pi3+pi4)      ########allocation probability of the experimental treatment 3 
p=c(pi11,pi22,pi33,pi44)

a1=(qnorm(1-alpha/(3))+qnorm(1-beta/3))^2*((sigma2^2)/pi22+(sigma1^2)/pi11)/delta^2
b1=(qnorm(1-alpha/(3))+qnorm(1-beta/3))^2*((sigma3^2)/pi33+(sigma1^2)/pi11)/delta^2
c1=(qnorm(1-alpha/(3))+qnorm(1-beta/3))^2*((sigma4^2)/pi44+(sigma1^2)/pi11)/delta^2
n_SRAB=a1      ##  oracle sample size for SRAB. We report it as the smallest integer greater than n_SRAB.Here a1,b1,c1 should be same. 

####SBB###

a2=(qnorm(1-alpha/(3))+qnorm(1-beta/3))^2*((sigma2^2)+(sigma1^2))*4/delta^2
b2=(qnorm(1-alpha/(3))+qnorm(1-beta/3))^2*((sigma3^2)+(sigma1^2))*4/delta^2
c2=(qnorm(1-alpha/(3))+qnorm(1-beta/3))^2*((sigma4^2)+(sigma1^2))*4/delta^2
n_SBB=max(a2,b2,c2)  ## oracle sample size for SBB. We report it as the smallest integer greater than n_SBB.


###Calculation of oracle sample size for SRAH Procedure###
####Define the objective function
obj=function(x){
  return(x[1]*psi1+x[2]*psi2
         +x[3]*psi3+x[4]*psi4)
}
###define the order of true means of experimental treatments 
d4=3
d3=2
d2=1
###define the constraints 
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
x0=c(10,10,10,10)   ###initial value for numerical optimization 
a=solnl(x0,obj,budgetconstraint,tolX = 1e-2,
        tolFun = 1e-5, tolCon = 1e-5, maxnFun = 1e+10,maxIter = 8000)$par       
pi=a/sum(a)       ###obtaining allocation proportion 

a3=(qnorm(1-alpha/(3-d2+1))+qnorm(1-beta/d2))^2*((sigma2^2)/pi[2]+(sigma1^2)/pi[1])/delta^2
b3=(qnorm(1-alpha/(3-d3+1))+qnorm(1-beta/d3))^2*((sigma3^2)/pi[3]+(sigma1^2)/pi[1])/delta^2
c3=(qnorm(1-alpha/(3-d4+1))+qnorm(1-beta/d4))^2*((sigma4^2)/pi[4]+(sigma1^2)/pi[1])/delta^2
 n_SRAH=a3    ##  Oracle sample size for SRAH. We report it as the smallest integer greater than n_SRAH. Here a3,b3,c3 should be same.

### Calculation for SBH 

a4=(qnorm(1-alpha/(3-d2+1))+qnorm(1-beta/d2))^2*(sigma2^2+sigma1^2)*4/delta^2
b4=(qnorm(1-alpha/(3-d3+1))+qnorm(1-beta/d3))^2*(sigma3^2+sigma1^2)*4/delta^2
c4=(qnorm(1-alpha/(3-d4+1))+qnorm(1-beta/d4))^2*(sigma4^2+sigma1^2)*4/delta^2

n_SBH=max(a4,b4,c4)  ## oracle sample size for SBH. We report it as the smallest integer greater than n_SBH.



