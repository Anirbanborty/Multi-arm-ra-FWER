####Oracle sample size calculation for t5 distribution##

alpha=.05            ###Fixing the value of alpha
beta=.1              ###Fixing the value of beta
v=5                  ###Fixing the degree of t distribution as 5
mu1=1                ###mean of the control
mu2=1                ##mean of the experimental treatment 1
mu3=2                ##mean of the experimental treatment 2
mu4=2                ##mean of the experimental treatment 3
sigma1=1.5          ##scale parameter of the control
sigma2=1.5          ##scale parameter of experimental treatment 1
sigma3=1.25          ##scale parameter of experimental treatment 2
sigma4=1.25          ##scale parameter of experimental treatment 3

h=1.5


delta0=0
delta1=1

delta=delta1-delta0





var1=v*sigma1^2/(v-2)       ##variance of the control
var2=v*sigma2^2/(v-2)       ##variance of experimental treatment 1
var3=v*sigma3^2/(v-2)       ##variance of experimental treatment 2
var4=v*sigma4^2/(v-2)       ##variance of experimental treatment 3

psi1=pt((h-mu1)/sigma1,df=v)
psi2=pt((h-mu2)/sigma2,df=v)
psi3=pt((h-mu3)/sigma3,df=v)
psi4=pt((h-mu4)/sigma4,df=v)


pi1=sqrt(var1/psi1)*
  sqrt(var2*psi2+
         var3*psi3+
         var4*psi4)
pi2=var2
pi3=var3
pi4=var4

pi11=pi1/(pi1+pi2+pi3+pi4)         ###allocation proportion of the control
pi22=pi2/(pi1+pi2+pi3+pi4)         ###allocation proportion of experimental treatment 1
pi33=pi3/(pi1+pi2+pi3+pi4)         ###allocation proportion of experimental treatment 1
pi44=pi4/(pi1+pi2+pi3+pi4)         ###allocation proportion of experimental treatment 1
p=c(pi11,pi22,pi33,pi44)

a1=(qnorm(1-alpha/(3))+qnorm(1-beta/3))^2*((var2)/pi22+(var1)/pi11)/delta^2
b1=(qnorm(1-alpha/(3))+qnorm(1-beta/3))^2*((var3)/pi33+(var1)/pi11)/delta^2
c1=(qnorm(1-alpha/(3))+qnorm(1-beta/3))^2*((var4)/pi44+(var1)/pi11)/delta^2

n_SRAB=a1      ##  oracle sample size for SRAB. We report it as the smallest integer greater than n_SRAB.Here a1,b1,c1 should be same. 

####SBB###



a2=(qnorm(1-alpha/(3))+qnorm(1-beta/3))^2*((var2)+(var1))*4/delta^2
b2=(qnorm(1-alpha/(3))+qnorm(1-beta/3))^2*((var1)+(var3))*4/delta^2
c2=(qnorm(1-alpha/(3))+qnorm(1-beta/3))^2*((var4)+(var1))*4/delta^2
n_SBB=max(a2,b2,c2)  ## oracle sample size for SBB. We report it as the smallest integer greater than n_SBB.

######For SRAH##########

##defining the objective function 
obj=function(x){
  return(x[1]*psi1+x[2]*psi2
         +x[3]*psi3+x[4]*psi4)
}
###defining the order of true means 
d4=1
d3=2
d2=3
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

x0=c(10,10,10,10)               ###initial value for numerical optimization
library(NlcOptim) 
a=solnl(x0,obj,budgetconstraint,tolX = 1e-2,
        tolFun = 1e-5, tolCon = 1e-5, maxnFun = 1e+10,maxIter = 8000)$par
pi=a/sum(a)       ###obtaining optimal allocation proportion 

a33=(qnorm(1-alpha/(3-d2+1))+qnorm(1-beta/d2))^2*((var2)/pi[2]+(var1)/pi[1])/delta^2
b33=(qnorm(1-alpha/(3-d3+1))+qnorm(1-beta/d3))^2*((var3)/pi[3]+(var1)/pi[1])/delta^2
c33=(qnorm(1-alpha/(3-d4+1))+qnorm(1-beta/d4))^2*((var4)/pi[4]+(var1)/pi[1])/delta^2

n_SRAH=a33    ##  Oracle sample size for SRAH. We report it as the smallest integer greater than n_SRAH. Here a33,b33,c33 should be same.

### Calculation for SBH 

a44=(qnorm(1-alpha/(3-d2+1))+qnorm(1-beta/d2))^2*(var2+var1)*4/delta^2
b44=(qnorm(1-alpha/(3-d3+1))+qnorm(1-beta/d3))^2*(var3+var1)*4/delta^2
c44=(qnorm(1-alpha/(3-d4+1))+qnorm(1-beta/d4))^2*(var4+var1)*4/delta^2

n_SBH=max(a44,b44,c44)  ## oracle sample size for SBH. We report it as the smallest integer greater than n_SBH.




