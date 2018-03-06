import scipy.special as sps
import numpy as np

def myrincbeta(x,a,b):
# compute the regularized incomplete beta function.
  if a < 0:
      cbf=(sps.gamma(a)*sps.gamma(b))/sps.gamma(a+b)
      res = (x**a * (1.0-x)**b) / (a * cbf)
      return myrincbeta(x,a+1.0,b) + res
  else:
#      cbf=(sps.gamma(a)*sps.gamma(b))/sps.gamma(a+b)
      cbf=1.0 # sps.betainc is the regularized inc. beta fun.
      res=(sps.betainc(a,b,x) / cbf)
      return res
    
def myredcosine(tmax,n):
# computes \int_0^tmax cos^n(x) dx

  if n < -2:
      res=np.cos(tmax)**(n+1)*np.sin(tmax)/(n+1) 
      return myredcosine(tmax,n+2)*(n+2)/(n+1) - res
  else:
      if n == 0:
          res=tmax
      if n == -1:
          res=np.log(np.absolute(1.0/np.cos(tmax) + np.tan(tmax)) )
      if n == -2:
          res=np.tan(tmax) 

      return res
