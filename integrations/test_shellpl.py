from shell_pl import shell_pl
import numpy as np
import matplotlib.pyplot as plt
import reduce_integrals as ri
import scipy.special as sps 
import num_los_int as nli

radarr=(np.arange(200)+1.0)/0.0005

#sindex=np.random.normal(0.0,2.0)
rmin=np.random.normal(2.0,0.5)
sindex=1.1
rmin=2.1
#rmax=np.random.normal(4.0,0.5)
rmax = -1.1
epsnot=0.001
print (rmin,rmax,sindex)
djr=nli.num_los_int(epsnot,sindex,rmin,rmax,radarr)
abc=shell_pl(epsnot,sindex,rmin,rmax,radarr)
plt.plot(radarr,abc)
plt.plot(radarr,djr)
plt.ylabel('integrated value')
plt.xlabel('radius')
plt.show()

mygi=(abc>0)
myratio=djr[mygi]/abc[mygi]

p = sindex/2.0 # e(r) = e_0 * (r^2)^(-p) for this notation / program
c1 = radarr<=rmax               # condition 1
sor=(radarr[c1]/rmax)           # scaled outer radii
tmax=np.arctan(np.sqrt(1.0 - sor**2)/sor)
sindex=0.0
abc=shell_pl(epsnot,sindex,rmin,rmax,radarr)
#abc=shell_pl(1.0,2.0,2.0,4.0,radarr)
#ri.myrincbeta(sor**2,p-0.5,0.5)
#ri.myredcosine(tmax,2.0*p-2.0)


n=-2
ri.myredcosine(tmin,tmax,n)
