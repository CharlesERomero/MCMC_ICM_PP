import numpy as np
import scipy as sp

def num_los_int(epsnot,sindex,rmin,rmax,radarr,c=1,ff=1e-3,epsatrmin=0):

  rrmm = (radarr==np.amin(radarr))
  if (radarr[rrmm] == 0) and (sindex > 0):
      radarr[rrmm]=ff

  if rmax < 0:
      if rmin == 0:
          scase=3
      else:
          scase=2
          epsatrmin=1
  else:
      if rmin == 0:
          scase=0
      else:
          scase=1
          epsatrmin=1

  rlen=len(radarr)
  prorad = np.outer(np.ones((rlen,)),radarr**2)
  zrads  = np.transpose(prorad)
#  import pdb; pdb.set_trace()
  radmat = prorad + zrads

  if epsatrmin > 0:
      scrad=radmat/rmin**2
      epsnorm=epsnot
      if scase == 1:
        epsnorm=epsnot*(rmax/rmin)**(sindex)
  else:
      epsnorm=epsnot
      scrad=radarr/rmax**2

  emvals=radmat*0.0
  if rmax > 0:
    c1 = radmat<=rmax**2         # condition 1
  else:
    c1 = radmat>=0                # condition 1
  c2 = radmat>=rmin**2           # condition 2
  c3 = (c1 & c2)  # using a bitwise AND here! Should work as I want.
  emvals[c3] = epsnorm*(scrad[c3])**(-sindex/2.0)

  dz=np.median(np.diff(radarr)) # !!! ASSUMES CONSTANT LINEAR SPACING!

  emint=emvals.sum(axis=0)*dz*2.0
#  import pdb; pdb.set_trace()

  return emint
