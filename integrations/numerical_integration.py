import numpy as np
import scipy as sp
from scipy.interpolate import interp1d
import astropy.units as u
import ellipsoidal_shells as es
import os
import retrieve_data_info as rdi

def int_profile(profrad, profile,radProjected,zmax=0):
    """
    This currently only integrates out to the max of profrad. If you want
    to give a *fixed z*, you should be sure it is LESS THAN the max of
    profrad, and then adjust the code below.

    You likely want profrad to be in kpc. In this way, you will integrate
    units of pressure over kpc, and the resultant units are comprehensible.

    radProjected is an array of z-coordinates along the line of sight.
    """
    
    nrP = len(radProjected); nPr=len(profrad)
    x = np.outer(profrad,np.zeros(nrP)+1.0)
    z = np.outer(np.zeros(nPr)+1.0,radProjected)
    rad = np.sqrt(x**2 + z**2)
    fint = interp1d(profrad, profile, bounds_error = False, fill_value = 0)
    radProfile = fint(rad.reshape(nrP*nPr))
    if zmax > 0:
        zre = z.reshape(nrP*nPr); settozero = (zre > zmax)
        radProfile[settozero] = 0.0
    foo =np.diff(z); bar =foo[:,-1];peloton=radProfile.reshape(nPr,nrP)
    diffz = np.insert(foo,-1,bar,axis=1)
    intProfile = 2.0*np.sum(radProfile.reshape(nPr,nrP)*diffz,axis=1)
    
    return intProfile

def create_profile_alla_cer(hk,dv,efv,nw=True,plot=False,finite=False):

    prmunit = (efv.prmunit*u.keV).to("cm**3")
    uless_r, edensity, etemperature, geoparams, inalphas = es.prep_a2146_binsky(hk.hk_ins,nw=nw)
    edensity = edensity*(prmunit.value) # Incorportate all the relevant factors. (~3.16)
    uless_rad = (uless_r/dv.cluster.d_a).value
    kpc_range = (dv.mapping.theta_range*dv.cluster.d_a).value;
    if finite == False:
        myrs=np.append(uless_r,500.0)
    else:
        myrs = uless_r
        
    if nw == True:
        ealp = hk.hk_ins.ealp_nw; talp = hk.hk_ins.talp_nw
    else:
        ealp = hk.hk_ins.ealp_se; talp = hk.hk_ins.talp_se

    myprof = kpc_range*0.0; sindex = ealp+talp
    for idx in range(len(myrs)-1):
        epsnot = edensity[idx]*etemperature[idx]; rmin = myrs[idx]; rmax = myrs[idx+1]
        tSZ,kSZ,int_factors = get_SZ_factors(etemperature[idx],dv,hk,efv,beta=0.0,betaz=0.0)
        if idx == 0:
            rmin = rmax/1000.0
        gi = (kpc_range > rmin) & (kpc_range < rmax)
        if idx == 0:
            myprof[gi] = epsnot*(kpc_range[gi]/rmax)**sindex[idx]
        else:
            myprof[gi] = (kpc_range[gi]/rmin)**sindex[idx]
        myprof[gi]*= int_factors*tSZ

    myind = np.where(myprof != 0.0)
    myprof = myprof[myind]; kpc_range = kpc_range[myind]
        
    if plot == True:
        import matplotlib.pyplot as plt
        plt.clf();  plt.figure(1); fig1,ax1 = plt.subplots()
        ax1.plot(kpc_range,-myprof);ax1.set_yscale("log")
        plt.axvline(uless_r[1],color="k", linestyle ="dashed")
        plt.axvline(uless_r[2],color="k", linestyle ="dashed")
        plt.xlabel("Radius (kpc)");  plt.title("Profile")
        filename = "Abell_2146_profile_checker.png"
        fullpath = os.path.join(hk.hk_outs.newpath,hk.hk_outs.prefilename+filename)
        plt.savefig(fullpath)

    return kpc_range,myprof

def get_SZ_factors(temp,dv,hk,efv,beta=0.0,betaz=0.0):

    int_factors = (efv.factors*u.cm)*u.kpc.to("cm"); int_factors = int_factors.value
    tSZ,kSZ = rdi.get_sz_bp_conversions(temp,hk.hk_ins.instrument,
                                        array="2",inter=False,beta=beta,
                                        betaz=betaz,rel=True)
    return tSZ,kSZ,int_factors
