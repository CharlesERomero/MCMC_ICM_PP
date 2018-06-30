import numpy as np
#import scipy.constants as spconst
#import astropy.constants as const
import astropy.units as u
import cosmolopy.distance as cd
from scipy.interpolate import interp1d
import scipy.signal
from numpy import sqrt
import scipy.ndimage
import image_filtering as imf
from astropy.io import fits
from cosmolopy import *
### 
from shell_pl import shell_pl
import retrieve_data_info as rdi

def log_profile(args,r_bins,radii,alphas=[],rintmax=[]):

    mybins=np.append([0],r_bins)
    #mybins=[0] + r_bins
    mybins[-1]=-1
    presprof=np.zeros(len(radii))
    #        ifht=('alphas' in dir())
    ifht=('alphas' in locals())
    if ifht: 
        if sum(alphas) != 0:
            aset=1
        else:
            aset=0
            alphas = np.zeros(len(r_bins))
    else:
        aset=0
        alphas = np.zeros(len(r_bins))
        
        #        print len(r_bins),len(alphas)
    ycyl=0.0 # Y_500,cyl integration.
    ### Actually, I am not calculating Ycyl here!!! Stupid me. I could 
    ### calculate Ysph here though...

    for idx, val in enumerate(r_bins):
        rin=mybins[idx]
        rout=mybins[idx+1]
        epsnot=args[idx]
        #            print idx,r_bins,len(r_bins),alphas,aset
        alpha=alphas[idx]
        if rin == 0:
            lr=np.log10(mybins[idx+2]/mybins[idx+1])
            lp=np.log10(args[idx+1]/args[idx])
            if aset == 0:
                alpha=-lp/lr
                #                myind=np.where((radii < rout) & (radii >= rin))
            myind=(radii < rout) & (radii >= rin)
            myrad=radii[myind]
            mypres=epsnot*(myrad/rout)**(-alpha)
            yint = 1.0 - (rin/rout)**(2-alpha)  # I could leave out the second term...
            # but this ensures an error if alpha >2 ...
            rnot=rout
        elif rout == -1:
            lr=np.log10(r_bins[idx]/r_bins[idx-1])
            lp=np.log10(args[idx]/args[idx-1])
            if aset == 0:
                alpha=-lp/lr
            epsnot=args[idx-1]
            #                myind=(radii >= rin)
            myind=np.where(radii >= rin)
            myrad=radii[myind]
            mypres=epsnot*(myrad/rin)**(-alpha)
            rnot=rin
            yint = -1.0
            if np.sum(rintmax) > 0:
                yint = (rintmax/rnot)**(2-alpha) - 1.0

        else:
            lr=np.log10(mybins[idx+1]/mybins[idx])
            lp=np.log10(args[idx]/args[idx-1])
            if aset == 0:
                alpha=-lp/lr
            myind=np.where((radii < rout) & (radii >= rin))
            #                myind=(radii < rout) & (radii >= rin)
            myrad=radii[myind]
            mypres=epsnot*(myrad/rout)**(-alpha)
            rnot=rin
            yint = (rout/rin)**(2-alpha) - 1.0

        presprof[myind]=mypres
        if aset == 0:
            alphas[idx]=alpha

        ypref = 2*np.pi*epsnot*(rnot**2)/(2-alpha)
        if np.sum(rintmax) > 0:
            if rin < rintmax:
                if (rout > 0) & (rout <= rintmax):
                    yint=(rintmax/rnot)**(2-alpha)-1.0
                ycyl=ycyl + ypref*yint

            return presprof,alphas,ycyl
        # back to this placent
        
    if aset == 0:
        return presprof,alphas
    else:
        return presprof

def binsky(args,r_bins,theta_range,theta,inalphas=[]):
    """
    Returns a surface brightness map for a binned profile, slopes, and radial integrals.
    
    Parameters
    __________
    args :  Pressure for each bin used
    Returns
    -------
    out: numpy.ndarray
    """
    Int_Pres,alphas,integrals = analytic_shells(r_bins,args,theta_range,alphas=inalphas)
    fint = interp1d(theta_range, Int_Pres, bounds_error = False, 
                    fill_value = 0)
    nx, ny = theta.shape
    map = np.float64(fint(theta.reshape(nx * ny))) # Type 17 = float? (Implicitly float 32?)
    map = map.reshape(nx,ny)

    return map,alphas,integrals

def prep_SZ_binsky(pressure, temp_iso, geoparams=None):
    """
    geoparams   :  [X_shift, Y_shift, Rotation, Ella*, Ellb*, Ellc*, Xi*, Opening Angle]
    """
    edensity = np.array(pressure) / temp_iso
    etemperature = np.array(pressure)*0 + temp_iso
    if geoparams == None:
        geoparams = [0,0,0,1,1,1,0,0] # Spherical Geometry

    return edensity, etemperature, geoparams

def integrate_profiles(density_proxy, etemperature, geoparams,r_bins,theta_range,hk,dv,inalphas=[],
                       beta=0.0,betaz=None,finint=False,narm=False,fixalpha=False,strad=False,
                       array="2",fullSZcorr=False,SZtot=False,columnDen=False,Comptony=True):
    """
    Returns a surface brightness map for a binned profile fit, with far more generality than previously done.
    
    Parameters
    __________
    density_proxy :  The electron density * boltzmann constant * kpc / m_e c**2
                     Its units are such that the integral of (density_proxy*etemperature) over theta_range
                     (itself in radians), results in the unitless Compton y parameter.
    etemperature  :  The electron temperature (k_B * T), but again, without units within Python
    geoparams     :  [X_shift, Y_shift, Rotation, Ella*, Ellb*, Ellc*, Xi*, Opening Angle]
    r_bins        :  The (elliptical) bins for the profile. 
    theta_range   :  The range of angles for which to create a 1D profile (which can then be interpolated)
    hk            :  The housekeeping class ("structure")
    dv            :  The data-varaibles class ("structure")
    inalphas      :  Nothing to see here. Move along.
    beta          :  Fraction of the speed of light of the cluster bulk (peculiar) motion.
    betaz         :  Fraction of the speed of light of the cluster along the line of sight.
    finint        :  Integrate out to last finite (defined) bin.
    narm          :  Normalized at R_Min. This is important for integrating shells.
    fixalpha      :  Fix alpha (to whatever inalpha is); useful for shocks.
    strad         :  STrict RADii; if the pressure model has to obey strict placements of radii, use this!

    Notes
    __________
    * Ella should be set to 1. Therefore, define Ellb relative to Ella (and likewise with Ellc)
    * Xi is a parameterization in a forthcoming memo (July 2017, CR)
    
    """
    if betaz == None:
        betaz = beta
### If geoparams[6] > 0, then we are modelling some non-ellipsoid...perhaps a shock. If the opening angle
### is not set, then this will create a bimodal (bipolar) model component, which we almost certainly don't
### want. If we do want a bimodal component, then I think a better override is to use geoparams[7]= 2 pi.

    #eff_pres = np.zeros(len(etemperature)); y_press= np.zeros(len(etemperature))
    if Comptony == True:
        vals  = density_proxy*etemperature
    if columnDen== True:
        vals  = density_proxy*dv.av.szcv["m_e_c2"]/dv.av.szcv["boltzmann"]
    if SZtot == True:
        vals = dv.tSZ*density_proxy*etemperature + dv.kSZ*density_proxy*dv.av.szcv["m_e_c2"]
                       
    if fullSZcorr == True:
        for i in range(len(etemperature)):
            tSZ,kSZ = rdi.get_sz_bp_conversions(etemperature[i],hk.hk_ins.instrument,array=array,inter=False,
                                                beta=beta,betaz=betaz,rel=True)
            vals[i] = tSZ*density_proxy[i]*etemperature[i] + kSZ*density_proxy[i]*dv.av.szcv["m_e_c2"]

    Int_Pres,alphas,integrals = analytic_shells(r_bins,vals,theta_range,alphas=inalphas,
                                                shockxi=geoparams[6],finite=finint,narm=narm,
                                                fixalpha=fixalpha,strad=strad)
    
    return Int_Pres,alphas,integrals
        
def general_gridding(xymap,theta_range,r_bins,geoparams,finite=False,taper='normal',
                     integrals=0,Int_Pres=0,ell_int=0,tap_int=0,oldvs=False):
    """
    Returns a surface brightness map for a binned profile fit, with far more generality than previously done.
    
    Parameters
    __________
    xymap         :  A tuple (x,y) where x and y are grids of their respective coordinates in << arceconds >>

    Notes
    __________
    * Ella should be set to 1. Therefore, define Ellb relative to Ella (and likewise with Ellc)
    * Xi is a parameterization in a forthcoming memo (July 2017, CR)
    
    """

    #print geoparams[6], (geoparams[6] > 0.0)
    #import pdb;pdb.set_trace()

    
    if geoparams[6] > 0.0:
        x,y = xymap;  mymap = np.zeros(x.shape); myrs = r_bins
        #taper = hk.hk_ins.taper      # Required field...
        if geoparams[7] == 0:
            geoparams[7] = np.pi 
        if finite == True:
            myrs = myrs[:-1]
        #for idx, val in enumerate(myrs):
        #    if val == 0: val=r_bins[idx+1] # Correct units? I think so.
        #    if taper == 'inverse':
        #        mymap += grid_profile(theta_range, ell_int[idx,:], xymap, geoparams=geoparams)
        #        mymap -= grid_profile(theta_range, tap_int[idx,:], xymap, geoparams=geoparams,myscale=val,axis='y')
        #    else:
        #        mymap += grid_profile(theta_range, integrals[idx,:], xymap, geoparams=geoparams,myscale=val,axis='x')
        ######################################################################################
        ### The following has been rewritten 30 Mar 2018, in hopes of being faster.
        if myrs[0] == 0: myrs[0]=myrs[1] # Correct units? I think so.
        if taper == 'inverse':
            for my_int_add, my_int_sub, val in zip(ell_int, tap_int,myrs):
                mymap += grid_profile(theta_range, my_int_add, xymap, geoparams=geoparams)
                mymap -= grid_profile(theta_range, my_int_sub, xymap, geoparams=geoparams,myscale=val,axis='y')
        else:
            if oldvs == True:
                for my_int, val in zip(integrals, myrs):
                    mymap += grid_profile(theta_range, my_int, xymap, geoparams=geoparams,myscale=val,axis='x')
            else:
                mymap=iter_grid_profile_v2(integrals, myrs, theta_range, xymap, geoparams=geoparams,axis='x')
         ######################################################################################
    else:  
        mymap = grid_profile(theta_range, Int_Pres, xymap, geoparams=geoparams,myscale=1.0)

### 03 August 2017 - WTF?????
    ###mymap = np.transpose(mymap)

    return mymap


#########################################################################################################
### + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ###
###+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +###
### + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ###
###+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +###
### + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ###
###+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +###
### + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ###
###                                                                                                   ###
###                         Let's try to do things in a more general way                              ###
###                                                                                                   ###
###+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +###
### + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ###
###+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +###
### + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ###
###+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +###
### + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ###
###+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +###
#########################################################################################################


def binsky_SZ_general(edensity, etemperature, geoparams,r_bins,theta_range,xymap,hk,dv,
                      inalphas=[],beta=0.0,betaz=None,finite=False,narm=False,fixalpha=False,
                      strad=False,array="2"):
    """
    Returns a surface brightness map for a binned profile fit, with far more generality than previously done.
    
    Parameters
    __________
    edensity    :  The electron density (no units in Python, but otherwise should be in cm**-3)
    etemperature:  The electron temperature (k_B * T), but again, without units within Python
    geoparams   :  [X_shift, Y_shift, Rotation, Ella*, Ellb*, Ellc*, Xi*, Opening Angle]
    r_bins      :  The (elliptical) bins for the profile. 
    theta_range :  The range of angles for which to create a 1D profile (which can then be interpolated)
    xymap       :  A tuple (x,y) where x and y are grids of their respective coordinates in << arceconds >>
    hk          :  The housekeeping class ("structure")
    dv          :  The data-varaibles class ("structure")
    inalphas    :  Nothing to see here. Move along.
    beta        :  Fraction of the speed of light of the cluster bulk (peculiar) motion.
    betaz       :  Fraction of the speed of light of the cluster along the line of sight.
    finite      :  Integrate out to last finite (defined) bin.
    narm        :  Normalized at R_Min. This is important for integrating shells.
    strad       :  STrict RADii; if the pressure model has to obey strict placements of radii, use this!

    Notes
    __________
    * Ella should be set to 1. Therefore, define Ellb relative to Ella (and likewise with Ellc)
    * Xi is a parameterization in a forthcoming memo (July 2017, CR)
    
    """
    if betaz == None:
        betaz = beta
### If geoparams[6] > 0, then we are modelling some non-ellipsoid...perhaps a shock. If the opening angle
### is not set, then this will create a bimodal (bipolar) model component, which we almost certainly don't
### want. If we do want a bimodal component, then I think a better override is to use geoparams[7]= 2 pi.

    taper='normal'
    if geoparams[6] > 0:
        taper = hk.hk_ins.taper      # Required field...
        if geoparams[7] == 0:
            geoparams[7] = np.pi 

    eff_pres = np.zeros(len(etemperature))
    for i in range(len(etemperature)):
        tSZ,kSZ = rdi.get_sz_bp_conversions(etemperature[i],hk.hk_ins.instrument,array=array,inter=False,
                                        beta=beta,betaz=betaz,rel=True)
        ### 05 July 2017 - I need to check the proper factors for the kSZ term. (CR)
#        print 'tSZ and kSZ values are: ',tSZ, kSZ
        eff_pres[i] = tSZ*edensity[i]*etemperature[i] + kSZ*edensity[i]*dv.av.szcv["m_e_c2"]
    
    map,alphas,integrals = binsky_general(eff_pres,geoparams,r_bins,theta_range,xymap,inalphas=inalphas,
                                          finite=finite,narm=narm,taper=taper,fixalpha=fixalpha,strad=strad)

    return map,alphas,integrals
        
def binsky_general(vals,geoparams,r_bins,theta_range,xymap,inalphas=[],
                   finite=False,narm=False,taper='normal',fixalpha=False,strad=False):
    """
    Returns a surface brightness map for a binned profile fit 
    
    Parameters
    __________
    vals      :  Pressure for each bin used
    geoparams :  [X_shift, Y_shift, Rotation, Ella*, Ellb*, Ellc*, Xi*, Opening Angle]
    
    Notes:
    __________
    --> We should consider Ella to be RESTRICTED to 1. That is, Ellb and Ellc should always be calculated
    relative to the x-axis parameter.
    
    Returns
    -------
    An map that accounts for a range of geometrical restrictions. The integrals may not be applicable.

    """

    if taper == 'inverse':
        Ell_Pres,alphas,ell_int = analytic_shells(r_bins,vals,theta_range,alphas=inalphas,
                                                    finite=finite,narm=narm,fixalpha=fixalpha,strad=strad)
        Tap_Pres,tap_alph,tap_int = analytic_shells(r_bins,vals,theta_range,alphas=inalphas,
                                                    shockxi=geoparams[6],finite=finite,narm=narm,
                                                    fixalpha=fixalpha,strad=strad)
        integrals = ell_int - tap_int
    else:
        Int_Pres,alphas,integrals = analytic_shells(r_bins,vals,theta_range,alphas=inalphas,
                                                    shockxi=geoparams[6],finite=finite,narm=narm,
                                                    fixalpha=fixalpha,strad=strad)

############################################################################

#    if geoparams[6] > 0:
#        x,y = xymap;  map = x*0.0; myrs = r_bins
#        if finite == True:
#            myrs = myrs[:-1]
#        for idx, val in enumerate(myrs):
#            if val == 0: val=r_bins[idx+1] # Correct units? I think so.
#            if taper == 'inverse':
#                map += grid_profile(theta_range, ell_int[idx,:], xymap, geoparams=geoparams)
#                map -= grid_profile(theta_range, tap_int[idx,:], xymap, geoparams=geoparams,myscale=val,axis='y')
#            else:
#                map += grid_profile(theta_range, integrals[idx,:], xymap, geoparams=geoparams,myscale=val,axis='x')
#    else:
#        map = grid_profile(theta_range, Int_Pres, xymap, geoparams=geoparams,myscale=1.0)

### 03 August 2017 - WTF?????
#    map = np.transpose(map)     
#    import pdb; pdb.set_trace()

    map = general_gridding(xymap,r_bins,geoparams,finite,narm,taper,strad,
                           integrals,Int_Pres,ell_int,tap_int)

    return map,alphas,integrals

def grid_profile(theta_range, profile, xymap, geoparams=[0,0,0,1,1,1,0,0],myscale=1.0,axis='z'):

    ### Get new grid:
    (x,y) = xymap
    x,y = rot_trans_grid(x,y,geoparams[0],geoparams[1],geoparams[2])
    x,y = get_ell_rads(x,y,geoparams[3],geoparams[4])
    radmap = np.sqrt(x**2 + y**2)
    theta = radmap*(u.arcsec).to("radian");  theta_min = np.min(theta_range)
#    import pdb;pdb.set_trace()
    bi=np.where(theta < theta_min);   theta[bi]=theta_min

    nx, ny = theta.shape
    fint = interp1d(theta_range, profile, bounds_error = False, fill_value = 0)
    mymap = np.float64(fint(theta.reshape(nx * ny))) # Type 17 = float? (Implicitly float 32?)
    mymap = mymap.reshape(nx,ny)
    ### And a couple more *necessary* modification:
    ### Where we want to scale it by a certain r_bin, given in radians. We also want to scale by "Ella", if axis='x':
    if axis == 'x':
        xell = (x/(geoparams[3]*myscale))*(u.arcsec).to("radian") # x is initially presented in arcseconds
        modmap = geoparams[5]*(xell**2)**(geoparams[6]) # Consistent with model creation??? (26 July 2017)
    if axis == 'y':
        yell = (y/(geoparams[4]*myscale))*(u.arcsec).to("radian") # x is initially presented in arcseconds
        modmap = geoparams[5]*(yell**2)**(geoparams[6]) # Consistent with model creation??? (26 July 2017)
    if axis == 'z':
        modmap = mymap*0.0 + geoparams[5]      # Just the plain old LOS elongation factor
    
    mymap = mymap*modmap   # Very important to be precise here.
#    print np.min(xell),np.max(xell),np.min(map),np.max(map)
#    import pdb;pdb.set_trace()
    if geoparams[7] > 0:
#        print 'Opening Angle Value Nonzero. Applying the cut.'
#        angmap = np.arctan2(x,y)
        angmap = np.arctan2(y,x)
        #        gi = np.where(abs(angmap) < geoparams[7]/2.0)
        bi = np.where(abs(angmap) > geoparams[7]/2.0)
        #        import pdb; pdb.set_trace()
        mymap[bi] = 0.0

    return mymap

def get_ell_rads(x,y,ella,ellb):

    xnew = x/ella ; ynew = y/ellb

    return xnew, ynew
    
def rot_trans_grid(x,y,xs,ys,rot_rad):

    xnew = (x - xs)*np.cos(rot_rad) + (y - ys)*np.sin(rot_rad)
    ynew = (y - ys)*np.cos(rot_rad) - (x - xs)*np.sin(rot_rad)

    return xnew,ynew

def ycylfromprof(Int_Pres,theta_range,theta_max):
    """
    Remember, theta_range is in radians
    """
    i=1 # Start at the second entry!
    Ycyl=0
    while i < len(theta_range):
        if theta_range[i] < theta_max:
            dtheta = theta_range[i]-theta_range[i-1]
            Yshell = 2.0*np.pi*theta_range[i]*dtheta*Int_Pres[i]
            Ycyl = Ycyl + Yshell
        i+=1

    return Ycyl

def analytic_shells(r_bins,vals,theta,correl=False,alphas=[],shockxi=0.0,fixalpha=False,
                    finite=False,narm=False,strad=False):
    """
    Returns an integrated map of some signal along the line of sight. This routine
    assumes that the pressure within a shell has a power law distribution.
    
    Parameters
    __________
    r_bins   : The radial bins (in radian, I believe)
    vals     :  Pressure for each bin used
    theta    : An array of radii (in radian) in the map, which will be used for gridding the model
    [correl] : Correlate?
    [alphas] : An array of power laws (indices) for 3d pressure distribution
    [shockxi]: Polar tapering, if used in a shock model.
    [finite] : Set this keyword if you do NOT want to integrate to infinity.
    [narm]   :  Normalize at R_min (within a bin)
    [strad]  : STrict RADii. When using a shock model (e.g. Abell 2146), where specific radii,
               ESPECIALLY inner radii are defined, this keyword SHOULD be set! Note that if
               the finite keyword is set, then this does not need to be set. 

    Returns
    -------
    out: numpy.ndarray
    Map convolved with the beam.          
    """
    if finite == False:
        if np.min(r_bins) != 0 and strad == False:
            mybins=np.append([0],r_bins)
            mybins[-1]=-1
        else:
            mybins=np.append(r_bins,-1)
    else:
        mybins = r_bins
        r_bins = r_bins[:-1]
            
#    import pdb; pdb.set_trace()
    nthetas = len(theta)
    integrals = np.zeros((len(r_bins),nthetas))
#########################################################################################
### TO REMOVE
    #        theta_val = (theta.reshape(nx * ny) * u.arcmin).to("radian").value
#    ifht=('alphas' in dir())
#
#    if ifht: 
#        if np.sum(alphas) != 0:
#            aset=1
#        else:
#            aset=0
#            alphas = np.zeros(len(r_bins))
#    else:
#        aset=0
#        alphas = np.zeros(len(r_bins))

#        print "Shape of bins: ", bins
#        print "Shape of vals: ", vals.shape
#########################################################################################
    if fixalpha == False:
        alphas = np.zeros(len(r_bins))

    for idx, myval in enumerate(r_bins):
        rin=mybins[idx]
        rout=mybins[idx+1]
        mypressure=vals[idx] # Gah, what a stupid way to do this.
        
        if fixalpha == False:                  
            if rin == 0:
                lr=np.log10(mybins[idx+2]/mybins[idx+1])
                lp=np.log10(vals[idx+1]/vals[idx])
                alphas[idx]=-lp/lr
            elif rout == -1:
                lr=np.log10(r_bins[idx]/r_bins[idx-1])
                lp=np.log10(vals[idx]/vals[idx-1])
#                    lr=np.log10(mybins[idx]/mybins[idx-1])
#                    lp=np.log10(vals[idx-1]/vals[idx-2])
                alphas[idx]=-lp/lr
                mypressure=vals[idx-1]
            else:
                lr=np.log10(mybins[idx+1]/mybins[idx])
                lp=np.log10(vals[idx]/vals[idx-1])
                alphas[idx]=-lp/lr
 
#        print '#-#--#---#----#-----#------#-------#------#-----#----#---#--#-#'
#        print alphas[idx], idx, rin, rout, shockxi
### Beware of 2.0*shockxi!!! (26 July 2017)
        #import pdb;pdb.set_trace()
        integrals[idx] = shell_pl(mypressure,alphas[idx]+2.0*shockxi,rin,rout,theta,narm=narm) #R had been in here.
#        print integrals[idx,0]
                
    totals = np.sum(integrals,axis=0)  # This should accurately produce Compton y values.

    return totals,alphas,integrals

##########################################################################################
##### I don't think I need the following module, but, I'll leave it for now          #####
##########################################################################################

def analytic_shock(r_bins,vals,alphas,theta,shockxi):
    """
    CURRENTLY UNUSED (September 2017)
    """
    
    mybins=np.append(r_bins,-1);    mybins[-1]=-1
    nthetas = len(theta)
    integrals = np.zeros((len(r_bins),nthetas))

    for idx, myval in enumerate(r_bins):
        rin=mybins[idx]
        rout=mybins[idx+1]
        mypressure=vals[idx] # Gah, what a stupid way to do this.
        integrals[idx] = shell_pl(mypressure,alphas[idx]+shockxi,rin,rout,theta) #R had been in here.
        
    totals = np.sum(integrals,axis=0)  # This should accurately produce Compton y values.

    return totals,integrals

##########################################################################################
##### End unecessary module. (19 July 2017)                                          #####
##########################################################################################

def shell_correl(integrals,r_bins,theta):

    nrad,nbin=integrals.shape
    mybins=np.append([0],r_bins)
#    mybins=[0] + r_bins
    mybins[-1]=-1
    avgs=np.zeros(len(r_bins),len(r_bins))
    for idx, val in enumerate(r_bins):
        for idy, val in enumerate(r_bins):
            rin=mybins[idy]
            rout=mybins[idy+1]
            myind=np.where((theta < rout) & (theta >= rin))
            avgs[idx,idy]=np.mean(integrals[idx,myind])      
            
    return avgs

def iter_grid_profile(integrals, myrs, theta_range, xymap, geoparams=[0,0,0,1,1,1,0,0],axis='z'):
    """
    This largely copies the functionality of grid_profile, but is designed to be much faster for
    iterative applications (same geoparams)
    """

    ### Get new grid:
    (x,y) = xymap
    x,y = rot_trans_grid(x,y,geoparams[0],geoparams[1],geoparams[2])
    x,y = get_ell_rads(x,y,geoparams[3],geoparams[4])
    radmap = np.sqrt(x**2 + y**2)
    theta = radmap*(u.arcsec).to("radian");  theta_min = np.min(theta_range)
    bi=np.where(theta < theta_min);   theta[bi]=theta_min
    nx, ny = theta.shape

    rsRads = myrs / np.min(myrs)   # This should be unitless from Python's perspective, but really in arcseconds.
    
    ### And a couple more *necessary* modification:
    ### Where we want to scale it by a certain r_bin, given in radians. We also want to scale by "Ella", if axis='x':
    if axis == 'x':
        xell = (x/(geoparams[3]*np.min(myrs)))*(u.arcsec).to("radian") # x is initially presented in arcseconds
        modmap = geoparams[5]*(xell**2)**(geoparams[6]) # Consistent with model creation??? (26 July 2017)
    if axis == 'y':
        yell = (y/(geoparams[4]*np.min(myrs)))*(u.arcsec).to("radian") # x is initially presented in arcseconds
        modmap = geoparams[5]*(yell**2)**(geoparams[6]) # Consistent with model creation??? (26 July 2017)
    if axis == 'z':
        modmap = map*0.0 + geoparams[5]      # Just the plain old LOS elongation factor

    mymap = np.zeros(x.shape)
    for profile, myscale in zip(integrals, rsRads):
        
        fint = interp1d(theta_range, profile, bounds_error = False, fill_value = 0)
        map = np.float64(fint(theta.reshape(nx * ny))) # Type 17 = float? (Implicitly float 32?)
        map = map.reshape(nx,ny)
        map = map * modmap * (myscale**(-2*geoparams[6]))
        mymap+=map
        
    if geoparams[7] > 0:
        angmap = np.arctan2(y,x)
        bi = np.where(abs(angmap) > geoparams[7]/2.0)
        mymap[bi] = 0.0

    return mymap

def iter_grid_profile_v2(integrals, myrs, theta_range, xymap, geoparams=[0,0,0,1,1,1,0,0],axis='z'):
    """
    This largely copies the functionality of grid_profile, but is designed to be much faster for
    iterative applications (same geoparams)
    """

    ### Get new grid:
    (x,y) = xymap
    x,y = rot_trans_grid(x,y,geoparams[0],geoparams[1],geoparams[2])
    x,y = get_ell_rads(x,y,geoparams[3],geoparams[4])
    radmap = np.sqrt(x**2 + y**2)
    theta = radmap*(u.arcsec).to("radian");  theta_min = np.min(theta_range)
    bi=np.where(theta < theta_min);   theta[bi]=theta_min
    nx, ny = theta.shape

    rsRads = myrs / np.min(myrs)   # This should be unitless from Python's perspective, but really in arcseconds.
    
    ### And a couple more *necessary* modification:
    ### Where we want to scale it by a certain r_bin, given in radians. We also want to scale by "Ella", if axis='x':
    if axis == 'x':
        xell = (x/(geoparams[3]*np.min(myrs)))*(u.arcsec).to("radian") # x is initially presented in arcseconds
        modmap = geoparams[5]*(xell**2)**(geoparams[6]) # Consistent with model creation??? (26 July 2017)
    if axis == 'y':
        yell = (y/(geoparams[4]*np.min(myrs)))*(u.arcsec).to("radian") # x is initially presented in arcseconds
        modmap = geoparams[5]*(yell**2)**(geoparams[6]) # Consistent with model creation??? (26 July 2017)
    if axis == 'z':
        modmap = map*0.0 + geoparams[5]      # Just the plain old LOS elongation factor

    Int_Prof = np.zeros((integrals.shape[1]))

    for profile, myscale in zip(integrals, rsRads):
        Int_Prof+=  (myscale**(-2*geoparams[6])) * profile
        
    fint = interp1d(theta_range, Int_Prof, bounds_error = False, fill_value = 0)
    mymap = np.float64(fint(theta.reshape(nx * ny))) # Type 17 = float? (Implicitly float 32?)
    mymap = mymap.reshape(nx,ny)
    mymap = mymap * modmap 
        
    if geoparams[7] > 0:
        angmap = np.arctan2(y,x)
        bi = np.where(abs(angmap) > geoparams[7]/2.0)
        mymap[bi] = 0.0

    return mymap
