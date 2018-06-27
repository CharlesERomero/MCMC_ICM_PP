import numpy as np
import retrieve_data_info as rdi
import astropy.units as u          # Install astropy

### June 20, 2017 - Charles Romero
### I want this to serve as a basic setup code for pressure profiles.
### Eventually, perhaps it can do some optimization of bin choice.

### Delete This on 01 August 2017!
def choose_radial_bins(instrument,nbins=6):
    
    fwhm1,norm1,fwhm2,norm2,fwhm,smfw,freq,Fov = rdi.inst_params(instrument)
    minstep=np.ceil(fwhm/5.0)*5.0 #Make steps be an integer of 5", by default
    binarr = (np.arange(nbins)+1) * minstep
#    inalphas=r_bins*0.0    # The variable needs to be defined...but if it is all zeros, I'll
#                           # overwrite the values.#

    ### I will want to add something more complicated - but let's
    ### Just start with this for now (20 June 2017)
    
    return binarr
    
def a10_gnfw(mycluster,mycosmo,radii):
    """
    Returns a pressure profile based on the "Universal Pressure Profile" (UPP)
    as presented in Arnaud+ 2010.
    
    Parameters
    __________
    my_cluster    - A structure of cluster properties
    mycosmo       - A dictionary of cosmological parameters
    r             - Radial bins (in kpc)
    
    Returns
    -------
    An array (the size of r) of pressure (e.g. a profile). The array will have
    units of keV cm**-3    !!!
    """
    
    P0 = mycluster.P_500 * 8.403 * mycosmo['h_70']**-1.5
    rP = mycluster.R_500 / 1.177 # C_500 = 1.177
#    ra = np.array(radii)              # In case it was a list.
    rf =  radii  / rP                # rP in kpc so rf dimensionless
    a=1.0510; b=5.4905; c=0.3081 # gNFW parameters found in Arnaud+ 2010

    result = (P0 / (((rf)**c)*((1 + (rf)**a))**((b - c)/a)))
    
    return result  ### HAS UNITS ASSOCIATED WITH IT!

### Can I condense some steps into a script??
def get_default_a10_profile(instrument,mycluster,mycosmo,nbins=6):

    # Get the radial bins (these are in arcseconds!):
    arcbins = choose_radial_bins(instrument,nbins=nbins)
### Delete This on 01 August 2017!
    if mycluster.name == 'abell_2146':
        arcbins*=1.3090075 # Perfectly matches SE model input...
#        arcbins*=1.5708075312038152 # Perfectly matches NW model input...
        
    r_bins = arcbins.to("rad")*mycluster.d_a/u.rad # r_bins in kpc
 
    pressure = a10_gnfw(mycluster,mycosmo,r_bins)

    return r_bins, pressure

def get_xymap(map,pixsize,xcentre=[],ycentre=[]):

    nx,ny=map.shape
    ypix = pixsize    # Generally pixel sizes are the same...
    xpix = pixsize    # ""
    if xcentre == []:
        xcentre = nx/2.0
    if ycentre == []:
        ycentre = ny/2.0

#    ny+=1  # Some maps require this...but generally we should not do it.
#    nx+=1  # Some maps require this...but generally we should not do it.
    x = np.outer(np.arange(0,xpix*(nx), xpix)- xpix* xcentre, np.zeros(ny)+1.0)   
    y = np.outer(np.zeros(nx) + 1.0, np.arange(0,ypix*(ny),ypix)- ypix * ycentre)
    return x,y

def get_radial_map(map,pixsize,xcentre=[],ycentre=[]):

    x,y = get_xymap(map,pixsize,xcentre=xcentre,ycentre=ycentre)
    r = np.sqrt(x*x +y*y)

    return r

def get_unitless_arrays(szcu,r_bins,pressure,mapping,mycluster,radmap):

    factors =  szcu['thom_cross']/(szcu['m_e_c2'])  # For SZ integrations (scalar terms)
    prmunit = factors*mycluster.d_a                 # Make a variable for the conversion.
    uless_p = (prmunit*pressure).decompose()        # Get "unitless pressure"
    uless_r = (r_bins/u.kpc).decompose()            # Get "unitless radius" (~kpc)
    t_range = mapping.theta_range*mycluster.d_a/u.kpc
    uless_t = t_range.decompose()                   # Unitless "theta" range (~kpc)
    uless_rm= (radmap*mycluster.d_a/(u.kpc*u.rad)).decompose() # Unitless "theta" map (~kpc)
    ### And make them just values (not "quantitites")
    uless_p=uless_p.value;uless_r=uless_r.value;uless_t=uless_t.value;uless_rm=uless_rm.value

    return uless_p,uless_r,uless_t,uless_rm,factors
