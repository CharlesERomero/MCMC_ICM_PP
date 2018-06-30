import numpy as np
import retrieve_data_info as rdi
import astropy.units as u          # Install astropy

def a10_gnfw(P500,R500,mycosmo,radii):
    """
    Returns a pressure profile based on the "Universal Pressure Profile" (UPP)
    as presented in Arnaud+ 2010.
    
    Parameters
    __________
    P500          - The representative pressure for a gNFW profile
    R500          - R_500 of the cluster, in units of your choosing.
    mycosmo       - A dictionary of cosmological parameters
    radii         - Radial bins (in the same units as R500)
    
    Returns
    -------
    An array (the size of r) of pressure (e.g. a profile). The array will have
    units of keV cm**-3    !!!
    """
    
    #P0 = P500 * 8.403 * mycosmo['h_70']**-1.5
    #rP = R500 / 1.177 # C_500 = 1.177
    #rf =  (radii/rP).decompose().value # rf must be dimensionless
    #a=1.0510; b=5.4905; c=0.3081 # gNFW parameters found in Arnaud+ 2010

    #result = (P0 / (((rf)**c)*((1 + (rf)**a))**((b - c)/a)))

    result = gnfw(mycosmo, radii, P500, R500)
    
    return result  ### HAS UNITS ASSOCIATED WITH IT!

def gnfw(mycosmo, radii, P500, R500, c500= 1.177, p=8.403, a=1.0510, b=5.4905, c=0.3081):
    
    P0 = P500 * p * mycosmo['h_70']**-1.5
    rP = R500 / c500 # C_500 = 1.177
    rf =  (radii/rP).decompose().value # rf must be dimensionless
    result = (P0 / (((rf)**c)*((1 + (rf)**a))**((b - c)/a)))

    return result

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
