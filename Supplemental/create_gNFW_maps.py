import numpy as np
import scipy.constants as spconst
import astropy.constants as const
import astropy.units as u
import numerical_integration as ni
import cluster_pressure_profiles as cpp  # Not C++ !!!
import my_astro_defaults as mad
import cosmolopy.distance as cd

#import create_gNFW_maps as cgm
cosmo,mycosmo = mad.get_cosmology("Concordance")
M500 = 7.8 * 10**14* const.M_sun; z = 0.89   # For CLJ 1226
#z=1.99, M500=0.63 * 10**14* const.M_sun      # For Stefano's cluster
d_a = cd.angular_diameter_distance(z, **cosmo) *u.Mpc

####################################################################################

def calculate_for_XLSSC():

    z=1.99; M500=0.63 * 10**14* const.M_sun      # For Stefano's cluster
    R500, P500 = R500_P500_from_M500_z(M500, z, mycosmo)
    yprof, inrad  = compton_y_profile_from_m500_z(M500, z, mycosmo)
    myprof, radii = pressure_profile_from_m500_z(M500, z, mycosmo)

    print myprof[1000]
    print inrad[1000]
    
    Rmax = R500
    Ycyl    = Y_cyl(yprof, inrad, Rmax)
    Ysphere = Y_sphere(myprof, radii, Rmax)
        
####################################################################################

def R500_P500_from_M500_z(M500, z, mycosmo):

    H = mycosmo['H0'] * (mycosmo['omega_m']*(1 + z)**3 +
                         mycosmo['omega_l'])**0.5
    E = H/mycosmo['H0']
    dens_crit = (3 * (H)**2)/(8 * np.pi * const.G)
    
    ### P500 will, by default, come back in units of keV/ cm3
    P500 =(1.65 * 10**-3) *((E)**(8./3)) *((
        M500 * mycosmo['h_70'] )/((3*10**14)  * const.M_sun)
        )**(2./3.)*(mycosmo['h_70'])**2  *u.keV /u.cm**3

    ### R500 will, by default, come back in units of meters
    R500 =(3 * M500/(4 * np.pi  * 500 * dens_crit))**(1/3.)
  
    return R500, P500
    
def pressure_profile_from_m500_z(M500, z, mycosmo, N_R500 = 5.0, inner_R = 10**(-4),
                                 Nsteps   = 5000):

    R500, P500 = R500_P500_from_M500_z(M500, z, mycosmo)
    
    R_scaled = np.logspace(np.log10(inner_R), np.log10(N_R500), Nsteps)
    radii    = R_scaled * R500
    
    myprof = cpp.a10_gnfw(P500,R500,mycosmo,radii)

    return myprof, radii

def compton_y_profile_from_m500_z(M500, z, mycosmo):
    """
    Computes the Compton y profile for an A10 (gNFW) profile, integrated from
    (z = +5*R500) to (z = -5*R500).

    Returns the Compton y profile alongside a radial profile in kpc.
    """

    myprof, radii = pressure_profile_from_m500_z(M500, z, mycosmo, N_R500 = 10.0)

    R500, P500 = R500_P500_from_M500_z(M500, z, mycosmo)
    R_scaled = np.logspace(np.log10(10.0**(-4)), np.log10(5.0), 9000)
    radProjected = R_scaled*R500

    m_e_c2 = (const.m_e *const.c**2).to("keV")
    thom_cross = const.sigma_T

    unitless_profile = (myprof * thom_cross * u.kpc / m_e_c2).decompose()

    inrad = radii.to("kpc"); zvals = radProjected.to("kpc")
    
    yprof = ni.int_profile(inrad.value, unitless_profile.value,zvals.value)

    return yprof, inrad

def Y_cyl(yprof, radii, Rmax, d_ang = d_a):

    print d_a
    try:
        angle = radii.to("rad")   # Try converting to radians.
        max_angle = Rmax.to("rad")
    except:
        print 'Radii are not given in angular units.'
        print 'Using angular diameter to convert radii to angular units.'
        angle = (radii/d_ang).decompose() * u.rad
        max_angle = (Rmax/d_ang).decompose() * u.rad
    
    goodR  = (angle < max_angle)
    goodprof = yprof[goodR]; goodangle = angle[goodR]
    
    prats = goodprof[:-1] / goodprof[1:]
    arats = goodangle[:-1] / goodangle[1:]
    alpha = np.log(prats) / np.log(arats)

    parint= ((goodangle[1:]/u.rad)**(2.0-alpha) - (goodangle[:-1]/u.rad)**(2.0-alpha) ) * \
            (goodprof[:-1]*(goodangle[:-1]/u.rad)**alpha) / (2.0 - alpha)
    tint  = 2.0*np.pi * np.sum(parint) * u.sr

    Ycyl  = tint.to("arcmin2")

    print 'Ycyl found to be ', Ycyl

    return Ycyl

def Y_sphere(myprof, radii, Rmax, d_ang = d_a):
    
    m_e_c2 = (const.m_e *const.c**2).to("keV")
    thom_cross = const.sigma_T
    unitless_profile = (myprof * thom_cross * u.kpc / m_e_c2).decompose()

    goodR  = (radii < Rmax)
    goodprof = unitless_profile[goodR]; goodradii = radii[goodR]

    prats = goodprof[:-1] / goodprof[1:]
    arats = goodradii[:-1] / goodradii[1:]
    alpha = np.log(prats) / np.log(arats)
     
    parint= ((goodradii[1:]/u.kpc)**(3.0-alpha) - (goodradii[:-1]/u.kpc)**(3.0-alpha) ) * \
            (goodprof[:-1]*(goodradii[:-1]/u.kpc)**alpha) / (3.0 - alpha)
    tint  = 4.0*np.pi * np.sum(parint) * (u.kpc)**2

    Ysphere = tint.to("Mpc2")

    Ysphere_rad = Ysphere / d_ang
    
    print 'Ysphere found to be ', Ysphere
    print 'or alternatively ', Ysphere_rad

    return Ysphere



    

