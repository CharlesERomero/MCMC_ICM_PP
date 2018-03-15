import scipy.constants as spconst
import astropy.constants as const
import astropy.units as u
import cosmolopy.distance as cd

class all_astro:

    def __init__(self,ref="Planck2016"):

        cosmo,mycosmo = get_cosmology(ref=ref)
        szcv,szcu = get_sz_values()
        ### These are all dicitonaries
        self.cosmo=cosmo      ### cosmo is formated for a cosmolopy routine
        self.mycosmo=mycosmo  ### mycosmo is formated for my use.
        self.szcv=szcv        ### SZ "constants" -> values
        self.szcu=szcu        ### SZ "constants" -> with units


def get_cosmology(ref = "Planck2016"):
    #Defines Universe Values and other useful Quantities
    
    if ref == "Planck2016":
        ### Cosmology from Planck+ 2016
        omega_m, omega_A = 0.3089, 0.6911 #Values from Remi's Thesis
        H0 = 67.74 *u.m /(const.kpc *u.s)
        #cosmo = {"omega_M_0":0.3089,"omega_lambda_0":0.6911, "h":0.6774}
    ########################################################
    
    if ref == "Hinshaw2013":
        ### Cosmology from Hinshaw+ 2013
        omega_m, omega_A = 0.2821, 0.7181 
        H0 = 69.7 *u.m /(const.kpc *u.s)
        #cosmo = {"omega_M_0":omega_m,"omega_lambda_0":omega_A, "h":0.697}

        
    ########################################################
    
    if ref == "Concordance":
        ### Cosmology from Hinshaw+ 2013
        omega_m, omega_A = 0.3, 0.7 
        H0 = 70.0 *u.m /(const.kpc *u.s)
        #cosmo = {"omega_M_0":omega_m,"omega_lambda_0":omega_A, "h":0.697}

        
    ########################################################
    H100 = 100.0 *u.m /(const.kpc *u.s)
    cosmo = {"omega_M_0":omega_m,"omega_lambda_0":omega_A, "h":(H0/H100).value}
    h_70 = H0/70. /u.m *(const.kpc *u.s)
    mycosmo = {"omega_m":omega_m,"omega_l":omega_A, "h":(H0/H100).value,
               "H0":H0,"h_70":h_70}
    
    cosmo = cd.set_omega_k_0(cosmo)

    return cosmo, mycosmo

def get_sz_values():
    ########################################################
    ### Astronomical value...
    tcmb = 2.72548*u.K # Kelvin (uncertainty = 0.00057)
    ### Reference:
    ### http://iopscience.iop.org/article/10.1088/0004-637X/707/2/916/meta
    
    ### Standard physical values.
    thom_cross = (spconst.value("Thomson cross section") *u.m**2).to("cm**2")
    m_e_c2 = (const.m_e *const.c**2).to("keV")
    kpctocm = 3.0856776 *10**21
    boltzmann = spconst.value("Boltzmann constant in eV/K")/1000.0 # keV/K  
    planck = spconst.value("Planck constant in eV s")/1000.0 # keV s
    c = const.c
    keVtoJ = (u.keV).to("J") # I think I need this...) 
    Icmb = 2.0 * (boltzmann*tcmb.value)**3 / (planck*c.value)**2
    Icmb *= keVtoJ*u.W *u.m**-2*u.Hz**-1*u.sr**-1 # I_{CMB} in W m^-2 Hz^-1 sr^-1
    JyConv = (u.Jy).to("W * m**-2 Hz**-1")
    Jycmb = Icmb.to("Jy sr**-1")  # I_{CMB} in Jy sr^-1
    MJycmb= Jycmb.to("MJy sr**-1")

    ### The following constants (and conversions) are just the values (in Python):
    sz_cons_values={"thom_cross":thom_cross.value,"m_e_c2":m_e_c2.value,
                    "kpctocm":kpctocm,"boltzmann":boltzmann,
                    "planck":planck,"tcmb":tcmb.value,"c":c.value,}
    ### The following "constants" have units attached (in Python)!
    sz_cons_units={"Icmb":Icmb,"Jycmb":Jycmb,"thom_cross":thom_cross,
                   "m_e_c2":m_e_c2}

    return sz_cons_values, sz_cons_units


#pconst = 6.62607004e-34
#c = 299792458.0
#kb = 1.38064852e-23
#tcmb = 2.725
#cp = c*pconst
#kt = kb*tcmb
#Inot = 2.0 * (kt/cp)**2 * kt
