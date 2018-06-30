import tSZ_spectrum as tsz
import kSZ_spectrum as ksz
import scipy.constants as spconst
import astropy.constants as const
import astropy.units as u
import cosmolopy.distance as cd
import numpy as np
import retrieve_data_info as rdi
import my_astro_defaults as mad

planck = spconst.value("Planck constant in eV s")/1000.0 # keV s
boltzmann = spconst.value("Boltzmann constant in eV/K")/1000.0 # keV/K  
thom_cross = (spconst.value("Thomson cross section") *u.m**2).to("cm**2")
m_e_c2 = (const.m_e *const.c**2).to("keV")
kpctocm = 3.0856776 *10**21
c = const.c
keVtoJ = (u.keV).to("J") # I think I need this...) 
tcmb = 2.72548*u.K # Kelvin (uncertainty = 0.00057)
Icmb = 2.0 * (boltzmann*tcmb.value)**3 / (planck*c.value)**2
Icmb *= keVtoJ*u.W *u.m**-2*u.Hz**-1*u.sr**-1 # I_{CMB} in W m^-2 Hz^-1 sr^-1
JyConv = (u.Jy).to("W * m**-2 Hz**-1")
Jycmb = Icmb.to("Jy sr**-1")  # I_{CMB} in Jy sr^-1
MJycmb= Jycmb.to("MJy sr**-1")

foo = np.log(6.8/1.9)
bar = np.log(260/150.0)
alpha = foo/bar
s353 = 6.8* (353/260.0)**alpha

################################################################

def old_single_freq(myfreq,temp=10.0):

    temp = 10.0 * u.keV
    myfreq = 353.0*u.K  # Really in GHz
    
    freq_conv = (planck*1.0e9)/(boltzmann*tcmb)
    temp_conv = 1.0/m_e_c2
#    bv = 100.0 # arcseconds^2
    fwhm = 10.0*60.0*u.arcsec
    sigma = fwhm / (2*np.sqrt(2*np.log(2)))
    bv = 2.0*np.pi*sigma**2
    bv = bv.to("steradian")
############################################
    llt  = temp*temp_conv  # Lower limit on temperature (well, theta)
    ult  = temp*temp_conv  # Upper limit on temperature
    st   = temp*temp_conv  # Temperature (theta) step
    flow = myfreq;    fhigh= myfreq*1.1; fstep = myfreq
    sx   = fstep*freq_conv
    llx  = flow*freq_conv
    nste = (fhigh-flow)/fstep
    ulx  = fhigh*freq_conv #+ sx/(2.0*nste)
    beta = 0.0
    betaz= 0.0

    tarr, xarr, T = tsz.tSZ_conv(llt,llx,ult,ulx,st,sx)
 #   tarr, xarr, K = ksz.kSZ_conv(beta,betaz,llt,llx,ult,ulx,st,sx,rel=rel)
    srtosas = (1.0*u.sr).to("arcsec**2")     # Square arcseconds in a steradian
    bpsr = srtosas/bv                        # Beams per steradiann
    JypB=Jycmb/bpsr                          # Intensity of the CMB in Jy/beam

    bT = T/tarr    # Divide by tarr (thetae) to get proper conversion units
    tSZ_JyBeam_per_y = JypB * bT  # Just multiply by Compton y to get Delta I (tSZ)

    return tSZ_JyBeam_per_y

def single_freq(myfreq,temp=10.0):

    temp = 10.0 * u.keV
    myfreq = 90.0 *u.K   #It's in GHz, but you put Kelvin so that things cancel...
    
    freq_conv = (planck*1.0e9)/(boltzmann*tcmb)
    temp_conv = 1.0/m_e_c2
    instrument='MUSTANG'
    bv = rdi.get_beamvolume(instrument)
    JypB = tsz.Jyperbeam_factors(bv)

############################################
    llt  = temp*temp_conv  # Lower limit on temperature (well, theta)
    ult  = temp*temp_conv  # Upper limit on temperature
    st   = temp*temp_conv  # Temperature (theta) step
    flow = myfreq;    fhigh= myfreq*1.1; fstep = myfreq
    sx   = fstep*freq_conv
    llx  = flow*freq_conv
    nste = (fhigh-flow)/fstep
    ulx  = fhigh*freq_conv #+ sx/(2.0*nste)
    beta = 0.0
    betaz= 0.0

    tarr, xarr, T = tsz.tSZ_conv(llt,llx,ult,ulx,st,sx)
 #   tarr, xarr, K = ksz.kSZ_conv(beta,betaz,llt,llx,ult,ulx,st,sx,rel=rel)
    srtosas = (1.0*u.sr).to("arcsec**2")     # Square arcseconds in a steradian
    bpsr = srtosas/bv                        # Beams per steradiann
    JypB=Jycmb/bpsr                          # Intensity of the CMB in Jy/beam

    bT = T/tarr    # Divide by tarr (thetae) to get proper conversion units
    tSZ_JyBeam_per_y = JypB * bT  # Just multiply by Compton y to get Delta I (tSZ)
    Kpy = tsz.TBright_factors(myfreq*freq_conv)
    tSZ_Kelvin_per_y = Kpy * bT  # Just multiply by Compton y to get Delta I (tSZ)

    import pdb;pdb.set_trace()
    
    print 'To go from Compton y to Jy/beam, multiply by: ', tSZ_JyBeam_per_y
    print 'To go from Compton y to Kelvin, multiply by: ', tSZ_Kelvin_per_y
    factor = Kpy/JypB
    print 'To go from Jy/Beam to Kelvin (CMB), multiply by: ',factor

def Jy2K(instrument):

    bv = rdi.get_beamvolume(instrument)
    JypB = tsz.Jyperbeam_factors(bv)
    fwhm1,norm1,fwhm2,norm2,fwhm,smfw,freq,FoV = rdi.inst_params(instrument)
    szcv,szcu = mad.get_sz_values()
    x = szcv["planck"]*(freq.to("Hz")).value / (szcv["boltzmann"]*szcv["tcmb"])
    Kpy = tsz.TBright_factors(x)

    print 'To go from Compton y to Jy/beam, multiply by: ', JypB
    print 'To go from Compton y to Kelvin, multiply by: ', Kpy
    
    factor = Kpy/JypB
    print 'To go from Jy/Beam to Kelvin (CMB), multiply by: ',factor
#    import pdb; pdb.set_trace()

def Carlstrom(instrument):

    bv = rdi.get_beamvolume(instrument)
    m2bvrdcs0910 = 122.047765545 * (u.arcsec)**2
    bv = m2bvrdcs0910
    #JypB = tsz.Jyperbeam_factors(bv)
    fwhm1,norm1,fwhm2,norm2,fwhm,smfw,freq,FoV = rdi.inst_params(instrument)
    freq = 87.0*u.GHz
    szcv,szcu = mad.get_sz_values()
    x = szcv["planck"]*(freq.to("Hz")).value / (szcv["boltzmann"]*szcv["tcmb"])

    Te = 6.4 # keV
    Icmb = szcu['Jycmb']
    Inot = (Icmb*(bv.to('sr')))
    term1 = (x * (np.exp(x) + 1.0)/(np.exp(x) - 1.0) -4.0)
    factors = (x**4 * np.exp(x))/(np.exp(x)-1.0)**2
    thetae  = Te/szcv["m_e_c2"]

    carlI  = term1*factors*Inot
    carlT  = term1*szcv["tcmb"]
    myboltz= spconst.value("Boltzmann constant") # J/K
    TfromI = (carlI/(bv.to('sr'))) * (const.c**2 /(2.0* myboltz*u.J/u.K * freq**2))
    TfromI = TfromI.decompose()

    Jy2K  = carlI / carlT

    xtilde = x * (np.exp(x) + 1.0)/(np.exp(x) - 1.0)
    stilde = x *(2.0*np.exp(x/2.0)) / (np.exp(x) - 1.0)
    Y0    = xtilde - 4.0
    Y1    = (-10.0 + 47.0*xtilde/2.0 - 42.0*xtilde**2/5.0 + 0.7*xtilde**3 +
             stilde**2 *(1.4*xtilde - 4.2))
    term2 = (1.0 )

    firstor = Inot*factors*(Y0 + thetae*Y1)
    print carlI
    print firstor

    
