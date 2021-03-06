import numpy as np
import scipy.constants as spconst
import astropy.constants as const
import astropy.units as u
import numerical_integration as ni
import cluster_pressure_profiles as cpp  # Not C++ !!!
import my_astro_defaults as mad
import cosmolopy.distance as cd
import mapping_modules as mm
import ellipsoidal_shells as es
import os
import retrieve_data_info as rdi
import instrument_processing as ip
import Azimuthal_Brightness_Profiles as ABP
import matplotlib.pyplot as plt
import retrieve_data_info as rdi
import image_filtering as imf

#import create_gNFW_maps as cgm
cosmo,mycosmo = mad.get_cosmology("Concordance")
#M500 = 7.8 * 10**14* const.M_sun; z = 0.89   # For CLJ 1226
#z=1.99, M500=0.63 * 10**14* const.M_sun      # For Stefano's cluster
defpath='/home/romero/Results_Python/MUSTANG2'

def get_d_a(z):

    d_a = cd.angular_diameter_distance(z, **cosmo) *u.Mpc

    return d_a

####################################################################################

def make_a10_grid():

    for i in range(5):
        mym500 = 2.0**i 
        for j in range(10)+1:
            z             = j*0.2
            M500          = mym500*const.M_sun*1.0e14
            R500, P500    = R500_P500_from_M500_z(M500, z, mycosmo)
            yprof, inrad  = compton_y_profile_from_m500_z(M500, z, mycosmo)
            ang_dist      = get_d_a(z)
            theta_range   = (inrad/ang_dist).decompose()
            theta_range   = theta_range.value

            zstr          = '%1.1f' % z
            mstr          = '%2.1f' % mym500
            
            create_gNFW_map(yprof, theta_range, xsize=512, ysize=512, mappixsize=2.0,savedir=None,
                            filename='A10_map_z'+zstr+'_m'+mstr+'e14.fits',
                            instrument='MUSTANG2',T_e = 5.0,units='Kelvin',
                            beta=0.0, betaz=0.0)

            
            
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

    R500, P500   = R500_P500_from_M500_z(M500, z, mycosmo)
    R_scaled     = np.logspace(np.log10(10.0**(-4)), np.log10(5.0), 9000)
    radProjected = R_scaled*R500

    m_e_c2 = (const.m_e *const.c**2).to("keV")
    thom_cross = const.sigma_T

    unitless_profile = (myprof * thom_cross * u.kpc / m_e_c2).decompose()

    inrad = radii.to("kpc"); zvals = radProjected.to("kpc")
    
    yprof = ni.int_profile(inrad.value, unitless_profile.value,zvals.value)

    return yprof, inrad

def Y_cyl(yprof, radii, Rmax, d_ang = 0, z=1.99):

    if d_ang == 0:
        d_a = cd.angular_diameter_distance(z, **cosmo) *u.Mpc
    
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

def Y_sphere(myprof, radii, Rmax, d_ang = 0, z=1.99):
    
    if d_ang == 0:
        d_ang = cd.angular_diameter_distance(z, **cosmo) *u.Mpc

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



def create_gNFW_map(yprof, theta_range, xsize=512, ysize=512, mappixsize=2.0,savedir=None,
                    filename='gNFW_map.fits', instrument='MUSTANG2',T_e = 5.0,units='Jy',
                    beta=0.0, betaz=0.0,prefile='ProfilePlot_',fig=None,obstarget='cluster'):
    """
    yprof         - A list/array of Compton y values
    theta_range   - The associated array of radii (on the sky) in radians.
    xsize         - map xsize in arcseconds
    ysize         - map ysize in arcseconds
    mappixsize    - pixel size, in arcseconds (on a side)
    savedir       - A string saying where the save directory is
    filename      - specifies the .fits file name
    instrument    - the instrument filtering to be used
    T_e           - assumed electron temperature
    units         - The sky brightness units - either 'Jy' (for Jy/beam) or 'Kelvin'.
    beta          - total speed of the cluster (as a fraction of c)
    betaz         - spped of the cluster along the line of sight (+ is away from us)
    fig           - Pass in a figure object to make overplots.

    """
    
    fwhm1,norm1,fwhm2,norm2,fwhm,smfw,freq,FoV = rdi.inst_params(instrument)
    bv = rdi.get_beamvolume(instrument)
    
    xymap = mm.create_xymap(xsize=xsize,ysize=ysize,mappixsize=mappixsize,
                            xcentre=[],ycentre=[])
    
    mymap = es.grid_profile(theta_range, yprof, xymap)
    w     = mm.create_astro(mymap,mappixsize=mappixsize)
    print fwhm

    xfer_class = mm.xfer_fxn(instrument)
    my_xfer    = rdi.get_xfer(xfer_class)
    beam_conv  = ip.conv_gauss_beam(mymap,mappixsize,fwhm)
    filt_img   = ip.apply_xfer(beam_conv, my_xfer,instrument)

    tSZ,kSZ = rdi.get_sz_bp_conversions(T_e,instrument,bv,units,array="2",
                                        inter=False,beta=beta,betaz=betaz,rel=True)
    Sky_Bright = filt_img * tSZ.value   # Map in measurable units (either Kelvin or Jy/beam)
    
    hdu   = mm.make_hdu(mymap, w, filt_img, Sky_Bright)

    myhome = os.path.expanduser("~")
    if type(savedir) == type(None):
        savedir = os.path.join(myhome,'Results_Python/gNFW_profiles')
    if not os.path.exists(savedir):
        os.makedirs(savedir)
    fullpath = os.path.join(savedir,filename)
    print fullpath
    mm.write_hdu_to_fits(hdu,fullpath)

    rmstarget=56.0
    rmsmap = make_rmsMap(xymap,theta_range,rmstarget,conv=tSZ)
    compyrms = np.abs(rmstarget/tSZ)
    rmsstr = "{:4.2f}".format(compyrms)
    
    if type(fig) == type(None):
        maxind=3
    else:
        maxind=2
        
    for thismap,thislab,index in zip([mymap, filt_img, rmsmap],['Sky Map','Filtered Map','center RMS: '+rmsstr+'e-6 Compton y'],
                                     [1,2,3]):
        if index > maxind: continue
        
        angmin = 0.0;  angmax = 2.0*np.pi
        mythetamask = np.zeros(thismap.shape)

        rads, prof, mythetamask = ABP.get_az_profile(thismap, xymap, angmin, angmax,
                                                 thetamask=mythetamask)
        if maxind == 2:
            prof*=tSZ.value
            mytarget = obstarget
        else:
            mytarget = obstarget+'_'+prefile
        mybins = np.arange(0.0,180.0,5.0)
        binres = ABP.radial_bin(rads, prof,10,rmax=180.0,bins=mybins,minangle=angmin,maxangle=angmax,
                                profunits='Compton y')
        if index == 2:
            filbinres = binres

            
        fig    = ABP.plot_one_slice(binres,myformat='png',fig = fig,target=mytarget,
                                     savedir=savedir,prefilename=prefile,mylabel=thislab)
        #import pdb;pdb.set_trace()
        
        if index == 3:
            bpixv = ((binres.npix*(mappixsize*u.arcsec)**2)/bv).decompose()
            print bpixv
            binres.profavg /= np.sqrt(bpixv.value)

            sig = filbinres.profavg / binres.profavg
            
            
            thislab = 'Standard error of the mean'
            #import pdb;pdb.set_trace()
            fig    = ABP.plot_one_slice(binres,myformat='png',fig = fig,target=prefile,
                                         savedir=savedir,prefilename=prefile,mylabel=thislab)
        
def make_rmsMap(xymap,theta_range,target=25.0,conv=1.0):
    """
    xymap      - you know!
    target     - target sensitivity in uJy/beam
    conv       - conversion factor (if you want it in Compton y)
    """

    arcrad   = theta_range * (u.rad).to('arcmin')
    rms_stg1 =  1.0 + np.sin(arcrad*np.pi/8.0)*2.5         # Outer part should equal 3.5
    rms_stg1 =  np.exp(1.0 - np.cos(arcrad*np.pi/8.0))     #
    rms_prof =  (rms_stg1**2)*target*1.0e-6 /np.abs(conv)    # Inner part ~ 56 uJy

    #import pdb;pdb.set_trace()
    mymap = es.grid_profile(theta_range, rms_prof, xymap)

    return mymap
    
def plot_y_from_z_m500(fig,z,mym500,myunits='Jy',instrument='MUSTANG2',temp=5,bv=120.0,
                       savedir=defpath,myformat='png',target='cluster',myfontsize=15,
                       pixsize=1.0):

    M500          = mym500 * 10**14* const.M_sun     
    R500, P500    = R500_P500_from_M500_z(M500, z, mycosmo)
    yprof, inrad  = compton_y_profile_from_m500_z(M500, z, mycosmo)
    ang_dist      = get_d_a(z)
    theta_range   = (inrad/ang_dist).decompose()
    theta_range   = theta_range.value * (u.rad).to('arcsec')

    
    ax = fig.add_subplot(111)
    ax.set_ylim(ax.get_ylim())
    ax.set_xlim(ax.get_xlim())

    tSZ,kSZ = rdi.get_sz_bp_conversions(temp,instrument,bv,units=myunits)
    mapprof = yprof*tSZ
    #import pdb;pdb.set_trace()

    linrads = (np.arange(int(np.max(theta_range)/(2.0*pixsize)))+1.0)*pixsize
    linprof = np.interp(linrads,theta_range,mapprof)

    xfer_class = mm.xfer_fxn(instrument)
    my_xfer    = rdi.get_xfer(xfer_class)
    myk        = my_xfer[0,0:]*pixsize

    filtprof = imf.fourier_filtering_1d(linprof,'tab',(myk,my_xfer[1,0:]))
    
    #import pdb;pdb.set_trace()
    ax.plot(theta_range,mapprof,label='Unfiltered A10 profile',color='g')
    ax.plot(linrads,filtprof,label='Filtered A10 profile',color='r')
    plt.legend(fontsize=myfontsize)    
    prefilename='Radial_profiles_'
    filename=target+'_v2.'
    fullpath = os.path.join(savedir,prefilename+filename+myformat)
    #if doleg == True: plt.legend()
    plt.savefig(fullpath,format=myformat)
