import my_astro_defaults as mad
import my_cluster_defaults as mcd
#import cosmolopy.distance as cd
#from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
from astropy.coordinates import Angle
import numpy as np
import astropy.constants as const
import retrieve_data_info as rdi   # Finds and loads data files (hard-coded)
import mapping_modules as mm       # Creates radius + pressure arrays.
import mfit_params as mfp          # Loads parameters for fitting procedure
import importlib,time,os
from astropy.wcs import WCS
from os.path import expanduser
myhome = expanduser("~")

today=time.strftime("%d%b%Y")
todaysp=time.strftime("%d %b %Y")
dandt = time.strftime("%d-%b-%Y %H:%M:%S")

###########################################################################
### VERSION HISTORY                                                     ###
###                                                                     ###
### Written 27 June 2017 by Charles Romero                              ###
###                                                                     ###
### SYNOPSIS:                                                           ###
### ------------------------------------------------------------------  ###
### While I had worked to divide up tasks into smaller routines, I      ###
### found many variables needing to be passed around. Sometimes I may   ###
### copy variables if I think it's appropriate to have in mutliple      ###
### structures (class) or dictionaries. Nevertheless, I realized I      ###
### needed to start nesting nearly all variables inside some kind.      ###
### Thus, this routine was born to collect all variables in a manner    ###
### that is hopefully logical.                                          ###
###                                                                     ###
###                                                                     ###
###########################################################################

def get_struct_of_variables(instruments,name,path=myhome+'/Results_Python/',
                            testmode=False, map_file='all',reduc='PCA',mymaxrad=0):

    """
    This module is the overarching module to collect variables.
    ------------------------------------
    INPUTS:
    ----
    instruments    - The instruments which took the data (e.g. "MUSTANG-2")
    name           - The name of the target (object) observed
    path           - The output directory
    testmode       - True or False; True uses few steps, so you test that the code works
                     and that the results seem reasonable (although with too few steps,
                     bugs may still persist).
    map_file       - Given that you may want to load noise, or a simulated model, you
                     can choose what "data" file you want to load.
    """
    
    ### Get data files and cluster/instrument specific variables:
    input_struct = importlib.import_module(name+'_info')
    shocks = input_struct.shocks()  # Needs to be defined...
    bulk   = input_struct.bulk()    # Needs to be defined...
    ptsrcs = input_struct.ptsrc()   # Locations of point sources; not instrument-specific.
    blobs  = input_struct.blob()   # Locations of point sources; not instrument-specific.
    priors = [input_struct.priors(cno+1) for cno in range(len(bulk.bins))]  # These are not fitting priors, but prior "known" (fixed) quantities.
    dv     = {}              # Data-related variables (e.g. maps)
    ifp    = {}              # Individual fitting parameters (components to be treated individually).
    # Minimum and maximum of angular scales "probed" by your instruments
    minmax     = np.array([3.0,60.0])*u.arcsec  # The minimum and maximum range you want your bins to cover
    mycluster  = [cluster_defaults(mypriors,number=i+1) for i,mypriors in enumerate(priors)]   # A class with cluster-specific attributes
    nifp       = 0                         # The number of individually-fitted parameters
    ninstptsrc = 0                         # The number of instruments for which to fit pt srcs.

    pie_frac = 1.0
    for cgeo in bulk.geoparams:
        if cgeo[7] > 0: pie_frac *= (2.0*np.pi)/cgeo[7]
    my_factor = pie_frac**(1.0/len(bulk.geoparams))
        
    for instrument in instruments:
        inputs = input_struct.files(instrument=instrument,map_file='all',reduction=reduc)
        ### Now, we can get the input data variables
        dv[instrument]  = data_vars(inputs,priors,mycluster,ptsrcs,bulk)
        ifp[instrument] = mfp.inst_fit_params(inputs,ptsrcs,blobs,instrument)
        ############## Let's think about what bins we can make #########################
        goodpix = dv[instrument].maps.masked_wts > 0             ## These are the pixels this instrument
        npix = len(dv[instrument].maps.masked_wts[goodpix])      ## is using.
        mask_rad = (np.sqrt(npix/np.pi))*dv[instrument].mapping.pixsize ## Which defines a radial limit.
        mask_rad *= my_factor                        ## accounts for fitting slices...imperfectly.
        if mymaxrad > 0:
            max_extent = mymaxrad*u.arcsec
        else:
            if mask_rad > dv[instrument].FoV/1.4:     ## And the FOV defines a transfer function radial limit.
                max_extent = dv[instrument].FoV/1.4   ## So, which one is limiting the range that we can fit?
            else:                                     ## 
                max_extent = mask_rad*0.9             ## With this figured out, we can look at other instruments too.

        print("Maximum radial extent: ",max_extent)
        #import pdb;pdb.set_trace()
        ### These minmax are for the bins. Here, I am taking the min and max across instruments
        minmax = angular_range(minmax,dv[instrument].fwhm/1.8,max_extent) #dv[instrument].fwhm
        
        #print dv[instrument].mapping.pixsize,dv[instrument].FoV
        ### Mean levels figure into instrument-specific parameters.
        nifp += ifp[instrument].n_add_params # Number of instrument-specific parameters 
        if ifp[instrument].pt_src: ninstptsrc += 1

    #xycen = (dv['MUSTANG2'].mapping.x0,dv['MUSTANG2'].mapping.y0)
    radec  = [dv['MUSTANG2'].mapping.ra,dv['MUSTANG2'].mapping.dec]

    ### Need to convert bins defined in kpc to arcseconds...
    ### Now, as the fitting parameters depends on the input file, let's get those:
    cfp = mfp.common_fit_params(bulk,shbins=shocks.bins,path=path,
                                shockgeo=shocks.geoparams,shockfin=shocks.shockfin,
                                ptsrcs=ptsrcs,testmode=testmode,
                                fstemps=shocks.fstemps,blobs=blobs.blobpars,
                                blobcens=[blobs.ra,blobs.dec],
                                minmax=minmax,cluster=mycluster,radec=radec)
    # A correction on the total number of dimensions:
    ### 20 Feb 2018 - This correction is true, but I had to make the common_fit_params
    ### routine in-line with this!
    cfp.ndim += nifp + ninstptsrc*len(cfp.ptsrc)

    print('Using the following bins for the bulk profile: ')
    for bulkbins in cfp.bulkarc:
        print(bulkbins*u.rad.to("arcsec"), ' arcseconds')
  
    hk = housekeeping(cfp,priors,instruments,name,mycluster)
    ### OK, now...can we collect other "defaults"
    ### Let's have a quick look at the data and see if we should be fitting
    ### a point source. If so, we'll update the fit parameters.
    ### ufp = mfp.update_fit_params(dv,hk)
    ### hk.fit_params = ufp # Maybe I want to update directly...
    
    return hk,dv,ifp
    
class housekeeping:

    def __init__(self,cfp,priors,instruments,name,mycluster):

        ### Here are the variables/classes that I put under the umbrella of "housekeeping":
        self.log         = logbook()                     # Ideally take notes, keep track of what was done.
        self.hk_ins      = priors                        # What is taken as known about this cluster?
        self.instruments = instruments                   # The instruments that we use for this fitting
        self.hk_outs     = hk_out(cfp,instruments,name)  # Some output variables (directories, filenames)
        self.cfp         = cfp                           # Common Fit Parameters
        self.cluster     = mycluster                     # A class with attributes of the cluster
        self.av          = mad.all_astro()               # A class with (general) astronomical values.
        
class data_vars:

    def __init__(self,inputs,priors,mycluster,ptsrcs,bulk,real=True,beta=0,betaz=0):
        """
        I need to update this so that it reads in all the data (files) and then I can discard the
        input file information.
        """

        self.maps  = rdi.maps(inputs,mycluster,ptsrcs,bulk)
        self.astro = rdi.astrometry(self.maps.header)
        self.xfer  = rdi.get_xfer(inputs)
        self.xferdims = inputs.tabdims
        
        try:
            Jy2K = self.maps.header['JY2K']
        except:
            try:
                Jy2K = 1.0 / self.maps.header['ANTGAIN']
            except:
                Jy2K = -1.0

        #print(Jy2K)
        #import pdb;pdb.set_trace()
                
        if Jy2K > 0:
            bv = rdi.get_bv_from_Jy2K(Jy2K,inputs.instrument)
        else:
            bv = rdi.get_beamvolume(inputs.instrument)

        tSZ = []; kSZ = []
        for myprior in priors:
            JytSZ,JykSZ = rdi.get_sz_bp_conversions(myprior.Tx,inputs.instrument,bv,'Jy/beam',array="2",
                                                inter=False,beta=beta,betaz=betaz,rel=True)
            K_tSZ,K_kSZ = rdi.get_sz_bp_conversions(myprior.Tx,inputs.instrument,bv,'Kelvin',array="2",
                                                inter=False,beta=beta,betaz=betaz,rel=True)
            mytSZ,mykSZ = rdi.get_sz_bp_conversions(myprior.Tx,inputs.instrument,bv,inputs.units,array="2",
                                                inter=False,beta=beta,betaz=betaz,rel=True)
            tSZ.extend([mytSZ])
            kSZ.extend([mykSZ])
            MyJy2K = JytSZ/K_tSZ  # Divide by this factor to go Jy -> K

        if Jy2K < 0: Jy2K = MyJy2K # Divide by this factor to go Jy -> K
            
        self.Jy2K = Jy2K # Divide by this factor to go Jy -> K
        self.bv   = bv
        
        av = mad.all_astro()     # av = astro vars
        self.mapping = mapping_info(mycluster,inputs,av,tSZ,kSZ,ptsrcs)
        
        self.tSZ = tSZ
        self.kSZ = kSZ
        fwhm1,norm1,fwhm2,norm2,fwhm,smfw,freq,FoV = rdi.inst_params(inputs.instrument)
        self.fwhm = fwhm if Jy2K < 0 else np.sqrt(bv/1.1330900354567985)
        self.freq = freq
        self.FoV  = FoV
        #self.rms_corr   = inputs.rmscorr
        self.conversion = rdi.get_conv_factor(inputs.instrument,bv=bv)
        
def print_attributes(myclass):

    attrs = vars(myclass)
    print(', '.join("%s" % item for item in attrs))

def angular_range(minmax,fwhm,fov):

    if fwhm < minmax[0]:
        minmax[0] = fwhm
    if fov > minmax[1]:
        minmax[1] = fov

    return minmax
    
class logbook:

    def __init__(self):

        self.today = today      # Might be useful for attaching to strings
        self.starttime = dandt  # For a full note of when it was done.
        self.endtime = dandt    # Will need to be updated
        self.notes = "None yet" # Might make this an array of strings
        self.fitresults = "None"# Same
        self.ptsrcresults=" NA" #
        self.shockresults=" NA" #
        self.blobresults =" NA" #
        self.miscresults =" NA" #

class hk_out:

    def __init__(self,fit_params,instruments,target):
    
        mlstr='ML-NO_'; ppstr='PP-NO_'; dtype='Virt_'; ubstr='POWER_'
        if len(instruments) == 1:
            instrument = instruments[0]
        else:
            instrument = "Combined"
            
        pre_filename = "{}_Virt_".format(instrument)
        if fit_params.real_data == True:   
            dtype='Real_'; pre_filename = "{}_Real_".format(instrument)
        if fit_params.pprior == True: ppstr='PP-YES_'
        if fit_params.ubmap  == True: ubstr='UNIFORM_'         
        #if fit_params.mn_lvl == True:
        #    fit_params.ndim=fit_params.ndim+1
        #    mlstr='ML-YES_'

        #print_attributes(fit_params)
        #import pdb;pdb.set_trace()

        nbstr=str(fit_params.ndim)+"_Dim_"
        nststr=str(fit_params.nsteps)+"S_"
        nburnstr=str(fit_params.burn_in)+"B_"
        walkstr=str(fit_params.nwalkers)+"W"
        newfolder = time.strftime("%d%m_")+format(instrument)+'_'  #change as wanted
        self.runstr=nbstr+dtype+nststr+nburnstr+mlstr+ppstr+ubstr+walkstr
        self.nnfolder=newfolder+self.runstr
        inst_path = r'{}/{}'.format(fit_params.path,instrument)
        if not os.path.exists(inst_path):
            os.makedirs(inst_path)     
        clus_path = r'{}/{}'.format(inst_path,target)
        if not os.path.exists(clus_path):
            os.makedirs(clus_path)

        self.newpath = clus_path
        self.prefilename=pre_filename

        
class cluster_defaults:

    def __init__(self,myinput=None,name="Unknown",number=1):
        """
        Returns a structure with important variables related to the cluster parameters 
        (which physically described the cluster), as well as parameters related to our
        viewing of the cluster (e.g. angular diameter) that depend on cosmology.

        Parameters
        __________
        name      - The name of the cluster
        M_500     - The mass enclosed within R_500
        z         - The redshift of the cluster
        ra        - The Right Ascenscion (in hours or degrees; degrees preferred)
        dec       - The Declination (in degrees)

        Returns
        -------
        The structure cluster
        """

    ### Get general cosmological parameters:
        
        cosmo,mycosmo=mad.get_cosmology()
        sz_constants=mad.get_sz_values()

        if not(myinput == None):
            z,ra,dec,M_500,Tx = myinput.z,myinput.ra,myinput.dec,\
                                myinput.M_500,myinput.Tx
            name = myinput.name
        else:
            z, ra, dec, M_500, Tx = clust_info(name) #myinput.name

        H0 = mycosmo['H0']
        h_70 = mycosmo['h_70']
        self.h70       = h_70
        self.H0        = H0
        self.name = name
        self.z         = z              # should be a list
        self.E         = (mycosmo['omega_m']*(1 + self.z)**3 + mycosmo['omega_l'])**0.5
        self.H         = H0 * self.E
        self.hofz      = self.E #* h_70        #* 0.7   Agh! This was incorrect. Fixed Aug. 30, 2018
        self.dens_crit = (3 * (self.H)**2)/(8 * np.pi * const.G)
        self.M_500     = M_500 * const.M_sun
        self.P_500     = (1.65 * 10**-3) *((self.E)**(8./3)) *((
            self.M_500 * h_70 )/((3*10**14)  * const.M_sun)
            )**(2./3.+0.11)*(h_70)**2  *u.keV /u.cm**3
        self.R_500     =(3 * self.M_500/(4 * np.pi  * 500 * self.dens_crit))**(1/3.)
        self.R_500     = self.R_500.to(u.kpc)
        self.R_max     = 5 * self.R_500
        self.d_a       = get_d_ang(self.z,cosmo)
        self.d_a       = self.d_a.to("kpc") #convert to kpc (per radian)
        self.scale     = self.d_a*(u.arcsec).to("rad") / u.arcsec
        self.theta_500 = (self.R_500 *u.rad/ self.d_a).to("arcmin")
        self.theta_max = (self.R_max *u.rad /self.d_a).to("arcmin")
        self.ra_deg    = ra.to("deg")
        self.dec_deg   = dec
        self.Tx        = Tx
        self.number    = number

def get_d_ang(z,cosmo):

    d_ang = cosmo.comoving_distance(z) / (1.0 + z)

    return d_ang

        
class mapping_info:

    def __init__(self,cluster,inputs,mycosmo,tSZ,kSZ,ptsrcs):
        """
        Returns a structure with important variables related to gridding a map.
        In this sense, I have called it "astrometry", even though every variable
        may not classically fit within astrometry.

        Parameters
        __________
        cluster   - A structure, from the class/routing above in this file.
        inputs    - Instrument-specific inputs
        priors    - 
        mycosmo   - A class with various cosmological parameters
        tSZ       - A scalar that converts Compton y to Jy/beam for the given frequency
        kSZ       - "" but for the kinetic SZ
        ptsrcs    - Point Source locations list of tuples: (RA, DEC)

        Returns
        -------
        The structure mapping
        """
        ############################################################################################
        ### First we do need to calculate a few variables
        
        #rintmax=cluster.R_500.value
        #theta_min=(0.2 * u.arcsec).to("radian").value
        #ltm=np.log10(theta_min/2.0)          # Log of Theta Min (radians)
        ### Theta_max is already 5* R_500
        #ltx=np.log10(cluster.theta_max.to("radian").value)  # Log of Theta Max (radians)
        w = WCS(inputs.fitsfile)
        #x0,y0=w.wcs_world2pix(cluster.ra_deg,cluster.dec_deg,0)
        image_data, ras, decs, hdr, pixs = rdi.get_astro(inputs.fitsfile)
        print('Pixel Size is: ',pixs)
        theta_min=(pixs/2.0).to("radian").value
        #import pdb; pdb.set_trace()
        #xymap  = mm.get_xymap(image_data,pixs,xcentre=x0.item(0),ycentre=y0.item(0))       
        #arcmap = mm.get_radial_map(image_data,pixs,xcentre=x0.item(0),ycentre=y0.item(0))  
        xymap,x0,y0  = mm.get_xymap(image_data,pixs,retcen=True)
        ra0,dec0 = w.wcs_pix2world(x0,y0,0)
        arcmap = mm.get_radial_map(image_data,pixs,xcentre=x0,ycentre=y0)  
        x,y    = xymap
        angmap = np.arctan2(y,x)                              # goes from -pi to +pi
        radmap = (arcmap*u.arcsec).to("rad").value            # In radians
        rmval  = radmap; bi=np.where(rmval < theta_min); rmval[bi]=theta_min
        radmap = rmval
        ### OK, the point source locations need to be defined here (Feb. 19, 2018)
        nx,ny  = image_data.shape
        ptxys=[]
        for myptsrc in ptsrcs.locs:
            x1,y1=w.wcs_world2pix(myptsrc[0].to('deg').value,myptsrc[1].to('deg').value,0)
            #x2,y2=w.all_world2pix(myptsrc[0].to('deg'),myptsrc[1].to('deg'),0)

            ### Calculate offsets relative to xymap
            xpt,ypt = (x1-x0)*(pixs.to('arcsec')).value, (y1-y0)*(pixs.to('arcsec')).value 
            #xp2,yp2 = (x2-x0)*(pixs.to('arcsec')).value, (y2-y0)*(pixs.to('arcsec')).value 
            print('Point source offsets are: ',xpt,ypt)
            #print('Point source v2 offsets are: ',xp2,yp2)
            ptxys.append((xpt,ypt))
            
        #import pdb;pdb.set_trace()

        ############################################################################################
        ### And now we can put variables into our structure.
        ######################################################### But first, another quick comment:
        ### ltrmax = Log(Theta_Range)_max - for LTRMAX*R_500
        ### Thus, the default is 5*R_500
        #self.theta_min=theta_min
        #self.theta_max= (15.0* u.arcmin).to("radian").value  # 15' in radians.
        #self.theta_range = np.logspace(ltm,ltx+np.log10(ltrmax), N)
        #self.theta_range = np.logspace(ltm,ltx, N)
        self.w          = w
        self.ra         = ra0.item(0)
        self.dec        = dec0.item(0)
        self.pixsize    = pixs
        self.x0         = x0
        self.y0         = y0
        self.tSZ        = tSZ
        self.kSZ        = kSZ
        self.tab        = rdi.get_xfer(inputs)
        self.tabdims    = inputs.tabdims     # Better to include this variable, if necessary
        #self.r_bins=r_bins           # Need just the value (do not want a "quantity")
        #self.a10_pressure = pressure # Need just the value (do not want a "quantity")        
        self.radmap     = radmap           # Need just the value (do not want a "quantity")
        self.arcmap     = arcmap
        self.angmap     = angmap
        self.xymap      = xymap
        self.nx         = nx
        self.ny         = ny
        self.instrument = inputs.instrument
        ### OK, how do I define this?
        self.ptsrclocs  = ptxys

        
