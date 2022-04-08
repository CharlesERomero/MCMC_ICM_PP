from importlib import reload
import numpy as np                      # A useful package...
import emcee, os, shelve, pickle, time  # A list of modules to import as-is
from astropy.io import fits             # To read/write fits
import ellipsoidal_shells as es         # Integrates ellipsoidal power law distributions
import instrument_processing as ip      # Determines instrument-specific values (e.g. xfer fxn)
import astropy.units as u               # Allows variables (quantities) to have units
import cluster_pressure_profiles as cpp # Creates radius + pressure arrays.
import rw_resultant_fits as rwrf        # Read/Wrte Resultant Fits
rwrf=reload(rwrf)
#from scipy import optimize              #
from scipy import stats                 # 
import scipy.special as scs             # For the Gamma function... could send this to ai...hmm
from os.path import expanduser          # A way to determine what the home directory is
import multiprocessing                  # A way to determine how many cores the computer has
import datetime                         # A more thorough module than time.
import Azimuthal_Brightness_Profiles as ABP # Modules for calculating what you think...
import numerical_integration as ni      # For the gNFW profile
#import emcee_stats
#import plot_mcmc_results as pmr         # Called here for plotting convergence tests.

#rwrf=reload(rwrf)                      # Reload - primarily of use during development of code.
myhome = expanduser("~")                # What is the user's home directory?
ncpus  = multiprocessing.cpu_count()    # One can decide how aggressively to parallel process.
#defaultYM = 'A10'

########################################################################
### March 6, 2018. I think the following bugs are fixed.
########################################################################
### Some *bugs* that had to be fixed:
### (1) We need to be sure that theta_range (dv.mapping.theta_range) encompasses the
###     range of radians subtended by the mapped region. Here, I needed to make sure
###     that the minimum of the range was still smaller than the "theta_min", which I
###     use an overwrite for any radii smaller than that value (so as to avoid errors
###     with a power law going to zero).
### (2) Fixing the transfer function, so as to accommodate very small scales (and thus
###     large wavenumber k)
###
########################################################################
### TO DO:
### (1) Document (comment) this and other codes.
########################################################################

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

class emcee_fitting_vars:

    def __init__(self,hk,dv,ifp,tag='',nthreads=1,ySph=True,YMrel='A10'):
        
        """
        The attributes of this class are attributes that are needed for fitting models
        to the data. In particular, we want to assemble attributes (variables) which
        are our initial guesses for the (all) parameters which are being modulated in
        the MCMC sampling.
        """
                
        ### cluster,mapping,fit_params,inalphas,p_init
        #geoparams = None # I think we don't need this at all....
        
        ### Get initial guess for the pressure of your bulk component(s):

        myval     = []   # 
        myalp     = []   #
        sbepos    = []   # Boolean array indicating whether myval should be positive or not.
        priors    = []   # An array of priors.
        priunc    = []   # If prior unc < 0, then there is no actual prior.
        compname  = []   #
        compcount = []   #
        blobsign  = []
        ptsrcsign = []
        ### Integral(P*dl) to Compton y (for dl given in radians):
        #Pdl2y = (hk.av.szcu['thom_cross']*hk.cluster.d_a/hk.av.szcu['m_e_c2']).to("cm**3 keV**-1")

        Pdl2y     = [] ;   R500 = [] ; P500 = [] ; largestR500 = 0.0 * u.rad/u.rad
        for mycluster in hk.cluster:
            Pdl2y.extend([(hk.av.szcu['thom_cross']*mycluster.d_a/hk.av.szcu['m_e_c2']).to("cm**3 keV**-1")])
            thisR500 = (mycluster.R_500/mycluster.d_a).decompose()
            R500.extend([thisR500])
            P500.extend([mycluster.P_500])
            if thisR500 > largestR500: largestR500 = thisR500

        nshockp = 0; nmnlvl = 0; nbulkp = 0; nblob = 0; ncent=0; nptsrc=0; ngeo=0
        minpixs = 10.0 #(arseconds); minimum pixel size among instruments (TBD)

        for myinst in hk.instruments:
            instpixs = (dv[myinst].mapping.pixsize).to("arcsec").value
            if instpixs < minpixs: minpixs = instpixs
            
            if ifp[myinst].mn_lvl:
                mnlvlguess = -5.0*np.mean(dv[myinst].maps.data)
                #import pdb;pdb.set_trace()
                myval.append(mnlvlguess) # Append better for scalars
                sbepos.extend([False])
                priors.extend([mnlvlguess])
                priunc.extend([-1.0])
                compname.extend(['mnlvl'])
                compcount.extend([nmnlvl])
                nmnlvl+=1

        minpixrad = (minpixs*u.arcsec).to("rad").value   # Profile calculated from array in radians.

        #############################################################################
        ### This is the min and max for creating the 1-d profile
        ### NOTE #1: The minimum given here is *the* minimum that will be gridded.
        ### That is, in a "radmap", any radius value less than the minimum here,
        ### takes on the minimum radial value.
        ### NOTE #2: While I should want to define the minimum value directly based on
        ### the minimum pixel size, I am assuming that most maps are gridded with a
        ### pixel size more than 3 times smaller than the FWHM, which is used to determine
        ### the minimum profile bin size. Thus, I scale by a factor of 10 from the minimum
        ### bin size. Note that this might not be a good choice!
        #############################################################################
        tnx = [minpixrad/2.0,10.0*largestR500.value] # In radians
        #tnx = [hk.cfp.bminmax[0].to("rad").value/20.0,10.0*R500.value] # In radians
        nb_theta_range = 150                   # The number of bins for theta_range
        theta_range = np.logspace(np.log10(tnx[0]),np.log10(tnx[1]), nb_theta_range)
        

        ### Array of fitting values; the order of components are:
        ### (1) mn lvl, (2) bulk, (3) shocks, (4) pt srcs, (5) "blobs"
        ### Within each component, each instrument is addressed. [2] and [3]
        ### are modelled jointly. Maybe I want to do the same with [5], depending
        ### on the supposed physical origin.

        bulkcount=0
        ylist=[]
        y2500=[]
        mlist=[]
        edens=[]
        for bulkbins,fit_cen,fit_geo,mybulkgeo,myPdl2y,mycluster in zip(hk.cfp.bulkarc,hk.cfp.bulk_centroid,
                                                                            hk.cfp.bulk_geometry,hk.cfp.bulkgeo,
                                                                            Pdl2y,hk.cluster):

            if mycluster.name == 'Zw3146':
                XMM_file = '/home/data/X-ray/XMM/ZW3146_density_v2.dat'
                try:
                    XMM_cols = np.loadtxt(XMM_file, comments='#')
                    xrads = XMM_cols[:,0]
                    xdens = XMM_cols[:,1]
                    xerrs = XMM_cols[:,2]
                    m_d_a = mycluster.d_a.to('kpc')
                    xdens = np.interp(theta_range*m_d_a.value,xrads,xdens)
                    dohse = True
                except:
                    xdens = np.zeros(len(theta_range))
                    dohse = False                    
            else:
                xdens = np.zeros(len(theta_range))
                dohse = False

            #############################################
            myR500 = (mycluster.R_500/mycluster.d_a).decompose()
            myP500 = mycluster.P_500
            
            arrys,arrms,arrps = yMP500_from_r500(theta_range,mycluster,ySZ=True,ySph=ySph,YMrel=YMrel)
            Y2500    = y_delta_from_mdelta(5*arrms.value,mycluster,delta=2500,ySph=ySph,YMrel=YMrel)
            #import pdb;pdb.set_trace()
            y2500.append(Y2500)
            ylist.append(arrys)
            mlist.append(arrms)
            edens.append(xdens)
            #print(mlist[0][0])
            
            if hk.cfp.model == 'NP':
                a10pres = cpp.a10_gnfw(myP500,myR500,hk.av.mycosmo,bulkbins)
                uless_p = (a10pres*myPdl2y).decompose().value
                sbepos.extend(np.ones((len(uless_p)),dtype=bool))
                #print(a10pres,bulkbins)
                #import pdb;pdb.set_trace()
                
            if hk.cfp.model == 'GNFW':
                uless_p = np.array([1.177,8.403,5.4905,0.3081])  # C500, P0,beta, gamma
                ### For cool-core clusters:
                #uless_p = np.array([1.128,3.249,5.4905,0.7736])  # C500, P0,beta, gamma
                #uless_p  = np.array([0.05,3.303,3.85,0.003])  # C500, P0,beta, gamma
                #import pdb;pdb.set_trace()
                if len(bulkbins) < len(uless_p):
                    #uless_p = uless_p[1:1+len(bulkbins)] # If I want to keep C500 fixed
                    uless_p = uless_p[0:len(bulkbins)]  # If I want to vary C500 and P0
                myextarr = np.ones((len(bulkbins)),dtype=bool)
                #myextarr[-1] = 0 # This doesn't make sense!?!
                sbepos.extend(myextarr)

            if hk.cfp.model == 'BETA':
                uless_p = np.array([myP500.to('keV / cm**3').value*10.0 ,myR500/10.0,1.0])
                sbepos.extend(np.ones((len(uless_p)),dtype=bool))
                
            ###################################################################################
            
            myval.extend(uless_p)
            myalp.extend(uless_p*0.0) # What if I want to use strict alphas??    
            ### I need to play with centroids (updated how I do the following):
            priors.extend(uless_p)
            priunc.extend(uless_p*0.0 -1.0)
            compname.extend(['bulk' for x in uless_p])
            compcount.extend([bulkcount for x in uless_p])
            nbulkp+=len(uless_p)
            ####################################################################################
                
            if fit_cen == True:
                myval.extend([1.0,1.0]) # In particular, how much should I allow this to vary?
                #myalp.extend([0.0,0.0]) # In particular, how much should I allow this to vary?
                # What units am I using for this?? (arcseconds? Pixel size? radians?)
                sbepos.extend([False,False])
                priors.extend([0.0,0.0])
                priunc.extend([-1.0,-1.0])
                compname.extend(['bulk','bulk']) # Important to keep this in the same compname!
                compcount.extend([bulkcount,bulkcount])
                ncent+=2
                print('Extending parameters to fit for a bulk centroid.')

            if fit_geo == True:
                #myval.extend([1.0,0.1]) # In particular, how much should I allow this to vary?
                myval.extend([1.0,0.8]) # In particular, how much should I allow this to vary?
                #myalp.extend([0.0,0.0]) # In particular, how much should I allow this to vary?
                # What units am I using for this?? (arcseconds? Pixel size? radians?)
                sbepos.extend([True,True])
                priors.extend([0.0,0.0])
                priunc.extend([-1.0,0.2])
                compname.extend(['bulk','bulk']) # Important to keep this in the same compname!
                compcount.extend([bulkcount,bulkcount])
                ngeo+=2
                print('Extending parameters to fit for a bulk geometry.')

            bulkcount+=1

        #import pdb; pdb.set_trace()

        for scount,shockbins in enumerate(hk.cfp.shockbin):
            myshockbins = shockbins
            if hk.cfp.shockfin[scount] == True:
                myshockbins = shockbins[:-1]
            a10pres = cpp.a10_gnfw(hk.cluster.P_500,R500,hk.av.mycosmo,myshockbins)
            uless_p = (a10pres*Pdl2y).decompose().value
            myval.extend(uless_p)
            myalp.extend(uless_p*0.0) # What if I want to use strict alphas??
            sbepos.extend(np.ones((len(uless_p)),dtype=bool))
            priors.extend(uless_p)
            priunc.extend(uless_p*0.0 -1.0)
            compname.extend(['shock' for x in uless_p])
            compcount.extend([scount for x in uless_p])
            nshockp+=len(uless_p)

        ### I need to think a bit more about whether this is best or not.
        for iptsrc, myptsrc in enumerate(hk.cfp.ptsrc):
            for myinst in hk.instruments:
                if ifp[myinst].pt_src == True:
                    print('Map units are: ',u.Unit(dv[myinst].maps.units))
                    print('Point source prior units are: ', ifp[myinst].prior[iptsrc].unit)
                    try:
                        if u.Unit(dv[myinst].maps.units) == ifp[myinst].prior[iptsrc].unit:
                            print('Setting Point Source Initial Guess via Prior.')
                            ptpr = ifp[myinst].prior[iptsrc]
                            ptun = ifp[myinst].priorunc[iptsrc]
                        else:
                            print('You need to convert the prior to the proper units.')
                            print('If it is just between Jy and K, we might be able to fix that.')
                            convunits = u.Unit(dv[myinst].maps.units)/ifp[myinst].prior[iptsrc].unit
                            factunits = u.Unit(dv[myinst].conversion.unit)
                            invfunits = u.Unit(1.0/dv[myinst].conversion.unit)
                            if convunits == factunits:
                                ptpr = ifp[myinst].prior[iptsrc] * dv[myinst].conversion
                                ptun = ifp[myinst].priorunc[iptsrc] * dv[myinst].conversion
                                print('Converting using ',dv[myinst].conversion)
                            elif convunits == invfunits:
                                ptpr = ifp[myinst].prior[iptsrc] / dv[myinst].conversion
                                ptun = ifp[myinst].priorunc[iptsrc] / dv[myinst].conversion
                                print('Converting using ', 1.0/dv[myinst].conversion)
                            else:
                                print('Looks like you have some weird units their bud.')
                                print('I will just let the errors occur. You have been warned.')


                        ptsign = -1.0 if ptpr.value < 0 else 1.0
                        ptpr *= ptsign
                        ptsrcsign.extend([ptsign])
                        myval.extend([ptpr.value]);
                        priors.extend([ptpr.value])
                        priunc.extend([ptun.value])
                            
                    except:
                        myval.extend([0.003]) # Estimate ~3 mK?? for most point sources...to start.
                        priors.extend([0.003])
                        priunc.extend([-1.0])

                    sbepos.extend([True]) # But really we found ~5 mJy for RXJ1347. WTF?
                    compname.extend(['ptsrc'])
                    compcount.extend([nptsrc])

                    #if hk.cfp.ptcens[iptsrc] == True:
                    if hk.cfp.psrc_centroid[iptsrc] == True:
                        myval.extend([1.0,1.0])          # Xoffset, Yoffset
                        priors.extend([0.0,0.0])         # (arcsec)  , (arcsec)  , (radians)
                        priunc.extend([1.0,1.0])         # 1" pointing error...
                        sbepos.extend([False,False])     # No need to be strictly positive 
                        ngeo+=2
                        compname.extend(['ptsrc','ptsrc'])
                        compcount.extend([nptsrc,nptsrc])
                        
                    if hk.cfp.ptshape[iptsrc] == True:
                        myval.extend([5.0,4.0,-0.5])    # FWHM_major, FWHM_minor, Rotation angle
                        priors.extend([5.0,4.0,1.0])   # (arcsec)  , (arcsec)  , (radians)
                        priunc.extend([-1.0,-1.0,-1.0])  # Priors do not apply...
                        sbepos.extend([True,True,False]) # No need to restrict rotation angle here. 
                        ngeo+=3
                        compname.extend(['ptsrc','ptsrc','ptsrc'])
                        compcount.extend([nptsrc,nptsrc,nptsrc])

                    nptsrc+=1

        blobcount=0
        for myblob,bra0,bdec0 in zip(hk.cfp.blob,hk.cfp.bras,hk.cfp.bdecs):

            for myinst in hk.instruments:
                if ifp[myinst].fitblob == True:
                    #import pdb;pdb.set_trace()
                    if bra0 != 0:
                        myra0 = bra0.to('deg'); mydec0 = bdec0.to('deg')
                        bx0 ,by0 = dv[myinst].mapping.w.wcs_world2pix(myra0,mydec0,0)
                        myblob[0]=np.asscalar(bx0 - dv[myinst].mapping.x0)
                        myblob[1]=np.asscalar(by0 - dv[myinst].mapping.y0)
                    print(myblob)
                    if myblob[-1] < 0:
                        blobsign.extend([-1.0])
                        myblob[-1] *= -1
                    else:
                        blobsign.extend([1.0])
                    myval.extend(myblob) 
                    priors.extend(myblob)
                    #priunc.extend([1.0,1.0,2.0,2.0,-1.0,-1.0])
                    priunc.extend([0.1,0.1,1.0,1.0,-1.0,-1.0])
                    sbepos.extend([False,False,True,True,True,True])
                    compname.extend(['blob','blob','blob','blob','blob','blob'])
                    compcount.extend([blobcount for x in range(6)])
                    nblob+=6
                    blobcount+=1

        ### Some attributes for starting conditions / model creation.
        self.alphas    = myalp
        self.pinit     = myval
        self.thetas    = theta_range
        self.thetamax  = tnx[1]    # Maximum angular scale in profile (radians)
        self.thetamin  = tnx[0]    # Minimum angular scale in profile (radians)
        self.Pdl2y     = Pdl2y     # Conversion of unitless pressure to y?
        self.sbepos    = sbepos    # Boolean list of whether something should be positive.
        self.priors    = priors    # List of prior "known" values
        self.priunc    = priunc    # List of uncertainties on priors.
        self.compname  = compname  # List of component names
        self.compcount = compcount # With identifying markers...
        self.blobsign  = blobsign  # Positive or Negative?
        self.ptsrcsign = ptsrcsign
        ### Some attributes for the results:
        self.t_mcmc    = 0.0       # Time that MCMC took.
        self.dohse     = dohse     # Calculate a mass based on HSE?
        self.ySph      = ySph      # Calculate Ysph or Ycyl?
        ####################################################################################
        ### The following variables pertain to numbers used in / calculated during MCMC
        self.samples   = None      # The normal return of parameters used in MCMC
        self.solns     = None      # Just the "solutions" (median + errors)
        self.psolns    = None      # The solutions, in pressure units
        self.MCblobs   = None      # This is for extra variables calculated during MCMC
        self.values    = None      # The best-fit values 
        self.errors    = None      # and their corresponding errors
        self.nthreads  = np.min([nthreads,ncpus])     # Number of threads to run over with emcee.
        ####################################################################################
        ### I'm not sure if this is a good idea (for being more compact) or not.
        ### I think it should be fine. (WRT ifp variable) 21 Feb 2018.
        self.ifp       = ifp     # Carry around IFP within this class! 
        self.paramdict = {}
        self.punits    = "keV cm**-3"
        self.runits    = "arcsec"
        self.tag       = tag
        self.r500      = R500
        self.y500s     = ylist
        self.y2500s    = y2500
        self.m500s     = mlist
        self.edens     = edens
        self.YMrel     = YMrel
        
        ######################################################################################
        #print bcolors.UNDERLINE + bcolors.HEADER + bcolors.BOLD + 'Time-estimator-NIKA2 '+version+' '+bcolors.ENDC
        #print bcolors.HEADER + bcolors.BOLD + ' See the "Guidelines for observing time estimates with the NIKA2 continuum camera' + bcolors.ENDC
        #print bcolors.HEADER + bcolors.BOLD + ' at the IRAM-30m Telescope" for details on used parameters and calculations.\n' + bcolors.ENDC
        ######################################################################################
        print(bcolors.OKGREEN + '#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-' + bcolors.ENDC)
        print(bcolors.OKGREEN + 'Found the following number of parameters to be fit for each type of component:' + bcolors.ENDC)
        print(bcolors.OKGREEN + 'Mean Levels: \t\t\t{}'.format(nmnlvl)   + bcolors.ENDC)
        print(bcolors.OKGREEN + 'Bulk pressure: \t\t\t{} over \t{} components'.format(nbulkp,len(hk.cfp.bulkarc)) + bcolors.ENDC)
        print(bcolors.OKGREEN + 'Shock pressure: \t\t{} over \t{} components'.format(nshockp,len(hk.cfp.shockbin)) + bcolors.ENDC)
        print(bcolors.OKGREEN + 'Point Source Amplitudes: \t{} over \t{} components'.format(nptsrc,len(hk.cfp.ptsrc)) + bcolors.ENDC)
        print(bcolors.OKGREEN + 'Centroids (any component): \t{}'.format(ncent) + bcolors.ENDC)
        print(bcolors.OKGREEN + 'Gaussian (blob) components: \t{} (this includes its own centroid)'.format(nblob) + bcolors.ENDC)
        print(bcolors.OKGREEN + 'Bulk Geometry parameters: \t{}'.format(ngeo) + bcolors.ENDC)
        print(bcolors.OKGREEN + '------------------------------------------------------------------' + bcolors.ENDC)
        ncomp = nmnlvl+nbulkp+nshockp+nptsrc+ncent+nblob+ngeo
        nptgrid = 1 if nptsrc > 0 else 0
        ngrid = (len(hk.cfp.bulkarc) + len(hk.cfp.shockbin) +nptgrid)*len(hk.instruments)
        print(bcolors.OKGREEN + 'Total parameters: \t\t{} with \t{} griddings over all instruments'.format(ncomp,ngrid) + bcolors.ENDC)
        print(bcolors.OKGREEN + '#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-' + bcolors.ENDC)
        #import pdb;pdb.set_trace()
        #print('#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-'
        #print('Found the following number of parameters to be fit for each type of component:'
        #print('Mean Levels: ', nmnlvl  
        #print('Bulk pressure: ', nbulkp,' over ',len(hk.cfp.bulkarc), ' components'
        #print('Shock pressure: ', nshockp,' over ',len(hk.cfp.shockbin), ' components'
        #print('Point Source Amplitudes: ', nptsrc,' over ',len(hk.cfp.ptsrc), ' components'
        #print('Centroids (any component): ', ncent
        #print('Gaussian (blob) components: ', nblob, ' (this includes its own centroid)'
        #print('Bulk Geometry parameters: ', ngeo
        #print('#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-'
        
def bulk_or_shock_component(pos,bins,hk,dv,efv,fit_cen,fit_geo,geom,alphas,n_at_rmin,maps={},posind=0,
                            fixalpha=False,fullSZcorr=False,SZtot=False,columnDen=False,Comptony=True,
                            finite=False,oldvs=False,cluster=None,bsind=0,retcurvs=False):

    nbins = len(bins)
    if finite == True:
        nbins-=1     # Important correction!!!
    ulesspres = pos[posind:posind+nbins]
    #myalphas  = alphas[posind:posind+nbins]
    myalphas = alphas   # I've updated how I pass alphas; indexing no longer necessary! (16 Nov 2017)
    ulessrad  = bins #.to("rad").value
    posind = posind+nbins
    if fit_cen == True:
        geom[0:2] = pos[posind:posind+2]  # I think this is what I want...
        posind = posind+2

    #import pdb;pdb.set_trace()
        
    if fit_geo == True:
        #geom[2:4] = pos[posind:posind+2]
        geom[2]     = pos[posind]           # Rotation angle
        #geom[3]     = pos[posind+1]+1.0     # Major axis (should be > 1.0, if minor is defined as 1).
        geom[4]     = pos[posind+1]     # Major axis (should be > 1.0, if minor is defined as 1).
        geom[5]     = np.sqrt(pos[posind+1])     # Assume that the axis along the l.o.s. is geom. mean of maj,min axes.
        #pos[posind] = pos[posind] % (2.0*np.pi) #if tdtheta > 2.0*np.pi
        pos[posind] = pos[posind] % (np.pi) #if tdtheta > np.pi b/c symmetry reduces it to just pi. (July 2018)
        posind      = posind+2

    Tx        = cluster.Tx
    R500      = cluster.R_500
    P500      = cluster.P_500
    ang_diam  = cluster.d_a
        
    density_proxy, etemperature, geoparams = es.prep_SZ_binsky(ulesspres,Tx,geoparams=geom)
    ### Can modify later to allow for X-ray images
    ### That is, I will want to do a loop over SZ images (reduce the number of integrations done),
    ### and then integrate over X-ray emissivity.
    #import pdb;pdb.set_trace()
    
    if fullSZcorr == False:
        #import pdb;pdb.set_trace()
        if hk.cfp.model == 'NP':
            Int_Pres,outalphas,integrals = es.integrate_profiles(density_proxy, etemperature, geom,bins,
                                           efv,hk,dv,myalphas,beta=0.0,betaz=None,finint=finite,
                                           narm=False,fixalpha=fixalpha,strad=False,array="2",SZtot=False,
                                            columnDen=False,Comptony=True)
        if hk.cfp.model == 'GNFW':
            radii    = (efv.thetas * ang_diam).to('kpc')     # This will be in kpc
            #import pdb;pdb.set_trace()
            myup     = np.array([1.177,8.403,5.4905,0.3081])
            if len(ulesspres) < len(myup):
                #myup[1:1+len(ulesspres)]=ulesspres
                myup[0:len(ulesspres)]=ulesspres  # Vary C500 and P0....
            else:
                myup = ulesspres

            A10all = 1.0510 # Alpha value in A10, all clusters
            A10cc  = 1.2223 # Alpha value in A10, cool-core clusters
            A10dis = 1.4063 # Alpha value in A10, disturbed clusters
            
            pprof    = cpp.gnfw(hk.av.mycosmo['h_70'], radii, P500, R500,
                                c500=myup[0], p=myup[1], a=A10all, b=myup[2], c=myup[3])
            #unitless_profile = pprof * efv.Pdl2y
            #unitless_profile = (pprof * hk.av.szcu['thom_cross'] * ang_diam / hk.av.szcu['m_e_c2']).decompose()
            unitless_profile = (pprof * hk.av.szcu['thom_cross'] * u.kpc / hk.av.szcu['m_e_c2']).decompose()
            #inrad = radii.to("kpc"); zvals = radProjected.to("kpc")

            Int_Pres = ni.int_profile(radii.value, unitless_profile.value,radii.value)
            outalphas = unitless_profile*0.0+2.0
            integrals = Int_Pres
            #print(pprof[50],unitless_profile[50])
            #import pdb;pdb.set_trace()

        if hk.cfp.model == 'BETA':
            d_a      = ang_diam.to('kpc').value
            radii    = efv.thetas 
            pprof     = ulesspres[0]*(1.0+(radii/ulesspres[1])**2)**(-1.5*ulesspres[2])        ### Beta model
            #scaling  = scs.gamma(1.5*ulesspres[2] - 0.5)/scs.gamma(1.5*ulesspres[2]) * ulesspres[1] * ulesspres[0]
            scaling  = scs.gamma(1.5*ulesspres[2]-0.5)/scs.gamma(1.5*ulesspres[2])*ulesspres[1] * d_a
            scaling *= ulesspres[0] * np.sqrt(np.pi) *(hk.av.szcu['thom_cross'] * u.kpc / hk.av.szcu['m_e_c2']).to('cm**3 / keV')
            Int_Pres    = scaling.value * (1.0+(radii/ulesspres[1])**2)**(0.5-1.5*ulesspres[2])
            outalphas = Int_Pres*0.0+2.0
            integrals = Int_Pres
            #yint = ni.Ycyl_from_yProf(yProf,ppbins,r500)
        
            #import pdb;pdb.set_trace()
        
        #yint=es.ycylfromprof(Int_Pres,efv.thetas,efv.thetamax) #
        #####################################
        #yint=es.ycylfromprof(Int_Pres,efv.thetas,efv.r500[bsind]) # As of July/August 2018
        #####################################
        #yint,newr500=consistent_Y_SZ(Int_Pres,efv.thetas,efv.r500[bsind],cluster.hofz,cluster.d_a.value/1000.0,
        #                             cluster.dens_crit,hk.av.mycosmo['h_70']) # As of Aug. 27th, 2018

        ### Both hse and ySph require an interpolated pressure profile...
        if efv.dohse or efv.ySph:
            if hk.cfp.model == 'NP':
                pprof,alphas   = es.log_profile(ulesspres,list(bins),efv.thetas) # Last pos is mn lvl
            elif hk.cfp.model == 'GNFW':
                pprof = unitless_profile *ang_diam.value
        #############################################################################################3
        if efv.dohse:
            palphas,norm  = es.ycyl_prep(pprof,efv.thetas)
            hsem500, hser500,miso,riso,hsem2500,hser2500,Mhse,Miso,M500 = hse_with_density_v2(pprof,palphas,efv.thetas,efv.edens[bsind],cluster.Tx,efv.m500s[bsind].value)
            #m500,r500,miso,riso,m2500,r2500,Mtot,Miso,m500v
        else:
            hsem500, hser500,miso,riso,hsem2500,hser2500 = 0,0,0,0,0,0
            Mhse,Miso,M500 = np.zeros(len(efv.thetas)),np.zeros(len(efv.thetas)),np.zeros(len(efv.thetas))
        
        if efv.ySph:
            yint ,newr500,y2500,r2500=Y_SZ_v2(pprof,efv.thetas,efv.r500[bsind],cluster,hk.av.mycosmo['h_70'],efv.y500s[bsind],geom,
                                              ySph=efv.ySph,yR2500=efv.y2500s[bsind],retcurvs=retcurvs) # As of Aug. 31, 2018
        else:
            yint ,newr500,y2500,r2500=Y_SZ_v2(Int_Pres,efv.thetas,efv.r500[bsind],cluster,hk.av.mycosmo['h_70'],efv.y500s[bsind],geom,
                                              ySph=efv.ySph,yR2500=efv.y2500s[bsind],retcurvs=retcurvs) # As of Aug. 31, 2018

        #print(yint,pprof[0],pprof[50],pprof[100])
        #import pdb;pdb.set_trace()
        ### I think I can just do "for myinst in hk.instruments:"
    for i in range(len(hk.instruments)):
        myinst = hk.instruments[i]; xymap=dv[myinst].mapping.xymap
        if fullSZcorr == True:
            IntProf,outalphas,integrals = es.integrate_profiles(density_proxy, etemperature, geom,bins,
                 efv,hk,dv,myalphas,beta=0.0,betaz=None,finint=finite,narm=False,fixalpha=fixalpha,
                 strad=False,array="2",SZtot=True,columnDen=False,Comptony=False)
            if efv.ySph:
                if hk.cfp.model == 'NP':
                    pprof,alphas   = es.log_profile(ulesspres,list(bins),efv.thetas) # Last pos is mn lvl
                elif hk.cfp.model == 'GNFW':
                    pprof = unitless_profile *ang_diam.value # Derp.
                yint ,newr500,y2500,r2500=Y_SZ_v2(pprof,efv.thetas,efv.r500[bsind],cluster,hk.av.mycosmo['h_70'],efv.y500s[bsind],geom,
                                                  ySph=efv.ySph,yR2500=efv.y2500s[bsind],retcurvs=retcurvs) # As of Aug. 31, 2018
            else:
                yint ,newr500,y2500,r2500=Y_SZ_v2(Int_Pres,efv.thetas,efv.r500[bsind],cluster,hk.av.mycosmo['h_70'],efv.y500s[bsind],geom
                                                  ,yR2500=efv.y2500s[bsind],ySph=efv.ySph,retcurvs=retcurvs) # As of Aug. 31, 2018
            yint=0 # A sure way to give something clearly wrong -> bring myself back here.
            import pdb;pdb.set_trace() # A better way...
        else:
            ### Convert to Column Density (for kSZ)....?
            #import pdb;pdb.set_trace() # A better way...
            ConvtoCD= hk.av.szcv["m_e_c2"]/(hk.av.szcv["boltzmann"]*cluster.Tx)
            IntProf = Int_Pres * (dv[myinst].tSZ[cluster.number-1] +
                                  dv[myinst].kSZ[cluster.number-1]*ConvtoCD)
            integrals = integrals * (dv[myinst].tSZ[cluster.number-1] +
                                     dv[myinst].kSZ[cluster.number-1]*ConvtoCD)

        ### Right...I need to have a zero-map to start with. Ughhh.
        #import pdb;pdb.set_trace()
        maps[myinst] += es.general_gridding(xymap,efv.thetas,bins,geom,finite=finite,integrals=integrals,
                                            Int_Pres=IntProf,oldvs=oldvs)

    if retcurvs:
        return maps,posind,yint,outalphas,hsem500,miso,y2500,hsem2500,Mhse,Miso,M500
    else:
        return maps,posind,yint,outalphas,hsem500,miso,y2500,hsem2500

def mk_twodgauss(pos,posind,mysign,hk,dv,maps={}):

    ### UNDER DEVELOPMENT!
    for myinst in hk.instruments:
        x,y     = dv[myinst].mapping.xymap            # Defined in arcseconds
        x0,y0,sx,sy,tdtheta,peak =  pos[posind:posind+6]
        ### Added 1 Dec 2021 (x0,y0 defined in pixel coords...why??)
        #import pdb;pdb.set_trace()
        pixs = dv[myinst].mapping.pixsize.to('arcsec').value
        ### <<<END ADDED LINE>>>
        xcoord  = (x - x0*pixs)/sx;  ycoord= (y - y0*pixs)/sy
        xrot = xcoord * np.cos (tdtheta) + ycoord * np.sin(tdtheta)
        yrot = ycoord * np.cos (tdtheta) - xcoord * np.sin(tdtheta)
        mygauss = np.exp(((-(xrot)**2 - (yrot)**2) )/(2.0))
        maps[myinst] += mygauss*peak*mysign
        #print(peak,min(maps[myinst])
        #import pdb;pdb.set_trace()
        pos[posind+4] = tdtheta % (2.0*np.pi) #if tdtheta > 2.0*np.pi
        posind += 6

    return maps,posind

### Need to rewrite this to correctly deal with point source.
def mk_ptsrc_v2(pos,posind,mysign,hk,dv,ifp,maps={}):

    for index in range(len(hk.cfp.ptsrc)):
        for myinst in hk.instruments:
            ### Centroid is initially defined by RA, Dec.
            if ifp[myinst].pt_src == True:
                myflux  = pos[posind]*mysign[index]
                x,y     = dv[myinst].mapping.xymap            # Defined in arcseconds
                fwhm    = dv[myinst].fwhm.to("arcsec").value  # Now this matches.
                if hk.cfp.psfwhm[index] > fwhm:
                    #print("Adjusted FWHM"
                    fwhm = hk.cfp.psfwhm[index]
                sigma    = fwhm/(2.0 * np.sqrt((2.0 * np.log(2.0))))
                #import pdb; pdb.set_trace()
                centroid = dv[myinst].mapping.ptsrclocs[index]
                xcoord   = x - centroid[0];  ycoord=y - centroid[1]
#                xcoord   = x - np.asscalar(centroid[0])
#                ycoord   = y - np.asscalar(centroid[1])
                mygauss  = np.exp(((-(xcoord)**2 - (ycoord)**2) )/( 2.0 * sigma**2))
                maps[myinst] += mygauss*myflux
                posind += 1   # posind must be incremented for *each* instrument!

                #print(np.max(maps[myinst]),np.min(maps[myinst])
                #print(centroid)
                #import pdb; pdb.set_trace()
        
    return maps,posind

def mk_ptsrc_v3(pos,posind,mysign,hk,dv,ifp,maps={}):
  
#    for index in range(len(hk.cfp.ptsrc)):
#        
#        extflux = hk.cfp.ptextint[index]
#        extfreq = hk.cfp.ptextfre[index]
#
#        if ifp['MUSTANG2'].pt_src == True:
#        for myinst in hk.instruments:
#            if myinst != 'MUSTANG2':
#                if ifp[myinst].link2M2[index] == True:
#                    
#                    
######################################################################3                

    linkind = posind+0
    M2flux=[]
    for index in range(len(hk.cfp.ptsrc)):
        M2flux.append(0.0)
        
        for myinst in hk.instruments:

            if ifp[myinst].pt_src == True:
                myflux  = pos[posind]*mysign[index]
                x,y     = dv[myinst].mapping.xymap            # Defined in arcseconds
                fwhm    = dv[myinst].fwhm.to("arcsec").value  # Now this matches.
                #import pdb;pdb.set_trace()
                if hk.cfp.psfwhm.ndim > 1:
                    ptshape = hk.cfp.psfwhm[index]
                    aa      = ptshape[0]/(2.0 * np.sqrt((2.0 * np.log(2.0))))
                    bb      = ptshape[1]/(2.0 * np.sqrt((2.0 * np.log(2.0))))
                    rot     = ptshape[2]
                else:
                    if hk.cfp.psfwhm[index] > fwhm:
                        fwhm = hk.cfp.psfwhm[index]
                    aa = fwhm/(2.0 * np.sqrt((2.0 * np.log(2.0))))
                    bb = fwhm/(2.0 * np.sqrt((2.0 * np.log(2.0))))
                    rot= 0.0
                #import pdb; pdb.set_trace()
                centroid = dv[myinst].mapping.ptsrclocs[index]
                xcoord   = x - centroid[0];  ycoord=y - centroid[1]
                #if hk.cfp.ptcens[index] == True:
                if hk.cfp.psrc_centroid[index] == True:
                    xcoord += pos[posind+1]
                    ycoord += pos[posind+2]
                    posind += 2   # posind must be incremented for *each* instrument!

                if hk.cfp.ptshape[index] == True:
                    #print("Adjusted FWHM"
                    ### Keep aa as the major axis.
                    if pos[posind+2] > pos[posind+1]:
                        majaxis = pos[posind+2]*1.0
                        pos[posind+2] = pos[posind+1]*1.0
                        pos[posind+1] = majaxis
                    aa  = pos[posind+1]
                    bb  = pos[posind+2]
                    pos[posind+3] = pos[posind+3] % (np.pi) #if tdtheta > np.pi b/c symmetry reduces it to just pi. (Nov 2018)
                    rot = pos[posind+3]
                    posind += 3   # posind must be incremented for *each* instrument!

                if rot != 0:
                    newx = xcoord*np.cos(rot) + ycoord*np.sin(rot)
                    newy = ycoord*np.cos(rot) - xcoord*np.sin(rot)
                    xcoord = newx
                    ycoord = newy
                
                mygX  = np.exp((-(xcoord)**2) /( 2.0 * aa**2))
                mygY  = np.exp((-(ycoord)**2) /( 2.0 * bb**2))
                mygauss = mygX * mygY
                maps[myinst] += mygauss*myflux
                posind += 1

                if myinst == 'MUSTANG2':
                    Nbeams =  2.0*np.pi*aa*bb/dv[myinst].bv.value
                    M2flux[index] = Nbeams*myflux
                    if dv[myinst].maps.units == 'Kelvin':
                        M2flux[index] *= dv[myinst].Jy2K
                
                #print(np.max(maps[myinst]),np.min(maps[myinst])
                #import pdb; pdb.set_trace()

    ### Now, if any other instruments are tied explicitly to MUSTANG-2, do that:
            
    #maps = mk_ptsrc_from_M2(hk,dv,ifp,M2flux,maps=maps)
        
    return maps,posind


def mk_ptsrc_from_M2(hk,dv,ifp,M2fluxes,maps={}):

    for index in range(len(hk.cfp.ptsrc)):
        if M2fluxes[index] == 0:
            continue
        else:
            ### In mJy; values:
            otherfluxes = hk.cfp.ptextint[index]
            otherfluxes.append(M2fluxes[index]*1.0e3)
                                                 
            ### In GHz; values:
            otherfreqs  = hk.cfp.ptextfre[index]
            otherfreqs.append(dv['MUSTANG2'].freq.value)

            #import pdb;pdb.set_trace()
            lx = np.log(np.array(otherfluxes)/1000.0) # Back to Jy
            ly = np.log(np.array(otherfreqs))
            #sindex
            slope, intercept, r_value, p_value, std_err = stats.linregress(lx,ly)
            S0 = np.exp(intercept)
            
            for myinst in hk.instruments:

                if ifp[myinst].link2M2[index] == True:
                    myfreq   = dv[myinst].freq.value              # in GHz
                    instflux = S0 * myfreq**slope                 # Voila
                    #import pdb;pdb.set_trace()
                    x,y     = dv[myinst].mapping.xymap            # Defined in arcseconds
                    fwhm    = dv[myinst].fwhm.to("arcsec").value  # Now this matches.
                    sigma    = fwhm/(2.0 * np.sqrt((2.0 * np.log(2.0))))
                    centroid = dv[myinst].mapping.ptsrclocs[index]
                    xcoord   = x - centroid[0];  ycoord=y - centroid[1]
                    mygauss  = np.exp(((-(xcoord)**2 - (ycoord)**2) )/( 2.0 * sigma**2))
                    #print(instflux,intercept,slope,myinst,dv[myinst].Jy2K,sigma
                    if instflux > 1.0e-2:
                        print(instflux,intercept,slope,myinst,dv[myinst].Jy2K,sigma)
                        import pdb;pdb.set_trace()
                    maps[myinst] += mygauss*instflux/dv[myinst].Jy2K
        
    return maps


### Need to rewrite this to correctly deal with point source.
def mk_ptsrc(pos,posind,centroid,hk,dv,ifp,maps={}):

    for myinst in hk.instruments:
        ### Centroid is initially defined by RA, Dec.
        if ifp[myinst].pt_src == True:
            myflux  = pos[posind]
            x,y     = dv[myinst].mapping.xymap            # Defined in arcseconds
            fwhm    = dv[myinst].fwhm.to("arcsec").value  # Now this matches.
            sigma   = fwhm/(2.0 * np.sqrt((2.0 * np.log(2.0))))
            xcoord  = x - centroid[0];  ycoord=y - centroid[1]
            mygauss = np.exp(((-(xcoord)**2 - (ycoord)**2) )/( 2.0 * sigma**2))
            maps[myinst] += mygauss*myflux
            posind += 1   # posind must be incremented for *each* instrument!
        
    return maps
             
def run_emcee(hk,dv,ifp,efv,init_guess=None,BSerr=False):

##################################################################################
### PSA / Reminder: EMCEE doesn't like it if the likelihood doesn't change much with respect
### to how many data points are being used. Or something like that.
    
    def lnlike(pos):                          ### emcee_fitting_vars
        #posind = 0
        #maps={}
        ### I need to have it return a "pressure profile parameter list only" (15 Feb 2018)
        mnlvlmaps,ptsrcmaps,maps,yint,outalphas,mhse,miso,y2500,h2500 = make_skymodel_maps(pos,hk,dv,ifp,efv)
        ### Previously, I allowed for multiple yints, thinking that ONE cluster might have subclusters.
        ### Now (Aug. 2018) I see that I may have more than 1 CLUSTER in a map (e.g. Lynx).
        damult = [(cluster.hofz**(-1./3)*cluster.d_a.value/1000.0)**2 for cluster in hk.cluster]
        #import pdb;pdb.set_trace()
        ytotint = [y*d for y,d in zip(yint,damult)] # I don't know how I would actually distinguish multiple yints...
        #ycyl    = np.sum(ytotint) # Now in Mpc**-2 (should be a value of ~10 e-5)
        ycyl    = ytotint         # Now in Mpc**-2 (should be a value of ~10 e-5)

        #print(ycyl)
        #import pdb;pdb.set_trace()
        models  = filter_maps(hk,dv,maps,ptsrcmaps,mnlvlmaps,mapin1d=True,BSerr=BSerr) # mapin1d added Aug. 2018
        loglike = 0.0
        #import pdb;pdb.set_trace()
        #import plot_mcmc_results as pmr
        #pmr.plot_best_sky(models,hk.hk_outs.newpath,hk.hk_outs.prefilename,dv,mycomp='test')
        #pmr.plot_best_sky(maps,hk.hk_outs.newpath,hk.hk_outs.prefilename,dv,mycomp='test1')
        #print 'LogLike initialized to zero'
        for myinst in hk.instruments:
            data    = dv[myinst].maps.data
            weights = dv[myinst].maps.masked_wts       # Best to mask weights (not use the entire map!)
            model   = models[myinst]
            loglike-= 0.5 * (np.sum(((model - data)**2) * weights))
            ### If I choose to apply a mask:
            mask=np.where(weights > 0)
            gim=model[mask]; gid=data[mask]


        #return loglike, outalphas,ycyl,mhse,miso,y2500,h2500
        return loglike, outalphas,yint,mhse,miso,y2500,h2500

##################################################################################    
    def fprior(pos):

        #prespos = pos[0:-1] # Define as separate variable here... 
        prespos = pos[efv.sbepos]
        #print(pos,efv.sbepos)
        #import pdb;pdb.set_trace()
        if all([param > 0.0 for param in prespos]):
            return True
        else:
            return False

    def lnprior(pos,outalphas,ycyl):        
        plike = 0.0
        if hk.cfp.pprior == True:
            plike = -0.5* (np.sum(((ycyl - yplanck)**2) * ywait)) #ywait = ycyl_weight(s)

        mypriind = (np.array(efv.priunc) > 0)
        pwithpri = np.array(pos)[mypriind]
        priorunc = np.array(efv.priunc)[mypriind]
        mypriors = np.array(efv.priors)[mypriind]
        addlike=0.0 
        if len(priorunc) > 0:
            addlike =  -0.5* (np.sum(((pwithpri - mypriors)**2) / priorunc**2))
        
        ### Pt src priors??
        #mnlvl_amp,ptsrc_amp,blob_amp = unpack_fit_params(pos,hk)
        #ulesspres = pos[0:hk.cfp.bins]
        #if hk.cfp.ptsrc == True:
        #    if len(hk.hk_ins.psfd) > 0:
### How do I want to deal with multiple point sources?
        #        pspr = hk.hk_ins.psfd[0]; psprunc = hk.hk_ins.psunc[0]
        #        plike += -0.5* (np.sum(((ptsrc_amp - pspr)**2)/(psprunc**2) ))

        prespos = pos[efv.sbepos]
        #import pdb;pdb.set_trace()
        for myouts in outalphas:
            if len(myouts) == 0:
                slopeok = True
            else:
                if myouts[-1] > 1.0: slopeok = True
                if myouts[-1] <= 1.0:
                    slopeok = False
                    break
            #print(slopeok)
            
        if all([param > 0.0 for param in prespos]) and (slopeok == True):
            #print 'Everything OK'
            return plike+addlike
        #print 'In LnPrior, LogLike set to infinity ', outalphas[-1]
        #import pdb;pdb.set_trace()
        return -np.inf

##################################################################################    
    def lnprob(pos):
        if fprior(pos):
            likel,outalphas,ycyl,mhse,miso,y2500,h2500 = lnlike(pos)
        else:
            likel = -np.inf
        #likel,outalphas,ycyl,mhse,miso = lnlike(pos)
        #yarr = np.array(ycyl)
        if not np.isfinite(likel):
            #print 'In LnProb, LogLike set to infinity'
            byc = [-1.0 for cluster in hk.cluster]
            bmh = [-1.0 for cluster in hk.cluster]
            bmi = [-1.0 for cluster in hk.cluster]
            return -np.inf, [byc,bmh,bmi,bmi,bmi]
            #return -np.inf, [-1.0 for ybad in ycyl]
        lp = lnprior(pos,outalphas,ycyl)
        #import pdb;pdb.set_trace()
        if not np.isfinite(lp):
            #print 'In LnProb, LogLike set to infinity'
            return -np.inf , [[-1.0 for ybad in ycyl] for i in range(5)]
        return lp + likel , [ycyl,mhse,miso,y2500,h2500]

#######################################################################################
###   ---   ###   ---   ###   ---   ###   ---   ###   ---   ###   ---   ###   ---   ###
#######################################################################################    
###   ---   ###   ---   ###   ---   ###   ---   ###   ---   ###   ---   ###   ---   ###
#######################################################################################    
###   ---   ###   ---   ###   ---   ###   ---   ###   ---   ###   ---   ###   ---   ###
#######################################################################################    

    t0 = time.time();    myargs = efv.pinit; dt0 = datetime.datetime.now()
    ndim = hk.cfp.ndim      #len(myargs)
    mnlvlmaps,ptsrcmaps,maps,yint,outalphas,mhse,miso,y2500,h2500 = make_skymodel_maps(myargs,hk,dv,ifp,efv)
    
    print(myargs)
    print("#########################################################################")
    ###timepergridding = 0.011 # seconds
    timepergridding = 0.0
    timeperpix = 4.0e-08
    for myinst in hk.instruments:
        myn = dv[myinst].mapping.nx*dv[myinst].mapping.ny
        timepergridding += (timeperpix*myn+0.0013) # 0.013
    
    numbergriddings = len(hk.cfp.bulkarc) + len(hk.cfp.shockbin) + len(hk.cfp.ptsrc)
    scaled_overhead = 1.2   # Other processing beyond gridding, e.g. FFTs
    flat_overhead   = 21.0  # The time to import modules??
    timepermodel    = timepergridding * numbergriddings * scaled_overhead
    esttime         = hk.cfp.nwalkers * hk.cfp.nsteps * timepermodel + flat_overhead
    print(bcolors.OKBLUE + 'Your estimated run duration is ',esttime/60.0,' minutes' + bcolors.ENDC)
    proccheck = 100  # Check in every 100 steps (iterations)
    fts = int(np.ceil(hk.cfp.nsteps*1.0/proccheck))
    #imo = int(hk.cfp.ndim+1) # Python 2 version
    imo = int(hk.cfp.ndim) # Python 3 version
    #import pdb;pdb.set_trace()
    gw2010 = np.zeros((fts,imo))
    newmet = np.zeros((fts,imo))
    import pdb;pdb.set_trace()
    pos = [(myargs + np.random.randn(ndim) * myargs/1e3) 
           for i in range(hk.cfp.nwalkers)]

#pos must be shape(nwalkers, ndim)
    t_premcmc = time.time()
    sampler = emcee.EnsembleSampler(hk.cfp.nwalkers, hk.cfp.ndim,lnprob, threads = efv.nthreads)
    print(bcolors.OKBLUE + "#########################################################################" + bcolors.ENDC)
    print(bcolors.OKBLUE + "Running emcee on ",efv.nthreads," threads." + bcolors.ENDC)
    print(bcolors.OKBLUE + "Start time is: ",dt0.strftime("%Y-%m-%d %H:%M:%S") + bcolors.ENDC)
    #sampler = emcee.EnsembleSampler(nwalkers, ndim, log_prob, args=[weights,gi],a=4)
    state = sampler.run_mcmc(pos, hk.cfp.nsteps,progress=True)
    """
    for i, result in enumerate(sampler.sample(pos, iterations=hk.cfp.nsteps)):
    #    print(i)
        if (i+1) % proccheck == 0:
            cind     = int((i+1) / proccheck)-1
            for jjj in range(hk.cfp.ndim):
                #dcheck = sampler.chain.shape
                #print(dcheck)
                gw2010[cind,:]   = sampler.get_autocorr_time
                #gw2010[cind,jjj]   = emcee_stats.autocorr_gw2010(sampler.chain[:,:i+1,jjj])
                #newmet[cind,jjj]   = emcee_stats.autocorr_new(sampler.chain[:,:i+1,jjj])
            #gw2010[cind,-1]   = emcee_stats.autocorr_gw2010(sampler._lnprob.T)
            #newmet[cind,-1]   = emcee_stats.autocorr_new(sampler._lnprob.T)
            
            t_so_far = time.time() - t_premcmc
            perdone  = float(i+1) / hk.cfp.nsteps
            t_total  = t_so_far / perdone
            t_remain = (t_total * (1.0 - perdone) * u.s).to("min")
            t_finish = dt0 + datetime.timedelta(seconds=t_total)
            print("{0:5.1%}".format(perdone)+' done; >>> Estimated Time Remaining: ',\
                t_remain.value,' minutes; for a finish at: ',t_finish.strftime("%Y-%m-%d %H:%M:%S"))
            #print "Average time per step so far: ", "{:5.1f}".format(t_so_far/(i+1.0))," seconds."

    ### sampler.chain should be in the dimensions of: (nwalkers,nsteps,nparams)...
    #import pdb;pdb.set_trace()
"""

    ### This will cause problems for the 'Test' case.
    #myabscissa = np.arange(np.ceil(hk.cfp.nsteps*1.0/ proccheck))*proccheck
    #ConvTests={'GoodmanWeare2010':gw2010,'Fardal_emcee':newmet,'Abscissa':myabscissa}
    ConvTests={'GoodmanWeare2010':sampler.get_autocorr_time(quiet=True)}

    #for jjj in range(hk.cfp.ndim):
    #    gw2010 = emcee_stats.autocorr_gw2010(sampler.chain[:,:,jjj])
    #    newmet = emcee_stats.autocorr_new(sampler.chain[:,:,jjj])
    #sampler.run_mcmc(pos,hk.cfp.nsteps)
    
    t_mcmc = time.time() - t_premcmc
    dt_end = datetime.datetime.now()
    
    print("MCMC time: ",t_mcmc/60.0,' minutes')
    print("Difference from predicted: ", (t_mcmc - esttime),' seconds')
    print("Finishing time: ", dt_end.strftime("%Y-%m-%d %H:%M:%S"))
    print("This is equivalent to ",t_mcmc/hk.cfp.nsteps," seconds per step.")
    print("Initial Guesses: ", efv.pinit)

    #import pdb;pdb.set_trace()

    return sampler, t_mcmc, ConvTests
    
def make_skymodel_maps(pos,hk,dv,ifp,efv):

    posind = 0
    mnlvlmaps={};ptsrcmaps={};maps={}; yint=[];miso=[];mhse=[];y2500=[];h2500=[]; outalphas=[]
    # pressure terms.
    for myinst in hk.instruments:
        nx                = dv[myinst].mapping.nx
        ny                = dv[myinst].mapping.ny
        maps[myinst]      = np.zeros((nx*ny))  # These maps will get convolved with a beam
        ptsrcmaps[myinst] = np.zeros((nx*ny))  # This comes "pre-convolved" with the beam
        # And this is not astronomical signal-related, so we don't want to filter it at all.
        if ifp[myinst].mn_lvl == True:
            mnlvlmaps[myinst] = pos[posind]       
            posind+=1
        else:
            mnlvlmaps[myinst] = 0
            
    testinst='MUSTANG2'
    ### Model the bulk pressure:
    for bins,fit_cen,fit_geo,geo,alp,narm,mycluster,cind in zip(hk.cfp.bulkarc,hk.cfp.bulk_centroid,hk.cfp.bulk_geometry,
                                                           hk.cfp.bulkgeo,hk.cfp.bulkalp,hk.cfp.bulknarm,hk.cluster,
                                                           range(len(hk.cluster))):
        #print 'Bulk map: ', maps[testinst].shape
        #import pdb;pdb.set_trace()
        #########################   NEED TO ADRESS MYCLUSTER (AUG. 11 2018)
        maps,posind,ynt,myalphas,mhe,mis,y25,h25 = bulk_or_shock_component(pos,bins,hk,dv,efv,fit_cen,fit_geo,geo,alp,narm,
                                                           maps,posind,fixalpha=hk.cfp.bulkfix,cluster=mycluster,
                                                                     bsind=cind)
        #print('I really hate this   ',ynt)
        yint.extend([np.real(ynt)]); outalphas.extend([np.real(myalphas)])
        mhse.extend([np.real(mhe)]); miso.extend([np.real(mis)]);y2500.extend([np.real(y25)]); h2500.extend([np.real(h25)])
        
    ### Model any shocks:
    for bins,fit_cen,geo,alp,narm,sfin,sind in zip(hk.cfp.shockbin,hk.cfp.shoc_centroid,hk.cfp.shockgeo,
                                                   hk.cfp.shockalp,hk.cfp.shocknarm,hk.cfp.shockfin,
                                                   range(len(hk.cfp.shockbin))):

        fit_geo=[False] # For now, I don't want to try to do this with shocks.
        maps,posind,shint,shout,shse,siso,s2500,s2500 = bulk_or_shock_component(pos,bins,hk,dv,efv,fit_cen,fit_geo,geo,alp,narm,
                                                          maps,posind,fixalpha=hk.cfp.shockfix,
                                                                    finite=sfin,oldvs=False,bsind=sind)
        #print 'Shock map: ', maps[testinst].shape

    ### Model any point sources (hk.cfp.ptsrc is a 2-tuple, the pt src. centroid):
    #print 'Position index used for point source: ',posind,' with flux ',pos[posind]
    #################################################################################
    #ptsrcmaps,posind = mk_ptsrc_v2(pos,posind,efv.ptsrcsign,hk,dv,ifp,ptsrcmaps)
    ptsrcmaps,posind = mk_ptsrc_v3(pos,posind,efv.ptsrcsign,hk,dv,ifp,ptsrcmaps)

    ### Model any "blobs" (2D Gaussians):
    ### This is currently not functional because I'm not sure exactly how I want to implement
    ### this feature.
    for myblob,mysign in zip(hk.cfp.blob,efv.blobsign):
        #import pdb;pdb.set_trace()
        maps,posind = mk_twodgauss(pos,posind,mysign,hk,dv,maps)

        #print 'Blob map: ', maps[testinst].shape

        
    ### Add any mean levels:
    ### (To be added)

    #print(mnlvlmaps[testinst].shape, ptsrcmaps[testinst].shape, maps[testinst].shape)
    #import pdb;pdb.set_trace()

    return mnlvlmaps,ptsrcmaps,maps,yint,outalphas,mhse,miso,y2500,h2500

def filter_one_comp(hk,dv,maps,myinst=None):

    filt_img={}
    if type(myinst) == type(None):
        for myinst in hk.instruments:
            saveimg = maps[myinst].reshape(dv[myinst].mapping.nx,dv[myinst].mapping.ny)
            pixs=dv[myinst].mapping.pixsize
            psv = pixs.to('arcsec').value
            #import pdb;pdb.set_trace()
            beam_conv = ip.conv_gauss_beam(saveimg,pixs,dv[myinst].fwhm)
            if myinst == 'ACT90' or myinst == 'ACT150':
                filt_img[myinst] = beam_conv
            else:
                filt_img[myinst]=ip.apply_xfer(beam_conv, dv[myinst].mapping.tab,dv[myinst].mapping.tabdims,psv)
    else:
        saveimg = maps[myinst].reshape(dv[myinst].mapping.nx,dv[myinst].mapping.ny)
        pixs=dv[myinst].mapping.pixsize
        psv = pixs.to('arcsec').value
        #import pdb;pdb.set_trace()
        beam_conv = ip.conv_gauss_beam(saveimg,pixs,dv[myinst].fwhm)
        if myinst == 'ACT90' or myinst == 'ACT150':
            filt_img[myinst] = beam_conv + mnlvlmaps[myinst]
        else:
            filt_img[myinst]=ip.apply_xfer(beam_conv, dv[myinst].mapping.tab,dv[myinst].mapping.tabdims,psv)
        
    return filt_img
    
def filter_maps(hk,dv,maps,ptsrcmaps,mnlvlmaps,mapin1d=True,BSerr=False):
    """
    This combines any sky signals (models) which have not yet been convolved with the appropriate telescope
    beam/PSF and subsequent transfer function, with any point source models (i.e. maps already beam-convolved),
    and finally adds in any mean level offsets. Even components are not being fitted, all variables should
    exist, even if only zero maps.


    ---- HISTORY:
    Last edit: 19 Mar 2018

    """
    
    models={}
    for myinst in hk.instruments:
        if mapin1d == True:
            ptsrcmaps[myinst] = ptsrcmaps[myinst].reshape((dv[myinst].mapping.nx,dv[myinst].mapping.ny))
            maps[myinst]      = maps[myinst].reshape((dv[myinst].mapping.nx,dv[myinst].mapping.ny))
            #mnlvlmaps[myinst] = mnlvlmaps[myinst].reshape(dv[myinst].mapping.nx,dv[myinst].mapping.ny)

        pixs=dv[myinst].mapping.pixsize
        psv = pixs.to('arcsec').value
        beam_conv = ip.conv_gauss_beam(maps[myinst],pixs,dv[myinst].fwhm) + ptsrcmaps[myinst]
        #models[myinst]=ip.apply_xfer(beam_conv, dv[myinst].mapping.tab,myinst,psv,BSerr=BSerr) + mnlvlmaps[myinst]
        if myinst == 'ACT90' or myinst == 'ACT150':
            models[myinst] = beam_conv + mnlvlmaps[myinst]
        else:
            models[myinst]=ip.apply_xfer(beam_conv, dv[myinst].mapping.tab,dv[myinst].mapping.tabdims,psv,BSerr=BSerr) + mnlvlmaps[myinst]

    return models

######################################################################################################
# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +#
#+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +#
#+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +#
#+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +#
#+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +#
######################################################################################################
# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +#
#+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +#
#+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +#
#+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +#
#+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +#
######################################################################################################
# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +#
#+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +#
#+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +#
#+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +#
#+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +#
######################################################################################################

def post_mcmc(sampler,t_mcmc,efv,hk,dv,ifp,BSerr=False):

### Note: I have excluded SAMPLER from the EFV class because it causes the following error:
### Can't pickle <function lnprob at 0x7fb2bd2c6848>: it's not found as max_like_fitting.lnprob

    txtout     = open(os.path.join(hk.hk_outs.newpath,'MCMC_key_results.txt'),'w')
    efv.t_mcmc = t_mcmc
    efv.samples = sampler.chain[:,hk.cfp.burn_in:, :].reshape((-1,hk.cfp.ndim))
    blobarr     = np.array(sampler.blobs)
    nbulk       = len(hk.cfp.bulkbins)
    newdim = 5 # if efv.dohse else 1
    efv.MCblobs = blobarr[hk.cfp.burn_in:,:,:].reshape((-1,newdim,nbulk))
    print(efv.samples.shape)
    print(efv.MCblobs.shape)
    
    #import pdb;pdb.set_trace()

    #myf = lambda v: (v[1], v[2]-v[1], v[1]-v[0])
    #myv = np.percentile(efv.samples, [16, 50, 84],axis=0)
    #yrv = myf(myv)
    efv.solns = np.array(list(map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                             zip(*np.percentile(efv.samples, [16, 50, 84],
                                                axis=0)))))
    mybulk     = np.array([(mcn == 'bulk') for mcn in efv.compname])
    mybcounts  = np.array(efv.compcount)[mybulk]
    myset      = np.unique(mybcounts)
    #import pdb;pdb.set_trace()
    psoln      = efv.solns[mybulk,:] * u.Unit(efv.punits)

    def lnlike(pos):                          ### emcee_fitting_vars
        #posind = 0
        #maps={}
        ### I need to have it return a "pressure profile parameter list only" (15 Feb 2018)
        mnlvlmaps,ptsrcmaps,maps,yint,outalphas,mhse,miso,y2500,h2500 = make_skymodel_maps(pos,hk,dv,ifp,efv)
        ### Previously, I allowed for multiple yints, thinking that ONE cluster might have subclusters.
        ### Now (Aug. 2018) I see that I may have more than 1 CLUSTER in a map (e.g. Lynx).
        damult = [(cluster.hofz**(-1./3)*cluster.d_a.value/1000.0)**2 for cluster in hk.cluster]
        #import pdb;pdb.set_trace()
        ytotint = [y*d for y,d in zip(yint,damult)] # I don't know how I would actually distinguish multiple yints...
        #ycyl    = np.sum(ytotint) # Now in Mpc**-2 (should be a value of ~10 e-5)
        ycyl    = ytotint         # Now in Mpc**-2 (should be a value of ~10 e-5)

        #print(ycyl)
        #import pdb;pdb.set_trace()
        models  = filter_maps(hk,dv,maps,ptsrcmaps,mnlvlmaps,mapin1d=True,BSerr=BSerr) # mapin1d added Aug. 2018
        loglike = 0.0
        #import pdb;pdb.set_trace()
        #import plot_mcmc_results as pmr
        #pmr.plot_best_sky(models,hk.hk_outs.newpath,hk.hk_outs.prefilename,dv,mycomp='test')
        #pmr.plot_best_sky(maps,hk.hk_outs.newpath,hk.hk_outs.prefilename,dv,mycomp='test1')
        #print 'LogLike initialized to zero'
        for myinst in hk.instruments:
            data    = dv[myinst].maps.data
            weights = dv[myinst].maps.masked_wts       # Best to mask weights (not use the entire map!)
            model   = models[myinst]
            loglike-= 0.5 * (np.sum(((model - data)**2) * weights))
            ### If I choose to apply a mask:
            mask=np.where(weights > 0)
            gim=model[mask]; gid=data[mask]


        #return loglike, outalphas,ycyl,mhse,miso,y2500,h2500
        return loglike, outalphas,yint,mhse,miso,y2500,h2500

##################################################################################    
    def fprior(pos):

        #prespos = pos[0:-1] # Define as separate variable here... 
        prespos = pos[efv.sbepos]
        #print(pos,efv.sbepos)
        #import pdb;pdb.set_trace()
        if all([param > 0.0 for param in prespos]):
            return True
        else:
            return False

    def lnprior(pos,outalphas,ycyl):        
        plike = 0.0
        if hk.cfp.pprior == True:
            plike = -0.5* (np.sum(((ycyl - yplanck)**2) * ywait)) #ywait = ycyl_weight(s)

        mypriind = (np.array(efv.priunc) > 0)
        pwithpri = np.array(pos)[mypriind]
        priorunc = np.array(efv.priunc)[mypriind]
        mypriors = np.array(efv.priors)[mypriind]
        addlike=0.0 
        if len(priorunc) > 0:
            addlike =  -0.5* (np.sum(((pwithpri - mypriors)**2) / priorunc**2))
        
        ### Pt src priors??
        #mnlvl_amp,ptsrc_amp,blob_amp = unpack_fit_params(pos,hk)
        #ulesspres = pos[0:hk.cfp.bins]
        #if hk.cfp.ptsrc == True:
        #    if len(hk.hk_ins.psfd) > 0:
### How do I want to deal with multiple point sources?
        #        pspr = hk.hk_ins.psfd[0]; psprunc = hk.hk_ins.psunc[0]
        #        plike += -0.5* (np.sum(((ptsrc_amp - pspr)**2)/(psprunc**2) ))

        prespos = pos[efv.sbepos]
        #import pdb;pdb.set_trace()
        for myouts in outalphas:
            if len(myouts) == 0:
                slopeok = True
            else:
                if myouts[-1] > 1.0: slopeok = True
                if myouts[-1] <= 1.0:
                    slopeok = False
                    break
            #print(slopeok)
            
        if all([param > 0.0 for param in prespos]) and (slopeok == True):
            #print 'Everything OK'
            return plike+addlike
        #print 'In LnPrior, LogLike set to infinity ', outalphas[-1]
        #import pdb;pdb.set_trace()
        return -np.inf

##################################################################################    
    def lnprob(pos):
        if fprior(pos):
            likel,outalphas,ycyl,mhse,miso,y2500,h2500 = lnlike(pos)
        else:
            likel = -np.inf
        #likel,outalphas,ycyl,mhse,miso = lnlike(pos)
        #yarr = np.array(ycyl)
        if not np.isfinite(likel):
            #print 'In LnProb, LogLike set to infinity'
            byc = [-1.0 for cluster in hk.cluster]
            bmh = [-1.0 for cluster in hk.cluster]
            bmi = [-1.0 for cluster in hk.cluster]
            return -np.inf, [byc,bmh,bmi,bmi,bmi]
            #return -np.inf, [-1.0 for ybad in ycyl]
        lp = lnprior(pos,outalphas,ycyl)
        #import pdb;pdb.set_trace()
        if not np.isfinite(lp):
            #print 'In LnProb, LogLike set to infinity'
            return -np.inf , [[-1.0 for ybad in ycyl] for i in range(5)]
        return lp + likel , [ycyl,mhse,miso,y2500,h2500]

    #import pdb;pdb.set_trace()
    logprob,(aaa,bbb,ccc,ddd,eee) = lnprob(efv.solns[:,0])
    
    
    for ii,mbc in enumerate(myset):
        wcm           = (mbc == mybcounts)
        psoln[wcm,:] /= efv.Pdl2y[ii].value

    efv.psolns  = psoln  #.to(efv.punits)
    #import pdb;pdb.set_trace()
    #efv.psolns = (efv.solns/(efv.Pdl2y)).to(efv.punits)
    #psolns = psolns.to("keV cm**-3")
    efv.values = efv.psolns.value[:,0]
    efv.errors = [efv.psolns.value[:,2],efv.psolns.value[:,1]] # Lower, then upper errors
    np.set_printoptions(precision=4)
    hdu = rwrf.saveimg(dv,hk,efv,ifp,component='Bulk')
#    rwrf.saveimg(dv,hk,efv,component='Residual')
    my_yint = post_check(sampler,t_mcmc,efv,hk,dv,ifp)
    ### Now including the E**(-2/3) term that is needed for cosmological relevancy
    damult = [(cluster.hofz**(-1./3)* cluster.d_a.value/1000.0)**2 for cluster in hk.cluster]
    ytotint = [y*d for y,d in zip(my_yint,damult)] # I don't know how I would actually distinguish multiple yints...
    ycyl    = ytotint         # Now in Mpc**-2 (should be a value of ~10 e-5)

    print("Unitless Pressure values from MCMC are:", efv.solns[:,0])
    print("Actual Pressure values from MCMC are:", efv.values)
    print("##################################################################")
    print("Integrated Y values are: ", ycyl)
    print("##################################################################")
    txtout.write(str(("Unitless Pressure values from MCMC are:", efv.solns[:,0])))
    txtout.write(str(("Actual Pressure values from MCMC are:", efv.values)))
    txtout.write("##################################################################")
    txtout.write(str(("Integrated Y values are: ", ycyl)))
    txtout.write("##################################################################")
    txtout.close()

    #print('Hi')
    tstr = hk.cfp.testmode+'_Run_'+efv.tag
    #savpik = os.path.join(hk.hk_outs.newpath,hk.hk_outs.nnfolder+'_pickle.sav')
    savpik = os.path.join(hk.hk_outs.newpath,tstr+'_pickle.sav')
    pfile = open(savpik,'wb')
    #print('Welcome')
    #import pdb;pdb.set_trace()
    #mstring = pickle.dumps(dv)
    #pfile.write(mstring)
    pickle.dump(dv,pfile,pickle.HIGHEST_PROTOCOL)
    pickle.dump(hk,pfile,pickle.HIGHEST_PROTOCOL)
    pickle.dump(efv,pfile,pickle.HIGHEST_PROTOCOL)
    #pfile.write(pickle.dumps(dv))
    #pfile.write(pickle.dumps(hk))
    #pfile.write(pickle.dumps(efv))
#    pfile.write(pickle.dumps(sampler))
    pfile.close
    npysavefile = os.path.join(hk.hk_outs.newpath,hk.hk_outs.nnfolder+'_solutions.npy')
    np.save(npysavefile, efv.solns)
    #print('Bye')
    npysavefile = os.path.join(hk.hk_outs.newpath,hk.hk_outs.nnfolder+'_logprob.npy')
    np.save(npysavefile, logprob)

    return hdu
    
#############################################################################################

def post_check(sampler,t_mcmc,efv,hk,dv,ifp):

    pos = efv.solns[:,0]
    #import pdb;pdb.set_trace()
    mnlvlmaps,ptsrcmaps,maps,yint,outalphas,mhse,miso,y2500,h2500 = make_skymodel_maps(pos,hk,dv,ifp,efv)

    return yint
    
def unpickle(file=None):

    #print 'THIS SECTION IS OLD! DO NOT USE!'
    #import pdb;pdb.set_trace()
    if file is None:
        sdir =myhome+'/Results_Python/MUSTANG2/a2146/'
        file=sdir+'2707_MUSTANG2_6_B_Real_2500S_500B_ML-NO_PP-NO_POWER_30W_pickle.sav'
    rdfl = open(file,'r')
    my_dv = pickle.load(rdfl)
    my_hk = pickle.load(rdfl)
    my_efv = pickle.load(rdfl)
    #    my_sampler = pickle.load(rdfl) ... Nope :(
    rdfl.close

    return my_dv,my_hk,my_efv

def get_initial_guess(pbins,uless_r,hk,conv=1.0):

    vals=pbins
    if hk.hk_ins[0].name == 'abell_2146':
        vals[0] /= 1.0  # Just spit-balling
        vals[1] /= 1.0  # Just spit-balling
    
    if hk.cfp.ndim > hk.cfp.bins:
        if hk.cfp.mn_lvl == True:
            vals=np.append(vals,0.01)
        if len(hk.cfp.ptsrc) > 0:      # I need to do this better (Feb 2018)
            vals=np.append(vals,0.01)
        if hk.cfp.blob == True:
            vals=np.append(vals,0.01)

#    hk.p_init=vals #??? Does it pass this back? --> Yes; yes it does.
            
    return vals

#def get_best_comp_maps(efv,hk,dv,mycomponent,hdu,returnwcs=False):
def get_best_comp_maps(efv,hk,dv,myinst,mycomponent,hdu,returnwcs=False):

    tstr = 'Test_Run_'
    tstr = hk.cfp.testmode+'_Run_'
    tstr = tstr+efv.tag

    modelsky={}
    fbase=tstr+myinst+"_"
    filename=tstr+myinst+"_Models_and_Residuals.fits"
    fullpath = os.path.join(hk.hk_outs.newpath,hk.hk_outs.prefilename+filename)
    hdulist = fits.open(fullpath)
    ########################### OVERRIDE 06 Feb 2018
    myhdu = hdu[myinst]
    hdulist = fits.HDUList(myhdu)

    #import pdb;pdb.set_trace()
    
    if returnwcs == False:
        modelsky[myinst]=get_best_inst_comp_map(efv,hk,dv,mycomponent,myinst,hdulist,returnwcs=False)
        return modelsky
    else:
        modelsky[myinst],w=get_best_inst_comp_map(efv,hk,dv,mycomponent,myinst,hdulist,returnwcs=True)
        return modelsky,w
        
def get_best_inst_comp_map(efv,hk,dv,mycomponent,myinst,hdulist,returnwcs=False):

    from astropy import wcs                 # Because WCS is useful, all things considered.

    comp_map=None
    for hdu in hdulist:
        myext = hdu.header['EXTNAME']
        #print myext.lower()
        if myext.lower() == myinst.lower()+'_'+mycomponent.lower():
            comp_map=hdu.data
            w = wcs.WCS(hdu.header)

    if type(comp_map) == type(None):
        print('hi', myext.lower(), myinst.lower()+'_'+mycomponent.lower())
        import pdb;pdb.set_trace()

    if returnwcs == True:
        return comp_map, w
    else:
        return comp_map
    


def r500_from_y500(yinteg,h,d_a,rho_crit,h70,ySZ=True,ySph=False):
    """
    Provide h and d_a as scalars:

    h        - little h (the one that changes with z)
    d_a      - angular distance, in Mpc, but just a value (not a quantity)
    rho_crit - has units of density!!!

    """

    Jofx    = 0.7398
    YMscale = 1.78
    if ySph:
        #ySZ  = False
        Jofx = 0.6145  # Actually I(x) in A10, but, it plays the same role, so use this variable

    ### So, when I do my initial Y500 calculation, it's just Y_SZ (without D_A**2, h(z)**-2/3).
    ### In this case, "yinteg" is Y_SZ, so set ySZ = True. However, when I come from the end of
    ### the MCMC run, then I've included these factors and thus ySZ = False.
    if ySZ == True:
        iv      = h**(-1./3)*d_a
        lside   = yinteg*iv**2
    else:
        lside   = yinteg
        
    Bofx    = 2.925e-5 * Jofx * h70**(-1)
    smu     = (lside/Bofx)**(1.0/YMscale)
    M500_i  = smu * 3e14 * u.Msun / h70   # In solar masses
    R500_i  = (3 * M500_i/(4 * np.pi  * 500.0 * rho_crit))**(1/3.)
    Mpc     = R500_i.decompose()
    Mpc     = Mpc.to('Mpc')
    r500    = (Mpc.value / d_a)

    return r500

def p500_m500_from_r500(r500,E,h_70,dens_crit):

    """ 
    Here, r500 and dens_crit need to have units!
    """

    M_500 = (4.0 * np.pi  * 500 * dens_crit)* r500**3 / 3.0
    m500  = M_500.decompose()
    m500  = m500.to('Msun')
    P_500  = (1.65 * 10**-3) *((E)**(8./3)) *((
            M_500 * h_70 )/((3*10**14)  * u.Msun)
            )**(2./3.)*(h_70)**2  *u.keV /u.cm**3
    p500 = P_500.to('keV cm**-3')

    return p500,m500


def Find_root(r_low,r_high,nsamp,rProf,alpha,norm,h70,h,d_a,rho_crit,r_thresh=1.0e-6):

    powspan     = r_high/r_low
    powers      = np.arange(nsamp) * np.log10(powspan) / (nsamp-1)
    radmax      = r_low*10**powers
    mydiff      = [Find_diff(r_in,rProf,alpha,norm,h70,h,d_a,rho_crit) for r_in in radmax]
    signs       = mydiff / mydiff[-1]
    negi        = np.where(signs < 0)[0]  # Interesting...it returns a tuple, which could have an empty array in it.
    #print 'Rocking this'
    if len(negi) == 0:
        print('Bad Compton y profile. Making bad assumptions now.')
        root  = np.sqrt(r_low*r_high)
        return root
    else:
        #print negi
        aind = np.max(negi)
        pt_a = radmax[aind]
        pt_b = radmax[aind+1] # By definition...
        if abs(mydiff[aind]) < r_thresh:
            root = pt_a
            return root
        elif abs(mydiff[aind+1]) < r_thresh:
            root = pt_b
            return root
        else:
            root = Find_root(pt_a,pt_b,nsamp,rProf,alpha,norm,h70,h,d_a,rho_crit,r_thresh=1.0e-6)
            return root

def Find_root_v2(r_low,r_high,y_low,y_high,rProf,alpha,norm,h70,h,d_a,rho_crit,r_thresh=1.0e-6):

    slope       = (y_high - y_low)/(r_high - r_low)
    delta       = -y_low / slope
    r_in        = r_low + delta
    mydiff      = Find_diff(r_in,rProf,alpha,norm,h70,h,d_a,rho_crit)
    #print 'Rocking this'
    if abs(mydiff) < r_thresh:
        root = r_in
        return root
    else:
        if mydiff < 0:
            r_low = r_in
            y_low = mydiff
        else:
            r_high= r_in
            y_high = mydiff
        root = Find_root_v2(r_low,r_high,y_low,y_high,rProf,alpha,norm,h70,h,d_a,rho_crit,r_thresh=1.0e-6)
        return root

def Find_diff(r_in,rProf,alpha,norm,h70,h,d_a,rho_crit):
    
    yinteg  = es.ycyl_integral(rProf,norm,alpha,r_in)
    r500    = r500_from_y500(yinteg,h,d_a,rho_crit,h70)

    return r_in - r500

def Y_SZ_v2(yProf,rProf,r_max,cluster,h70,ylist,geom,debug=False,yR2500=None,ySph=False,retcurvs=False):
    """
    yProf        - an integrated profile, i.e. in Compton y that matches theta_range
    rProf        - the radial profile on the sky (in radians)
    r_max        - the maximum (e.g. R500)
    cluster      - a class containing redshift/cosmological information about the cluster.
    h70          - H_0 scaled to 70 km/s/Mpc.
    ylist        - reference list for Y500 (at all radii).
    nsamp        - number of trial points
    
    I'm adopting equations 25-27 in Arnaud+ 2010, which makes use of Y_SZ, or Y_cyl and the
    Universal Pressure Profile (UPP). I tried to find just a straight empirical Y_cyl(R500)-M_500,
    but that doesn't seem to exist?!?

    (hofz,ang_dist,dens_crit) are defined within (cluster)

    """

    hofz      = cluster.hofz             # 
    ang_dist  = cluster.d_a.value/1000.0 # kpc
    dens_crit = cluster.dens_crit        # critical density. No units??
                                     #r_low       = r_max * 10**(-0.8)
    #r_high      = r_max * 10**(0.8)
    alpha,norm  = es.ycyl_prep(yProf,rProf)
    #if not (yR2500 is None):
    #    retcurvs=True
        
    #import pdb;pdb.set_trace()
    #yinteg, root = ycyl_simul_r500(rProf,yProf,alpha,r_max,h70,hofz,ang_dist,dens_crit,r_thresh=3.0e-2)
    if not ySph:
        yinteg, root,yR,yM = ycyl_simul_v2(rProf,yProf,alpha,r_max,h70,hofz,ang_dist,dens_crit,ylist,geom,r_thresh=3e-2,retcurvs=True)
        #import pdb;pdb.set_trace()
    else:
        ### yProf is ***NOT*** a Compton y profile!!! IT IS A PRESSURE PROFILE! -- May 2019
        yinteg, root,yR,yM = ysph_simul(ylist,rProf,yProf,alpha,geom,cluster,h70,ythresh=3e-8,retcurvs=True)

    if not (yR2500 is None):
        y2500,r2500 = my_root_finder(rProf,yM,yR2500[:-1].value,rProf[51],yR2500[50].value,ythresh=1.0e16)
        #import pdb;pdb.set_trace()
    else:
        y2500,r2500 = 0,0

    #import pdb;pdb.set_trace()

    if retcurvs:
        return yinteg, root,y2500,r2500,yR,yM
    else:
        return yinteg, root,y2500,r2500

    
def ycyl_simul_v2(rads,yProf,alpha,maxrad,h70,hofz,d_a,rho_crit,ylist,geom,r_thresh=3e-2,retcurvs=False):

    #rads,yProf,alpha):
    """
    Remember, theta_range is in radians
    """

    ### In accordance with how Arnaud+ 2010 defines these terms...
    #h70      = (cosmo.H(0) / (70.0*u.km / u.s / u.Mpc))
    #rho_crit = cosmo.critical_density(map_vars["z"])
    #hofz     = cosmo.H(map_vars["z"])/cosmo.H(0)                    

    fgeo          = geom[3]*geom[4] # Scale by ellipsoidal radii scalings.
    Ycyl          = 0
    if alpha[0]  <= -2: alpha[0]=-1.9
    badalp        = (alpha == -2)
    alpha[badalp] = -2.01 # Va fanculo.
    rolledrad     = np.roll(rads,-1)
    intupper      = rolledrad**2 * (rolledrad/rads)**(alpha) #* myrads
    intlower      = rads**2
    intlower[0]   = 0.0
    integrand     = intupper - intlower
    Yshell        = 2.0*np.pi*yProf[:-1]*integrand[:-1]/(alpha[:-1]+2.0)
    Ycyl          = np.cumsum(Yshell)*fgeo
    #import pdb;pdb.set_trace()
    Yref          = ylist[:-1].value
    #Yref          = ylist[:-1]
    mydiff        = Yref - Ycyl

    nrbins        = len(rads)
    #mydiff = rads[1:] - r500
    #absdif = np.abs(mydiff)
    posdiffs = (mydiff > 0)
    turnover = mydiff[posdiffs]
    #import pdb;pdb.set_trace()
    bestr = -rads[1]
    bestY = -Ycyl[0]
    bisca = 0
    if len(turnover) > 1:
        besti  = np.where(mydiff == np.min(turnover))
        bisca  = np.asscalar(besti[0])
    if bisca < nrbins -3 and bisca > 10: 
        myinds = bisca + np.asarray([-2,-1,0,1,2],dtype='int')
        #myinds = np.intersect1d(naind,
        myrs   = rads[myinds+1]
        myYs   = Ycyl[myinds]
        myds   = mydiff[myinds]
        myp2   = np.polyfit(myrs,myds,2)
        myY2   = np.polyfit(myrs,myYs,2)
        myroot = np.roots(myp2)
        rdiff  = np.abs(myroot - myrs[2])

        bestr  = myroot[0] if rdiff[0] < rdiff[1] else myroot[1]
        Y2fxn  = np.poly1d(myY2)
        bestY  = Y2fxn(bestr)

        
    #import pdb;pdb.set_trace()
    #bestr  = r500[besti]
    #bestY  = Ycyl[besti]
    if retcurvs:
        return bestY,bestr, Yref,Ycyl
    else:
        return bestY,bestr
   
def ysph_simul(ylist,rads,pProf,alpha,geom,cluster,h70,ythresh=3e-8,retcurvs=False):
    """
    To be developped.
    """

    ### In accordance with how Arnaud+ 2010 defines these terms...
    #h70      = (cosmo.H(0) / (70.0*u.km / u.s / u.Mpc))
    #rho_crit = cosmo.critical_density(map_vars["z"])
    #hofz     = cosmo.H(map_vars["z"])/cosmo.H(0)                    
    h   = cluster.hofz
    d_a = cluster.d_a.to('Mpc')
    rho_crit = cluster.dens_crit
    E   = cluster.hofz
    #h70 = cluster.h70

    fgeo          = geom[3]*geom[4]*geom[5] # Scale by ellipsoidal radii scalings.
    Ysph          = 0
    if alpha[0]  <= -3: alpha[0]=-2.9
    badalp        = (alpha == -3)
    alpha[badalp] = -3.01 # Va fanculo.
    rolledrad     = np.roll(rads,-1)
    intupper      = rolledrad**3 * (rolledrad/rads)**(alpha) #* myrads
    intlower      = rads**3
    intlower[0]   = 0.0
    integrand     = intupper - intlower
    Yshell        = 4.0*np.pi*pProf[:-1]*integrand[:-1]/(alpha[:-1]+3.0) 
    Ysph          = np.cumsum(Yshell) *fgeo
    Yref          = ylist[:-1].value
    #Yref          = ylist[:-1]
    mydiff        = Yref - Ysph

    ### Look here (April 13, 2019)
    #import matplotlib.pyplot as plt
    #plt.plot(rads[1:]*3600.0*180/np.pi,Yref)
    #plt.plot(rads[1:]*3600.0*180/np.pi,Ysph)
    #plt.yscale('log');plt.xscale('log')
    #plt.show()
    #import pdb;pdb.set_trace()

    #mydiff = rads[1:] - r500
    #absdif = np.abs(mydiff)
    posdiffs = (mydiff > 0)
    turnover = mydiff[posdiffs]
    #import pdb;pdb.set_trace()
    bestr = rads[51]
    bestY = Ysph[50]
    bisca = 0
    nrbins        = len(rads)
    if len(turnover) > 1:
        besti  = np.where(mydiff == np.min(turnover))
        bisca  = np.asscalar(besti[0])
    if bisca < nrbins-3 and bisca > 10: 
        myinds = bisca + np.asarray([-2,-1,0,1,2],dtype='int')
        #myinds = np.intersect1d(naind,
        myrs   = rads[myinds+1]
        myYs   = Ysph[myinds]
        myds   = mydiff[myinds]
        myp2   = np.polyfit(myrs,myds,2)
        myY2   = np.polyfit(myrs,myYs,2)
        myroot = np.roots(myp2)
        rdiff  = np.abs(myroot - myrs[2])

        bestr  = myroot[0] if rdiff[0] < rdiff[1] else myroot[1]
        Y2fxn  = np.poly1d(myY2)
        bestY  = Y2fxn(bestr)

        if bestY > ythresh:
            print(bestr, bestY, np.max(rads))
            stupid = np.random.normal(0,1)
            if stupid > 5: import pdb;pdb.set_trace()
            
        
    #import pdb;pdb.set_trace()
    #bestr  = r500[besti]
    #bestY  = Ysph[besti]
    if retcurvs:
        return bestY,bestr, Yref,Ysph
    else:
        return bestY,bestr

def yMP500_from_r500(radian500,cluster,ySZ=True,ySph=True,delta=500,YMrel='A10'):

    h   = cluster.hofz
    d_a = cluster.d_a.to('Mpc')
    rho_crit = cluster.dens_crit
    E   = cluster.E
    h70 = cluster.h70

    m500    = 4*np.pi*delta*rho_crit/3 *(radian500*d_a)**3
    m500    = m500.to('M_sun')
    Y500    = y_delta_from_mdelta(m500.value,cluster,delta=delta,ySph=ySph,YMrel=YMrel)
    #print('And checking again... ',Y500)
    #import pdb;pdb.set_trace()
    
    P500 = (1.65 * 10**-3) * ((E)**(8./3)) * ((
        m500 * h70)/ ((3*10**14) * u.M_sun)
        )**(2./3+0.11) * h70**2 * u.keV / u.cm**3

    return Y500, m500, P500

def rMP500_from_y500(yinteg,cluster,ySZ=True,ySph=False,delta=500,YMrel='A10'):
    """
    Provide h and d_a as scalars:

    h        - little h (the one that changes with z)
    d_a      - angular distance, in Mpc, but just a value (not a quantity)
    rho_crit - has units of density!!!

    """
    h   = cluster.hofz
    d_a = cluster.d_a.to('Mpc')
    rho_crit = cluster.dens_crit
    E   = cluster.E
    h70 = cluster.h70
    #print("In rMP500",d_a,E,h,yinteg)

    M500_i  = m_delta_from_ydelta(yinteg,cluster,delta=delta,ySph=ySph,YMrel=YMrel)*u.M_sun
    #M500_i  = smu * 3e14 * u.Msun / h70   # In solar masses
    R500_i  = (3 * M500_i/(4 * np.pi  * delta * rho_crit))**(1/3.)
    Mpc     = R500_i.decompose()
    Mpc     = Mpc.to('Mpc')
    r500    = (Mpc / d_a).value

    P500 = (1.65 * 10**-3) * ((h)**(8./3)) * ((
        M500_i * h70)/ ((3*10**14) * u.M_sun)
        )**(2./3 +0.11) * h70**2 * u.keV / u.cm**3

    if ySZ == True:
        iv      = h**(-1./3)*d_a
        lside   = iv**2
    else:
        lside   = 1.0

    logy  = np.log10(lside.value*yinteg)
    msys = get_YM_sys_err(logy,YMrel,delta=500,ySph=ySph,h70=h70)
  
    return r500, M500_i, P500, msys

def hse_with_density_v2(pprof,alphas,rads,edens,Tx,m500v,pause=False):

    mu    = 0.61
    Gmp   = 1.116e-37 # m**3 s**-2
    myrho = 5.5427366e15 # M_sun / cm**3
    Mtot = -pprof*alphas * rads * myrho / (edens * mu)
    Miso = -Tx *alphas * rads * myrho / mu
    m2500v = m500v*5
    
    m500,r500 = my_root_finder(rads,Mtot,m500v,rads[51],m500v[50],ythresh=3.0e16,minty=65,pause=pause)
    miso,riso = my_root_finder(rads,Miso,m500v,rads[51],m500v[50],ythresh=3.0e16,minty=65,pause=pause)
    m2500,r2500 = my_root_finder(rads,Mtot,m2500v,rads[41],m2500v[40],ythresh=3.0e16,minty=63,pause=pause)

    if pause:
        import pdb;pdb.set_trace()

    
    #if retcurvs:
    return m500,r500,miso,riso,m2500,r2500,Mtot,Miso,m500v
    #else:
    #    return m500,r500,miso,riso

def my_root_finder(xarr,y1,y2,bestx,besty,ythresh=3e-8,minty=63,scale=1e14,pause=False):

    my1      = y1[minty:]/scale
    my2      = y2[minty:]/scale
    myx      = xarr[minty:]
    mydiff   = (my1 - my2)   # y1 = "reference array"
    if mydiff[0] > 0:
        posdiffs = (mydiff < 0)       # Correct condition is
        #print('hey')
    else:
        posdiffs = (mydiff > 0)       # Correct condition is
        #print('yo')

    turnover = mydiff[posdiffs]
    #print('Length of turnover is ',len(turnover),map_vars["nrbins"])

    bisca = 0
    if len(turnover) > 1:
        #besti  = np.where(mydiff == np.min(turnover))
        besti  = np.where(mydiff == turnover[0])
        bisca  = np.asscalar(besti[0])
        #print(bisca,len(myx))
    if bisca < len(myx)-3 and bisca > 1: 
        myinds = bisca + np.asarray([-2,-1,0,1,2],dtype='int')
        #myinds = np.intersect1d(naind,
        myrs   = myx[myinds+1]
        myYs   = my1[myinds]
        myds   = mydiff[myinds]
        myp2   = np.polyfit(myrs,myds,2)
        myY2   = np.polyfit(myrs,myYs,2)
        #if any(np.isnan(myp2)) or any(np.isinf(myp2)):
        #print(myinds,myrs,myds)
        #print('######################################')
        #print(myYs,my2[myinds])
        #print('--------------------------',myp2,myY2)

        myroot = np.roots(myp2)
        #print(myroot)
        rdiff  = np.abs(myroot - myrs[2])
        #print(rdiff)

        bestx  = myroot[0] if rdiff[0] < rdiff[1] else myroot[1]
        #print('########################################################')
        #print(myroot,rdiff,bestx)
        Y2fxn  = np.poly1d(myY2)
        besty  = Y2fxn(bestx)*scale
        #print(besty/1e14,bestx)
        #print(Y2fxn(myrs))
        
        if besty > ythresh:
            print('besty beyond threshold: ',bestx, besty, np.max(myx))
            stupid = np.random.normal(0,1)
            if stupid > 5: import pdb;pdb.set_trace()

    if pause:
        import pdb;pdb.set_trace()

    return besty,bestx

###################################################################################################################

###################################################################################################################

###################################################################################################################

def virial_fits_v2(efv,hk,doR2500=True):
    """

    Need to  replace map_vars['thetas'] with efv.thetas
    and map_vars['d_ang'] with...?   mycluster.R_500/mycluster.d_a

    yprof is in units of rads **-2 ??
    pprof is in units of keV cm**-3
    """

    solns=efv.solns
    model=hk.cfp.model
    mdists=[]; popts=[]; pcovs=[];m2dists=[]
    posind=1
    print('You are in virial fits')
    for mycluster,Pdl2y,bins,m500v in zip(hk.cluster,efv.Pdl2y,hk.cfp.bulkarc,efv.m500s):
        d_ang = mycluster.d_a
        mysolns = solns[posind:posind+len(bins)]
        posind+=len(bins)
        kpcprof = (d_ang * efv.thetas).to('kpc').value    
        if model == 'NP':
            bprof,sprof,Yprof,SYprof = np_pprof_yprof_werrs(mysolns,Pdl2y,bins)
            rprof = (d_ang * bins).to('kpc').value
        elif model == 'GNFW':
            bprof,sprof,Yprof,SYprof = gnfw_pprof_yprof_werrs(mysolns,efv,hk,mycluster,Pdl2y)
            rprof = kpcprof[:-1]
        elif model == 'BETA':
            bprof,sprof,Yprof,SYprof = beta_pprof_yprof_werrs(mysolns,efv,hk,mycluster,Pdl2y)
            rprof = kpcprof[:-1]
        else:
            print('You made a mistake. Look at your life. Look at your choices.')
        
        import scipy.optimize as opt    
        mu_e  = 1.17       # mean electron density (Mroczkowski+ 2011)
        mu    = 0.61
        m_e   = 9.11e-31   # electron mass (kg)
        f_gas = 0.13       # gas fraction
        TXS   = 6.65e-29   # Thomson cross section (m**2)
        G     = 6.67e-11   # Gravitational constant (m**3 kg**-1 s**-2)
        c     = 3.00e8     # speed of light (m s**-1)

        ls = 3.0 * m_e * c**2 / TXS
        mypprof = bprof * u.keV.to('kg m**2 s**-2') * 1.0e6 # J / m**3
        mySpprof = sprof * u.keV.to('kg m**2 s**-2') * 1.0e6 # J / m**3
        els = - 4.0 * np.pi * rprof**3 * mypprof * u.kpc.to('m')
        Sigels = - 4.0 * np.pi * rprof**3 * mySpprof * u.kpc.to('m')
        #cf  = (1 + 1/mu_e) / (16 * np.pi**2 * G * f_gas)
        cf  = (mu_e / mu) / (16 * np.pi**2 * G * f_gas)

        #h = cosmo.H(mycluster.z)/cosmo.H(0)
        h  = mycluster.hofz
        myyprof = Yprof * (d_ang.to('kpc').value)**2
        mySyprof = SYprof * (d_ang.to('kpc').value)**2
        #import pdb;pdb.set_trace()
        fls = cf * (ls*myyprof+els) #
        Sigfls = cf * (ls*mySyprof + Sigels)
        gi1 = (fls > 0)       # 
        gi2 = (rprof > 100.0) # Cool core is usually ~100 kpc
        gi  = [gii1 and gii2 for gii1,gii2 in zip(gi1,gi2)]
    
        initial_guess = [0.01,1/500.0] # (m_p) /cm**3 , kpc
        mybounds = (np.array([0.0,0.0]),np.array([np.inf,np.inf]))
        #import pdb;pdb.set_trace()

        try:
            popt2, pcov2 = opt.curve_fit(NFW_int_v3, rprof[gi], fls[gi], p0=initial_guess,sigma=Sigfls[gi],bounds=mybounds)
            print(popt2,pcov2)
        except:
            import pdb;pdb.set_trace()
            popt2 = initial_guess
            pcov2 = np.array([[0.005**2,-1e-6],[-1e-6,2e-6]])
        mymasses = Mtot_NFW(kpcprof, popt2[0], 1.0/popt2[1])
        m500,r500 = my_root_finder(efv.thetas,mymasses,m500v.value,efv.thetas[51],m500v.value[50],ythresh=1.0e16)
        lowmasses = Mtot_NFW(kpcprof, popt2[0]*np.sqrt(0.9), 1.0/popt2[1])
        highmasses = Mtot_NFW(kpcprof, popt2[0]/np.sqrt(0.9), 1.0/popt2[1])
        m500l,r500l = my_root_finder(efv.thetas,lowmasses,m500v.value,efv.thetas[51],m500v.value[50],ythresh=1.0e16)
        m500h,r500h = my_root_finder(efv.thetas,highmasses,m500v.value,efv.thetas[51],m500v.value[50],ythresh=1.0e16)
        mcunc = [m500-m500l,m500h-m500]

        #if doR2500:
        m2500,r2500   = my_root_finder(efv.thetas,mymasses , 5*m500v.value,efv.thetas[51],5*m500v.value[50],ythresh=1.0e16)
        m2500l,r2500l = my_root_finder(efv.thetas,lowmasses, 5*m500v.value,efv.thetas[51],5*m500v.value[50],ythresh=1.0e16)
        m2500h,r2500h = my_root_finder(efv.thetas,highmasses,5*m500v.value,efv.thetas[51],5*m500v.value[50],ythresh=1.0e16)
        m2cunc = [m2500-m2500l,m2500h-m2500]
 
        from scipy import linalg as sla
        #cho = sla.cholesky(pcov,lower=True)
        #rotmat = np.transpose(cho)
        cho2 = sla.cholesky(pcov2,lower=True)
        rotmat2 = np.transpose(cho2)

        nstep=1000
        arrm500 = np.zeros(nstep);        arrr500  = np.zeros(nstep)
        arrm2500 = np.zeros(nstep);       arrr2500 = np.zeros(nstep)
        for i in range(nstep):
            mypars = np.array([-1.0,-1.0])
            badpars = any(mypars <= 0)
            iters=0
            while badpars:
                iters  += 1
                norstep = np.random.normal(0,1,2)
                #my_step = np.matmul(rotmat,norstep)
                #my_step = np.matmul(norstep,rotmat)
                my_step = np.matmul(norstep,rotmat2)
                mypars  = popt2 + my_step
                badpars = any(mypars <= 0) and (iters < 100)
            #print(iters)
            mymasses = Mtot_NFW(kpcprof, mypars[0], 1.0/mypars[1])
            arrm500[i],arrr500[i] = my_root_finder(efv.thetas,mymasses,m500v.value,efv.thetas[51],m500v.value[50],ythresh=1.0e16)
            arrm2500[i],arrr2500[i] = my_root_finder(efv.thetas,mymasses,5*m500v.value,efv.thetas[51],5*m500v.value[50],ythresh=1.0e16)

        mperc  = np.percentile(arrm500, [16, 50, 84],axis=0)
        m2perc  = np.percentile(arrm2500, [16, 50, 84],axis=0)
        #mdist  = np.array([mperc[1],mperc[2]-mperc[1],mperc[1]-mperc[0]])
        mdist  = np.array([m500,mperc[2]-mperc[1],mperc[1]-mperc[0]])
        m2dist  = np.array([m2500,m2perc[2]-m2perc[1],m2perc[1]-m2perc[0]])
        print(mdist)
        mdists.append(mdist); popts.append(popt2); pcovs.append(pcov2); m2dists.append(m2dist)

    #import pdb;pdb.set_trace()

    if doR2500:
        return mdists,popts,pcovs,mcunc,m2dists,m2cunc
    else:
        return mdists,popts,pcovs,mcunc
    
def NFW_int_v2(rads,rho_not,Rs):

    x = (1 + rads/Rs)

    integral = (x - 1 - np.log(x))*Rs/x
    parens   = integral - Rs/(2.*(1.+Rs/rads)**2) 
    rs       = (rho_not * Rs**2)**2 * parens

    frs      = rs * 10**12 * (u.kpc.to('m') * u.M_p.to('kg'))**2 # I had 10**6 before
    frs     *= u.kpc.to('m')
    
    return frs

def NFW_int_v3(rads,rho_not,invRs):

    Rs = 1/invRs
    x = (1 + rads/Rs)

    integral = (x - 1 - np.log(x))*Rs/x
    parens   = integral - Rs/(2.*(1.+Rs/rads)**2) 
    rs       = (rho_not * Rs**2)**2 * parens

    frs      = rs * 10**12 * (u.kpc.to('m') * u.M_p.to('kg'))**2 # I had 10**6 before
    frs     *= u.kpc.to('m')
    
    return frs

def NFW_int(rads,rho_not,Rs):

    myvals = np.log(1 + rads/Rs) / (1 + rads/Rs)**2
    alpha, pnorm = es.ycyl_prep(myvals,rads)
    #print(pnorm[50:60])
    #print(myvals[50:60])
    #import pdb;pdb.set_trace()
    #print(alpha[0],alpha[49],alpha[99],alpha[149])
    
    rs         = 0
    if alpha[0]  <= -1: alpha[0]=-0.9
    badalp        = (alpha == -1)
    if any(badalp):
        import pdb;pdb.set_trace()
    alpha[badalp] = -1.01 # Va fanculo.
    rolledrad     = np.roll(rads,-1)
    intupper      = rolledrad * (rolledrad/rads)**(alpha) #* myrads
    intlower      = rads * 1.0
    intlower[0]   = 0.0
    integrand     = intupper - intlower
    ginteg        = myvals[:-1] * integrand[:-1]/(alpha[:-1]+1.0) 
    parens        = np.cumsum(ginteg) - Rs/(2.*(1.+Rs/rads[:-1])**2) 
    rs            = (rho_not * Rs**2)**2 * parens

    frs           = rs * 10**12 * (u.kpc.to('m') * u.M_p.to('kg'))**2 # I had 10**6 before
    frs          *= u.kpc.to('m')
    
    return frs

def Mtot_NFW(rads,rho_not,Rs):

    myvals  = np.log(1 + rads/Rs) - (1 + Rs/rads)**(-1)
    mscale  = 4*np.pi*rho_not * Rs**3 * u.M_p.to('kg')
    mts     = mscale * (u.kpc.to('cm'))**3 * u.kg.to('M_sun')
    mtot    = mts * myvals
    
    return mtot

def plfit(y,sy,x,xp=0.0):

    lny  = np.log(y)
    slny = sy/y
    lnx  = np.log(x)

    N    = len(y)
    #Del  = N * np.sum(lnx**2) - (np.sum(lnx)**2)
    #A    = (np.sum(lnx**2)*np.sum(lny) - np.sum(lnx)*np.sum(lnx*lny))/Del
    #B    = (N*np.sum(lnx*lny) - np.sum(lnx)*np.sum(lny))/Del
    #sA   = slny * np.sqrt(np.sum(lnx**2)/Del)
    #sB   = slny * np.sqrt(N/Del)

    wts = 1.0/slny**2
    w   = np.sum(wts)
    wxx = np.sum(wts*lnx**2)
    wy  = np.sum(wts*lny)
    wx  = np.sum(wts*lnx)
    wxy = np.sum(wts*lnx*lny)
    Del = w*wxx - wx**2
    A   = (wxx*wy - wx*wxy)/Del
    B   = (w*wxy - wx*wy)/Del
    sA  = np.sqrt( wxx / Del)
    sB  = np.sqrt( w   / Del)
    
    if xp > 0.0:
        lnxp = np.log(xp)

        lnS = A + B*lnxp
        S   = np.exp(lnS)
        sS  = (sA + sB*lnxp)*S

    return A,B, sA, sB

def np_pprof_yprof_werrs(solns,Pdl2y,bins):

    bprof = (solns[:,0]/Pdl2y).value
    sprof = (np.sqrt(solns[:,1]*solns[:,2])/Pdl2y).value

    ulprof  = solns[:,0]
    sulprof = np.sqrt(solns[:,1]*solns[:,2])
    
    Yprof  = np.zeros(len(ulprof))
    SYprof = np.zeros(len(ulprof))

    #import pdb;pdb.set_trace()

    for jj in range(len(ulprof)):
        ii = jj if jj == 0 else jj-1
        y  = ulprof[ii:ii+2];  sy = sulprof[ii:ii+2]; x = bins[ii:ii+2]
        A,B,sA,sB = plfit(y,sy,x,xp=0.0)
        if jj == 0:
            intupper      = x[1]**3  #* myrads
            intlower      = 0
            supper        = x[1]**3 / (B+3.0)
            slower        = 0
        else:
            intupper      = x[1]**3 * (x[1]/x[0])**(B)#* myrads
            intlower      = x[0]**3
            supper        = x[1]**3*(np.log(x[1]/x[0])*(x[1]/x[0])**(B) - 1/(B+3.0))
            slower        = x[0]**3 / (B+3.0)
                                     
            
        integrand     = intupper - intlower
        sintegrand    = (supper - slower)*sB
        Yshell        = 4.0*np.pi*ulprof[ii]*integrand/(B+3.0) 
        SaYshell      = 4.0*np.pi*ulprof[ii]*sintegrand/(B+3.0)
        SpYshell      = sulprof[ii]*Yshell/ulprof[ii]
        VarYshell     = SpYshell**2 + SaYshell**2
        Yprof[jj]     = Yshell
        SYprof[jj]    = VarYshell

    Yprof  = np.cumsum(Yprof)
    SYprof = np.sqrt(np.cumsum(SYprof))

    return bprof,sprof,Yprof,SYprof
    
def gnfw_pprof_yprof_werrs(solns,efv,hk,mycluster,Pdl2y,ySph=True,geom=[0,0,0,1,1,1,0,0],bsind=0):

    radii    = (efv.thetas* mycluster.d_a).to('kpc')
    myup     = np.array([1.177,8.403,5.4905,0.3081,1.0510])
    mypos    = solns[:,0]
    spos     = np.sqrt(solns[:,1]*solns[:,2])
    if len(mypos) < len(myup):
        #myup[1:1+len(pos)]=pos
        #myup[1:1+len(mypos)]=mypos
        myup[0:len(mypos)]=mypos
    else:
        myup = mypos    # If len(pos) > 5, that's OK...we won't use those!
           
    A10all = 1.0510 # Alpha value in A10, all clusters
    A10cc  = 1.2223 # Alpha value in A10, cool-core clusters
    A10dis = 1.4063 # Alpha value in A10, disturbed clusters

    R500 = mycluster.R_500.to('kpc')
    P500 = mycluster.P_500.to('keV cm**-3')
    r500 = (R500 / mycluster.d_a).decompose().value
        
    #pprof    = gdi.gnfw(R500, P500, radii, c500=myup[0], p=myup[1], a=myup[4], b=myup[2], c=myup[3])
    pprof    = cpp.gnfw(hk.av.mycosmo['h_70'], radii, P500, R500,
                                c500=myup[0], p=myup[1], a=myup[4], b=myup[2], c=myup[3])
    unitless_profile = (pprof * Pdl2y).decompose().value
    yProf = ni.int_profile(efv.thetas, unitless_profile,efv.thetas)
    outalphas = unitless_profile*0.0+2.0
    integrals = yProf
    #yint = ni.Ycyl_from_yProf(yProf,ppbins,r500)
    palphas, pnorm = es.ycyl_prep(unitless_profile,efv.thetas)

    #print(pprof,efv.thetas)
    #import pdb;pdb.set_trace()
    
    if efv.ySph:
        yint ,newr500,y2500,r2500,yR,yM=Y_SZ_v2(unitless_profile,efv.thetas,efv.r500[bsind],mycluster,hk.av.mycosmo['h_70'],efv.y500s[bsind],geom,
                                          ySph=efv.ySph,yR2500=efv.y2500s[bsind],retcurvs=True) # As of Aug. 31, 2018
    else:
        yint ,newr500,y2500,r2500,yR,yM=Y_SZ_v2(Int_Pres,efv.thetas,efv.r500[bsind],mycluster,hk.av.mycosmo['h_70'],efv.y500s[bsind],geom,
                                          ySph=efv.ySph,yR2500=efv.y2500s[bsind],retcurvs=True) # As of Aug. 31, 2018

    srad = efv.thetas/r500
    pf1  = -myup[3] / myup[0]
    pf2a = (myup[4]-myup[2])/(myup[0]*myup[3])
    pf2b = myup[4]*(myup[0]*srad)**myup[4] / (1 + (myup[0]*srad)**myup[4])
    sprof = np.sqrt( (pprof*spos[0]/myup[1])**2 + (pprof*spos[1]*(pf1+pf2a*pf2b))**2 )

    SYprof = yM * sprof[:-1]/pprof[:-1]
    
    return pprof[:-1].value,sprof[:-1].value,yM,SYprof

def beta_pprof_yprof_werrs(solns,efv,hk,mycluster,Pdl2y,ySph=True,geom=[0,0,0,1,1,1,0,0],atR2500=False,noerrs=False,bsind=0):

    #radii    = (efv.thetas* mycluster.d_a).to('kpc').value
    radii    = efv.thetas
    myup     = np.array([1.177,8.403,5.4905,0.3081,1.0510])
    mypos    = solns[:,0]
    spos     = np.sqrt(solns[:,1]*solns[:,2])
    R500 = mycluster.R_500.to('kpc')
    P500 = mycluster.P_500.to('keV cm**-3')
    r500 = (R500 / mycluster.d_a).decompose().value

    pprof     = mypos[0]*(1.0+(radii/mypos[1])**2)**(-1.5*mypos[2])*u.keV / (u.cm**3)        ### Beta model
    scaling  = scs.gamma(1.5*mypos[2]-0.5)/scs.gamma(1.5*mypos[2])*(mypos[1]/mycluster.d_a.to('kpc').value)
    scaling *= mypos[0] * np.sqrt(np.pi)
    unitless_profile = (pprof * Pdl2y).decompose().value
    yProf    = scaling * Pdl2y.value *  (1.0+(radii/mypos[1])**2)**(0.5-1.5*mypos[2])
    outalphas = unitless_profile*0.0+2.0

    dpdrc    = pprof * (3*mypos[2]*radii**2/mypos[1]**3) / (1.0+(radii/mypos[1])**2)
    dpdbeta  = pprof * -1.5 * np.log(1.0+(radii/mypos[1])**2)
    dpdpnot  = pprof / mypos[0]

    ###########################################################################################################
    if efv.ySph:
        yint ,newr500,y2500,r2500,yR,yM=Y_SZ_v2(unitless_profile,efv.thetas,efv.r500[bsind],mycluster,hk.av.mycosmo['h_70'],efv.y500s[bsind],geom,
                                          ySph=efv.ySph,yR2500=efv.y2500s[bsind],retcurvs=True) # As of Aug. 31, 2018
    else:
        yint ,newr500,y2500,r2500,yR,yM=Y_SZ_v2(Int_Pres,efv.thetas,efv.r500[bsind],mycluster,hk.av.mycosmo['h_70'],efv.y500s[bsind],geom,
                                          ySph=efv.ySph,yR2500=efv.y2500s[bsind],retcurvs=True) # As of Aug. 31, 2018

    sprof    = np.sqrt( (dpdpnot*spos[0])**2 + (dpdrc*spos[1])**2 + (dpdbeta*spos[2])**2 )
    SYprof = yM * sprof[:-1]/pprof[:-1]

    #import pdb;pdb.set_trace()
    if noerrs:
        return pprof.value
    else:
        print('Making bogus values!')
        return pprof[:-1].value,sprof[:-1].value,yM,SYprof.value


##########################################################################################################################
##########################################################################################################################


def m2r_delta(m_delta,cluster,delta=500):

    rho_crit = cluster.dens_crit
    d_a = cluster.d_a.to('Mpc')

    R_delta  = (3 * m_delta/(4 * np.pi  * delta * rho_crit))**(1/3.)
    Mpc     = R_delta.to('Mpc')
    r500    = (Mpc.value / d_a)

    return r500.value

def m_delta_from_ydelta(y_delta,cluster,delta=500,ySph=True,YMrel='A10'):
    """
    Basically just a repository of Y-M relations.
    
    """
    h   = cluster.hofz
    d_a = cluster.d_a.to('Mpc')
    rho_crit = cluster.dens_crit
    E   = cluster.E
    h70 = cluster.h70
    iv      = h**(-1./3)*d_a

    myYdelta = y_delta * (iv**2)

    AAA,BBB = get_AAA_BBB(YMrel,delta,ySph=ySph,h70=h70)

    m_delta = ( myYdelta.value / 10**BBB )**(1./AAA)

    return m_delta

def y_delta_from_mdelta(m_delta,cluster,delta=500,ySph=True,YMrel='A10'):
    """
    Finds A,B (scaling law terms, in get_AAA_BBB()) and applies them.
    """

    h   = cluster.hofz
    d_a = cluster.d_a.to('Mpc')
    rho_crit = cluster.dens_crit
    E   = cluster.E
    h70 = cluster.h70
    iv      = h**(-1./3)*d_a

    AAA,BBB = get_AAA_BBB(YMrel,delta,ySph=ySph,h70=h70)

    y_delta = m_delta**AAA * 10**BBB / (iv**2)

    return y_delta

def get_AAA_BBB(YMrel,delta,ySph=True,h70=1.0):
    """
    Basically just a repository of Y-M relations.
    YMrel must be either:
       (1) 'A10'
       (2) 'M12', or
       (3) 'P17'

    All are converted to Y = 10^A * M^B; mass (M) is in units of solar masses; Y is in Mpc^2 (i.e. with D_A^2 * E(z)^-2/3)
    
    """

    ycyl = not ySph
    if delta == 2500:
        if YMrel == 'A10':
            AAA    = 1.637;   BBB = -28.13  # From Comis+ 2011
        #if ycyl:
        #    AAA = 1.60;   BBB = -27.4   # From Comis+ 2011
        elif YMrel == 'M12':
            #BBB = -29.66909090909 # How did I get this wrong value?
            BBB = -30.66909090909 # recalculated... July 10, 2019
            AAA = 1.0 / 0.55
        elif YMrel == 'M12-SS':
            BBB = -28.501666666666667
            AAA = 5.0/3.0
        elif YMrel == 'P17':
            AAA = 1.755
            BBB = -29.6833076    # -4.585
        else:
            print('using Comis+ 2011 values')
            
    elif delta == 500:
        if YMrel == 'A10':
            AAA   = 1.78
            Jofx  = 0.7398 if ycyl else 0.6145  # Actually I(x) in A10, but, it plays the same role, so use this variable
            print(Jofx,' ycyl: ',ycyl)
            Bofx  = 2.925e-5 * Jofx * h70**(-1) / (3e14/h70)**AAA
            BBB = np.log10(Bofx)
        elif YMrel == 'M12':
            #BBB = -30.66909090909 # BBB = -16.567
            BBB = -37.65227272727
            AAA = 1.0 / 0.44
        elif YMrel == 'M12-SS':
            BBB = -28.735
            AAA = 5.0/3.0
        elif YMrel == 'P17':
            AAA = 1.685
            BBB = -29.0727644    # -4.585
        else:
            print('Woops')
    else:
        import pdb;pdb.set_trace()

    return AAA,BBB

def get_YM_sys_err(logy,YMrel,delta=500,ySph=True,h70=1.0):

    if delta == 500:
        if YMrel == 'A10':
            pivot = 3e14; Jofx  = 0.6145 if ySph else 0.7398
            Norm  = 2.925e-5 * Jofx * h70**(-1); PL = 1.78
            #t1   = ((logy - 1)/PL )*0.024
            t1   = 0.024 / PL
            t2   = ((np.log10(Norm) - logy)/PL**2)*0.08
            xer  = np.sqrt(t1**2 + t2**2) * np.log(10)
            #import pdb;pdb.set_trace()
        elif YMrel == 'M12':
            t1   = np.array([1.0,logy+5])
            #t1   = np.array([0.367,0.44])
            tcov = np.array([[0.098**2,-0.012],[-0.012,0.12**2]])
            #tcov = np.array([[0.098**2,-(0.012**2)],[-(0.012**2),0.12**2]])
            #tcov = np.array([[0.098**2,0],[0,0.12**2]])
            t2   = np.abs(np.matmul(t1,np.matmul(tcov,t1)))
            xer  = np.sqrt(t2) * np.log(10)
            print(xer)
            #import pdb;pdb.set_trace()
        elif YMrel == 'M12-SS':
            t1   = 0.0 # Fixed slope
            t2   = 0.036
            xer  = np.sqrt(t1**2 + t2**2) * np.log(10)
            print(xer)
        elif YMrel == 'P17':
            Norm = -4.305; PL = 1.685
            t1   = 0.009 / PL
            t2   = ((Norm - logy)/PL**2)*0.013
            xer  = np.sqrt(t1**2 + t2**2) * np.log(10)
            #xer = 0.104
        else:
            print('No match!')
            import pdb;pdb.set_trace()
    elif delta == 2500:
        if YMrel == 'A10':
            #LogNorm = -28.13; PL = 1.637
            #t1   = 0.88 / PL 
            #t2   = ((logy - LogNorm)/PL**2)*0.062
            #xer  = np.sqrt(t1**2 + t2**2) * np.log(10)
            xer  = np.log(1 + 0.23)
        elif YMrel == 'M12':
            t1   = np.array([1.0,logy+5])
            #t1   = np.array([0.367,0.44])
            #tcov = np.array([[0.063**2,-0.008],[-0.008,0.14**2]])
            tcov = np.array([[0.063**2,-(0.008**2)],[-(0.008**2),0.14**2]])
            #tcov = np.array([[0.098**2,0],[0,0.12**2]])
            t2   = np.abs(np.matmul(t1,np.matmul(tcov,t1)))
            xer  = np.sqrt(t2) * np.log(10)
            print(xer)
        elif YMrel == 'M12-SS':
            t1   = 0.0 # Fixed slope
            t2   = 0.033
            xer  = np.sqrt(t1**2 + t2**2) * np.log(10)
            print(xer)
        elif YMrel == 'P17':
            Norm = -4.5855; PL = 1.755
            t1   = 0.014 / PL
            t2   = ((Norm - logy)/PL**2)*0.020
            xer  = np.sqrt(t1**2 + t2**2) * np.log(10)
            #xer = 0.104
        else:
            print('No match!')
            import pdb;pdb.set_trace()
    else:
        print('No match for delta!')
        import pdb;pdb.set_trace()

    return xer
        
