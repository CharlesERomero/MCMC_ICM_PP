from astropy.coordinates import Angle
from astropy.coordinates import SkyCoord
#from astropy.coordinates import separation
import numpy as np
import astropy.units as u          # Install astropy
from os.path import expanduser
myhome = expanduser("~")

###############################################################################
### This is only for MODELLING OF ABELL 2146!!! 
###############################################################################

### CLASSES WHICH ARE USED BY OTHER ROUTINES:
#
# (1) files
# (2) priors
# (3) shocks
#
# That is... private_vars is only called within this routine, so that I (or you)
# can organize have some presets if you need to look back on them later, but you'll
# have to reformulate how they get used in the actual fitting code.

######################################################
### Let's list all directories by instrument:
mdir   = "/home/data/MUSTANG/"
m2dir  = "/home/data/MUSTANG2/"
bodir  = "/home/data/Bolocam/"
nk2dir = "/home/data/NIKA2/"
nkdir  = "/home/data/NIKA/"

class files:

    def __init__(self,instrument="MUSTANG2",map_file='all',reduction='PCA'):

        ###############################################################################
        ### Some fundamental files (for data), as well as their formatting.

        if instrument == "MUSTANG2": 
        
            self.instrument="MUSTANG2"
            self.name='idcs'
            #############################################################
            if reduction == 'PCA':
                if map_file == 'noise':
                    self.indir= m2dir+"AGBT17_Products/IDCS1426+35/"
                    #self.fitsfile=self.indir+"Kelvin_idcs1426+35_2aspix_pca3_0f05_bugfixed_noise_iter1.fits"
                    self.fitsfile=self.indir+"Kelvin_idcs1426+35_2aspix_pca4_0f045_fullcov_noise_iter1.fits"
                if map_file == 'all':
                    self.indir= m2dir+"AGBT17_Products/IDCS1426+35/"
                    #self.fitsfile=self.indir+"Kelvin_idcs1426+35_2aspix_pca3_0f05_bugfixed_map_iter1.fits"
                    self.fitsfile=self.indir+"Kelvin_idcs1426+35_2aspix_pca4_0f045_fullcov_map_iter1.fits"
                #############################################################
                self.wtfile=self.fitsfile # It's in the same file; just a different extension.
                self.wtext=1         # The extension of the weight (or RMS) array
                self.wtisrms=False   # The "weight" file is actual the RMS of pixels
                self.units='Kelvin'  # A fundamental, critical, and wholly important variable!!
                #            self.units='Kelvin'  # A fundamental, critical, and wholly important variable!!
            else:
                raise Exception
                
            ###############################################################################
            ### Here's what I need to know about the transfer function format:
            ### If there are issues, please see code XYZ

            #self.tabfile = m2dir+"Template_M2_xfer_fxn_srcsz_210.txt"
            self.tabfile = m2dir+"IDL_Xfer_Fxns/2XMMJ0830+5241_highres_xfer.txt"
            self.tabcomments='#'
            self.tabformat = 'ascii'
            self.tabdims = '1D'
            self.tabextend = True    # Do we need to extent to higher k numbers?
            self.calunc = 0.1      # 10% calibration accuracy.
            self.fitptsrcs = True
            self.fitmnlvl  = True
            self.rmscorr   = 1.00

        if instrument == "NIKA2":

            self.calunc = 0.07      # 10% calibration accuracy.
            self.fitptsrcs = True
            self.fitmnlvl  = True
            print 'This section not developed yet!'
            self.units='Jy/beam'
            import pdb; pdb.set_trace()
            self.rmscorr   = 1.00
            
        if instrument == "BOLOCAM":

            self.indir= bodir+"MACS_J1347.5-1144/"
            if map_file == 'noise':
                self.fitsfile=self.indir+"filtered_image_noise_realizations.fits"
            if map_file == 'all':
                self.fitsfile=self.indir+"filtered_image.fits"
            self.wtfile=self.indir+"filtered_image_rms.fits"
            
            self.calunc = 0.05      # 10% calibration accuracy.
            #fitsfile=boldir+"filtered_image.fits"
            self.wtext=1         # The extension of the weight (or RMS) array
            self.wtisrms=True   # The "weight" file is actual the RMS of pixels
            self.units='Kelvin'  # A fundamental, critical, and wholly important variable!!
        
            ###############################################################################
            ### Here's what I need to know about the transfer function format:
            ### If there are issues, please see code XYZ

            #self.tabfile = m2dir+"Template_M2_xfer_fxn_srcsz_210.txt"
            self.tabfile = self.indir+"filtered_image_signal_transfer_function.fits"
            self.tabcomments='#'
            self.tabformat = 'fits'
            self.tabdims = '2D'
            self.tabextend = True    # Do we need to extent to higher k numbers?
            self.calunc = 0.1      # 10% calibration accuracy.

            self.fitptsrcs = False
            self.fitmnlvl  = False
            print 'This section not developed yet!'
            self.units='Kelvin'
            self.rmscorr   = 1.00
            import pdb; pdb.set_trace()
            
class priors:
        
    def __init__(self):
        
        ###############################################################################
        ### Prior known values regarding the RXJ1053. Redshift, ra, and dec *MUST* be
        ### known / accurate. M_500 and Tx are useful for creating initial guesses.
        ### Tx is still important if relativistic corrections may be severe.
        
        self.z=1.75                        # Redshift
        #self.ra = Angle('02h21m45.184s')  # Right Ascencion, in hours
        #self.dec= Angle('-03d46m14.94s')  # Declination, in degrees
        self.ra = Angle('14h26m33.0s')   # Right Ascencion, in hours
        self.dec= Angle('+35d08m02.6s')   # Declination, in degrees
        self.M_500 = 2.6e14                # Solar masses
        self.Tx    = 5.0                   # keV
        self.name  = 'idcs'
        
        ###  For when the time comes to use the *actual* coordinates for Abell 2146,
        ###  Here they are. Even now, it's useful to calculate the offsets of the centroids
        ###  for the radius of curvature of the shocks.

#class private_vars:
#
#    def __init__(self):
#        
#        RXJ1347_ra    = Angle('14h26m33.089s'); RXJ1347_dec    = Angle('+35d08m34.01s') 
#        RXJ1347_ra    = Angle('14h26m33.0s'); RXJ1347_dec    = Angle('+35d08m02.6s')
 
class shocks:

    def __init__(self):

        ### GEOPARAMS = [X-offset (arcsec), Y-offset (arcsec), Major-Axis rotation,
        ###              X-axis scaling, Y-axis scaling, Z-axis scaling,
        ###              Taper scaling (power law), Opening Angle]
        ### NOTE: If Taper scaling (geoparams[6]) is 0, then no tapering is applied
        geoparams = [0,0,0,1,1,1,0,0] # Spherical Geometry
        ### These parameters were used in March. I think I want to change things a bit...
        #rxjshock  = [0,0,4.06,1,1,1,2.0,3.1415]   # Angles specified in radians!!! (3.4 to 4.2 rad)
        ### Let's change 
        #rxjshock  = [0,0,4.06,1,1,1,0.0,0.8]   # Angles specified in radians!!! (3.4 to 4.2 rad)
        ### 3.75 radians = angle of shock  (SW)
        ### 0.7 radians  = opening angle       ; 1.31

        ### For RXJ1347, 1" ~ 6 kpc, so 60 kpc bins ~ 10" is well matched to our beam.
        mybins = np.array([60.0,120.0,180.0,240.0,300.0])*u.kpc
        
### If you don't want to fit for shock components:
        self.geoparams = []                       # Array of geometric parameters
        self.bins      = []                       # Bins, specified in kpc
        self.fstemps   = []                       # Fit for shock temperatures?
        self.shockalp  = []                       # set
        self.taper     = ['normal']               # Type of taper. 'normal' is recommended.
        self.narm      = [True]                   # Normalize at R_min
        self.shockfin  = [True]
### Otherwise, for RXJ1347:
#        self.geoparams = [rxjshock]               # Array of geometric parameters
#        self.bins      = [mybins]                 # Array OF arrays. Units necessary.
#        self.fstemps   = [False]                  # Fit for shock temperatures?
#        self.shockalp  = [np.zeros(len(mybins))]  # set
#        self.taper     = ['normal']               # Type of taper. 'normal' is recommended.
#        self.narm      = [True]                   # Normalize at R_min
#        self.shockfin  = [True]                   # Finite integration (out to last bin)
#############################################################################
### LEGACY / REFERENCE CODE. TO BE DELETED WITH ENOUGH FAMILIARITY WITH CODE
### Reference commands from Abell 2146:
#        self.geoparams=[geoparams]
#        self.bins=[kpc_rad1,kpc_rad2]
#        self.fstemps=[False,False]
#        self.shockalp  = [np.zeros(len(kpc_rad1)),np.zeros(len(kpc_rad2))] # set
#        self.taper = ['normal','normal']
#        self.narm  = [True, True]  # Normalize at R_min

class bulk:

    def __init__(self):

        geoparams = [0,0,0,1,1,1,0,0] # Spherical Geometry

        self.geoparams=[geoparams]
        ### You can specify the number of bins (as a LIST, as below):
        self.bins = [4]      
        #self.bins = []   # No bulk component!      
        ### Or, you can specify an array, which *MUST* then have units
        ### attached to it.
        #self.minarc = 2.0*u.arcsec
        #self.maxarc = 
        self.fbtemps = [False]
        self.bulkalp  = np.zeros(self.bins) # set to zeros -> these will be "fit for".
        self.narm  = [True]   # Normalize at R_min
        self.taper = ['normal']   # A bit silly to have, but it's better...
        self.fit_cen = [True]
        self.fit_geo = [False]

class blob:

    def __init__(self):

        self.ra = [Angle('09h10m45.663s')]     # Right Ascencion, in hours
        self.dec= [Angle('+54d22m04.28s')]       # Declination, in degrees
        #blobparams = [0.5,0.5,2.0,0.5,0.2,-1.0e-3]
        ### Try without the blob...
        blobparams = []
        self.blobpars = [blobparams]
        self.dofit    = {'MUSTANG2':False}
        
class ptsrc:

    def __init__(self):

        ### The NEGATIVE POINT SOURCE:
        neglocs  = [Angle('14h26m33.089s'),Angle('+35d08m34.01s')]
        negfwhm  = 9.0
        m2negpr  = -0.0004*u.K
        m2negunc = 0.1*u.K
        #self.priorunc = {'MUSTANG2':[0.0004*u.K]}

        ### THE POSITIVE POINT SOURCE
        poslocs  = [Angle('14h26m32.252s'),Angle('+35d08m14.59s')]
        posfwhm  = 9.0
        m2pospr  = 0.0012*u.K
        m2posunc = 0.1*u.K

        self.locs = [neglocs,poslocs]
        self.fwhm = [negfwhm,posfwhm]
        self.prior    = {'MUSTANG2':[m2negpr,m2pospr]}
        self.priorunc = {'MUSTANG2':[m2negunc,m2posunc]}

        #self.locs = [poslocs]
        #self.fwhm = [posfwhm]
        #self.prior    = {'MUSTANG2':[m2pospr]}
        #self.priorunc = {'MUSTANG2':[m2posunc]}
    
