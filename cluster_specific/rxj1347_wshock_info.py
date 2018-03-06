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
m2dir = "/home/data/MUSTANG2/"


class files:

    def __init__(self,instrument="MUSTANG2",map_file='all'):

        ###############################################################################
        ### Some fundamental files (for data), as well as their formatting.

        if instrument == "MUSTANG2": 
        
            self.instrument="MUSTANG2"
            self.name='rxj1347'
#############################################################
            if map_file == 'noise':
                self.indir= m2dir+"AGBT17_Products/RXJ1347/"
                self.fitsfile=self.indir+"pca7_f0.09_noise.fits"
            if map_file == 'all':
                self.indir= m2dir+"AGBT17_Products/RXJ1347/"
                self.fitsfile=self.indir+"pca7_f0.09_map.fits"
#############################################################
            self.wtfile=self.fitsfile # It's in the same file; just a different extension.
            self.wtext=1         # The extension of the weight (or RMS) array
            self.wtisrms=False   # The "weight" file is actual the RMS of pixels
        
            ###############################################################################
            ### Here's what I need to know about the transfer function format:
            ### If there are issues, please see code XYZ

            #self.tabfile = m2dir+"Template_M2_xfer_fxn_srcsz_210.txt"
            self.tabfile = m2dir+"AGBT17_Products/pca7_f0.09_onHSC_2.txt"
            self.tabcomments='#'
            self.tabformat = 'ascii'
            self.tabdims = '1D'
            self.tabextend = True    # Do we need to extent to higher k numbers?
            self.calunc = 0.1      # 10% calibration accuracy.
            self.fitptsrcs = True
            self.fitmnlvl  = True

        if instrument == "NIKA2":

            self.calunc = 0.07      # 10% calibration accuracy.
            self.fitptsrcs = True
            self.fitmnlvl  = True
            print 'This section not developed yet!'
            import pdb; pdb.set_trace()
            
        if instrument == "BOLOCAM":

            self.calunc = 0.05      # 10% calibration accuracy.
            self.fitptsrcs = False
            self.fitmnlvl  = False
            print 'This section not developed yet!'
            import pdb; pdb.set_trace()
            
class priors:
        
    def __init__(self):
        
        ###############################################################################
        ### Prior known values regarding the RXJ1053. Redshift, ra, and dec *MUST* be
        ### known / accurate. M_500 and Tx are useful for creating initial guesses.
        ### Tx is still important if relativistic corrections may be severe.
        
        self.z=0.4510                      # Redshift
        self.ra = Angle('13h47m30.5s')     # Right Ascencion, in hours
        self.dec= Angle('-11d45m9s')       # Declination, in degrees
        self.M_500 = 2.2e15                # Solar masses
        self.Tx    = 10.8                  # keV
        self.name  = 'rxj1347'
        
        ###  For when the time comes to use the *actual* coordinates for Abell 2146,
        ###  Here they are. Even now, it's useful to calculate the offsets of the centroids
        ###  for the radius of curvature of the shocks.

class private_vars:

    def __init__(self):
        
        RXJ1347_ra    = Angle('13h47m30.5s'); RXJ1347_dec    = Angle('-11d45m9s')
 
class shocks:

    def __init__(self):

        ### GEOPARAMS = [X-offset (arcsec), Y-offset (arcsec), Major-Axis rotation,
        ###              X-axis scaling, Y-axis scaling, Z-axis scaling,
        ###              Taper scaling (power law), Opening Angle]
        ### NOTE: If Taper scaling (geoparams[6]) is 0, then no tapering is applied
        geoparams = [0,0,0,1,1,1,0,0] # Spherical Geometry
        rxjshock  = [0,0,3.75,1,1,1,0,0.7]   # Angles specified in radians!!!
        ### 3.75 radians = angle of shock  (SW)
        ### 0.7 radians  = opening angle

        ### For RXJ1347, 1" ~ 6 kpc, so 60 kpc bins ~ 10" is well matched to our beam.
        mybins = np.array([60.0,120.0,180.0,240.0])*u.kpc
        
### If you don't want to fit for shock components:
#        self.geoparams = []                       # Array of geometric parameters
#        self.bins      = []                       # Bins, specified in kpc
#        self.fstemps   = []                       # Fit for shock temperatures?
#        self.shockalp  = []                       # set
#        self.taper     = ['normal']               # Type of taper. 'normal' is recommended.
#        self.narm      = [True]                   # Normalize at R_min
### Otherwise, for RXJ1347:
        self.geoparams = [rxjshock]               # Array of geometric parameters
        self.bins      = [mybins]                 # Array OF arrays. Units necessary.
        self.fstemps   = [False]                  # Fit for shock temperatures?
        self.shockalp  = [np.zeros(len(mybins))]  # set
        self.taper     = ['normal']               # Type of taper. 'normal' is recommended.
        self.narm      = [True]                   # Normalize at R_min
        self.shockfin  = [True]                   # Finite integration (out to last bin)
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
        self.bins = [6]      
        ### Or, you can specify an array, which *MUST* then have units
        ### attached to it.
        #self.minarc = 2.0*u.arcsec
        #self.maxarc = 
        self.fbtemps = [False]
        self.bulkalp  = np.zeros(self.bins) # set to zeros -> these will be "fit for".
        self.narm  = [True]   # Normalize at R_min
        self.taper = ['normal']   # A bit silly to have, but it's better...

class ptsrc:

    def __init__(self):

        self.locs = [(Angle('13h47m30.6s'),Angle('-11d45m10.24s'))]
        ### Enter 0 if a point source is truly point-like.
        self.fwhm = [9.0]
