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
            self.name='rdcs0910'
            #############################################################
            if reduction == 'CMCORR':
                if map_file == 'noise':
                    self.indir= m2dir+"AGBT17_Products/RDCS0910/"
                    #self.fitsfile=self.indir+"Kelvin_RDCS0910_2aspix_cmcorr_v3_noise_iter1.fits"
                    self.fitsfile=self.indir+"Kelvin_RDCS0910_2aspix_pca3_0f05_bugfixed_noise_iter1.fits"
                if map_file == 'all':
                    self.indir= m2dir+"AGBT17_Products/RDCS0910/"
                    #self.fitsfile=self.indir+"Kelvin_RDCS0910_2aspix_cmcorr_v3_map_iter1.fits"
                    self.fitsfile=self.indir+"Kelvin_RDCS0910_2aspix_pca3_0f05_bugfixed_map_iter1.fits"
                self.wtfile=self.fitsfile # It's in the same file; just a different extension.
                self.wtext=1         # The extension of the weight (or RMS) array
                self.wtisrms=False   # The "weight" file is actual the RMS of pixels
                self.units='Kelvin'  # A fundamental, critical, and wholly important variable!!

            if reduction == 'PCA':
                if map_file == 'noise':
                    self.indir= m2dir+"AGBT17_Products/RDCS0910/"
                    #self.fitsfile=self.indir+"pca3_f_0.09_roi1.5_half_kelvin_noise.fits"
                    self.fitsfile=self.indir+"Kelvin_RDCS0910_2aspix_pca2_0f055_bugfixed_noise_iter1.fits"
                if map_file == 'all':
                    self.indir= m2dir+"AGBT17_Products/RDCS0910/"
                    #self.fitsfile=self.indir+"RDCS0910_PCA3_map.fits"
                    self.fitsfile=self.indir+"Kelvin_RDCS0910_2aspix_pca2_0f055_bugfixed_map_iter1.fits"
                self.wtfile=self.fitsfile # It's in the same file; just a different extension.
                self.wtext=1         # The extension of the weight (or RMS) array
                self.wtisrms=False   # The "weight" file is actual the RMS of pixels
                self.units='Kelvin'  # A fundamental, critical, and wholly important variable!!

            if reduction == 'MINKASI':
                if map_file == 'noise':
                    self.indir= m2dir+"MINKASI/RDCS0910/"
                    self.fitsfile=self.indir+"Kelvin_RDCS0910_2aspix_cmcorr_v3_noise_iter1.fits"
                if map_file == 'all':
                    self.indir= m2dir+"MINKASI/RDCS0910/"
                    self.fitsfile=self.indir+"Kelvin_RDCS0910_2aspix_cmcorr_v3_map_iter1.fits"
                self.wtfile=None     # It's in the same file; just a different extension.
                self.wtext=0         # The extension of the weight (or RMS) array
                self.wtisrms=False   # The "weight" file is actual the RMS of pixels
                self.units='Kelvin'  # A fundamental, critical, and wholly important variable!!
                
            ###############################################################################
            ### Here's what I need to know about the transfer function format:
            ### If there are issues, please see code XYZ

            #self.tabfile = m2dir+"Template_M2_xfer_fxn_srcsz_210.txt"
            #self.tabfile = m2dir+"AGBT17_Products/pca7_f0.09_onHSC_2.txt"
            longname='Kelvin_2XMMJ0830+5241_2aspix_xfer_all_v2_transfer_function_all_runs.txt'
            shortname='Xfer_2XMMJ0830+5241.txt'
            self.tabfile = m2dir+"IDL_Xfer_Fxns/"+shortname
            self.tabcomments='#'
            self.tabformat = 'ascii'
            self.tabdims = '1D'
            self.tabextend = True    # Do we need to extent to higher k numbers?
            self.calunc = 0.1      # 10% calibration accuracy.
            self.fitptsrcs = True
            self.fitmnlvl  = True
            self.rmscorr   = 1.17
            
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
        
        self.z=1.10                      # Redshift
        #self.ra = Angle('09h10m45.761s')     # Right Ascencion, in hours
        #self.dec= Angle('+54d22m3.97s')       # Declination, in degrees
        #self.ra = Angle('09h10m45.086s')     # Right Ascencion, in hours
        #self.dec= Angle('+54d22m04.17s')       # Declination, in degrees
        #self.ra = Angle('09h10m43.919s')     # Right Ascencion, in hours
        #self.dec= Angle('+54d22m14.46s')       # Declination, in degrees
        self.ra = Angle('09h10m44.294s')     # Right Ascencion, in hours
        self.dec= Angle('+54d22m14.28s')       # Declination, in degrees
        self.M_500 = 7.0e14                # Solar masses
        self.Tx    = 6.4                  # keV
        self.name  = 'rdcs0910'
        
        ###  For when the time comes to use the *actual* coordinates for Abell 2146,
        ###  Here they are. Even now, it's useful to calculate the offsets of the centroids
        ###  for the radius of curvature of the shocks.

class private_vars:

    def __init__(self):
        
        #RXJ1347_ra    = Angle('09h10m45.086s') ; RXJ1347_dec    = Angle('+54d22m04.17s')
        RXJ1347_ra    = Angle('09h10m43.919s') ; RXJ1347_dec    = Angle('+54d22m14.46s')

class shocks:

    def __init__(self):

        ### GEOPARAMS = [X-offset (arcsec), Y-offset (arcsec), Major-Axis rotation,
        ###              X-axis scaling, Y-axis scaling, Z-axis scaling,
        ###              Taper scaling (power law), Opening Angle]
        ###              X-axis is taken to be the major axis...
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
        self.bins = [6]      
        ### Or, you can specify an array, which *MUST* then have units
        ### attached to it.
        #self.minarc = 2.0*u.arcsec
        #self.maxarc = 
        self.fbtemps = [False]
        self.bulkalp  = np.zeros(self.bins) # set to zeros -> these will be "fit for".
        self.narm  = [True]   # Normalize at R_min
        self.taper = ['normal']   # A bit silly to have, but it's better...
        self.fit_cen = [False]
        self.fit_geo = [True]

class blob:

    def __init__(self):

        #self.ra = [Angle('09h10m45.663s')]     # Right Ascencion, in hours
        #self.dec= [Angle('+54d22m04.28s')]       # Declination, in degrees
        self.ra = [Angle('09h10m45.502s')]     # Right Ascencion, in hours
        self.dec= [Angle('+54d22m04.04s')]       # Declination, in degrees
        blobparams = [0.5,0.5,2.0,0.5,0.2,-1.0e-3]
        ### A second blob...
        #self.ra = [Angle('09h10m42.875s')]     # Right Ascencion, in hours
        #self.dec= [Angle('+54d22m19.72s')]       # Declination, in degrees
        #blobparams = [0.5,0.5,2.0,0.5,0.2,-1.0e-3]
        ### Try without the blob...
        #blobparams = []
        self.blobpars = [blobparams]
        self.dofit    = {'MUSTANG2':True}
        
class ptsrc:

    def __init__(self):

        ### What I found with mpfitfun or whatever:
        #self.locs  = [(Angle('13h47m30.6s'),Angle('-11d45m10.24s'))]
        ### From ALMA paper:
        #self.locs  = [(Angle('13h47m30.62s'),Angle('-11d45m09.5s'))]
        ### Let me just try something...:
        #self.locs  = [(Angle('02h21m45.184s'),Angle('-02d46m14.94s'))]
        self.locs = []
        
        ### Enter 0 if a point source is truly point-like.
        self.fwhm  = [9.0]
        ### From Kitayama et al. 2016: ~4.00 +/- 0.03 +/- 0.25 mJy at 90 GHz
        #self.prior    = {'MUSTANG2':[0.003*u.K]}
        #self.priorunc = {'MUSTANG2':[0.0002*u.K]}
        self.prior    = {'MUSTANG2':[0.000*u.K]}
        self.priorunc = {'MUSTANG2':[0.02*u.K]}
