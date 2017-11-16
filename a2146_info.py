from astropy.coordinates import Angle
from astropy.coordinates import SkyCoord
#from astropy.coordinates import separation
import numpy as np
import astropy.units as u          # Install astropy

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

class files:

    def __init__(self,instrument="MUSTANG2",map_file='nw_model'):

        ###############################################################################
        ### Some fundamental files (for data), as well as their formatting.

        if instrument == "MUSTANG2": 
        
            self.instrument="MUSTANG2"
            self.name='abell_2146'
#############################################################
            if map_file == 'noise':
                self.indir= "/home/data/MUSTANG2/"
                self.fitsfile=self.indir+"Noise_Map_M2_27_Jan_2017.fits"
            if map_file == 'nw_model':
                self.indir= "/home/romero/Results_Python/MUSTANG2/a2146/"
                #self.fitsfile=self.indir+"MUSTANG2_Real_Mock_Observation_Map.fits"
                self.fitsfile=self.indir+"MUSTANG2_Real_Full_Run_Mock_Observation_Map_wNoise.fits"
            if map_file == 'se_model':
                self.indir= "/home/romero/Results_Python/MUSTANG2/a2146/"
                #self.fitsfile=self.indir+"MUSTANG2_Real_Mock_Observation_Map.fits"
                self.fitsfile=self.indir+"MUSTANG2_Real_Full_Run_Mock_Observation_Map_se_v5.fits"
#############################################################
            self.wtfile=self.fitsfile # It's in the same file; just a different extension.
            self.wtext=1         # The extension of the weight (or RMS) array
            self.wtisrms=False   # The "weight" file is actual the RMS of pixels
        
            ###############################################################################
            ### Here's what I need to know about the transfer function format:
            ### If there are issues, please see code XYZ

            m2dir = "/home/data/MUSTANG2/"
            self.tabfile = m2dir+"Template_M2_xfer_fxn_srcsz_210.txt"
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
        
        self.z=0.2323                      # Redshift
        self.ra = Angle('10h53m41.8369s')  # Hour angle
        self.dec= Angle('57d35m21.0754s')  # Degrees
        self.M_500 = 5.4e14                # Solar masses
        self.Tx = 5.6                      # keV
        self.name='abell_2146'
        
        ###  For when the time comes to use the *actual* coordinates for Abell 2146,
        ###  Here they are. Even now, it's useful to calculate the offsets of the centroids
        ###  for the radius of curvature of the shocks.

class private_vars:

    def __init__(self):
        
        A2146_ra    = Angle('15h56m08.9s'); A2146_dec    = Angle('66d21m21s')
        A2146_se_ra = Angle('15h56m07.1s'); A2146_se_dec = Angle('66d22m45s')
        A2146_nw_ra = Angle('15h56m08.5s'); A2146_nw_dec = Angle('66d21m35s')
        se_ra_off = A2146_se_ra-A2146_ra  ; se_dec_off = A2146_se_dec-A2146_dec 
        nw_ra_off = A2146_nw_ra-A2146_ra  ; nw_dec_off = A2146_nw_dec-A2146_dec 
        nw_x_off = -nw_ra_off.to('arcsec') * np.cos(A2146_dec.to('rad'))
        se_x_off = -se_ra_off.to('arcsec') * np.cos(A2146_dec.to('rad'))
        nw_y_off = nw_dec_off.to('arcsec'); se_y_off = se_dec_off.to('arcsec')

        self.nw_x_off = nw_x_off
        self.nw_y_off = nw_y_off
        self.se_x_off = se_x_off
        self.se_y_off = se_y_off

        ###############################################################################
        ### Prior known values regarding the cluster, point sources, or other values
        ### relevant to fitting parameters (i.e putting priors on parameters).
        ### For point sources, I will list their fluxes and uncertaints (*in a list*).

#        psfd, psfdunc = calc_ptsrc_flux() # Take prior measurements, and apply MUSTANG bandwidth.
#        self.psfd = np.array([psfd])  # For a band-averaged flux
#        self.psunc= np.array([psfdunc])

        ###############################################################################
        ### Supplementary files (maps, for now)
        
        self.psfile = None    # If you have a known point source (model) fits file
        self.blfile = None    # If you have a known "blob" (model) fits file
        self.shfile = None    # If you have a known shock (model) fits file
        self.miscfile1 = None # I'm trying allow for extra files for future needs
        self.miscfile2 = None # I think...it's not so necessary. But surely it's not
        self.miscfile3 = None # draining memory to do it like this.
        
        ###############################################################################
        ### Bulk Profile (just for a a really simple model)
        
        self.rads = np.arange(11)*100.0 # 0 through 1 Mpc
        eden = np.arange(11); eden[0]=eden[1]; eden = 10**(-((eden-1.0)/5.0)**1.7)
        self.eden = eden*0.015        # cm**-3
        self.temp = np.zeros(11)+5.5  # keV        
        # I should already have the infrastructure to find the alphas??

        
        ###############################################################################
        ### Shock Profiles /// THESE ARE FOR YOUR OWN USE!!!
        
        self.rads_se = np.array([450.0,498.0,650.0])*u.kpc  # in kpc
        self.eden_se = np.array([0.01,0.003,0.0008])  # cm**-3
        self.temp_se = np.array([2.8,10.5,7.0])       # keV
#        self.ealp_se = np.array([0.0,-0.29,-1.26])  # Actual values inferred from Russell+ 2012
        self.ealp_se = np.array([0.0,-0.29,-2.2]) # Force to have a steeper outer profile.
        self.talp_se = np.array([0.0,-0.1,-0.3])
        self.angl_se = (240.0*u.deg).to("radian").value
        self.opan_se = np.pi
        self.shxi_se = 2.0
        self.lban_se = -np.pi
        self.uban_se = -np.pi/2.0
        
        self.rads_nw = np.array([60.0,225.0,300.0])*u.kpc   # 0.0,165.0,240.0
        self.eden_nw = np.array([0.0038,0.0038,0.0018])
        self.temp_nw = np.array([10.0, 16.5, 7.5])
#        self.ealp_nw = np.array([-0.1,-0.2,-2.0]) # Actual values inferred from Russell+ 2012
        self.ealp_nw = np.array([-0.1,-0.2,-2.5]) # Force to have a steeper outer profile.
        self.talp_nw = np.array([0.3,-0.75,0.0])
        self.angl_nw = (55.0*u.deg).to("radian").value
        self.opan_nw = np.pi
        self.shxi_nw = 1.0
        self.lban_nw = 0
        self.uban_nw = np.pi/2.0
        
        self.narm = True # Normalize at R_min
        self.taper = 'normal'  #'inverse'

class shocks:

    def __init__(self):

        geoparams1 = [0,0,0,1,1,1,0,0] # Spherical Geometry
        mypriors = private_vars()
        etemperature1 = mypriors.temp_nw
        geoparams1[0] = mypriors.nw_x_off.to("arcsec").value
        geoparams1[1] = mypriors.nw_y_off.to("arcsec").value
        geoparams1[2] = mypriors.angl_nw;  geoparams1[6] = mypriors.shxi_nw
        geoparams1[7] = mypriors.opan_nw;  kpc_rad1 = mypriors.rads_nw

        geoparams2 = [0,0,0,1,1,1,0,0] # Spherical Geometry
        etemperature2 = mypriors.temp_se;
        geoparams2[0] = mypriors.se_x_off.to("arcsec").value
        geoparams2[1] = mypriors.se_y_off.to("arcsec").value
        geoparams2[2] = mypriors.angl_se;  geoparams2[6] = mypriors.shxi_se
        geoparams2[7] = mypriors.opan_se;  kpc_rad2 = mypriors.rads_se

### If you don't want to fit for shock components:
        #self.geoparams=[]
### Otherwise:
        self.geoparams=[geoparams1,geoparams2]

### If you don't want to use the prior-determined radii:
        #self.bins=None
### Else:
        self.bins=[kpc_rad1,kpc_rad2]
        self.fstemps=[False,False]
        self.shockalp  = [np.zeros(len(kpc_rad1)),np.zeros(len(kpc_rad2))] # set
        self.taper = ['normal','normal']
        self.narm  = [True, True]  # Normalize at R_min


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

        self.locs = []
