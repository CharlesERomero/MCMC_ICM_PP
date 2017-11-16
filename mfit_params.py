import time
import os
import find_ptsrc as fp
import retrieve_data_info as rdi
import astropy.units as u
import numpy as np

class common_fit_params:

    def __init__(self,bins=[6],shbins=[6],path='/home/romero/Results_Python/',
                 bulkgeo=[],bulknarm=[False],bulkcen=[],bulkalp=[],
                 shockgeo=[],shocknarm=[True],shockalp=[],
                 ptsrcs=[],blobs=[],
                 fbtemps=[False],fstemps=[False],minmax=np.array([20.0,50.0])*u.arcsec,
                 cluster=None,testmode=False,autodetect=False):

    ##################################################################################
    #####      Let's first prepare the bins for our bulk and shock components    #####
    ##################################################################################

        bulkarc=[]
        mygeo=[]
        geoparams=[0,0,0,1,1,1,0,0]
        totbins=0
        for mybins in bins:
            if hasattr(mybins,"__len__"):
                myarc = mybins.to("rad").value
            else:
                nbins = mybins
                radnx = minmax.to("rad").value
                nzbins= np.logspace(np.log10(radnx[0]),np.log10(radnx[1]), nbins)
                myarc = nzbins

            ### Create lists of the bin positions and power-law slope (default=0).
            bulkarc.append(myarc)
            if len(bulkalp) == 0:
                bulkalp.append(myarc*0)
            totbins  += len(myarc)
            mygeo.append(geoparams)

        if len(bulkgeo) == 0:
            bulkgeo = mygeo

        myshbins=[]
        anotset= (len(shockalp) == 0)
        for mybins in shbins:
            if mybins.unit.is_equivalent("kpc"):
                mybins = convert_kpc_to_rad(mybins,cluster.d_a)
            elif mybins.unit.is_equivalent("rad"):
                mybins = mybins.to("rad").value
            else:
                raise Exception("Your shock bins have no units; Add units in your INFO file.")
            
            ### Create lists of the bin positions and power-law slope (default=0).
            myshbins.append(mybins)
            if anotset:
                shockalp.append(np.zeros(len(mybins)))
            totbins+=len(mybins)
                                    
    #########################################################################################
    #####  Now let's define some attributes for the common fitting parameter (class)    #####
    #########################################################################################

        self.bulk     = bulkcen    # Provide a list of centroids.
        self.bulkbins = bins       # The number of (radial) bins
        self.bminmax  = minmax     #
        self.bulkarc  = bulkarc    # It will be set, based on all instruments used.
                                   # We want this to be an array of bins (in arcseconds).
        self.bulkalp  = bulkalp    # set the power law index (alpha)
        self.bulkgeo  = bulkgeo    # set the geometry for the bulk model
        self.bulknarm = bulknarm   # set the normalization method
        self.bulkfix  = False
        self.fbtemps  = fbtemps    # Fit for bulk profile temperatures (if X-ray data is present)
        
        self.ptsrc    = ptsrcs     # Provide a list of centroids.
        self.shockalp = shockalp   # Provide a list of shock log pressure slope
        self.shockgeo = shockgeo   # Provide a list of shock geometries
        self.shockbin = myshbins   # How many bins to use; OR, an array of bins positions.
        self.shocknarm= shocknarm  # set
        self.shockfix = True       # Fix alphas for shock? By default, yes.
        self.fstemps  = fstemps    # Fit for shock profile temperatures (if X-ray data is present)
        self.trimOut  = True       # Trim outer (shock) shell.
        
        ### What additional components do we care to fit for?
        #### IMPORTANT: ###
        # If you want to specify the bin locations, you should do so within the file
        # that spits out "priors"
        
        self.blob     = blobs      # Do not fit for a "blob" component

### I think I will want to add masks for bulk and shock components. Or at least
### the option to have them... I'm not sure how I want to do this yet, so I guess I will
### hold off on this for now (08 Nov 2017).
        
### I CURRENTLY HAVE A BUG WITH mn_lvl. (20 June 2017)
        self.pprior = False      # Use a Planck Prior on Y_int
        self.ubmap = False       # Uniform bins (pressure is flat)
        self.real_data = True    # 
        self.joint = False       # I never made it work jointly...
        self.path = path         # directory where figures are saved
        #### Num. of free dimensions
        #import pdb;pdb.set_trace()
        self.ndim    = totbins + len(self.ptsrc) + len(self.blob)     

### Some "advanced" features, which I hope to implement at some point
        self.bulk_centroid = [False]*len(bulkarc)    # Fit for a galaxy cluster centroid
        self.shoc_centroid = [False]*len(myshbins)   # Fit for shock centroid(s)
        self.psrc_centroid = [False]*len(self.ptsrc) # Fit for point source centroid(s) 
        self.blob_centroid = [False]*len(self.blob)  # Fit for the blob centroid
        self.testmode = testmode
### "Testing" values:
        if testmode == False:
            self.nwalkers= 30
            self.nsteps  = 2500
            self.burn_in = 500
        else:
            self.nwalkers= 20
            self.nsteps  = 200
            self.burn_in = 40


#########################################################################################
### + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ###
###+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +###
### + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ###
###+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +###
#########################################################################################
### / / / / / / / / / / / / / /                           \ \ \ \ \ \ \ \ \ \ \ \ \ \ ###
### / / / / /                                                               \ \ \ \ \ ###
###                                                                                   ###
###                         Done with Common Fitting Parameters.                      ###
###                  Below: fitting parameters for individual datasets                ###
###                  (Namely: (1) a mean level, if necessary                          ###
###                       and (2) a calibration offset                                ###
###                                                                                   ###
### \ \ \ \ \                                                               / / / / / ###
### \ \ \ \ \ \ \ \ \ \ \ \ \ \                           / / / / / / / / / / / / / / ###
#########################################################################################
### + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ###
###+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +###
### + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ###
###+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +###
#########################################################################################
            
class inst_fit_params:

    def __init__(self,inputs,maskrad=0):

        fwhm1,norm1,fwhm2,norm2,fwhm,smfw,freq,FoV = rdi.inst_params(inputs.instrument)
        if maskrad == 0: maskrad = FoV*1.5

        self.mn_lvl  = inputs.fitmnlvl     # Do not fit for a mean level
        self.pt_src  = inputs.fitptsrcs
        
        self.maskrad = maskrad
            
        n_add_params=0
        if self.mn_lvl == True : n_add_params+=1
        self.n_add_params = n_add_params
        self.calunc = inputs.calunc
        #self.caloff = 1.0
    
def convert_kpc_to_rad(bins,ang_dist):

    rads=[]
    for mybin in bins:
        myrad= (mybin/ang_dist).decompose().value # This will be in radians.
        rads.append(myrad)

    return rads
