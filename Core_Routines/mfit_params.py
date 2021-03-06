import time
import os
import find_ptsrc as fp
import retrieve_data_info as rdi
import astropy.units as u
import numpy as np
from os.path import expanduser
myhome = expanduser("~")

class common_fit_params:

    def __init__(self,bulk,shbins=[6],path=myhome+'/Results_Python/',
                 bulkgeo=[],bulknarm=[False],bulkcen=[],bulkalp=[],
                 shockgeo=[],shocknarm=[True],shockalp=[],shockfin=[True],
                 ptsrcs=[],psfwhm=[],blobs=[],fstemps=[False],
                 minmax=np.array([2.0,100.0])*u.arcsec,blobcens=[[0],[0]],
                 cluster=None,testmode='Test',autodetect=False,
                 fitblobcen=[False],fitptcen=[False],fitshockcen=[False]):

    ### Pre-August 9 2018.
        #def __init__(self,bins=[6],shbins=[6],path=myhome+'/Results_Python/',
        #         bulkgeo=[],bulknarm=[False],bulkcen=[],bulkalp=[],
        #         shockgeo=[],shocknarm=[True],shockalp=[],shockfin=[True],
        #         ptsrcs=[],psfwhm=[],blobs=[],fbtemps=[False],fstemps=[False],
        #         minmax=np.array([2.0,100.0])*u.arcsec,blobcens=[[0],[0]],
        #         cluster=None,testmode='Test',autodetect=False,fitbulkcen=[False],
        #         fitblobcen=[False],fitptcen=[False],fitshockcen=[False],
        #         fitbulkgeo=[False]):

    ##################################################################################
    #####      Let's first prepare the bins for our bulk and shock components    #####
    ##################################################################################

        #if type(bulk.fbtemps) == type(None):
        #    fbtemps=[False]
        fbtemps     = bulk.fbtemps
        bins        = bulk.bins
        fitbulkgeo  = bulk.fit_geo
        fitbulkcen  = bulk.fit_cen
        model       = bulk.model
        
        bulkarc = []
        mygeo   = []
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
        for scount,mybins in enumerate(shbins):
            
            if shockfin[scount] == True:
                thesebins = mybins[:-1]
            else:
                thesebins = mybins
                
            if mybins.unit.is_equivalent("kpc"):
                mybins = convert_kpc_to_rad(mybins,cluster.d_a)
            elif mybins.unit.is_equivalent("rad"):
                mybins = mybins.to("rad").value
            else:
                raise Exception("Your shock bins have no units; Add units in your INFO file.")
            
            ### Create lists of the bin positions and power-law slope (default=0).
            myshbins.append(mybins)
            if anotset:
                shockalp.append(np.zeros(len(thesebins)))
            totbins+=len(thesebins)
        myshbins = np.array(myshbins)    # This should work correctly now (05 Mar 2018)
                                    
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
        self.model    = model      # one of ['NP','GNFW','BETA'] pressure profile models

        self.ptsrc    = ptsrcs     # Provide a list of centroids.
        self.psfwhm   = psfwhm     # list of FWHM...if not truly point-like
        self.shockalp = shockalp   # Provide a list of shock log pressure slope
        self.shockgeo = shockgeo   # Provide a list of shock geometries
        self.shockbin = myshbins   # How many bins to use; OR, an array of bins positions.
        self.shocknarm= shocknarm  # set
        self.shockfin = shockfin   # Do we do a finite integration (out to last bin)? (vs. infinite)
        self.shockfix = True       # Fix alphas for shock? By default, yes.
        self.fstemps  = fstemps    # Fit for shock profile temperatures (if X-ray data is present)
        self.trimOut  = True       # Trim outer (shock) shell.
        
        ### What additional components do we care to fit for?
        #### IMPORTANT: ###
        # If you want to specify the bin locations, you should do so within the file
        # that spits out "priors"
        
        self.blob     = blobs       # List of blob components
        self.bras     = blobcens[0] # List of Right Ascensions
        self.bdecs    = blobcens[1] # List of Declinations

### I think I will want to add masks for bulk and shock components. Or at least
### the option to have them... I'm not sure how I want to do this yet, so I guess I will
### hold off on this for now (08 Nov 2017).
        
        self.pprior = False      # Use a Planck Prior on Y_int
        self.ubmap = False       # Uniform bins (pressure is flat)
        self.real_data = True    # 
        self.joint = False       # I never made it work jointly...
        self.path = path         # directory where figures are saved
        #### Num. of free dimensions
        #import pdb;pdb.set_trace()
        print totbins,len(self.ptsrc), len(self.blob)
        ### I will need to add the blob stuff outside of this routine (20 Feb 2018)
        nblobbins = 0
        for blobpars in self.blob: nblobbins = nblobbins + len(blobpars)
        self.ndim    = totbins + nblobbins #+ len(self.ptsrc)

### Some "advanced" features, which I hope to implement at some point
        self.bulk_centroid = fitbulkcen  * len(bulkarc)    # Fit for a galaxy cluster centroid
        self.shoc_centroid = fitshockcen * len(myshbins)   # Fit for shock centroid(s)
        self.psrc_centroid = fitptcen    * len(self.ptsrc) # Fit for point source centroid(s) 
        self.blob_centroid = fitblobcen  * len(self.blob)  # Fit for the blob centroid
        self.bulk_geometry = fitbulkgeo  * len(bulkarc)    # Fit for a galaxy cluster centroid

        ncen  = 0   # Number of fitted centroids
        ncen = ncen+2*len(bulkarc)    if fitbulkcen  == [True] else ncen+0
        ncen = ncen+2*len(myshbins)   if fitshockcen == [True] else ncen+0
        ncen = ncen+2*len(self.ptsrc) if fitptcen    == [True] else ncen+0
        ### By default, I'll fit for the blob center...
        #ncen = ncen+2*len(self.blob)  if fitblobcen  == [True] else ncen+0

        self.ndim += ncen

        nextra = 0
        nextra = nextra+2*len(bulkarc)    if fitbulkgeo  == [True] else nextra+0
        
        self.ndim += nextra
        totaldim = self.ndim + len(self.ptsrc) + 1
        
        self.testmode = testmode      
### "Testing" values:
        ### Here is the longest I would think to do:
        if testmode == 'XLong':
            self.nwalkers= int(totaldim *2)*2
            self.nsteps  = 20000
            self.burn_in = 3000
        elif testmode == 'Long':
            self.nwalkers= int(totaldim *1.5)*2
            self.nsteps  = 5000
            self.burn_in = 1500
        ### And here is a test mode which really just verifies that the code will run.
        elif testmode == 'Test':
            self.nwalkers= totaldim *2 + 10   # Whatever...as it needs to be.
            self.nsteps  = 25                  # Seriously need to push to small numbers
            self.burn_in = 5                   # Really small numbers
        ### The following is designed to show the burn-in steps:
        elif testmode == 'Burn':             
            self.nwalkers= totaldim *2 + 10
            self.nsteps  = 500
            self.burn_in = 100
        ### This can be the "standard" ("Full") run:
        else:
            self.nwalkers= int(totaldim *1.5)*2
            self.nsteps  = 2500
            self.burn_in = 500

        print('You are set to use the following MCMC parameters under mode '+testmode+':')
        print('N_Walkers: ',self.nwalkers,' N_steps: ',self.nsteps,' Burn_in: ',self.burn_in)


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

    def __init__(self,inputs,ptsrcs,blobs,instrument,maskrad=0):

        fwhm1,norm1,fwhm2,norm2,fwhm,smfw,freq,FoV = rdi.inst_params(inputs.instrument)
        if maskrad == 0: maskrad = FoV*1.5
        self.mn_lvl  = inputs.fitmnlvl          # Do not fit for a mean level
        self.pt_src  = inputs.fitptsrcs         # Should be True or False...
        self.fitblob = blobs.dofit[instrument]  # Should be True or False...
        self.maskrad = maskrad                  #
        
        self.prior=[myprior for myprior in ptsrcs.prior[instrument]]
        self.priorunc=[mypriorunc for mypriorunc in ptsrcs.priorunc[instrument]]
            
        n_add_params=0
        if self.mn_lvl == True : n_add_params+=1
        self.n_add_params = n_add_params
        self.calunc = inputs.calunc

        print('Your mean level and point source fitting values are: ')
        print('mn_lvl: ',self.mn_lvl)
        print('pt_src: ',self.pt_src)
        print('<<<<{{{{(((([[[[-_-_-_-_-_-_-_-_-]]]]))))}}}}>>>>')
        #self.caloff = 1.0
    
def convert_kpc_to_rad(bins,ang_dist):

    rads=[]
    for mybin in bins:
        myrad= (mybin/ang_dist).decompose().value # This will be in radians.
        rads.append(myrad)

    return rads
