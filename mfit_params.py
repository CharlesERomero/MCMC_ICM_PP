import time
import os
import find_ptsrc as fp
import retrieve_data_info as rdi
import astropy.units as u
import numpy as np

class fit_params:

    def __init__(self,dv,inputs,bins=6,path='/home/romero/Results_Python/',
                 testmode=False,autodetect=False):

        hk_inputs = inputs.files_and_priors(map_type=map_type)
        data_map = dv.maps.data; wtmap = dv.maps.wts; pixs = dv.mapping.pixsize
        instrument = hk_inputs.instrument; name=hk_inputs.name):

### What additional components do we care to fit for?
        self.mn_lvl = False      # Do not fit for a mean level
        self.ptsrc  = False      # Do not fit for a point source
        self.psfile = None       # No point source file
        self.shock = False       # Do not fit for a shock component
        self.shfile = None       # No shock file
        self.blob  = False       # Do not fit for a "blob" component
        self.blfile = None       # No blob file (user-defined)
### I CURRENTLY HAVE A BUG WITH mn_lvl. (20 June 2017)
        self.pprior = False      # Use a Planck Prior on Y_int
        self.ubmap = False       # Uniform bins (pressure is flat)
        self.bins = bins         # The number of (radial) bins
        self.real_data = True    #
        self.joint = False       # I never made it work jointly...
        self.path = path         #directory where figures are saved
        self.ndim    = bins      # Num. of free dimensions
### A few more fitting parameters:
        self.mask = True         # Use a mask? (Use only certain pixels?)
        self.maskrad = 60.0      # Radius within which to use data
### Some "advanced" features, which I hope to implement at some point
        self.gc_centroid = False # Fit for the galaxy cluster centroid
        self.sh_centroid = False # Fit for the shock centroid 
        self.ps_centroid = False # Fit for the point source centroid 
        self.bl_centroid = False # Fit for the blob centroid
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

### Gah...I should do this in a better way; but for a quick fix:
### For Abell 2146:
        self.maskrad = 150
            
### "Standard" values:
#        self.nwalkers= 50
#        self.nsteps  = 3000
#        self.burn_in = 500

        if autodetect == True:
            indloc,maxsnr,snr,peak = fp.find_ptsrc_fixed_shape(data_map,pixs,issmo=False,
                                  issnr=False,wtmap=wtmap, instrument=instrument)
            if maxsnr > 2.0:
                ### 14 July 2017 - switch to fixing the shape
                ptsrc,fpsfn=fp.make_ptsrc_model(instrument=instrument,name=name,
                             writefile=True,fixshape=True,normalize=True)
                self.ptsrc=True
                self.psfile=fpsfn
            else:
                print 'No updates performed'
        #raise Exception("This section is under development!")

        fwhm1,norm1,fwhm2,norm2,fwhm,smfw,freq,Fov = rdi.inst_params(instrument)

        if hk.fit_params.mask == True:
            arcbins = (dv.mapping.r_bins/dv.cluster.d_a)*(u.rad).to("arcsec")
            arcbins = arcbins.value
            print 'Max Arcbins: ',np.max(arcbins)
            self.maskrad = np.max(arcbins) + fwhm.value # We want at least another FWHM beyond the last edge.
        
#    import pdb; pdb.set_trace()

    n_add_params=0
    if fit_params.mn_lvl == True : n_add_params+=1
    if fit_params.ptsrc  == True : n_add_params+=1
    if fit_params.shock  == True : n_add_params+=1
    if fit_params.blob   == True : n_add_params+=1
    fit_params.ndim = fit_params.bins + n_add_params
    
class hk_out:

    def __init__(self,fit_params,instrument,target):
    
        mlstr='ML-NO_'; ppstr='PP-NO_'; dtype='Virt_'; ubstr='POWER_'
        pre_filename = "{}_Virt_".format(instrument)
        if fit_params.real_data == True:   
            dtype='Real_'; pre_filename = "{}_Real_".format(instrument)
        if fit_params.pprior == True: ppstr='PP-YES_'
        if fit_params.ubmap  == True: ubstr='UNIFORM_'         
#        if fit_params.mn_lvl == True:
#            fit_params.ndim=fit_params.ndim+1
#            mlstr='ML-YES_'
        
        nbstr=str(fit_params.bins)+"_B_"
        nststr=str(fit_params.nsteps)+"S_"
        nburnstr=str(fit_params.burn_in)+"B_"
        walkstr=str(fit_params.nwalkers)+"W"
        newfolder = time.strftime("%d%m_")+format(instrument)+'_'  #change as wanted
        self.runstr=nbstr+dtype+nststr+nburnstr+mlstr+ppstr+ubstr+walkstr
        self.nnfolder=newfolder+self.runstr
        inst_path = r'{}/{}'.format(fit_params.path,instrument)
        if not os.path.exists(inst_path):
            os.makedirs(inst_path)
        else:
            new_dir_query()        
        clus_path = r'{}/{}'.format(inst_path,target)
        if not os.path.exists(clus_path):
            os.makedirs(clus_path)
        else:
            new_dir_query()        

        self.newpath = clus_path
        self.prefilename=pre_filename


def choose_radial_bins(hk,dv,sector='nw',preset=False):

    fwhm1,norm1,fwhm2,norm2,fwhm,smfw,freq,Fov = rdi.inst_params(instrument)
    minstep=np.ceil(fwhm/5.0)*5.0 #Make steps be an integer of 5", by default
    binarr = (np.arange(nbins)+1) * minstep
#    inalphas=r_bins*0.0    # The variable needs to be defined...but if it is all zeros, I'll
#                           # overwrite the values.#
    if dv.cluster.name == 'abell_2146':
        arcbins*=1.3090075 # Perfectly matches SE model input...
#        arcbins*=1.5708075312038152 # Perfectly matches NW model input...

    ### I will want to add something more complicated - but let's
    ### Just start with this for now (20 June 2017)
    
    return binarr
    
