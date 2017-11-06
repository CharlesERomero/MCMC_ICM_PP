import my_astro_defaults as mad
import my_cluster_defaults as mcd
import cosmolopy.distance as cd
import astropy.units as u
from astropy.coordinates import Angle
import numpy as np
import astropy.constants as const
import retrieve_data_info as rdi   # Finds and loads data files (hard-coded)
import pressure_profiles as pp     # Creates radius + pressure arrays.
import mfit_params as mfp          # Loads parameters for fitting procedure
import importlib
import time

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

def get_struct_of_variables(instrument,name,path='/home/romero/Results_Python/',
                            testmode=False,sanche=False,map_type='nw_model'):

    ### Get my ducks in a row first!
    hk = housekeeping(instrument,name,path=path,testmode=testmode,sanche=sanche,map_type=map_type)
    ### OK, now...can we collect other "defaults"
    dv = data_vars(hk,sanche=sanche)
    ### Let's have a quick look at the data and see if we should be fitting
    ### a point source. If so, we'll update the fit parameters.
    ufp = mfp.update_fit_params(dv,hk)
    hk.fit_params = ufp # Maybe I want to update directly...
    
    return hk,dv
    
class housekeeping:

    def __init__(self,instrument,name,path=None,testmode=None,sanche=False,map_type='nw_model'):

#        if name =="a2146":
#            bins=3
#        else:
#            bins=6
            
            
        inputs = importlib.import_module(name+'_'+instrument+'_info')
        self.log = logbook()
        self.hk_ins = inputs.files_and_priors(map_type=map_type)
        # South and North CHEck. 
        if sanche == True:
            bins = len(self.hk_ins.rads_nw)
            fit_params = mfp.fit_params(bins=bins,path=path,testmode=testmode)
        else:
            fit_params = mfp.fit_params(path=path,testmode=testmode)
        self.hk_outs = mfp.hk_out(fit_params,instrument,name)
        self.fit_params = fit_params
            
        ### I should add variables regarding logs here

class data_vars:

    def __init__(self,hk,real=True,sanche=False,beta=0,betaz=0):

        self.maps = rdi.maps(hk)
        tSZ,kSZ = rdi.get_sz_bp_conversions(hk.hk_ins.Tx,hk.hk_ins.instrument,array="2",inter=False,
                                        beta=beta,betaz=betaz,rel=True)
        self.av = mad.all_astro()  # av = astro vars
        self.cluster=mcd.cluster(hk)
### The optional input N (defaults to 120) may be something to change.../make more visible in "mapping": 
        self.mapping = mcd.mapping(self.cluster,hk,hk.fit_params,self.av.mycosmo,
                                   tSZ,kSZ,hk.hk_ins.instrument,sanche=sanche)
        
def print_attributes(myclass):

    attrs = vars(myclass)
    print ', '.join("%s" % item for item in attrs)

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
