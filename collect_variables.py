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

def get_struct_of_variables(instruments,name,path='/home/romero/Results_Python/',
                            testmode=False,map_file='nw_model'):

    """
    This module is the overarching module to collect variables.
    ------------------------------------
    INPUTS:
    ----
    instrument     - The instrument which took the data (e.g. "MUSTANG-2")
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
    priors = input_struct.priors() # These are not fitting priors, but prior "known" (fixed) quantities.
    dv = {}
    for instrument in instruments:
        inputs = input_struct.files(instrument=instrument,map_file='nw_model')
        ### Now, we can get the input data variables
        dv[instrument] = data_vars(inputs,priors)
        
    ### Now, as the fitting parameters depends on the input file, let's get those:
    fit_params = mfp(dv,inputs,path=path,testmode=testmode)
  
    hk = housekeeping(dv,inputs,map_file=map_file)
    ### OK, now...can we collect other "defaults"
    ### Let's have a quick look at the data and see if we should be fitting
    ### a point source. If so, we'll update the fit parameters.
    ufp = mfp.update_fit_params(dv,hk)
    hk.fit_params = ufp # Maybe I want to update directly...
    
    return hk,dv
    
class housekeeping:

    def __init__(self,priors,fit_params,instruments):
          
        self.log = logbook()
        self.hk_ins = inputs
        self.fit_params = fit_params
            
        ### I should add variables regarding logs here

class data_vars:

    def __init__(self,inputs,priors,real=True,beta=0,betaz=0):
        """
        I need to update this so that it reads in all the data (files) and then I can discard the
        input file information.
        """

        self.maps  = rdi.maps(inputs)
        self.astro = rdi.astrometry(self.maps.header)
        self.xfer  = rdi.get_xfer(inputs)
        self.xferdims = inputs.tabdims
        tSZ,kSZ = rdi.get_sz_bp_conversions(priors.Tx,inputs.instrument,array="2",inter=False,
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
