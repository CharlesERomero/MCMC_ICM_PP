###########################################################################
### VERSION HISTORY                                                     ###
###                                                                     ###
### Written 20 June 2017 by Charles Romero                              ###
### Large revisions November 2017.                                      ###
### ---> Recast many variables/classes                                  ###
### ---> Allow multiple components to be fit                            ###
### ---> Allow multiple data sets to be fit (simultaneaously)           ###
### This revisions enter testing phase 11 Nov. 2017                     ###
###                                                                     ###
### PROJECTED REVISIONS:                                                ###
###                                                                     ###
### (1) Allow for point source CENTROID fitting.                        ###
### (2) Get PICKLE to work properly                                     ###
###                                                                     ###
###########################################################################
###                                                                     ###
### PRIVATE CODE DEPENDENCIES:                                          ###
### --> This has moved to the README.                                   ###
###                                                                     ###
###########################################################################
### For 2600 steps, 24 walkers (and a burn-in of 600), and 8 parameters ###
### the current runtime is very close to 12 hours. Conversely, for 300  ###
### steps, 26 walkers, and 8 parameters, the code takes somewhere close ###
### to 1.5 hours.                                                       ###
###                                                                     ###
###########################################################################
import numpy as np
############################### CER codes:  ###############################
import save_image as si
import plot_mcmc_results as pmr    # Each plot is its own function.     ###
import max_like_fitting as mlf     # Some setup & running of emcee      ###
import collect_variables as cv     # Calls many other modules           ###
import Azimuthal_Brightness_Profiles as ABP 
from os.path import expanduser
myhome = expanduser("~")

########## Allow a few defaults / parameters to be set here: ##############
instrument='MUSTANG2'; name='rxj1347_wshock'; savedir=myhome+'/Results_Python'
tag='Re_9FWHM_v1_'; testmode='Test'; nthreads=1
# Available testmodes: 'Test', 'Burn', 'Long', and 'Full' (the default).

################ Get parameters for fitting procedure: ####################
hk,dv,ifp= cv.get_struct_of_variables([instrument],name,savedir,testmode=testmode)

####### Create another class with variables for running the fits: #########
efv = mlf.emcee_fitting_vars(hk,dv,ifp,tag=tag,nthreads=nthreads)

#######          Do simple analyses on just the maps as-is        #########
mylist, aminarr, amaxarr = ABP.iter_two_slices(dv,hk,myformat='eps')

angmap=dv['MUSTANG2'].mapping.angmap
snrmap,hdr  = si.get_snrmap()
mywcs = si.get_wcs(hdr)
maskang = si.get_slice_mask(angmap,np.pi,6.0*np.pi/4.0)
si.plot_image_wcs(snrmap, maskang, mywcs,dpi=200,myfontsize=15,zoom=True,filename='My_Map',
                   plotmask=True,savedir=savedir,format='eps')

profsl = ABP.get_slices(dv,hk)
ABP.plot_all_slices(profsl,myformat='eps')
import pdb; pdb.set_trace()
