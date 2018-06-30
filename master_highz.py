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
############################### CER codes:  ###############################
import plot_mcmc_results as pmr    # Each plot is its own function.     ###
import max_like_fitting as mlf     # Some setup & running of emcee      ###
import collect_variables as cv     # Calls many other modules           ### 
import cProfile, sys, pstats       # I may not need pstats
#import re
from os.path import expanduser
myhome = expanduser("~")

########## Allow a few defaults / parameters to be set here: ##############
instrument='MUSTANG2'; savedir=myhome+'/Results_Python'
testmode='Full'; nthreads=1

#name='2xmm'; reduc='PCA'
name='idcs'; reduc='PCA'
#name='rdcs0910'; reduc='PCA'
#name='rxj1053'; reduc='CMCORR'
# Available testmodes: 'Test', 'Burn', 'Long', and 'Full' (the default).

tag=reduc+'_ptsrcs_'
################ Get parameters for fitting procedure: ####################
hk,dv,ifp= cv.get_struct_of_variables([instrument],name,savedir,testmode=testmode,
                                      reduc=reduc)

####### Create another class with variables for running the fits: #########
efv = mlf.emcee_fitting_vars(hk,dv,ifp,tag=tag,nthreads=nthreads)

####################### And now  run emcee! ###############################
pr = cProfile.Profile()
pr.enable()
sampler,t_mcmc = mlf.run_emcee(hk,dv,ifp,efv)
pr.disable()

####### Compile some results, and save (i.e. shelve) the results! #########
hdu = mlf.post_mcmc(sampler,t_mcmc,efv,hk,dv,ifp)    

######################### Plot the results ################################
pmr.plot_results(hk,dv,efv,sampler,hdu,ifp)
print 'type: mlf = reload(mlf)'
import pdb; pdb.set_trace()
#mlf.post_check(sampler,t_mcmc,efv,hk,dv,ifp)

#p = pstats.Stats('pr')
prof_out = myhome+'/PythonProfilerResults.txt'
sys.stdout = open(prof_out, 'w')
pr.print_stats()
#close(prof_out)
sys.stdout = sys.__stdout__

#pr.sort_stats('time').print_stats(10)

##################### Load and Plot the results ###########################
#my_dv,my_hk,my_efv = mlf.unpickle()
#pmr.plot_res_no_sampler(my_hk,my_dv,my_efv,overlay='input')
