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

########## Allow a few defaults / parameters to be set here: ##############
instrument='MUSTANG2'; name='rxj1347_wshock'; savedir='/home/romero/Results_Python'
tag='Re_9FWHM_v0_'; testmode='Burn'
# Available testmodes: 'Test', 'Burn', 'Long', and 'Full' (the default).

################ Get parameters for fitting procedure: ####################
hk,dv,ifp= cv.get_struct_of_variables([instrument],name,savedir,testmode=testmode)

####### Create another class with variables for running the fits: #########
efv = mlf.emcee_fitting_vars(hk,dv,ifp,tag=tag)

####################### And now  run emcee! ###############################
sampler,t_mcmc = mlf.run_emcee(hk,dv,ifp,efv)
#cv.print_attributes(sampler)

####### Compile some results, and save (i.e. shelve) the results! #########
hdu = mlf.post_mcmc(sampler,t_mcmc,efv,hk,dv)    

######################### Plot the results ################################
pmr.plot_results(hk,dv,efv,sampler,hdu)
import pdb; pdb.set_trace()

##################### Load and Plot the results ###########################
#my_dv,my_hk,my_efv = mlf.unpickle()
#pmr.plot_res_no_sampler(my_hk,my_dv,my_efv,overlay='input')
