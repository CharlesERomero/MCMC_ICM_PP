###########################################################################
### VERSION HISTORY                                                     ###
###                                                                     ###
### Written 20 June 2017 by Charles Romero                              ###
###                                                                     ###
### PROJECTED REVISIONS:                                                ###
###                                                                     ###
### (1) Allow for point source CENTROID fitting.                        ###
### (2) Have the PICKLE variable save working.                          ###
###                                                                     ###
###########################################################################
###                                                                     ###
### PRIVATE CODE DEPENDENCIES:                                          ###
### --> This has moved to the README.                                   ###
###                                                                     ###
###########################################################################
###                                                                     ###
###########################################################################
############################### CER codes:  ###############################
import collect_variables as cv     # Calls many other modules           ### 
import max_like_fitting as mlf     # Some setup & running of emcee      ###
import plot_mcmc_results as pmr    # Each plot is its own function.     ###

########## Allow a few defaults / parameters to be set here: ##############
instrument='MUSTANG2'; name='a2146'; savedir='/home/romero/Results_Python'
tag='Re_'; sanche=False; testmode=False; map_type = 'se_model'

if map_type == 'se_model':
    nw=False; sector='se'
if map_type == 'nw_model':
    nw=True; sector='nw'

################ Get parameters for fitting procedure: ####################
hk,dv= cv.get_struct_of_variables(instrument,name,savedir,testmode=testmode,
                                  sanche=sanche,map_type=map_type)

####### Create another class with variables for running the fits: #########
efv = mlf.emcee_fitting_vars(hk,dv,tag=tag,sanche=sanche,sector=sector)
import pdb; pdb.set_trace()

####################### And now  run emcee! ###############################
sampler,t_mcmc = mlf.run_emcee(hk,dv,efv)
cv.print_attributes(sampler)

####### Compile some results, and save (i.e. shelve) the results! #########
mlf.post_mcmc(sampler,t_mcmc,efv,hk,dv)    

######################### Plot the results ################################
pmr.plot_results(hk,dv,efv,sampler,overlay='input')
import pdb; pdb.set_trace()

##################### Load and Plot the results ###########################
#my_dv,my_hk,my_efv = mlf.unpickle()
#pmr.plot_res_no_sampler(my_hk,my_dv,my_efv,overlay='input')
