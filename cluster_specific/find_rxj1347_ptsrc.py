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
import find_ptsrc as fp            # Will fit for location and shape    ###


########## Allow a few defaults / parameters to be set here: ##############
instrument='MUSTANG2'; name='rxj1347_wshock'; savedir='/home/romero/Results_Python'
tag='Re_12FWHM_v1_'; testmode='Full'
# Available testmodes: 'Test', 'Burn', 'Long', and 'Full' (the default).

################ Get parameters for fitting procedure: ####################
hk,dv,ifp= cv.get_struct_of_variables([instrument],name,savedir,testmode=testmode)

#for myinst in hk.instruments:
myinst = hk.instruments[0]
map = dv[myinst].maps.data
weights = dv[myinst].maps.masked_wts
pixs=dv[myinst].mapping.pixsize
ptsrc = fp.find_ptsrc_loc_shape(map,pixs,issmo=False,issnr=False,
                                wtmap=weights,instrument=myinst)

