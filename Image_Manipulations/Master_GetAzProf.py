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
import Azimuthal_Brightness_Profiles as ABP 
from os.path import expanduser
from astropy.io import fits
import mapping_modules as mm       # Creates radius + pressure arrays.
import retrieve_data_info as rdi
from astropy.wcs import WCS
from astropy.coordinates import Angle
import file_presets as fp

myhome = expanduser("~")

########## Allow a few defaults / parameters to be set here: ##############
#datadir='/home/data/MUSTANG2/MINKASI/HighZ/'
#fitsfile=datadir+'2XMMJ0830+5241_precon_2_arcsec_pass_4.fits'
#ra_deg  = Angle('8h30m25.875s').to('deg')
#dec_deg = Angle('52d41m33.96s').to('deg')

src              = '2XMM'
datadir,fitsfile,ra_deg,dec_deg = fp.get_info(src)

data_map, header = fits.getdata(fitsfile, header=True)
w                = WCS(fitsfile)
x0,y0            = w.wcs_world2pix(ra_deg,dec_deg,0)
ras,decs,pixs    = rdi.astro_from_hdr(header)
xymap            = mm.get_xymap(data_map,pixs,xcentre=x0,ycentre=y0) # In arcseconds
rads, prof       = ABP.get_az_profile(data_map, xymap, 0.0, 2.0*np.pi-1.0e5, geoparams=[0,0,0,1,1,1,0,0])
slices           = []
mybins           = np.arange(0.0,180.0,5.0)
binres           = ABP.radial_bin(rads, prof,10,rmax=60.0,bins=mybins,minangle=0.0,maxangle=2.0*np.pi-1.0e5)
ABP.plot_one_slice(binres,myformat='eps',savedir=datadir,target=src)
slices.append(binres)


#######          Do simple analyses on just the maps as-is        #########
#mylist, aminarr, amaxarr = ABP.iter_two_slices(dv,hk,myformat='eps')

#angmap=dv['MUSTANG2'].mapping.angmap
#snrmap,hdr  = si.get_snrmap()
#mywcs = si.get_wcs(hdr)
#maskang = si.get_slice_mask(angmap,np.pi,6.0*np.pi/4.0)
#si.plot_image_wcs(snrmap, maskang, mywcs,dpi=200,myfontsize=15,zoom=True,filename='My_Map',
#                   plotmask=True,savedir=savedir,format='eps')

#profsl = ABP.get_slices(dv,hk)
#ABP.plot_all_slices(profsl,myformat='eps')
#import pdb; pdb.set_trace()
