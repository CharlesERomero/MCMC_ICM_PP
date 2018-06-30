import create_gNFW_maps as cgm
import astropy.constants as const
import numpy as np
import my_astro_defaults as mad


def calculate_for_z_M(z=1.2,mym500=3.0,myunits='Jy'):

    ### Defaults are for Stefano's cluster
    cosmo,mycosmo = mad.get_cosmology("Concordance")
    M500          = mym500 * 10**14* const.M_sun     
    R500, P500    = cgm.R500_P500_from_M500_z(M500, z, mycosmo)
    yprof, inrad  = cgm.compton_y_profile_from_m500_z(M500, z, mycosmo)
    myprof, radii = cgm.pressure_profile_from_m500_z(M500, z, mycosmo)
    Rmax          = R500
    #Ycyl    = cgm.Y_cyl(yprof, inrad, Rmax,z=z)
    #Ysphere = cgm.Y_sphere(myprof, radii, Rmax,z=z)
    d_a     = cgm.get_d_a(z)
    
    theta_range = ((inrad / d_a).decompose()).value
    msolare14 = (M500 / (10**14* const.M_sun)).value
    myfilename = 'gNFW_map_z'+"{:4.2f}".format(z)+ \
                 '_M{:.2f}'.format(msolare14) + 'E14.fits'
    prefile='gNFW_profile_z'+"{:4.2f}".format(z)+ \
                 '_M{:.2f}'.format(msolare14) + 'E14'
    
    cgm.create_gNFW_map(yprof, theta_range, xsize=512, ysize=512,
                    mappixsize=2.0,savedir=None,
                        filename=myfilename,instrument='MUSTANG2',
                        prefile=prefile,units=myunits)

    
def make_grid():

    for i in range(5):
        mym500 = 2.0**i 
        for j in range(10):
            z             = (j+1.0)*0.2
            calculate_for_z_M(z,mym500)


make_grid()
#calculate_for_z_M()
