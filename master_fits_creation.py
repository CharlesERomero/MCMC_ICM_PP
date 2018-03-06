from astropy.io import fits
import max_like_fitting as mlf
import astropy.units as u
import ellipsoidal_shells as es
import instrument_processing as ip
import numpy as np
import os

mydir='/home/romero/Results_Python/Combined/rxj1347/'
fullpath=mydir+'Combined_Real_Test_Run_Residual.fits'
hdulist = fits.open(fullpath)

mydata= hdulist[0].data
myhdr = hdulist[0].header
hdu0  = fits.PrimaryHDU(mydata,header=myhdr)

hdu1 = fits.ImageHDU(mydata)
hdu1.header = myhdr
hdu1.name='Test'

myhdu = [hdu0,hdu1]
mylist = fits.HDUList(myhdu)

myhdu.info()

################################################

hdu3 = fits.ImageHDU(modelsky[mykey])
hdu3.header = dv[mykey].mapping.w.to_header()
hdu3.name = title+str(count)
hdu3.header.append(("Title",title))
hdu3.header.append(("Target",hk.hk_ins.name))
hdu3.header.append(("XTENSION","What Mate"))
hdu3.header.append(("SIMPLE","T")) 
hdu3.verify('fix')
myhdu.append(hdu2)

nfilename=tstr+"TEST_TEST_TEST.fits"
nfullpath = os.path.join(hk.hk_outs.newpath,hk.hk_outs.prefilename+nfilename)
newhdulist.writeto(fullpath,overwrite=True)

testread = fits.open(fullpath)
