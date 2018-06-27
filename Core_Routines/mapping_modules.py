import numpy as np
import astropy.units as u          # Install astropy
from astropy.coordinates import Angle
from astropy.coordinates import SkyCoord
from astropy import wcs
from astropy.io import fits

def get_xymap(map,pixsize,xcentre=[],ycentre=[]):
    """
    Returns a map of X and Y offsets (from the center) in arcseconds.

    INPUTS:
    -------
    map      - a 2D array for which you want to construct the xymap
    pixsize  - a quantity (with units of an angle)
    xcentre  - The number of the pixel that marks the X-centre of the map
    ycentre  - The number of the pixel that marks the Y-centre of the map

    """

    ny,nx=map.shape
    ypix = pixsize.to("arcsec").value # Generally pixel sizes are the same...
    xpix = pixsize.to("arcsec").value # ""
    if xcentre == []:
        xcentre = nx/2.0
    if ycentre == []:
        ycentre = ny/2.0

    x = np.outer(np.zeros(ny)+1.0,np.arange(0,xpix*(nx), xpix)- xpix* xcentre)   
    y = np.outer(np.arange(0,ypix*(ny),ypix)- ypix * ycentre, np.zeros(nx) + 1.0)
    return x,y

def get_radial_map(map,pixsize,xcentre=[],ycentre=[]):

    x,y = get_xymap(map,pixsize,xcentre=xcentre,ycentre=ycentre)
    r = np.sqrt(x*x +y*y)

    return r

def create_xymap(xsize=512,ysize=512,mappixsize=2.0,xcentre=[],ycentre=[]):

    nx = int(xsize/mappixsize); ny = int(ysize/mappixsize)
    ypix = mappixsize; xpix = mappixsize
    
    if xcentre == []:
        xcentre = nx/2.0
    if ycentre == []:
        ycentre = ny/2.0

    x = np.outer(np.zeros(ny)+1.0,np.arange(0,xpix*(nx), xpix)- xpix* xcentre)   
    y = np.outer(np.arange(0,ypix*(ny),ypix)- ypix * ycentre, np.zeros(nx) + 1.0)
    return x,y

def create_astro(mymap,crpix=[], crval=[Angle('13h00m0.0s'), Angle('20d00m0.0s')],
                 cdelt=np.array([1.0,1.0]),mappixsize=2.0):

    w = wcs.WCS(naxis=2)

    if crpix == []:
        crpix = [mymap.shape[0]/2.0, mymap.shape[0]/2.0]

    degcrval = [ (crval[0].to('deg')).value, (crval[1].to('deg')).value ]
        
    w.wcs.crpix = crpix
    w.wcs.cdelt = cdelt
    w.wcs.crval = degcrval
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.cd    = [[-mappixsize/3600.0, -0.0],[0.0,mappixsize/3600.0]]

    return w
    
def make_hdu(mymap, w,img2=[], img3=[]):

    myhdu = []
    
    hdr = w.to_header()
    newshape = mymap.shape
    hdr.append(('NAXIS1',newshape[0],'Number of pixels along dimension1'))
    hdr.append(('NAXIS2',newshape[1],'Number of pixels along dimension2'))
    hdu = fits.PrimaryHDU(mymap,header=hdr)
    hdu.name = 'Compton y'
    hdu.verify('fix')
    myhdu.append(hdu)

    if img2 != []:
        #newshape = mymap.shape
        #hdr.append(('NAXIS1',newshape[0],'Number of pixels along dimension1'))
        #hdr.append(('NAXIS2',newshape[1],'Number of pixels along dimension2'))
        #hdu = fits.PrimaryHDU(img2,header=hdr)
        
        hdu2 = fits.ImageHDU(img2)
        hdu2.header = hdr
        hdu2.name = 'Filtered y map'
        hdu2.header.append(("Title",'Filtered y map'))
        hdu2.header.append(("Target",'Unkown'))
        hdu2.header.append(("XTENSION","What Mate"))
        hdu2.header.append(("SIMPLE","T")) 
        hdu2.verify('fix')
        myhdu.append(hdu2)
       
    if img3 != []:
        #newshape = mymap.shape
        #hdr.append(('NAXIS1',newshape[0],'Number of pixels along dimension1'))
        #hdr.append(('NAXIS2',newshape[1],'Number of pixels along dimension2'))
        #hdu = fits.PrimaryHDU(img2,header=hdr)
        
        hdu3 = fits.ImageHDU(img3)
        hdu3.header = hdr
        hdu3.name = 'Filtered unit map'
        hdu3.header.append(("Title",'Filtered unit map'))
        hdu3.header.append(("Target",'Unkown'))
        hdu3.header.append(("XTENSION","What Mate"))
        hdu3.header.append(("SIMPLE","T")) 
        hdu3.verify('fix')
        myhdu.append(hdu3)

    return myhdu

def write_hdu_to_fits(hdu,fullpath):

    hdulist = fits.HDUList(hdu)
    hdulist.writeto(fullpath,overwrite=True)

class xfer_fxn:

    def __init__(self,instrument="MUSTANG2"):
    
        if instrument == "MUSTANG2":
            m2dir            = '/home/data/MUSTANG2/'
            self.tabfile     = m2dir+"AGBT17_Products/pca7_f0.09_onHSC_2.txt"
            self.tabcomments = '#'
            self.tabformat   = 'ascii'
            self.tabdims     = '1D'
            self.tabextend   = True    # Do we need to extent to higher k numbers?
            self.instrument  = instrument

        
