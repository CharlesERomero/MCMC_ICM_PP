from astropy.io import fits
import astropy.units as u
import numpy as np

def create_xymap(xsize=512;ysize=512;mappixsize=2.0,xcentre=[],ycentre=[]):

    nx = xsize/mappixsize; ny = ysize/mappixsize
    ypix = mappixsize; npix = mappixsize
    
    if xcentre == []:
        xcentre = nx/2.0
    if ycentre == []:
        ycentre = ny/2.0

    x = np.outer(np.zeros(ny)+1.0,np.arange(0,xpix*(nx), xpix)- xpix* xcentre)   
    y = np.outer(np.arange(0,ypix*(ny),ypix)- ypix * ycentre, np.zeros(nx) + 1.0)
    return x,y

