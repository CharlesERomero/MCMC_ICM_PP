import numpy as np
import astropy.units as u          # Install astropy

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

    
    nx,ny=map.shape
    ypix = pixsize.to("arcsec").value # Generally pixel sizes are the same...
    xpix = pixsize.to("arcsec").value # ""
    if xcentre == []:
        xcentre = nx/2.0
    if ycentre == []:
        ycentre = ny/2.0

    x = np.outer(np.arange(0,xpix*(nx), xpix)- xpix* xcentre, np.zeros(ny)+1.0)   
    y = np.outer(np.zeros(nx) + 1.0, np.arange(0,ypix*(ny),ypix)- ypix * ycentre)
    return x,y

def get_radial_map(map,pixsize,xcentre=[],ycentre=[]):

    x,y = get_xymap(map,pixsize,xcentre=xcentre,ycentre=ycentre)
    r = np.sqrt(x*x +y*y)

    return r
