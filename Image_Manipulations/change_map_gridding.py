import FITS_tools.hcongrid as fth
import numpy as np
import astropy.coordinates as apc  # http://www.astropy.org/
from astropy import wcs
from astropy.io import fits
from astropy import units as u
from astropy.time import Time as APT
from astropy.visualization.wcsaxes import WCSAxes
import os

def convert_image_with_hdr(infile):

    hdulist = fits.open(infits)
    header = Coverage.w.to_header()
    wgm = wcs.WCS(hdulist[0].header)
    wgm.naxes = 2
    hgm = wgm.to_header()
    hgm['wcsaxes']=2
    hgm.remove('crpix3'); hgm.remove('crpix4'); hgm.remove('cdelt3'); hgm.remove('cdelt4')
    hgm.remove('ctype3'); hgm.remove('ctype4'); hgm.remove('crval3'); hgm.remove('crval4')
    hgm.remove('cunit4')
    #hgm.update(('WCSAXES',2,'Number of coordinate axes'))
    inmap = hdulist[0].data; goodmap=inmap.squeeze(); values=goodmap.flatten()
    sordid = np.sort(values); mylen=len(sordid)
    ContLvls = [ sordid[int(mylen*c)] for c in inContours ]
    hgm.append(('NAXIS1',goodmap.shape[0],'Number of pixels along dimension1'))
    hgm.append(('NAXIS2',goodmap.shape[1],'Number of pixels along dimension2'))
    header.append(('NAXIS1',Coverage.weight1mm.shape[0],'Number of pixels along dimension1'))
    header.append(('NAXIS2',Coverage.weight1mm.shape[1],'Number of pixels along dimension2'))
    indata = fth.hcongrid(goodmap, hgm, header, preserve_bad_pixels=False)      
    #pdb.set_trace()

def rescale_image(infile,pixsize=2.0*u.arcsec,powersoftwo=True):

    hdulist = fits.open(infile); myhdu=[]
    
    for i,hdu in enumerate(hdulist):

        header  = hdu.header
        wgm = wcs.WCS(hdu.header)
        inmap = hdu.data; goodmap=inmap.squeeze(); values=goodmap.flatten()
        wgmcopy = wgm
        
        pixs_orig = np.sqrt(np.abs(np.linalg.det(np.array(wgm.wcs.cd)))) # Original pixel size
        pix_factor = pixsize.to(wgmcopy.wcs.cunit[0]) / (pixs_orig*u.Unit(wgm.wcs.cunit[0]))

        ### Somehow this is wrong? Why?
        #if hdu.name == 'wtmap':
        #    goodmap*= pix_factor**2   # Scale the weightmap appropriately
        
        wgmcopy.wcs.cd = [[-pixsize.to(wgmcopy.wcs.cunit[0]).value,-0.0],
                          [0.0,pixsize.to(wgmcopy.wcs.cunit[1]).value]]
        wgmcopy.wcs.cdelt = [pixsize.to(wgmcopy.wcs.cunit[0]).value,
                             pixsize.to(wgmcopy.wcs.cunit[0]).value]
        
        origshape = goodmap.shape
        newshape=tuple([round(length / pix_factor.value) for length in origshape])
        
        if powersoftwo == True:
            nx,ny = newshape
            lnx,lny = np.log(nx)/np.log(2), np.log(ny)/np.log(2)
            ilnx,ilny = int(lnx),int(lny)
            nnx,nny = 2**ilnx , 2**ilny
            newshape = (nnx,nny)
            
        crpix = [newshape[0]/2.0,newshape[1]/2.0]
        wgmcopy.wcs.crpix = crpix
        hgm = wgmcopy.to_header()
        #hgm.append(('PC1_2',0,'Coordinate transformation matrix element'))
        #hgm.append(('PC2_1',0,'Coordinate transformation matrix element'))
        hgm.append(('NAXIS1',newshape[0],'Number of pixels along dimension1'))
        hgm.append(('NAXIS2',newshape[1],'Number of pixels along dimension2'))
        copy_missing_header_keys(header,hgm)       # Update the header to include desired fields

        #for myline in hgm:
        #    print myline, hgm[myline] 
        #import pdb;pdb.set_trace()
        
        outmap = fth.hcongrid(goodmap, header, hgm, preserve_bad_pixels=False)

        if i == 0:
            thishdu = fits.PrimaryHDU(outmap,header=hgm)
            myhdu.append(thishdu)
        else:
            hdu1 = fits.ImageHDU(outmap)
            hdu1.header = hgm
            hdu1.verify('fix')            
            myhdu.append(hdu1)
            
    return myhdu

def write_hdu_to_fits(hdu,fullpath):

    hdulist = fits.HDUList(hdu)
    hdulist.writeto(fullpath,overwrite=True)

def example_rewrite():

    ############################################################
    ### COMPLETE SCRATCH:
    #mydir     = '/home/data/MUSTANG2/AGBT17_Products/RXJ1347/post2018/'
    #myfile    = 'grid_pca7_f_Low0.080__map.fits'
    #myfile    = "grid_pca7_f_Low0.080__noise.fits"
    
    ############################################################
    #mydir     = '/home/data/MUSTANG2/AGBT17_Products/HSC/'
    #myfile    = 'pca3_lowf0.8_roi1.5_Kelvin_map.fits'
    
    ############################################################
    #mydir     = '/home/data/MUSTANG2/AGBT17_Products/RDCS0910/'
    #myfile    = 'Kelvin_RDCS0910_2aspix_pca3_0f05_bugfixed_map_iter1.fits'
    #filename  = 'RDCS0910_PCA3_fixed_map.fits'
    
    ############################################################
    mydir     = '/home/data/MUSTANG2/AGBT17_Products/2XMM/'
    myfile    = 'Kelvin_2XMMJ0830+5241_2aspix_pca3_0f05_v1_map_iter1.fits'
    filename  = '2XMM_PCA3_fixed_map.fits'
     
    infile    = mydir+myfile
    myhdu     = rescale_image(infile,pixsize=2.0*u.arcsec)
    fullpath  = os.path.join(mydir,filename)
    write_hdu_to_fits(myhdu,fullpath)

def copy_missing_header_keys(full_hdr,inc_hdr):
    """
    Takes a full header (full_hdr) and copies any keys not found in the incomplete 
    header (inc_hdr).
    """

    ign_keys = ['CD1_1','CD1_2','CD2_1','CD2_2']           # Set of keys to ignore
    exc_keys = set(inc_hdr.keys()).union(ign_keys)         # Set of keys to exclude
    diffkeys = set(full_hdr.keys()).difference(exc_keys)   # The difference of the sets
    print '-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#'
    for key in diffkeys:
        #print type(full_hdr[key]), full_hdr[key]
        if hasattr(full_hdr[key],'__len__') and not(type(full_hdr[key]) is str):
            for mycomment in full_hdr[key]:
                inc_hdr.append((key,mycomment))
        else:
            print key,full_hdr[key]
            inc_hdr.append((key,full_hdr[key]))


def add_snr_ext():

    print 'howdy do'
