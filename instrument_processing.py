import numpy as np
import scipy.ndimage
import scipy.signal
import image_filtering as imf
import retrieve_data_info as rdi
from numpy import sqrt

def convsky(map, pix_sigma,tab,instrument="MUSTANG", convolve=True,
            joint=False,jbeam=True,beammap=[]):
    
    convolved = scipy.ndimage.filters.gaussian_filter(map, pix_sigma)
    if instrument == "NIKA":
        if 'beammap' in vars():
            if len(beammap) > 0:
                convolved =  scipy.signal.fftconvolve(map,beammap,mode='same')
                    
    if jbeam != True:
        return map
    if convolve == True:
        if instrument == "MUSTANG" or instrument == "MUSTANG2" or instrument == "NIKA2":
            mapfilt = imf.fourier_filtering_2d(convolved,'tab',(tab[0,0:],tab[1,0:]))
        if instrument == "BOLOCAM":
            centre = convolved.shape[0]/2
            convolved = convolved[centre-21:centre+21,centre-21:centre+21]
            mapfilt_ft = tab*np.fft.fft2(convolved)
            mapfilt = np.real(np.fft.ifft2(mapfilt_ft))
        if instrument == "NIKA":
            mapfilt = imf.fourier_filtering_2d(convolved,'tab',(tab[0][0],tab[0][1]))
        return mapfilt
    else:
        return convolved

def conv_inst_beam(map,pixs,instrument="MUSTANG"):

    fwhm1,norm1,fwhm2,norm2,fwhm,smfw,freq,FoV = rdi.inst_params(instrument)
    pix_sigma1 = (fwhm1 /(pixs * 2 * sqrt((2 * np.log(2))))).decompose()
    pix_sigma2 = (fwhm2 /(pixs * 2 * sqrt((2 * np.log(2))))).decompose()
    pix_sigma1 = pix_sigma1.value
    pix_sigma2 = pix_sigma2.value

    conv1 = scipy.ndimage.filters.gaussian_filter(map, pix_sigma1)*norm1
    conv2 = scipy.ndimage.filters.gaussian_filter(map, pix_sigma2)*norm2
    convolved = conv1+conv2

    return convolved

def conv_gauss_beam(map,pixs,fwhm):

    pix_sigma = (fwhm /(pixs * 2 * sqrt((2 * np.log(2))))).decompose()
    pix_sigma = pix_sigma.value
    conv = scipy.ndimage.filters.gaussian_filter(map, pix_sigma)

    return conv

def apply_xfer(mymap, tab,instrument="MUSTANG"):
    """
    Applies a transfer function based on the instrument for which you want to
    simulate an observation for. This is namely dependent on the format of the
    transfer function. It would be great if the format were standardized...but
    for now, we'll stick with this. 

    Caution
    -------
    The transfer function for MUSTANG2 and NIKA2 are *not* accurate as of 
    June 2017. (They use the old transfer function...just to have something
    defined.)
    
    Parameters
    ----------
    mymap       -  float 2D numpy array
    tab         -  A tabulated (or 2D array) of the transfer function
    instrument  -  A string that indicates which instrument (transfer function)
                   is being used.
    Returns
    -------
    mapfilt     -  The filtered mymap (2D array)
    """

    if instrument == "MUSTANG" or instrument == "MUSTANG2" or instrument == "NIKA2":
        mapfilt = imf.fourier_filtering_2d(mymap,'tab',(tab[0,0:],tab[1,0:]))
    if instrument == "BOLOCAM":
        centre = mymap.shape[0]/2
        mymap = mymap[centre-21:centre+21,centre-21:centre+21]
        mapfilt_ft = tab*np.fft.fft2(mymap)
        mapfilt = np.real(np.fft.ifft2(mapfilt_ft))
    if instrument == "NIKA":
        mapfilt = imf.fourier_filtering_2d(mymap,'tab',(tab[0][0],tab[0][1]))
    
    return mapfilt
