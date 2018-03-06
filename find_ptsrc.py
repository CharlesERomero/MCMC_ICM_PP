from mpfit import mpfit # https://github.com/scottransom/presto/blob/master/lib/python/mpfit.py
import scipy.ndimage
import retrieve_data_info as rdi
import instrument_processing as ip
import numpy as np
import astropy.units as u
import time
from astropy.io import fits

today=time.strftime("%d%b%Y")
todaysp=time.strftime("%d %b %Y")

def find_multiple_ptsrcs(map,pixs,issmo=False,issnr=False,wtmap=[],instrument=[],
                         nptsrcs=3):

    mymap=map
    loclist=[]
    smooth = True
    snrmap = make_snr(map,wtmap,pixs,instrument=instrument,smooth=True)
    issnr=True
    
    
    for i in range(nptsrcs):

        indloc,maxsnr,snr,peak = find_ptsrc_fixed_shape(snrmap,pixs,issmo=smooth,issnr=issnr,
                                                   wtmap=wtmap,instrument=instrument)
        print maxsnr
        ptsrc = create_ptsrc(mymap,pixs,instrument=instrument,indloc=indloc)
        snrmap= snrmap-ptsrc*maxsnr
        loclist.append(indloc)

    return loclist

def find_ptsrc_fixed_shape(map,pixs,issmo=False,issnr=False,wtmap=[],instrument=[]):
    """
    Returns the location of a pre-defined beam shape, at the highest
    significant peak.
    
    Parameters
    __________
    map       - 
    
    
    Returns
    -------
    The structure mapping
    """

    fwhm1,norm1,fwhm2,norm2,fwhm,smfw,freq,FoV = rdi.inst_params(instrument)
    if issmo == False:
        smoothed = ip.conv_inst_beam(map,pixs,instrument=instrument)
    else:
        smoothed=map

    if issnr == False:
        mywarning="The input map is not an SNR map, but a weight map is"+\
            "needed to create an SNR map."
        if wtmap == []:
            raise Exception(mywarning)
        else:
            smwt = ip.conv_inst_beam(wtmap,pixs,instrument=instrument)
            snr  = smoothed * smwt**0.5
    else:
        snr = smoothed

    maxsnr=np.max(snr)
    indloc=np.where(snr == maxsnr)
    peak = smoothed[indloc]

    return indloc,maxsnr,snr,peak

def create_ptsrc(map,pixs,instrument=[],indloc=[]):

    ptsrcmap = map*0.0   # Clear the map
    fwhm1,norm1,fwhm2,norm2,fwhm,smfw,freq,FoV = rdi.inst_params(instrument)
    sigma1 = (fwhm1/(pixs * 2.0 * np.sqrt((2.0 * np.log(2.0))))).decompose()  
    sigma2 = (fwhm2/(pixs * 2.0 * np.sqrt((2.0 * np.log(2.0))))).decompose()
    sigma1 = sigma1.value ;
    sigma2 = sigma2.value  # Don't want "quantities"
    x = np.arange(0,map.shape[0],1.0)
    y = np.arange(0,map.shape[1],1.0)
    xcen = indloc[0].item(0)
    ycen = indloc[1].item(0)
    
    g1 = np.exp(((-(x -xcen)**2 - (y -ycen)**2) )/( 2.0 * sigma1**2))
    g2 = np.exp(((-(x -xcen)**2 - (y -ycen)**2) )/( 2.0 * sigma2**2))

    ptsrc = norm1*g1 + norm2*g2

    return ptsrc

def tdg_sigmas(p,fjac=None, x=None, y=None, data=None, err=None):

    status = 0
    model = twoD_Gauss(p, x=x, y=y)
    return [status, (data-model)/err]
    
def twoD_Gauss(p, x=None, y=None):
    """
    Returns a 2D Gaussian based on the input parameters, p
    
    Parameters
    __________
    p         - p[0] = mean level
                p[1] = x-centroid
                p[2] = y-centroid
                p[3] = sigma_a ("major axis" - but not necessarily)
                p[4] = sigma_b ("minor axis" - but not necessarily)
                p[5] = rotation
                p[6] = peak
    
    Returns
    -------
    The structure mapping
    """
    xrot = (x -p[1]) * np.cos(p[5]) - (y -p[2]) * np.sin(p[5])
    yrot = (y -p[2]) * np.cos(p[5]) + (x -p[1]) * np.sin(p[5])
    xgau = np.exp((-(xrot)**2) / ( 2.0 * p[3]**2))
    ygau = np.exp((-(yrot)**2) / ( 2.0 * p[4]**2))
    tdg  = xgau*ygau*p[6] + p[0]
    
    return tdg

def make_snr(map,wtmap,pixs,instrument=[],smooth=True):

    if smooth == True:
        mymap = ip.conv_inst_beam(map,pixs,instrument=instrument)
        mywtm = ip.conv_inst_beam(wtmap,pixs,instrument=instrument)
    else:
        mymap = map
        mywtm = wtmap

    snr  = mymap * mywtm**0.5

    return snr

def find_ptsrc_loc_shape(map,pixs,issmo=False,issnr=False,wtmap=[],instrument=[]):
    """
    Returns the location and fitted 2D Gaussian.
    
    Parameters
    __________
    map       - 
    
    
    Returns
    -------
    The structure mapping
    """

    indloc,maxsnr,snr,peak = find_ptsrc_fixed_shape(map,pixs,issmo=issmo,issnr=issnr,
                                               wtmap=wtmap,instrument=instrument)
    fwhm1,norm1,fwhm2,norm2,fwhm,smfw,freq,FoV = rdi.inst_params(instrument)
    mean_guess=0.0;xguess=indloc[0].item(0);yguess=indloc[1].item(0)
    sigma_maj = (fwhm/(pixs * 2.0 * np.sqrt((2.0 * np.log(2.0))))).decompose()  
    sigma_min = sigma_maj; rot_guess=1.0; peak_height=peak
    p0=np.array([mean_guess,xguess,yguess,sigma_maj,sigma_min,rot_guess,peak_height])
    x = np.outer(np.arange(0,map.shape[0],1.0),np.zeros(map.shape[1])+1.0)
    y = np.outer(np.zeros(map.shape[0])+1.0,np.arange(0,map.shape[1],1.0))
    nzi=np.where(wtmap > 0)
    err=wtmap*0.0 + 1.0
    err[nzi]=wtmap[nzi]**(-0.5)
    xf=x.reshape(map.shape[0]*map.shape[1])
    yf=y.reshape(map.shape[0]*map.shape[1])
    ef=err.reshape(map.shape[0]*map.shape[1])
    df=map.reshape(map.shape[0]*map.shape[1])    
    ptsrc = twoD_Gauss(p0,x=xf, y=yf)
    ### OK, it looks like flattening the arrays also helps / is necessary.
    fa = {'x':xf, 'y':yf,'data':df,'err':ef} # Flattened
 #   fa = {'x':x, 'y':y,'data':map,'err':err} # Not flattened??
    ### You have to specify parinfo, otherwise @!$@# go crazy. WTF. Just WTF.
    parinfo = [{'value':0., 'fixed':0, 'limited':[0,0], 'limits':[0.,0.]}]*len(p0)
#    parinfo[5]['limited'][0]=1  ; parinfo[5]['limited'][1]=1
#    parinfo[5]['limits'][0]=0.0 ; parinfo[5]['limits'][1] =2.0*np.pi    
### Also, don't call with a string; call with the FUNCTION (as known by Python)
    myp = mpfit(tdg_sigmas, p0, functkw=fa,parinfo=parinfo)
#    myp = mpfit('twoD_Gauss', p0, functkw=fa)

    return myp.params

def make_ptsrc_model(instrument="MUSTANG",name="rxj1347",writefile=True,fixshape=False,
                     normalize=False):

    resdir='/home/romero/MUSTANG2/Python_Products/'
    fpsfn = resdir+'fitted_ptsrc_for_'+name+'.fits'
    data_map, wt_map, header, ras, decs, pixs, w, tab = \
        rdi.get_maps_and_info(instrument,name,real=True)
    pixs*=u.arcsec       # Good to have unit cancellation (ahead)
    indloc,maxsnr,snr,peak = find_ptsrc_fixed_shape(data_map,pixs,issmo=False,issnr=False,
                                                       wtmap=wt_map,instrument=instrument)
    fwhm1,norm1,fwhm2,norm2,fwhm,smfw,freq,FoV = rdi.inst_params(instrument)
    mean_guess=0.0;xguess=indloc[0].item(0);yguess=indloc[1].item(0)
    sigma_maj = (fwhm/(pixs * 2.0 * np.sqrt((2.0 * np.log(2.0))))).decompose()  
    sigma_min = sigma_maj; rot_guess=0.0; peak_height=peak
    p=[mean_guess,xguess,yguess,sigma_maj,sigma_min,rot_guess,peak_height]
    xsz=data_map.shape[0] ; ysz=data_map.shape[1]
    xar = np.outer(np.arange(xsz),np.zeros(ysz)+1.0)
    yar = np.outer(np.zeros(xsz)+1.0,np.arange(ysz))
    nzi=np.where(wt_map > 0)        # We don't want zero-weight pixels
    err=wt_map*0.0 + 1.0            # Let the default error values be 1 (for now)
    err[nzi]=wt_map[nzi]**(-0.5)    # Calculate appropriate error values where valid.

    fittedp = find_ptsrc_loc_shape(data_map,pixs,issmo=False,issnr=False,wtmap=wt_map,
                                      instrument=instrument)
    foo = np.floor(fittedp[5]/(2.0*np.pi));bar = fittedp[5] - (foo*2.0*np.pi);
    fittedp[5] = bar  ### This just does a "mod" over 2*pi
    if normalize==True:
        fittedp[6]=1.0
    mnlvl = fittedp[0]; fittedp[0]=0.0
    ptsrc = twoD_Gauss(fittedp,x=xar, y=yar); ptsrchdr = header
    fittedp[0]=mnlvl
    hdu1 = fits.PrimaryHDU(ptsrc,header=ptsrchdr)    
    pfields=['mean_lvl','x0','y0','major','minor','angle','height']
    pcomments=['This is the fitted value, but it is NOT applied!',
               'X coordinate of centroid','Y coordinate of centroid',
               'Gaussian Sigma (pixels)','Gaussian Sigma (pixels)',
               'Rotation angle','Fitted amplitude of the Gaussian']
    for i in range(len(fittedp)):
        hdu1.header.append((pfields[i],fittedp[i],pcomments[i]))
        
#    ptsrchdr.add_history("Fits were performed on "+todaysp+ "using MPFIT.")
#    if normalize==True:
#        ptsrchdr.add_history("Peak was forced (not fitted) to 1.")

    if fixshape==True:
        sigma1 = (fwhm1/(pixs * 2.0 * np.sqrt((2.0 * np.log(2.0))))).decompose()  
        sigma2 = (fwhm2/(pixs * 2.0 * np.sqrt((2.0 * np.log(2.0))))).decompose()  
        p1=[0.0,fittedp[1],fittedp[2],sigma1,sigma1,0.0,norm1*fittedp[6]]
        p2=[0.0,fittedp[1],fittedp[2],sigma2,sigma2,0.0,norm2*fittedp[6]]
        ptsrc = twoD_Gauss(p1,x=xar,y=yar) + twoD_Gauss(p2,x=xar,y=yar)
        ptsrchdr.append(('Sigma1',sigma1.value,'Used this fixed Gaussian Sigma1'))
        ptsrchdr.append(('Sigma2',sigma2.value,'Used this fixed Gaussian Sigma2'))
        ptsrchdr.append(('Norm1',norm1,'Used this fixed Gaussian Norm1'))
        ptsrchdr.append(('Norm1',norm2,'Used this fixed Gaussian Norm2'))
        ptsrchdr.add_history("The location was found with MPFIT, but we have then opted to")
        ptsrchdr.add_history("use a pre-determined double Gaussian, with the above parameters.")
                        
    if writefile == True:
        hdulist = fits.HDUList([hdu1])
        hdulist.writeto(fpsfn,overwrite=True)
#        fits.writeto(fpsfn,ptsrc,ptsrchdr,overwrite=True)
        return ptsrc,fpsfn
    else:
        return ptsrc
