import importlib,time,os
import numpy as np
from astropy.io import fits
import scipy.ndimage
import matplotlib.pyplot as plt
import mpfit
from scipy.optimize import curve_fit

def get_maps(name='rxj1347_wshock',instrument='MUSTANG2',smfw=9.0,pixs=1.0):
    
    input_struct = importlib.import_module(name+'_info')
    inputdata    = input_struct.files(instrument=instrument,map_file='all')
    inputnoise   = input_struct.files(instrument=instrument,map_file='noise')
    
    data_map, header = fits.getdata(inputdata.fitsfile, header=True)
    wt_map = fits.getdata(inputdata.wtfile,inputdata.wtext)

    noise_map, header = fits.getdata(inputnoise.fitsfile, header=True)

    #noise_smooth = np.zeros(noise_map.shape)
    sigma = (smfw/pixs)/(2.0*np.sqrt(2.0*np.log(2.0)))
    bv = 2.0*np.pi*(sigma**2)
    data_smooth=scipy.ndimage.filters.gaussian_filter(data_map,sigma)
    noise_smooth=scipy.ndimage.filters.gaussian_filter(noise_map,sigma)
    wt_smooth=scipy.ndimage.filters.gaussian_filter(wt_map,sigma)
    wt_smooth*= bv

    nzi   = (wt_map > 0)
    wtcut = 0.1 * np.median(wt_map[nzi])
    gi    = (wt_map > wtcut)
    nsnr  = noise_map[gi] * 1.0e3 * np.sqrt(wt_map[gi])
    snsnr = noise_smooth[gi] * 1.0e3 * np.sqrt(wt_smooth[gi])

    nhist,nedge = np.histogram(nsnr,'auto')
    shist,sedge = np.histogram(snsnr,'auto')

    nmids = (nedge[:-1]+ nedge[1:])/2.0
    smids = (sedge[:-1]+ sedge[1:])/2.0

    noisenorm,noisesig = crgauss(nmids**2,nhist)
    smoothnorm,smoothsig = crgauss(smids**2,shist)

    nopt, ncov = curve_fit(gauss_function, nmids, nhist, p0 = [max(nhist), 0.0, 1.0])
    sopt, scov = curve_fit(gauss_function, smids, shist, p0 = [max(shist), 0.0, 1.0])
    #import pdb;pdb.set_trace()
    
    #sfit = smoothnorm* np.exp(-smids**2 / (2.0 * smoothsig**2))
    #nfit = noisenorm * np.exp(-nmids**2 / (2.0 * noisesig**2))
    nfit = gauss_function(nmids,nopt[0],nopt[1],nopt[2])
    sfit = gauss_function(smids,sopt[0],sopt[1],sopt[2])
    
    #print noisesig, smoothsig
    print nopt[2],sopt[2]
    
    plt.figure(2,figsize=(20,12));    plt.clf()
    #plt.axvline(rin,color=axcol, linestyle ="dashed")
    #plt.axvline(rout,color=axcol, linestyle ="dashed")
    plt.plot(nmids,nhist,label = "Noise SNR")
    plt.plot(smids,shist,label = "Smoothed noise SNR")

    plt.plot(nmids,nfit,label  = "Fit to Noise SNR "+str(nopt[2]))
    plt.plot(smids,sfit,label  = "Fit to Smoothed SNR "+str(sopt[2]))
   
    plt.yscale("log")
    #plt.xscale("log")
    plt.xlabel("SNR")
    plt.ylabel("Count")
    plt.title("Testing")
    plt.grid()
    plt.ylim((1.0,np.max(nhist)))

    plt.legend()
    
    filename = "example_plot.png"
    #savedir  = "/home/romero/Python/StandAlone/Comptony_Modelling/"
    #fullpath = savedir+filename
    plt.savefig(filename)
    #plt.close()

   
def crgauss(x,y,yerr=1.0):

    lny = np.log(y)
    gi  = (lny > 0)
    
    if hasattr(yerr,'__len__'):
        wts = 1.0/yerr**2
    else:
        wts = y*0.0 + yerr

    sw = np.sum(wts[gi]); swx = np.sum(wts[gi]*x[gi])
    swy = np.sum(wts[gi]*lny[gi])
    swxy = np.sum(wts[gi]*x[gi]*lny[gi])
    swxx = np.sum(wts[gi]*x[gi]*x[gi])
    delta = sw*swxx - swx**2

    A = (swxx*swy - swx*swxy)/delta
    B = (sw*swxy - swx*swy)/delta

    plt.figure(1,figsize=(20,12));    plt.clf()
    plt.plot(x,lny,label  = "Data")
    fit = A + B*x
    plt.plot(x,fit,label  = "Fit")
    filename = "my_linear_fit.png"
    plt.savefig(filename)

    #import pdb;pdb.set_trace()
    sig  = np.sqrt(-1.0/2*B)
    norm = np.exp(A)

    return norm, sig

def gauss_function(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))
