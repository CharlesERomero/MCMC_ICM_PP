import numpy as np
import astropy.coordinates as apc  # http://www.astropy.org/
from astropy import wcs
from astropy.io import fits
from astropy import units as u
from astropy.time import Time as APT
from astropy.visualization.wcsaxes import WCSAxes       
import matplotlib.pyplot as plt
import matplotlib.dates
from matplotlib import colors
from matplotlib import colorbar as cb
import os,re,copy

def plot_image_wcs(mymap, mapmask, wcs,dpi=200,myfontsize=15,zoom=True,filename='My_Map',
                   plotmask=True,savedir='',cblim=False,mymin=-10,mymax=10,mytitle='SNR map of ',
                   target='RXJ1347',format='png',islog=False,cbtit='SNR',cbar=True,
                   subtitle=None,cmap='bwr',flipsign=False):

    n_title_lines = 1                  # TBD what you want
    yaxoff = n_title_lines*0.025  #
    fig = plt.figure(dpi=dpi,figsize=(8,8)); axpos=[0.2, 0.20-yaxoff, 0.7, 0.7]
    ax = WCSAxes(fig, axpos, wcs=wcs); fig.add_axes(ax)
    fig.add_axes(ax)  # note that the axes have to be added to the figure
    if cblim == False:
        cax = ax.imshow(mymap,interpolation='none',cmap=cmap,origin='lower')
        #plt.contour(mymap, [0],colors=('black'),linewidths=3)        
    else:
        if islog == True:
            if flipsign == True:
                mymap*=-1.0
                mycmap = cmap+'_r'
                cax = ax.imshow(mymap,interpolation='none',origin='lower',
                                norm=colors.LogNorm(vmin=mymin,vmax=mymax),cmap=mycmap)
                mymap*=-1.0
            else:
                cax = ax.imshow(mymap,interpolation='none',origin='lower',
                            norm=colors.LogNorm(vmin=mymin,vmax=mymax),cmap=cmap)
        else:
            cax = ax.imshow(mymap,interpolation='none',origin='lower',
                            vmin=mymin,vmax=mymax,cmap=cmap)

    if cbar == True:
        #print mymin,mymax,np.min(mymap),np.max(mymap),flipsign
        mycb = fig.colorbar(cax)
        mycb.set_label(cbtit,fontsize=myfontsize)
        #mycb.ax.set_label("Jy/beam")
        if flipsign == True:
            myytick_labels = mycb.ax.get_yticklabels()  # horizontal colorbar
            ### Define a regular expression; here, the first digit:
            pattern = re.compile(r'\d')
            #print 'Checking for the labels: ',type(myytick_labels)
            #for ytl in myytick_labels:
                #print 'Checking for individual labels: ',type(ytl)
                #mytext = ytl.get_text()
                #print mytext, '   And it is TypeNone: ',type(mytext) == type(None)
                #fuckthisshit = pattern.search(mytext)
                #imout = fuckthisshit.span()
                #pythonsucksballs = insertChar(mytext, imout[0], u'-')
                #imout = fuckthisshit[0].start()
                ### Find the first occurence of this regular expression:
                #significand = pattern.search(mytext).group()
                #print significand
                #myreplacement = ''.join((u'-',significand))
                #myfc = copy.deepcopy(myreplacement)
                #import pdb;pdb.set_trace()
                ### Last number = number of occurences:
                #new_label = mytext.replace(significand,myreplacement,1)
                #new_label = pythonsucksballs
                #ytl.set_text(new_label)

            #cbar.ax.set_xticklabels(['Low', 'Medium', 'High'])  # horizontal colorbar
            
    plt.title(mytitle+target,fontsize=myfontsize*1.2,y=1.05+2*yaxoff)
    #plt.suptitle('Slice angle on sky measured as indicated by arrow.',x=0.45,y=0.84,fontsize=14,color='blue')
    if type(subtitle) != type(None):
        plt.suptitle(subtitle,x=0.45,y=0.84,fontsize=14,color='blue')

    if zoom == True:
        xsz,ysz = mymap.shape
        ax.set_xlim(xsz/4,3*xsz/4)
        ax.set_ylim(ysz/4,3*ysz/4)
        filename=filename+'_zoom'

    ### I hope this works....
    #my_annotations(ax)
    ax.set_xlabel('Right Ascension (J2000)')
    ax.set_ylabel('Declination (J2000)')
        
    if plotmask == True:
        plot_mask(ax,mapmask)

    cwd = os.getcwd();
    if savedir == '': savedir=cwd
    fullbase = os.path.join(savedir,filename)
    fulleps = fullbase+'.eps'; fullpng = fullbase+'.png'
    if format == 'png':
        plt.savefig(fullpng,format='png')
    else:
        plt.savefig(fulleps,format='eps')

    #plt.close(fig)
    plt.close()

        
def get_slice_mask(angmap,minangle,maxangle):

    if np.min(angmap) < 0:
        ltz = angmap < 0
        angmap[ltz]+= 2.0*np.pi

    
    nx, ny = angmap.shape
    maskang = np.zeros((nx,ny))
    if (minangle < maxangle):
        cond1 = angmap < maxangle;        cond2 = angmap >= minangle
        odc1  = cond1.reshape(nx*ny);     odc2  = cond2.reshape(nx*ny)
        mask = [all(tup) for tup in zip(odc1,odc2)]
        avgangle = (minangle+maxangle)/2.0
    else:
        cond1 = angmap >= minangle;       cond2 = angmap < maxangle
        odc1  = cond1.reshape(nx*ny);     odc2  = cond2.reshape(nx*ny)
        mask = [any(tup) for tup in zip(odc1,odc2)]
        avgangle = (minangle+maxangle+2.0*np.pi)/2.0 % (2.0*np.pi)

    odma=maskang.reshape(nx*ny)
    odma[mask]=1.0
    maskang = odma.reshape(nx,ny)
    
    return maskang
    
def plot_mask(ax,mask,color='k'):

    cset = ax.contour(mask, [0.5], colors=color,linewidths=3)
    plt.contour(mask, [0.5], colors=('black'),linewidths=3)

    
def plot_maskrad(ax,image,hdr,mrad,pixs,color='k'):

    dxa,dya = get_xymap(hdr)
    contmap = np.zeros(image.shape)
    drr = (dxa**2 + dya**2)**0.5
    masked = (drr < mrad/pixs)
    contmap[masked]=1
    cset = plt.contour(contmap, [0.5], colors=color)
#    cset = ax.contour(contmap, [0.5], colors=color)

def get_xymap(hdr):

    xsz = hdr['naxis1']
    ysz = hdr['naxis2']
    xar = np.outer(np.arange(xsz),np.zeros(ysz)+1.0)
    yar = np.outer(np.zeros(xsz)+1.0,np.arange(ysz))
    ####################
    
    xcen = hdr['CRPIX1']
    ycen = hdr['CRPIX2']
    dxa = xar - xcen
    dya = yar - ycen
    
    return dxa,dya

def get_snrmap():

    mydir='/home/data/MUSTANG2/AGBT17_Products/RXJ1347/post2018/'
    myfile='grid_pca7_f_Low0.080__snr.fits'
    snrmap, ras, decs, hdr, pixs = get_astro(mydir+myfile,ext=0)

    return snrmap,hdr

def get_astro(file,ext=1):

    hdu = fits.open(file)
    hdr = hdu[ext].header
    image_data = hdu[ext].data
    dxa,dya = get_xymap(hdr)
    ### RA and DEC in degrees:
    if 'CD1_1' in hdr.keys():
        ras = dxa*hdr['CD1_1'] + dya*hdr['CD2_1'] + hdr['CRVAL1']
        decs= dxa*hdr['CD1_2'] + dya*hdr['CD2_2'] + hdr['CRVAL2']
        pixs= abs(hdr['CD1_1'] * hdr['CD2_2'])**0.5 * 3600.0
    if 'PC1_1'  in hdr.keys():
        ras = dxa*hdr['PC1_1']*hdr['CDELT1'] + \
              dya*hdr['PC2_1']*hdr['CDELT2'] + hdr['CRVAL1']
        decs= dxa*hdr['PC1_2']*hdr['CDELT1'] + \
              dya*hdr['PC2_2']*hdr['CDELT2'] + hdr['CRVAL2']
        pixs= abs(hdr['PC1_1']*hdr['CDELT1'] * \
                  hdr['PC2_2']*hdr['CDELT2'])**0.5 * 3600.0
       

    return image_data, ras, decs, hdr, pixs

def get_wcs(hdr):

    mywcs = wcs.WCS(hdr)

    return mywcs

def my_annotations(ax):

    ax.annotate(r'$\theta$',
                xy=(225, 325), xycoords='data',
                xytext=(200, 0), textcoords='offset points',
                arrowprops=dict(arrowstyle="->",
                                connectionstyle="arc3,rad=1"))



    #connectionstyle="angle3,angleA=100,angleB=190"
    #connectionstyle="arc3,rad=1"

def insertChar(mystring, position, chartoinsert):
    longi = len(mystring)
    mystring   =  mystring[:position] + chartoinsert + mystring[position:] 
    return mystring  
