import numpy as np
import astropy.units as u
import cosmolopy.distance as cd
from scipy.interpolate import interp1d
import scipy.signal, scipy.ndimage
import os
import image_filtering as imf
from astropy.io import fits
import ellipsoidal_shells as es
import matplotlib.pyplot as plt
import plot_mcmc_results as pmr
import rw_resultant_fits as rwrf


newpath='/home/romero/Results_Python/MUSTANG2/rxj1347_wshock'


def get_az_profile(mymap, xymap, minangle, maxangle, geoparams=[0,0,0,1,1,1,0,0],
                   wtmap=[],thetamask=[]):
    """
    This program is designed to take a radial profile within a given azimuthal slice.

    ----
    INPUTS:
    mymap         - A 2D (numpy) array of which the profile should be taken.
    xymap         - a 2-tuple, each element is a 2D array of either X or Y coordinates.
    minangle      - The minimum angle of the slice, in radians.
    maxangle      - The maximum angle of the slice, in radians.
    ***NOTE:        You can have minangle > maxangle if, for example, you want to go
                    between 7*pi/4 and pi/4 (a slice subtending only pi/2). Of course,
                    if you want to go from pi/4 to 7*pi/4, that slice would subtend 3*pi/2
                    radians.
    geoparams     - This follows the geoparams used in the rest of my code. The entries
                    are as follows:
                    geoparams[0] = x0 (offset from 0 in the x-coordinate map)
                    geoparams[1] = y0 (offset from 0 in the y-coordinate map)
                    geoparams[2] = theta_rot (rotation after  translation)
                    geoparams[3] = elliptical a - scaling for the major axis (in the plane)
                    geoparams[4] = ellpitical b - scaling for the minor axis (in the plane)
                    geoparams[5] = elliptical c - scaling for the axis along the line of sight
                    geoparams[6] = Shock "xi" - a power law for tapering.
                    geoparams[7] = Opening angle for a shock
                    ** Note that [6] and [7] are not so relevant here.
    wtmap         - If specified, it goes with the input map and specifies the weights
                    (inverse variance) associated with each pixel.

    ----
    OUTPUTS:


    """
    
    (x,y) = xymap
    angmap = np.arctan2(y,x)        # Calculate the angles for all pixels.
    negind = (angmap < 0)           # Arctan2 goes -pi < theta < pi
    angmap[negind] += 2.0*np.pi     # I prefer to go 0 < theta < 2*pi
    x,y = es.rot_trans_grid(x,y,geoparams[0],geoparams[1],geoparams[2])
    x,y = es.get_ell_rads(x,y,geoparams[3],geoparams[4])
    radmap = np.sqrt(x**2 + y**2)
    nx, ny = radmap.shape

    #theta = radmap*(u.arcsec).to("radian");  theta_min = np.min(theta_range)
    ##########################################################################
    ### I need to do some troubleshooting (12 March 2018)
    nx, ny = radmap.shape
    xx     = np.outer(np.ones(nx),np.arange(0,(ny), 1))
    yy     = np.outer(np.arange(0,(nx), 1),np.ones(ny))
    trad   = radmap < 3.0      # Test Rad = less than 3 arcseconds
    tval   = mymap[trad]       # Test Val = those values within 3 arcseconds
    xcoord = xx[trad]
    ycoord = yy[trad]

    foo = mymap > 0.00111263; foo1 = mymap[foo]; foo2 = xx[foo]; foo3 = yy[foo]
    bar = foo1 < 0.00111265;  bar1 = foo1[bar]; bar2 = foo2[bar]; bar3 = foo3[bar]

    #print 'found a pixel of value ',bar1,' located at x = ',bar2,' and y = ',bar3
    #print tval, xcoord, ycoord
    
    #import pdb;pdb.set_trace()
    #bi=np.where(theta < theta_min);   theta[bi]=theta_min
    ##########################################################################

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

    #mask    = np.array(mask);   mask     = mask.reshape(nx,ny)
    maskrad = radmap.reshape(nx*ny)[mask]
    maskprof= mymap.reshape(nx*ny)[mask]
    ### Ah, I was sorting each array separately! Stupid me!
    ### I need to sort by radius for both arrays. Here is my solution:
    dtype   = [('radius',float),('values',float)]          # Not the most efficient
    profstr = np.array(zip(maskrad,maskprof),dtype=dtype)  # But, it's explicit
    sortbyr = np.sort(profstr,order='radius')              # and it gets the job done
    sortrad, sortprof= zip(*sortbyr)                       # Voila.
    sortrad=np.array(sortrad); sortprof=np.array(sortprof) # OK, now voila.
    
    if len(thetamask) > 0:
        if thetamask.shape == angmap.shape:
            onedarr = thetamask.reshape(nx*ny)
            onedarr[mask] = avgangle
            thetamask = onedarr.reshape(nx,ny)
        else:
            raise Exception

        return sortrad, sortprof, thetamask

    else:
        return sortrad, sortprof
    
def get_slices(dv,hk,angstep=np.pi/6.0,format='png'):

    angarr = np.arange(0,2.0*np.pi,angstep)
    profsl = {}
    
    for myinst in hk.instruments:
        data    = dv[myinst].maps.data
        weights = dv[myinst].maps.masked_wts  # Best to mask weights (not use the entire map!)
        raw_wts = dv[myinst].maps.wts  
        xymap   = dv[myinst].mapping.xymap
        hdr     = dv[myinst].maps.header
        slices  = []
        mythetamask = np.zeros(data.shape)
        for angmin in angarr:
            angmax = angmin + angstep
            print angmin, angmax
            rads, prof, mythetamask = get_az_profile(data, xymap, angmin, angmax,
                                                     thetamask=mythetamask)
            mybins=np.arange(0.0,60.0,4.0)
            binres = radial_bin(rads, prof,10,rmax=60.0,bins=mybins,minangle=angmin,maxangle=angmax)
            #plot_one_slice(binres)
            slices.append(binres)

        profsl[myinst]=slices
        prefilename='Mask_of_various_slices';        filename='v0.'
        mapaxisunits='arcseconds';                   mapunits=dv[myinst].maps.units
        title='Mask defining the various slices I used for RXJ 1347'
        fullpath = os.path.join(newpath,prefilename+filename+format)
        pmr.plot_sky_map(fullpath,mythetamask,title,mapaxisunits,mapunits,format=format)
        filename='_no_coordinates.fits'
        fullpath = os.path.join(newpath,prefilename+filename)
        rwrf.savemap(mythetamask,fullpath,header=hdr)

    return profsl

def get_two_slices(dv,hk,angmin,angmax):

    profsl = {}
    
    for myinst in hk.instruments:
        data    = dv[myinst].maps.data
        weights = dv[myinst].maps.masked_wts  # Best to mask weights (not use the entire map!)
        raw_wts = dv[myinst].maps.wts  
        xymap   = dv[myinst].mapping.xymap
        hdr     = dv[myinst].maps.header
        slices  = []
        mythetamask = np.zeros(data.shape)

        ###############################################################################
        ### Get the first slice (minangle = angmin)
        
        rads, prof, mythetamask = get_az_profile(data, xymap, angmin, angmax,
                                                 thetamask=mythetamask)
        mybins=np.arange(0.0,60.0,4.0)
        binres = radial_bin(rads, prof,10,rmax=60.0,bins=mybins,minangle=angmin,maxangle=angmax)
        slices.append(binres)
        
        ###############################################################################
        ### Get the second slice (minangle = angmax)
        
        rads, prof, mythetamask = get_az_profile(data, xymap, angmax, angmin,
                                                 thetamask=mythetamask)
        mybins=np.arange(0.0,60.0,4.0)
        binres = radial_bin(rads, prof,10,rmax=60.0,bins=mybins,minangle=angmax,maxangle=angmin)
        slices.append(binres)
        
        ###############################################################################

        profsl[myinst]=slices
 #       prefilename='Mask_of_various_slices';        filename='v0.png'
 #       mapaxisunits='arcseconds';                   mapunits='radians'
 #       title='Mask defining the various slices I used for RXJ 1347'
 #       fullpath = os.path.join(newpath,prefilename+filename)
 #       pmr.plot_sky_map(fullpath,mythetamask,title,mapaxisunits,mapunits)
 #       filename='_two_slices.fits'
 #       fullpath = os.path.join(newpath,prefilename+filename)
 #       rwrf.savemap(mythetamask,fullpath,header=hdr)

    return profsl

def iter_two_slices(dv,hk,myformat='png'):

    angmin = np.arange(3.0*np.pi/4.0, 5.0*np.pi/4.0, np.pi/8.0)  # steps of 22.5 degrees
    angmax = np.arange(5.0*np.pi/4.0, 7.0*np.pi/4.0, np.pi/8.0)  # steps of 22.5 degrees
    mylist=[]; aminarr = []; amaxarr=[]
    radunits = 'arcseconds' # This is in fact, hard-coded into my xymaps
    myinst='MUSTANG2'
    profunits= dv[myinst].maps.units
    
    fig, ax, fullpath = myplotsetup(radunits=radunits,profunits=profunits,format=myformat)
    linestyles = ['-','--','-.',':']
    markers    = ['o','*','+','D']
    
    for i,thisamin in enumerate(angmin):
        linestyle = linestyles[i]
        for j,thisamax in enumerate(angmax):
            marker = markers[j]
            profsl = get_two_slices(dv,hk,thisamin,thisamax)
            mylist.append(profsl)
            aminarr.append(thisamin)
            amaxarr.append(thisamax)
            #kwargs=[linestyle,marker]
            kwargs={'linestyle':linestyle,'marker':marker,'offset':j}
            quick_plot_two_slices(fig,ax,profsl,**kwargs)

    plt.legend()
    plt.savefig(fullpath,format=myformat)
    return mylist, aminarr, amaxarr

            
def quick_plot_two_slices(fig,ax,profsl,**kwargs):

    myls = kwargs['linestyle']; myma=kwargs['marker']; offset=(kwargs['offset']-2)*0.5
    for myinst in profsl:
        myslices = profsl[myinst]
        name0    = "{:5.2f}".format(myslices[0].minangle)
        name1    = "{:5.2f}".format(myslices[0].maxangle)
        mylabel  = name0+'_to_'+name1
        ax.plot(myslices[0].rads+offset,myslices[0].profavg,label=mylabel,color='g',
                 linestyle=myls,marker=myma,markersize=10.0)
        name0    = "{:5.2f}".format(myslices[1].minangle)
        name1    = "{:5.2f}".format(myslices[1].maxangle)
        mylabel  = name0+'_to_'+name1
        #print myslices[1].npix[10]
        ax.plot(myslices[1].rads+offset,myslices[1].profavg,label=mylabel,color='r',
                 linestyle=myls,marker=myma,markersize=10.0)


def myplotsetup(thispath=newpath,prefilename='Testing_2slice_profiles_',filename='v0.',
                radunits='arcseconds',profunits='Kelvin',format='png'):

    fig = plt.figure(2,figsize=(20,12));    plt.clf()
    ax = fig.add_subplot(111)
    ax.set_xlabel("Radius ("+radunits+")")
    ax.set_ylabel("Map Intensity ("+profunits+")")
    ax.set_title('RXJ1347')
    plt.grid()
    fullpath = os.path.join(thispath,prefilename+filename+format)

    return fig, ax, fullpath
            
class radial_bin:
    
    def __init__(self,radii, profile, nbins=0, rmax=0,bins=[],radunits='arcseconds',
                 profunits='Kelvin',minangle=0,maxangle=2.0*np.pi):
        """
        Note that nbins = len(bins)-1, because bins defines the upper and lower bound for the
        bins.

        """
        if len(bins) == 0:
            if len(radii) == 0:
                print 'shit just got weird.'
                import pdb;pdb.set_trace()
            
            rmin = np.min(radii)
            if rmax == 0:  rmax  = np.max(radii)
            if nbins == 0: nbins = np.ceil(len(radii)**(0.333))
            rstep = (rmax-rmin)/(nbins-1)
            bins = np.arange(rmin,rmax+rstep,rstep)

        radavg=[]; radmin=[]; radmax=[]
        profavg=[]; profrms=[]; npix=[]
        for i,binmin in enumerate(bins[:-1]):
            cond1   = radii > bins[i];     cond2 = radii < bins[i+1]
            radmin.append(bins[i]);        radmax.append(bins[i+1])
            mask    = [all(tup) for tup in zip(cond1,cond2)]
            radavg.append(np.mean(radii[mask]))
            profavg.append(np.mean(profile[mask]))
            profrms.append(np.std(profile[mask]))
            npix.append(len(profile[mask]))

        if minangle < maxangle:
            avgangle = (minangle+maxangle)/2.0
        else:
            avgangle = (minangle+maxangle+2.0*np.pi)/2.0 % (2.0*np.pi)

        self.rads      = np.array(radavg)
        self.radmin    = np.array(radmin)
        self.radmax    = np.array(radmax)
        self.radunits  = radunits
        self.profavg   = np.array(profavg)
        self.profrms   = np.array(profrms)
        self.npix      = np.array(npix)
        self.profunits = profunits
        self.minangle  = minangle
        self.maxangle  = maxangle
        self.name      = "{:5.2f}".format(avgangle)


        
def plot_one_slice(myslice,myformat='png'):

    plt.figure(2,figsize=(20,12));    plt.clf()
    xerr = [myslice.rads - myslice.radmin, myslice.radmax-myslice.rads]
    yerr = [myslice.profrms/np.sqrt(myslice.npix),myslice.profrms/np.sqrt(myslice.npix)]
    plt.errorbar(myslice.rads,myslice.profavg,xerr=xerr,yerr=yerr,fmt='.',
                 label="Radial profile, azimuthal slice",capsize=5)
#    plt.yscale("log")
#    plt.xscale("log")
    plt.xlabel("Radius ("+myslice.radunits+")")
    plt.ylabel("Map Intensity ("+myslice.profunits+")")
    plt.title('RXJ1347')
    plt.grid()
    #newpath='/home/romero/Results_Python/MUSTANG2/rxj1347_wshock'
    prefilename='Testing_slice_profiles_'
    filename='v0.'
    fullpath = os.path.join(newpath,prefilename+filename+myformat)
    plt.savefig(fullpath,format=myformat)

def plot_slices_one_inst(myslices,myformat='png'):

    plt.figure(2,figsize=(20,12));    plt.clf()
    for myslice in myslices:
        xerr = [myslice.rads - myslice.radmin, myslice.radmax-myslice.rads]
        yerr = [myslice.profrms/np.sqrt(myslice.npix),myslice.profrms/np.sqrt(myslice.npix)]
        plt.errorbar(myslice.rads,myslice.profavg,xerr=0,yerr=yerr,
                     label=myslice.name,capsize=5)
        
    plt.xlabel("Radius ("+myslice.radunits+")")
    plt.ylabel("Map Intensity ("+myslice.profunits+")")
    plt.title('RXJ1347')
    plt.grid()
    plt.legend()
    #newpath='/home/romero/Results_Python/MUSTANG2/rxj1347_wshock'
    prefilename='Testing_slice_profiles_'
    filename='v0.'
    fullpath = os.path.join(newpath,prefilename+filename+myformat)
    plt.savefig(fullpath,format=myformat)
        
def plot_all_slices(profsl,myformat='png'):

    for myinst in profsl:
        myslices = profsl[myinst]
        plot_slices_one_inst(myslices,myformat=myformat)

        
