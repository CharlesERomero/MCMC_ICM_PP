from astropy.io import fits
import max_like_fitting as mlf
import astropy.units as u
import ellipsoidal_shells as es
import instrument_processing as ip
import numpy as np
import os

def savemap(map,filename,wtmap=None,header=None):

    hdu1 = fits.PrimaryHDU(map,header=header)
    hdulist = fits.HDUList([hdu1])
    hdulist.info()
    hdulist.writeto(filename,overwrite=True,output_verify="exception")

def saveimg(dv,hk,efv,modelmap=None,weightmap=None,
            component='Bulk',filename="Map.fits"):

    hdu = make_and_save_model_maps(hk,dv,efv)

    return hdu

######################################################################################################
# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +#
#+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +#
#+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +#
#+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +#
#+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +#
######################################################################################################

def clear_maps(maps,hk,dv):

    for myinst in hk.instruments:
        nx,ny = dv[myinst].mapping.radmap.shape
        maps[myinst]=np.zeros((nx,ny))

def make_and_save_model_maps(hk,dv,efv):

    pos = efv.values # Given as is...
    posind = 0;  parbycomp={}  #; indbycomp={}
    zeromaps={}; component={}; maps=[]; yint=[]; outalphas=[]; hdu={}; mapind=0
   
    
    tstr = 'Test_Run_'
    if hk.cfp.testmode == False:
        tstr = 'Full_Run_'

    for myinst in hk.instruments:
        ### Set up the necessary variables for the model map(s) of each instrument
        nx,ny = dv[myinst].mapping.radmap.shape
        zeromaps[myinst]=np.zeros((nx,ny))
        component[myinst]=zeromaps[myinst]+ pos[posind]
        parbycomp[myinst]={'mnlvl':pos[posind]}
        parbycomp[myinst]={'mnlvl_ind':np.array([posind])}
        posind+=1
        hdu[myinst] = []

        ### And go ahead and setup the primary header (data unit)

    maps=[component] 
    conc_hdu(dv,hk,maps[mapind],title='Mean_Level',hdu=hdu); mapind+=1
    clear_maps(zeromaps,hk,dv)

    ### Model the bulk pressure:
    count=1
    for bins,fit_cen,geo,alp,narm in zip(hk.cfp.bulkarc,hk.cfp.bulk_centroid,hk.cfp.bulkgeo,
                                         hk.cfp.bulkalp,hk.cfp.bulknarm):
        compname = 'Bulk'+str(count)
        nbins = len(bins)
        parbycomp[myinst][compname]=pos[posind:posind+nbins]
        parbycomp[myinst][compname+'_bins']=bins # bins are in radians!
        parbycomp[myinst][compname+'_ind']=np.arange(posind,posind+nbins)
        outmaps,posind,ynt,myalphas = mlf.bulk_or_shock_component(pos,bins,hk,dv,efv,fit_cen,geo,alp,narm,zeromaps,posind,
                                              fixalpha=hk.cfp.bulkfix)
        maps.append(outmaps)
        yint.append(ynt); outalphas.extend(myalphas)
        conc_hdu(dv,hk,maps[mapind],title='Bulk',hdu=hdu,count=count); mapind+=1
        print 'hi'
        clear_maps(zeromaps,hk,dv)
        
    ### Model any shocks:
    count=1
    for bins,fit_cen,geo,alp,narm in zip(hk.cfp.shockbin,hk.cfp.shoc_centroid,hk.cfp.shockgeo,
                                         hk.cfp.shockalp,hk.cfp.shocknarm):
        compname = 'Shock'+str(count)
        nbins = len(bins)
        parbycomp[myinst][compname]=pos[posind:posind+nbins]
        parbycomp[myinst][compname+'_bins']=bins
        parbycomp[myinst][compname+'_ind']=np.arange(posind,posind+nbins)
        outmaps,posind,shint,shout = mlf.bulk_or_shock_component(pos,bins,hk,dv,efv,fit_cen,geo,alp,narm,zeromaps,posind,
                                                                 fixalpha=hk.cfp.shockfix,finite=hk.cfp.shockfin)
        maps.append(outmaps)
        conc_hdu(dv,hk,maps[mapind],title='Shock',hdu=hdu,count=count)
        clear_maps(zeromaps,hk,dv)

    ### Model any point sources (hk.cfp.ptsrc is a 2-tuple, the pt src. centroid):
    count=1
    for myptsrc in hk.cfp.ptsrc:
        compname = 'PtSrc'+str(count)
        parbycomp[myinst][compname]=pos[posind]
        parbycomp[myinst][compname+'_ind']=np.array([posind])
        #outmaps,posind = mlf.mk_ptsrc(pos,posind,myptsrc,hk,dv,zeromaps)
        outmaps,posind = mlf.mk_ptsrc_v2(pos,posind,hk,dv,efv.ifp,zeromaps)
        maps.append(outmaps)
        conc_hdu(dv,hk,maps[mapind],title='PtSrc',hdu=hdu,count=count); mapind+=1
        clear_maps(zeromaps,hk,dv)

    ### Model any "blobs" (2D Gaussians):
    ### This is currently not functional because I'm not sure exactly how I want to implement
    ### this feature.
    count=1
    for myblob in hk.cfp.blob:
        compname = 'Blob'+str(count)
        parbycomp[myinst][compname]=pos[posind:posind+3]
        parbycomp[myinst][compname+'_ind']=np.arange(posind,posind+3)
        outmaps,posind = mlf.mk_twodgauss(pos,posind,myptsrc,hk,dv,zeromaps)
        maps.append(outmaps)
        conc_hdu(dv,hk,maps[mapind],title='Blob',hdu=hdu,count=count); mapind+=1
        clear_maps(zeromaps,hk,dv)
        
    write_hdu_to_fits(hk,tstr,hdu)
    efv.paramdict = parbycomp

    return hdu
    
################################################################################################

def conc_hdu(dv,hk,modelsky,title='Unknown',hdu={},count=1):
   
    for mykey in hdu:
        myhdu = hdu[mykey]
        myhdr=dv[mykey].maps.header
    
        if len(myhdu) == 0:
            weights = dv[mykey].maps.wts
            hdu0 = fits.PrimaryHDU(weights,header=myhdr)
            hdu0.header.append(("Title","Weight Map"))
            hdu0.header.append(("Target",hk.hk_ins.name))          
            hdu0.name = 'Weights'

            print modelsky[mykey].shape, len(modelsky),np.min(modelsky[mykey]),np.max(modelsky[mykey])
            
            hdu1 = fits.ImageHDU(modelsky[mykey])
            hdu1.header = dv[mykey].mapping.w.to_header()
            hdu1.name = title+str(count)
            #hdu1.header = dv[key].mapping.w.to_header()
            hdu1.header.append(("Title","Jansky/beam Sky Map"))
            hdu1.header.append(("Target",hk.hk_ins.name))          
            hdu1.header.append(("XTENSION","What Mate"))
            hdu1.header.append(("SIMPLE","T")) 
            hdu1.verify('fix')
            
            myhdu=[hdu0,hdu1]
        else:
            hdu2 = fits.ImageHDU(modelsky[mykey])
            hdu2.header = dv[mykey].mapping.w.to_header()
            hdu2.name = title+str(count)
            hdu2.header.append(("Title",title))
            hdu2.header.append(("Target",hk.hk_ins.name))
            hdu2.header.append(("XTENSION","What Mate"))
            hdu2.header.append(("SIMPLE","T")) 
            hdu2.verify('fix')
            myhdu.append(hdu2)

        count+=1
        hdu[mykey]=myhdu
        print len(myhdu)

def write_hdu_to_fits(hk,tstr,hdu):

    for myinst in hk.instruments:
        myhdu = hdu[myinst]
        hdulist = fits.HDUList(myhdu)
        fbase=tstr+myinst+"_"
        filename=tstr+"Residual.fits"
        fullpath = os.path.join(hk.hk_outs.newpath,hk.hk_outs.prefilename+filename)
        hdulist.writeto(fullpath,overwrite=True)
