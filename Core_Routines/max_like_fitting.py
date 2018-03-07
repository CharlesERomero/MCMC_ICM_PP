import numpy as np                      # A useful package...
import emcee, os, shelve, pickle, time  # A list of modules to import as-is
from astropy.io import fits             # To read/write fits
import ellipsoidal_shells as es         # Integrates ellipsoidal power law distributions
import instrument_processing as ip      # Determines instrument-specific values (e.g. xfer fxn)
import astropy.units as u               # Allows variables (quantities) to have units
import cluster_pressure_profiles as cpp # Creates radius + pressure arrays.
import rw_resultant_fits as rwrf        # Read/Wrte Resultant Fits
from scipy.interpolate import interp1d  # 1-dimensional interpolation (for the integrated profiles)
from os.path import expanduser          # A way to determine what the home directory is
import multiprocessing                  # A way to determine how many cores the computer has
rwrf=reload(rwrf)                       # Reload - primarily of use during development of code.
myhome = expanduser("~")                # What is the user's home directory?
ncpus  = multiprocessing.cpu_count()    # One can decide how aggressively to parallel process.

########################################################################
### March 6, 2018. I think the following bugs are fixed.
########################################################################
### Some *bugs* that had to be fixed:
### (1) We need to be sure that theta_range (dv.mapping.theta_range) encompasses the
###     range of radians subtended by the mapped region. Here, I needed to make sure
###     that the minimum of the range was still smaller than the "theta_min", which I
###     use an overwrite for any radii smaller than that value (so as to avoid errors
###     with a power law going to zero).
### (2) Fixing the transfer function, so as to accommodate very small scales (and thus
###     large wavenumber k)
###
########################################################################
### TO DO:
### (1) Document (comment) this and other codes.
########################################################################

class emcee_fitting_vars:

    def __init__(self,hk,dv,ifp,tag='',nthreads=1):
        """
        The attributes of this class are attributes that are needed for fitting models
        to the data. In particular, we want to assemble attributes (variables) which
        are our initial guesses for the (all) parameters which are being modulated in
        the MCMC sampling.
        """
                
        ### cluster,mapping,fit_params,inalphas,p_init
        #geoparams = None # I think we don't need this at all....
        
        ### Get initial guess for the pressure of your bulk component(s):

        myval=[]
        myalp=[]
        sbepos=[]  # Boolean array indicating whether myval should be positive or not.
        ### Integral(P*dl) to Compton y (for dl given in radians):
        Pdl2y = (hk.av.szcu['thom_cross']*hk.cluster.d_a/hk.av.szcu['m_e_c2']).to("cm**3 keV**-1")
        R500 = (hk.cluster.R_500/hk.cluster.d_a).decompose()

        tnx = [hk.cfp.bminmax[0].to("rad").value/20.0,10.0*R500.value] # In radians
        nb_theta_range = 150                   # The number of bins for theta_range
        theta_range = np.logspace(np.log10(tnx[0]),np.log10(tnx[1]), nb_theta_range)

        ### Array of fitting values; the order of components are:
        ### (1) mn lvl, (2) bulk, (3) shocks, (4) pt srcs, (5) "blobs"
        ### Within each component, each instrument is addressed. [2] and [3]
        ### are modelled jointly. Maybe I want to do the same with [5], depending
        ### on the supposed physical origin.
        
        for myinst in hk.instruments:
            if ifp[myinst].mn_lvl:
                mnlvlguess = np.mean(dv[myinst].maps.data)
                #import pdb;pdb.set_trace()
                myval.append(mnlvlguess) # Append better for scalars
                sbepos.extend([False])
        
        for bulkbins,fit_cen in zip(hk.cfp.bulkarc,hk.cfp.bulk_centroid):
            a10pres = cpp.a10_gnfw(hk.cluster.P_500,R500,hk.av.mycosmo,bulkbins)
            uless_p = (a10pres*Pdl2y).decompose().value
            myval.extend(uless_p)
            myalp.extend(uless_p*0.0) # What if I want to use strict alphas??
            ### I need to play with centroids (updated how I do the following):
            sbepos.extend(np.ones((len(uless_p)),dtype=bool))
            if fit_cen == True:
                myval.extend([1.0,1.0]) # In particular, how much should I allow this to vary?
                myalp.extend([0.0,0.0]) # In particular, how much should I allow this to vary?
                # What units am I using for this?? (arcseconds? Pixel size? radians?)
                sbepos.extend([False,False])

        for shockbins in hk.cfp.shockbin:
            a10pres = cpp.a10_gnfw(hk.cluster.P_500,R500,hk.av.mycosmo,shockbins)
            uless_p = (a10pres*Pdl2y).decompose().value
            myval.extend(uless_p)
            myalp.extend(uless_p*0.0) # What if I want to use strict alphas??
            sbepos.extend(np.ones((len(uless_p)),dtype=bool))
 
        ### I need to think a bit more about whether this is best or not.
        for myptsrc in hk.cfp.ptsrc:
            for myinst in hk.instruments:
                if ifp[myinst].pt_src == True:
                    myval.extend([0.002]) # Estimate ~2 mJy?? for most point sources...to start.
                    sbepos.extend([True]) # But really we found ~5 mJy for RXJ1347. WTF?

        for myblob in hk.cfp.blob:
            for myinst in hk.instruments:
                if ifp[myinst].fitblob == True:
                    myval.extend([1.0,1.0,1.0]) 
                    sbepos.extend([True,True,True])

        ### Some attributes for starting conditions / model creation.
        self.alphas  = myalp
        self.pinit   = myval
        self.thetas  = theta_range
        self.thetamax= tnx[1]
        self.thetamin= tnx[0]
        self.Pdl2y   = Pdl2y
        self.sbepos  = sbepos
        ### Some attributes for the results:
        self.t_mcmc  = 0.0       # Time that MCMC took.
        self.samples = None
        self.solns   = None
        self.psolns  = None
        self.values  = None
        self.errors  = None
        self.nthreads= np.min([nthreads,ncpus])     # Number of threads to run over with emcee.
        ### I'm not sure if this is a good idea (for being more compact) or not.
        ### I think it should be fine. (WRT ifp variable) 21 Feb 2018.
        self.ifp     = ifp     # Carry around IFP within this class! 
        self.paramdict ={}
        self.punits  = "keV cm**-3"
        self.runits  = "arcsec"
        self.tag = tag
        ######################################################################################

def bulk_or_shock_component(pos,bins,hk,dv,efv,fit_cen,geom,alphas,n_at_rmin,maps={},posind=0,
                            fixalpha=False,fullSZcorr=False,SZtot=False,columnDen=False,Comptony=True,
                            finite=False):

    nbins = len(bins)        
    ulesspres = pos[posind:posind+nbins]
    #myalphas  = alphas[posind:posind+nbins]
    myalphas = alphas   # I've updated how I pass alphas; indexing no longer necessary! (16 Nov 2017)
    ulessrad  = bins #.to("rad").value
    posind = posind+nbins
    if fit_cen == True:
        geom[0:2] = pos[posind:posind+2]  # I think this is what I want...
        posind = posind+2

    density_proxy, etemperature, geoparams = es.prep_SZ_binsky(ulesspres,hk.hk_ins.Tx,geoparams=geom)
    ### Can modify later to allow for X-ray images
    ### That is, I will want to do a loop over SZ images (reduce the number of integrations done),
    ### and then integrate over X-ray emissivity.
    #import pdb;pdb.set_trace()
    
    if fullSZcorr == False:
        #import pdb;pdb.set_trace()
        Int_Pres,outalphas,integrals = es.integrate_profiles(density_proxy, etemperature, geom,bins,
                 efv.thetas,hk,dv,myalphas,beta=0.0,betaz=None,finint=finite,narm=False,fixalpha=fixalpha,
                 strad=False,array="2",SZtot=False,columnDen=False,Comptony=True)
        yint=es.ycylfromprof(Int_Pres,efv.thetas,efv.thetamax) #

    ### I think I can just do "for myinst in hk.instruments:"
    for i in range(len(hk.instruments)):
        myinst = hk.instruments[i]; xymap=dv[myinst].mapping.xymap
        if fullSZcorr == True:
            IntProf,outalphas,integrals = es.integrate_profiles(density_proxy, etemperature, geom,bins,
                 efv.thetas,hk,dv,myalphas,beta=0.0,betaz=None,finint=finite,narm=False,fixalpha=fixalpha,
                 strad=False,array="2",SZtot=True,columnDen=False,Comptony=False)
            ### The following is not really correct. As this is under development, I'll leave it for later
            ### to solve.
            yint=es.ycylfromprof(IntProf,efv.thetas,efv.thetamax) #
            yint=0 # A sure way to give something clearly wrong -> bring myself back here.
            import pdb;pdb.set_trace() # A better way...
        else:
            ConvtoCD= hk.av.szcv["m_e_c2"]/(hk.av.szcv["boltzmann"]*hk.hk_ins.Tx)
            IntProf = Int_Pres * (dv[myinst].tSZ + dv[myinst].kSZ*ConvtoCD)

        ### Right...I need to have a zero-map to start with. Ughhh.
        maps[myinst] += es.general_gridding(xymap,efv.thetas,bins,geom,integrals=integrals,Int_Pres=IntProf)
            
    return maps,posind,yint,outalphas

def mk_twodgauss(pos,posind,centroid,hk,dv,maps={}):

    ### UNDER DEVELOPMENT!
    for myinst in hk.instruments:
        x,y     = dv[myinst].mapping.xymap            # Defined in arcseconds
        s1,s2,peak =  pos[posind:posind+3]
        xcoord  = x - centroid[0];  ycoord=y - centroid[1]
        mygauss = np.exp(((-(xcoord)**2 - (ycoord)**2) )/( 2.0 * sigma**2))
        maps[myinst] += mygauss
        posind += 3

    return maps

### Need to rewrite this to correctly deal with point source.
def mk_ptsrc_v2(pos,posind,hk,dv,ifp,maps={}):

    for index in range(len(hk.cfp.ptsrc)):
        for myinst in hk.instruments:
            ### Centroid is initially defined by RA, Dec.
            if ifp[myinst].pt_src == True:
                myflux  = pos[posind]
                x,y     = dv[myinst].mapping.xymap            # Defined in arcseconds
                fwhm    = dv[myinst].fwhm.to("arcsec").value  # Now this matches.
                if hk.cfp.psfwhm[index] > fwhm:
                    #print "Adjusted FWHM"
                    fwhm = hk.cfp.psfwhm[index]
                sigma   = fwhm/(2.0 * np.sqrt((2.0 * np.log(2.0))))
                #import pdb; pdb.set_trace()
                centroid= dv[myinst].mapping.ptsrclocs[index]
                xcoord  = x - centroid[0];  ycoord=y - centroid[1]
                mygauss = np.exp(((-(xcoord)**2 - (ycoord)**2) )/( 2.0 * sigma**2))
                maps[myinst] += mygauss
                posind += 1   # posind must be incremented for *each* instrument!
        
    return maps,posind

### Need to rewrite this to correctly deal with point source.
def mk_ptsrc(pos,posind,centroid,hk,dv,ifp,maps={}):

    for myinst in hk.instruments:
        ### Centroid is initially defined by RA, Dec.
        if ifp[myinst].pt_src == True:
            myflux  = pos[posind]
            x,y     = dv[myinst].mapping.xymap            # Defined in arcseconds
            fwhm    = dv[myinst].fwhm.to("arcsec").value  # Now this matches.
            sigma   = fwhm/(2.0 * np.sqrt((2.0 * np.log(2.0))))
            xcoord  = x - centroid[0];  ycoord=y - centroid[1]
            mygauss = np.exp(((-(xcoord)**2 - (ycoord)**2) )/( 2.0 * sigma**2))
            maps[myinst] += mygauss
            posind += 1   # posind must be incremented for *each* instrument!
        
    return maps
             
def run_emcee(hk,dv,ifp,efv,init_guess=None):

##################################################################################
### PSA / Reminder: EMCEE doesn't like it if the likelihood doesn't change much with respect
### To how many data points are being used. Or something like that.
    
    def lnlike(pos):                          ### emcee_fitting_vars
        posind = 0
        maps={}
        ### I need to have it return a "pressure profile parameter list only" (15 Feb 2018)
        maps,yint,outalphas = make_skymodel_maps(pos,hk,dv,ifp,efv)
        ytotint = np.sum(yint) # I don't know how I would actually distinguish multiple yints...
        ycyl=ytotint*((hk.cluster.d_a.value/1000.0)**2) # Now in Mpc**-2 (should be a value of ~10 e-5)
        models = filter_maps(hk,dv,maps)
        loglike=0.0
        #print 'LogLike initialized to zero'
        for myinst in hk.instruments:
            data    = dv[myinst].maps.data
            weights = dv[myinst].maps.masked_wts       # Best to mask weights (not use the entire map!)
            model   = models[myinst]
            loglike-= 0.5 * (np.sum(((model - data)**2) * weights))
            ### If I choose to apply a mask:
            mask=np.where(weights > 0)
            gim=model[mask]; gid=data[mask]


        return loglike, outalphas,ycyl

##################################################################################    
    def lnprior(pos,outalphas,ycyl):        
        plike = 0.0
        if hk.cfp.pprior == True:
            plike = -0.5* (np.sum(((ycyl - yplanck)**2) * ywait)) #ywait = ycyl_weight(s)

        ### Pt src priors??
        #mnlvl_amp,ptsrc_amp,blob_amp = unpack_fit_params(pos,hk)
        #ulesspres = pos[0:hk.cfp.bins]
        #if hk.cfp.ptsrc == True:
        #    if len(hk.hk_ins.psfd) > 0:
### How do I want to deal with multiple point sources?
        #        pspr = hk.hk_ins.psfd[0]; psprunc = hk.hk_ins.psunc[0]
        #        plike += -0.5* (np.sum(((ptsrc_amp - pspr)**2)/(psprunc**2) ))

        addlike=0.0; #prespos = pos[1:]
        prespos = pos[efv.sbepos]
        if all([param > 0.0 for param in prespos]) and (outalphas[-1] > 1.0):
            #print 'Everything OK'
            return plike+addlike
        print 'In LnPrior, LogLike set to infinity'
        #import pdb;pdb.set_trace()
        return -np.inf

##################################################################################    
    def lnprob(pos):
        likel,outalphas,ycyl = lnlike(pos)
        if not np.isfinite(likel):
            #print 'In LnProb, LogLike set to infinity'
            return -np.inf
        lp = lnprior(pos,outalphas,ycyl)
        if not np.isfinite(lp):
            #print 'In LnProb, LogLike set to infinity'
            return -np.inf
        return lp + likel

##################################################################################    
    t0 = time.time();    myargs = efv.pinit
    ndim = len(myargs)
    maps,yint,outalphas = make_skymodel_maps(myargs,hk,dv,ifp,efv)
    print myargs
    import pdb;pdb.set_trace()
    pos = [(myargs + np.random.randn(ndim) * myargs/1e3) 
           for i in range(hk.cfp.nwalkers)]

#pos must be shape(nwalkers, ndim)
    t_premcmc = time.time()
    sampler = emcee.EnsembleSampler(hk.cfp.nwalkers, hk.cfp.ndim,
                                    lnprob, threads = efv.nthreads)
    print "Running emcee on ",efv.nthreads," threads."
    for i, result in enumerate(sampler.sample(pos, iterations=hk.cfp.nsteps)):
    #    print i
        if (i+1) % 10 == 0:
            print "{0:5.1%}".format(float(i) / hk.cfp.nsteps)
            
    #import pdb;pdb.set_trace()
    #sampler.run_mcmc(pos,hk.cfp.nsteps)
    
    t_mcmc = time.time() - t_premcmc
    print "MCMC time: ",t_mcmc
    print "Initial Guesses: ", efv.pinit
    
    return sampler, t_mcmc
    
def make_skymodel_maps(pos,hk,dv,ifp,efv):

    posind = 0
    maps={}; yint=[]; outalphas=[]
    # pressure terms.
    for myinst in hk.instruments:
        nx,ny = dv[myinst].mapping.radmap.shape
        maps[myinst]=np.zeros((nx,ny))+ pos[posind]
        posind+=1
    
    ### Model the bulk pressure:
    for bins,fit_cen,geo,alp,narm in zip(hk.cfp.bulkarc,hk.cfp.bulk_centroid,hk.cfp.bulkgeo,
                                         hk.cfp.bulkalp,hk.cfp.bulknarm):
        #print 'Bulk bins: ', bins
        maps,posind,ynt,myalphas = bulk_or_shock_component(pos,bins,hk,dv,efv,fit_cen,geo,alp,narm,
                                                           maps,posind,fixalpha=hk.cfp.bulkfix)
        yint.append(ynt); outalphas.extend(myalphas)
        
    ### Model any shocks:
    for bins,fit_cen,geo,alp,narm in zip(hk.cfp.shockbin,hk.cfp.shoc_centroid,hk.cfp.shockgeo,
                                         hk.cfp.shockalp,hk.cfp.shocknarm):
        #print 'Shock Geometry is: ',geo
        #print 'Shock bins: ',bins
        maps,posind,shint,shout = bulk_or_shock_component(pos,bins,hk,dv,efv,fit_cen,geo,alp,narm,
                                                          maps,posind,fixalpha=hk.cfp.shockfix,
                                                          finite=hk.cfp.shockfin)
        #import pdb;pdb.set_trace()
        
    ### Model any point sources (hk.cfp.ptsrc is a 2-tuple, the pt src. centroid):
    #for myptsrc in hk.cfp.ptsrc:
    #    maps,posind = mk_ptsrc(pos,posind,myptsrc,hk,dv,ifp,maps)

    ### Let's try a second version
    maps,posind = mk_ptsrc_v2(pos,posind,hk,dv,ifp,maps)

    ### Model any "blobs" (2D Gaussians):
    ### This is currently not functional because I'm not sure exactly how I want to implement
    ### this feature.
    for myblob in hk.cfp.blob:
        maps,posind = mk_twodgauss(pos,posind,myptsrc,hk,dv,maps)

    ### Add any mean levels:
    ### (To be added)
    
    return maps,yint,outalphas

def filter_maps(hk,dv,maps):

    models={}
    for myinst in hk.instruments:
        pixs=dv[myinst].mapping.pixsize
        beam_conv = ip.conv_gauss_beam(maps[myinst],pixs,dv[myinst].fwhm)
        models[myinst]=ip.apply_xfer(beam_conv, dv[myinst].mapping.tab,myinst)

    return models

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

def post_mcmc(sampler,t_mcmc,efv,hk,dv):

### Note: I have excluded SAMPLER from the EFV class because it causes the following error:
### Can't pickle <function lnprob at 0x7fb2bd2c6848>: it's not found as max_like_fitting.lnprob

    
    efv.t_mcmc = t_mcmc
    efv.samples = sampler.chain[:,hk.cfp.burn_in:, :].reshape((-1,hk.cfp.ndim))
    efv.solns = np.array(map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                             zip(*np.percentile(efv.samples, [16, 50, 84],
                                                axis=0))))
    efv.psolns = (efv.solns/(efv.Pdl2y)).to(efv.punits)
    #psolns = psolns.to("keV cm**-3")
    efv.values = efv.psolns.value[:,0]
    efv.errors = [efv.psolns.value[:,1],efv.psolns.value[:,2]]
    np.set_printoptions(precision=4)
    hdu = rwrf.saveimg(dv,hk,efv,component='Bulk')
#    rwrf.saveimg(dv,hk,efv,component='Residual')
    
    print "Unitless Pressure values from MCMC are:", efv.solns[:,0]
    print "Actual Pressure values from MCMC are:", efv.values
    savpik = os.path.join(hk.hk_outs.newpath,hk.hk_outs.nnfolder+'_pickle.sav')
    file = open(savpik,'w')
    file.write(pickle.dumps(dv))
    file.write(pickle.dumps(hk))
    file.write(pickle.dumps(efv))
#    file.write(pickle.dumps(sampler))
    file.close

    return hdu
    
#############################################################################################

def unpickle(file=None):

    sdir =myhome+'/Results_Python/MUSTANG2/a2146/'
    file=sdir+'2707_MUSTANG2_6_B_Real_2500S_500B_ML-NO_PP-NO_POWER_30W_pickle.sav'
    rdfl = open(file,'r')
    my_dv = pickle.load(rdfl)
    my_hk = pickle.load(rdfl)
    my_efv = pickle.load(rdfl)
#    my_sampler = pickle.load(rdfl) ... Nope :(

    return my_dv,my_hk,my_efv

def get_initial_guess(pbins,uless_r,hk,conv=1.0):

    vals=pbins
    if hk.hk_ins.name == 'abell_2146':
        vals[0] /= 1.0  # Just spit-balling
        vals[1] /= 1.0  # Just spit-balling
    
    if hk.cfp.ndim > hk.cfp.bins:
        if hk.cfp.mn_lvl == True:
            vals=np.append(vals,0.01)
        if len(hk.cfp.ptsrc) > 0:      # I need to do this better (Feb 2018)
            vals=np.append(vals,0.01)
        if hk.cfp.blob == True:
            vals=np.append(vals,0.01)

#    hk.p_init=vals #??? Does it pass this back? --> Yes; yes it does.
            
    return vals

def get_best_comp_maps(efv,hk,dv,myinst,mycomponent,hdu):

    tstr = 'Test_Run_'
    #if hk.cfp.testmode == False:
    #    tstr = 'Full_Run_'
    tstr = hk.cfp.testmode+'_Run_'
    
    #import pdb;pdb.set_trace()
    modelsky={}
    #for myinst in hk.instruments:
    fbase=tstr+myinst+"_"
    filename=tstr+"Residual.fits"
    fullpath = os.path.join(hk.hk_outs.newpath,hk.hk_outs.prefilename+filename)
    hdulist = fits.open(fullpath)
    ########################### OVERRIDE 06 Feb 2018
    myhdu = hdu[myinst]
    hdulist = fits.HDUList(myhdu)

    modelsky[myinst]=get_best_inst_comp_map(efv,hk,dv,mycomponent,myinst,hdulist)

    return modelsky
        
def get_best_inst_comp_map(efv,hk,dv,mycomponent,myinst,hdulist):

    for hdu in hdulist:
        myext = hdu.header['EXTNAME']
        #print 'hi', myext.lower(), mycomponent.lower()
        if myext.lower() == mycomponent.lower():
            comp_map=hdu.data
        #else:
            #print 'hi', myext.lower(), mycomponent.lower()
            #import pdb;pdb.set_trace()
            
    return comp_map
    
def old_get_best_map(efv,hk,dv):

    uless_fact = (efv.prmunit).to("cm**3 keV**-1").value
    uless_p = efv.values*uless_fact
#    import pdb; pdb.set_trace()
    edensity, etemperature, geoparams = es.prep_SZ_binsky(uless_p,hk.hk_ins.Tx,geoparams=efv.geoparams)
    
    modelsky,outalphas,integrals = es.binsky_SZ_general(edensity, etemperature, geoparams,
                                                        efv.uless_rad,dv.mapping.theta_range,
                                                        dv.mapping.xymap,hk,dv,inalphas=efv.inalphas,
                                                        beta=0.0,betaz=0.0,narm=hk.hk_ins.narm)

    return modelsky

def pack_fit_params(uless_p,hk):

    edensity, etemperature, geoparams = es.prep_SZ_binsky(uless_p,hk.hk_ins.Tx)
    vals=pbins
    vals.extend(etemperature)
    vals.extend(geoparams)
    pos = -1 # Not valid for the moment!
    
    return pos
    
def unpack_fit_params(pos,hk):

    mnlvl_amp=0.0; ptsrc_amp = 0.0; blob_amp = 0.0
    if hk.cfp.ndim > hk.cfp.bins:
        addind=1
        if hk.cfp.mn_lvl == True:
            mnlvl_amp = pos[hk.cfp.bins+addind-1]
            addind+=1
        if len(hk.cfp.ptsrc) > 0:
            ptsrc_amp = pos[hk.cfp.bins+addind-1]
            addind+=1
        if hk.cfp.blob == True:
            blob_amp = pos[hk.cfp.bins+addind-1]
            addind+=1               

    return mnlvl_amp,ptsrc_amp,blob_amp 
