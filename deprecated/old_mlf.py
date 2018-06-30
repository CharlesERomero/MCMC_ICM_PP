
        
def run_emcee_old(hk,dv,efv,init_guess=None):

    #data    = dv.maps.data
    #weights = dv.maps.wts

    if hk.cfp.mask == True:
        bi=np.where(dv.mapping.arcmap > hk.cfp.maskrad)
        weights[bi] = 0.0
        if hk.hk_ins.name == 'abell_2146':
            bi = np.where(dv.mapping.angmap < hk.hk_ins.lban_nw)
            weights[bi] = 0.0
            bi = np.where(dv.mapping.angmap > hk.hk_ins.uban_nw)
            weights[bi] = 0.0
       
    if init_guess == None:
        myargs = efv.p_init # Just use them all
    else:
        myargs = init_guess
    yplanck = 1.0        # Bogus value for now
    ywait = 1.0         # Bogus value for now

    #if hk.cfp.ptsrc == True:
    #    ptsrc_map, ptsrc_hdr = fits.getdata(hk.cfp.psfile, header=True)

    def lnlike(pos):                          ### emcee_fitting_vars
##################################################################################
### Sequencing of parameters passed into LNLIKE:
### [Pressure bins, [mean level],[pt_src],[blob],[shock]]
        ulesspres = pos[0:hk.cfp.bins]
        add_models = data*0.0
        if hk.cfp.ndim > hk.cfp.bins:
            addind=1
            if hk.cfp.mn_lvl == True:
                mnlvl_model=(data*0.0+1.0)*pos[hk.cfp.bins+addind-1]
                addind+=1; add_models+= mnlvl_model
            if hk.cfp.ptsrc == True:
                ptsrc_model=ptsrc_map*pos[hk.cfp.bins+addind-1]
                addind+=1; add_models+= ptsrc_model
            if hk.cfp.blob == True:
                blob_model=ptsrc_map*pos[hk.cfp.bins+addind-1]
                addind+=1; add_models+= blob_model                
            
    ##################################################################################
    ### The "old" way, which is rather restrictive on the geometry.
#        modelsky,outalphas,integrals = es.binsky(ulesspres,efv.uless_rad,dv.mapping.theta_range,
#                                                 dv.mapping.radmap,inalphas=efv.inalphas)
#        modelsky *= dv.mapping.tSZ;
    ##################################################################################
    ### The "new" way: allows for ellispoids, and centroid shift(s)
        edensity, etemperature, geoparams = es.prep_SZ_binsky(ulesspres,hk.hk_ins.Tx,geoparams=efv.geoparams)
        modelsky,outalphas,integrals = es.binsky_SZ_general(edensity, etemperature, efv.geoparams,
                                                            efv.uless_rad,efv.theta_range,
                                                            dv.mapping.xymap,hk,dv,inalphas=efv.inalphas,
                                                            beta=0.0,betaz=0.0,narm=hk.hk_ins.narm)
    ##################################################################################
        Int_Pres=np.sum(integrals,axis=0)
        yint=es.ycylfromprof(Int_Pres,dv.mapping.theta_range,dv.mapping.theta_max) #
        ycyl=yint*((hk.cluster.d_a.value/1000.0)**2) # Now in Mpc**-2 (should be a value of ~10 e-5)
        beam_conv = ip.conv_inst_beam(modelsky,dv.mapping.pixsize*u.arcsec,instrument=dv.mapping.instrument)
 #       gc_model=ip.apply_xfer(beam_conv, dv.mapping.tab,instrument=dv.mapping.instrument)
 ### Wait! I've made a mistake. The point source map is not filtered. Thus I want to add that in
 ### first, and then filter the combined model.
        unfilt_model = beam_conv + add_models # Unfiltered models
        model=ip.apply_xfer(unfilt_model, dv.mapping.tab,instrument=dv.mapping.instrument)

#        import pdb; pdb.set_trace()
#        model = gc_model+add_models  ### (The previous way; see above comments)
        loglike= - 0.5 * (np.sum(((model - data)**2) * weights))
        mask=np.where(weights > 0)
        gim=model[mask]; gid=data[mask]
#        print '###########################################################################################'
#        print np.min(gim),np.max(gim),np.min(gid),np.max(gid)
#        print np.min(modelsky[mask]),np.max(modelsky[mask]),np.min(beam_conv[mask]),np.max(beam_conv[mask])
#        print '-------------------------------------------------------------------------------------------'

        return loglike, outalphas,ycyl,ulesspres

    def lnprior(pos,outalphas,ycyl):
        
        plike = 0.0
        if hk.cfp.pprior == True:
            plike = -0.5* (np.sum(((ycyl - yplanck)**2) * ywait)) #ywait = ycyl_weight(s)

        mnlvl_amp,ptsrc_amp,blob_amp = unpack_fit_params(pos,hk)
        ulesspres = pos[0:hk.cfp.bins]
        if hk.cfp.ptsrc == True:
            if len(hk.hk_ins.psfd) > 0:
### How do I want to deal with multiple point sources?
#                import pdb; pdb.set_trace()
                pspr = hk.hk_ins.psfd[0]; psprunc = hk.hk_ins.psunc[0]
                plike += -0.5* (np.sum(((ptsrc_amp - pspr)**2)/(psprunc**2) ))

        addlike=0.0    
        
        if all([param > 0.0 for param in ulesspres]) and (outalphas[-1] > 1.0):        
            return plike+addlike
        return -np.inf

    def lnprob(pos):
        likel,outalphas,ycyl,pb = lnlike(pos)
        if not np.isfinite(likel):
            return -np.inf
        lp = lnprior(pos,outalphas,ycyl)
        if not np.isfinite(lp):
            return -np.inf
        return lp + likel


    t0 = time.time(); ulesspres = myargs[0:hk.cfp.bins]
#    import pdb; pdb.set_trace()
#    modelsky,outalphas,integrals = es.binsky(ulesspres,efv.uless_rad,dv.mapping.theta_range,
#                                             dv.mapping.radmap,inalphas=efv.inalphas)
#    modelsky *= dv.mapping.tSZ
#    geoparams = None
#    if hk.hk_ins.name == 'abell_2146':
#        uless_r, edensity, etemperature, geoparams, inalphas = es.prep_a2146_binsky(hk.hk_ins,nw=True)
    edensity, etemperature, geoparams = es.prep_SZ_binsky(ulesspres,hk.hk_ins.Tx,geoparams=efv.geoparams)
    modelsky,outalphas,integrals = es.binsky_SZ_general(edensity, etemperature, efv.geoparams,
                                                        efv.uless_rad,dv.mapping.theta_range,
                                                        dv.mapping.xymap,hk,dv,inalphas=efv.inalphas,
                                                        beta=0.0,betaz=0.0,narm=hk.hk_ins.narm)
    if hk.cfp.ubmap == True:
        startmap = 0.0 # .... Need to write something here (20 June 2017)

    Int_Pres=np.sum(integrals,axis=0)
    yint=es.ycylfromprof(Int_Pres,dv.mapping.theta_range,dv.mapping.theta_max) # will be in radians^-2
    ycyl=yint*((hk.cluster.d_a.value/1000.0)**2) # Now in Mpc**-2 (should be a value of ~10 e-5)
    beam_conv = ip.conv_inst_beam(modelsky,dv.mapping.pixsize,instrument=dv.mapping.instrument)
    model=ip.apply_xfer(beam_conv, dv.mapping.tab,instrument=dv.mapping.instrument)

    rwrf.saveimg(dv,hk,efv,modelmap=model,weightmap=weights,
                 filename="Inital_Guess_Map.fits")
    
#    model=ip.convsky(modelsky,pix_sigma,dv.mapping.tab,map_type=map_type,convolve=True,
#                     joint=joint,jbeam=True,beammap=beammap)       

    print "time to make map: ", time.time()-t0

#################################################################################

#    import pdb; pdb.set_trace()

    pos = [(myargs + np.random.randn(hk.cfp.ndim) * myargs/1e3) 
           for i in range(hk.cfp.nwalkers)]

#pos must be shape(nwalkers, ndim)
    t_premcmc = time.time()
    print "Input Pressures:", ulesspres

    sampler = emcee.EnsembleSampler(hk.cfp.nwalkers, hk.cfp.ndim,
                                    lnprob, threads = 1)
    sampler.run_mcmc(pos,hk.cfp.nsteps)
    
    t_mcmc = time.time() - t_premcmc
    print "MCMC time: ",t_mcmc
    print "Input Pressures:", ulesspres

    
    return sampler, t_mcmc
