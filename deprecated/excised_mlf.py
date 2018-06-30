        ##############################################################################################
        ### Let's do some double checks. (02 April 2018)
        #
        #for myinst in hk.instruments:
        #    maps[myinst]      = np.zeros((nx,ny))
        #origposind = posind
        #print posind
        #maps,posind,shint,shout = bulk_or_shock_component(pos,bins,hk,dv,efv,fit_cen,geo,alp,narm,
        #                                                  maps,posind,fixalpha=hk.cfp.shockfix,
        #                                                  finite=sfin,oldvs=True)
        #minangle = geo[2]-0.1; maxangle=geo[2]+0.1
        #mybins=np.arange(0.0,120.0,4.0)
        #for myinst in hk.instruments:
        #    rads, prof = ABP.get_az_profile(maps[myinst], dv[myinst].mapping.xymap, minangle, maxangle)
        #    binres = ABP.radial_bin(rads, prof+0.0002,10,rmax=120.0,bins=mybins,minangle=minangle,maxangle=maxangle)
        #    fig = ABP.plot_one_slice(binres,myformat='png',fig = None,target='RXJ1347_Shock',savedir=hk.hk_outs.newpath,
        #                             prefilename='Testing_tapers_', mylabel="Old_version")
        #
        #for myinst in hk.instruments:
        #    maps[myinst]      = np.zeros((nx,ny))
        #posind = origposind
        #print posind

        # Here I had bulk_or_shock_component()
        
        #for myinst in hk.instruments:
        #    rads, prof = ABP.get_az_profile(maps[myinst], dv[myinst].mapping.xymap, minangle, maxangle)
        #    binres = ABP.radial_bin(rads, prof,10,rmax=120.0,bins=mybins,minangle=minangle,maxangle=maxangle)
        #    fig = ABP.plot_one_slice(binres,myformat='png',fig = fig,target='RXJ1347_Shock',savedir=hk.hk_outs.newpath,
        #                             prefilename='Testing_tapers_', mylabel="New_version")
        #
        #import pdb;pdb.set_trace()
        ##############################################################################################

