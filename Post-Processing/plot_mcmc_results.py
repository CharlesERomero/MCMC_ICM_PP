import matplotlib.pyplot as plt
import retrieve_data_info as rdi
import max_like_fitting as mlf
import instrument_processing as ip
import numpy as np
import os, re
import corner
import astropy.units as u          # Install astropy
import collect_variables as cv
import save_image as si
import ellipsoidal_shells as es
import matplotlib.ticker as ticker

myfontsize=20

def plot_res_no_sampler(hk,dv,efv,overlay=None):
##################################### I'm not using this generally!

    tstr = 'Test_Run_'
    #if hk.cfp.testmode == False:
    #    tstr = 'Full_Run_'
    tstr = hk.cfp.testmode
    tstr = tstr+efv.tag

    r_bins= hk.cfp.bulkarc # Bulkarc is defined in radians.
    arcbins = r_bins * (u.rad).to("arcsec")  # added this line 1 Mar 2018
    #arcbins = (r_bins*u.kpc/hk.cluster.d_a)*(u.rad).to(efv.runits)
    #arcbins = arcbins.value #...*sigh*
#    plot_steps(sampler,hk.cfp,hk.hk_outs.newpath,hk.hk_outs.prefilename+tstr,efv)
### I need to develop a way to sift through multiple bulk profiles, multiple shocks.
    plot_pres_bins(arcbins,efv,hk,hk.cluster,tstr,overlay=overlay)
    pres_samples = ((efv.samples/(efv.Pdl2y)).to(efv.punits)).value # In keV cm**-3
    plot_correlations(pres_samples,hk.hk_outs.newpath, hk.hk_outs.prefilename+tstr)
    maxlikesky=mlf.get_best_comp_maps(efv,hk,dv)
    plot_best_sky(maxlikesky,hk.hk_outs.newpath, hk.hk_outs.prefilename+tstr,dv,
                  mycomp='data')
##########################################################################

def plot_results(hk,dv,efv,sampler,hdu,ifp,overlay=None):

    tstr = 'Test_Run_'
    #if hk.cfp.testmode == False:
    #    tstr = 'Full_Run_'
    tstr = hk.cfp.testmode
    tstr = tstr+'_'+efv.tag

    ibulk=0; ishock=0; iptsrc=0   # I'm not sure what to do here.
    #mydatamaps={}
    plot_comps_and_resids(hk,dv,efv,hdu,tstr,ifp,overlay=None)
    
    ### I am 99.99% sure I don't want to loop over instruments here.
    for myinst in hk.instruments:
 
        r_bins= hk.cfp.bulkarc[ibulk] # Bulkarc is defined in radians.
        arcbins = r_bins * (u.rad).to("arcsec")
        plot_steps(sampler,hk.cfp,hk.hk_outs.newpath,hk.hk_outs.prefilename+tstr,efv)
        
        ###############################################################################
        ### Need to list keys in efv.paramdict[inst]
        mycomponent = 'Bulk'; mycount=1
        plot_pres_bins(arcbins,efv,hk,hk.cluster,tstr+'bulk_',center=False,overlay=overlay,
                       inst=myinst,count=mycount,component=mycomponent)
        ###############################################################################
        #import pdb;pdb.set_trace()
        if len(hk.cfp.shockbin) > 0:
            mycomponent = 'Shock'; mycount=1
            plot_pres_bins(arcbins,efv,hk,hk.cluster,tstr+'shock_',overlay=overlay,inst=myinst,count=mycount,
                           component=mycomponent)
        ###############################################################################
        ###
        pres_samples = ((efv.samples/(efv.Pdl2y)).to(efv.punits)).value # In keV cm**-3
        plot_correlations(pres_samples,hk.hk_outs.newpath, hk.hk_outs.prefilename+tstr)
        ###############################################################################
        ### plot_best_sky already loops over instrument keys - so I need to reorganize this.
        ### (March 7, 2018)
        ###

def plot_comps_and_resids(hk,dv,efv,hdu,tstr,ifp,overlay=None):
    
    ibulk=0; ishock=0; iptsrc=0   # I'm not sure what to do here.
    pos_comps  = set(efv.compname) #.intersection(['bulk','ptsrc','shock','blob'])
    comp_maps={}               # Here, dictionary by component and instrument
    mydatamaps={}; residual={}; filtimg={}; wbi={}; doublecheck={} # And here, by just instrument
    pos = efv.solns[:,0]       # Correct units...the unitless kind, mostly

    mnlvlmaps,ptsrcmaps,maps,yint,outalphas = mlf.make_skymodel_maps(pos,hk,dv,ifp,efv)
    ytotint = np.sum(yint) # I don't know how I would actually distinguish multiple yints...
    ycyl=ytotint*((hk.cluster.d_a.value/1000.0)**2) # Now in Mpc**-2 (should be a value of ~10 e-5)
    models = mlf.filter_maps(hk,dv,maps,ptsrcmaps,mnlvlmaps)

    for myinst in hk.instruments:
        mydatamaps[myinst]=dv[myinst].maps.data
        residual[myinst] = np.zeros(mydatamaps[myinst].shape) + mydatamaps[myinst]

        for component in pos_comps:
            ### REMINDER: Each component map will have maps *PER INSTRUMENT* too!

            #print component
            ### How should I figure out if there are multiple of a component?
            mycount=1
            
            comp_maps[component],wbi[myinst] = mlf.get_best_comp_maps(efv,hk,dv,myinst,component+str(mycount),
                                                          hdu,returnwcs=True)
            filtimg[myinst] = plot_best_sky(comp_maps[component],hk.hk_outs.newpath, hk.hk_outs.prefilename+tstr,
                                            dv,mycomp=component,count=mycount,w=wbi)
            residual[myinst] = mydatamaps[myinst] - filtimg[myinst]
            plot_best_sky(residual,hk.hk_outs.newpath, hk.hk_outs.prefilename+tstr,dv,
                          mycomp='residual_after_'+component,w=wbi)

        doublecheck[myinst] = mydatamaps[myinst]-models[myinst]

    ############################################################################
    plot_best_sky(mydatamaps,hk.hk_outs.newpath, hk.hk_outs.prefilename+tstr,
                  dv,mycomp='data',w=wbi)
    plot_best_sky(residual,hk.hk_outs.newpath, hk.hk_outs.prefilename+tstr,dv,
                  mycomp='residual',w=wbi)
    plot_best_sky(doublecheck,hk.hk_outs.newpath, hk.hk_outs.prefilename+tstr,dv,
                  mycomp='residual_dbl_chk',w=wbi)
    #import pdb;pdb.set_trace()

def plot_steps(sampler,fit_params,newpath,pre_filename,efv):
    stepmap    = plt.figure(1,figsize=(20,len(efv.compname))); plt.clf()
    pos_comps  = set(efv.compname)
    pos_colors = ['k','g','r','c','m','y']
    assoc_colo = {}
    
    for i,x in enumerate(pos_comps):
        assoc_colo[x] = pos_colors[i]

    
    for i in range(fit_params.ndim):
        ax = stepmap.add_subplot(fit_params.ndim+1,1,i+1)
        color = assoc_colo[efv.compname[i]]
        #import pdb;pdb.set_trace()
        ax.plot(np.array([sampler.chain[:,j,i] for j in range(fit_params.nsteps)]),color)
        isgtz = (sampler.chain[:,:,i] > 0)
        isgtz1d = isgtz.reshape(np.product(isgtz.shape))
        ax.get_xaxis().set_visible(False) # Maybe works?
        #import pdb;pdb.set_trace()
        if i == 0:
            ylims = ax.get_ylim()
            yval  = (ylims[1] - ylims[0])*0.5 + ylims[1]
            for j,comp in enumerate(pos_comps):
                xval = float(fit_params.nsteps*j)/(len(pos_comps)+0.5)
                ax.text(xval, yval, comp, color=assoc_colo[comp],fontsize=20)
            xval = float(fit_params.nsteps*len(pos_comps))/(len(pos_comps)+0.5)
            ax.text(xval, yval, "likelihood", color="b",fontsize=20)

        
        if all(isgtz1d):
            ax.set_yscale("log", nonposy='clip')
            ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y,pos: ('{{:.{:1d}f}}'.format(int(np.maximum(-np.log10(y),0)))).format(y)))
            #ax.get_yaxis().set_visible(False) # Maybe works?
        ax.set_ylabel((r'$P_{0}$',r'$P_{1}$',r'$P_{2}$',r'$P_{3}$',r'$P_{4}$',r'$P_{5}$',r'$P_{6}$',r'$P_{7}$',
                       r'$P_{8}$',r'$P_{9}$',r'$P_{10}$',r'$P_{11}$',r'$P_{12}$',r'$P_{13}$',r'$P_{14}$',
                       r'$P_{15}$',r'$P_{16}$',r'$P_{17}$',r'$P_{18}$',r'$P_{19}$',r'$P_{20}$',r'$P_{21}$',
                       r'$P_{22}$',r'$P_{23}$',r'$P_{24}$')[i],fontsize=myfontsize)
    ax = stepmap.add_subplot(fit_params.ndim+1,1,fit_params.ndim+1)
    ax.plot(np.array([sampler._lnprob[:,j] for j in range(fit_params.nsteps)]),"b")
    myyr = [np.min(sampler._lnprob[:,fit_params.burn_in*2:])*1.1,np.max(sampler._lnprob[:,fit_params.burn_in*2:])*0.9]
    ax.set_ylim(myyr)
    ax.set_ylabel(r'$ln(\mathcal{L}$)',fontsize=myfontsize)
  
    plt.xlabel('Steps',fontsize=myfontsize)
    filename = "step.png"
    fullpath = os.path.join(newpath, pre_filename+filename)
#2106_MUSTANG_6_B_Real_200S_40B_ML-NO_PP-NO_POWER_20W/
#    fullpath='/home/romero/Results_Python/plots_to_show/MUSTANG_Real_step.png'
    plt.savefig(fullpath)
    #plt.close(stepmap)
    plt.close()

### Now, radarr should be able to be removed (1 Mar 2018)
def plot_pres_bins(radarr,efv,hk,cluster,tstr,center=True,overlay=None,inst=None,count=1,
                   component='Bulk'):

    fit_params=hk.cfp; plt.figure(2,figsize=(20,12));    plt.clf()
    if inst == None:
        inst=hk.instruments[0]

    punits=efv.punits; runits=efv.runits
    fwhm1,norm1,fwhm2,norm2,fwhm,smfw,freq,FoV = rdi.inst_params(inst)
    valstoplot=efv.values   # Could probably delete this line
    errstoplot=efv.errors   # Could probably delete this line
    #psolns = efv.psolns
    compname = component+str(count)
    psolns   = (efv.solns[efv.paramdict[inst][compname+'_ind'],:]/(efv.Pdl2y)).to(efv.punits)
    ### But I don't use this???
    #import pdb;pdb.set_trace()
    binsolns = efv.paramdict[inst][compname+'_bins'] * (u.rad).to("arcsec") # Now correct for plotting.
    
    nbins = fit_params.bulkbins[count-1] # How many bins did I have in this model?
    errstoplot= [psolns[:,1].value,psolns[:,2].value]
    valstoplot= psolns[:,0].value #already just a value (not quantity)
    xerr=None
    #cents=radarr * (u.rad).to("arcsec") # Now correct for plotting.
    cents = binsolns
    rin = fwhm.to("arcsec").value / 2.0; rout = FoV.to("arcsec").value / 2.0; axcol = 'r'
    if center == True:
        if hk.cfp.shockfin[count-1] == True and component =='Shock':
            #print 'hey dude'
            cents   = (binsolns[:-1] + binsolns[1:])/2.0
            diffs   = np.diff(binsolns)
            xerr    = [diffs/2.0,diffs/2.0]
        else:
            #print 'normal stuff'
            lradarr = np.append(0,binsolns[:-1])
            cents   = (lradarr + binsolns)/2.0
            offlow  = cents - lradarr
            offhig  = binsolns - cents
            xerr    = [offlow,offhig]
            
        plt.errorbar(cents,valstoplot,xerr=xerr,yerr=errstoplot,fmt='.',label="Deprojected profile",capsize=5)
    else:
        plt.errorbar(cents,valstoplot,xerr=0,yerr=errstoplot,fmt='.',label="Deprojected profile",capsize=5)
        radii =  np.logspace(np.log10(rin / 5.0), np.log10(rout*3.0),150)
        pprof,alphas = es.log_profile(valstoplot,binsolns,radii)
        #import pdb; pdb.set_trace()
        plt.plot(radii, pprof, label='Power law interpolation')
        
        #plt.plot(range_r, clj1226.pprofile(range_r, P0true, rptrue, 
#         atrue,ctrue,btrue), "b", label = "Nominal Profile")
#if map_type != "BOLOCAM":
#    plt.axvline(FoV.value,color="r", linestyle ="dashed")
#    rin = fwhm.to("arcsec").value / 2.0; rout = FoV.to("arcsec").value / 2.0; axcol = 'r'
    if hk.hk_ins.name == "abell_2146":
        rin = (hk.hk_ins.rads_nw[1]*u.kpc/cluster.scale).value
        rout= (hk.hk_ins.rads_nw[2]*u.kpc/cluster.scale).value

    if overlay == 'input':
        overplot_input(hk,efv,cluster)

    if overlay == 'russell':
        overplot_russell(efv,cluster)
        
    #print cents
    #print 'Inner and outer radii... supposed bin limits: ',rin, rout
    ### I want to change this....
    plt.axvline(rin,color=axcol, linestyle ="dashed")
    plt.axvline(rout,color=axcol, linestyle ="dashed")
    plt.yscale("log")
    plt.xscale("log")
    plt.xlabel("Radius ("+runits+")",fontsize=myfontsize)
    plt.ylabel("Pressure ("+punits+")",fontsize=myfontsize)
    plt.title(cluster.name)
    plt.grid()
    filename = tstr+"pressure.png"
    fullpath = os.path.join(hk.hk_outs.newpath,hk.hk_outs.prefilename+filename)
    plt.savefig(fullpath)
    filename = tstr+"pressure.eps"
    fullpath = os.path.join(hk.hk_outs.newpath,hk.hk_outs.prefilename+filename)
    plt.savefig(fullpath,format='eps')
    #import pdb;pdb.set_trace()
    plt.close()

def plot_correlations(samples,newpath, pre_filename):

    plt.figure(2,figsize=(20,12))
    plt.clf()
    fig = corner.corner(samples, bins = 45,quantiles=[0.16,0.50,0.84],
                        labels=["$P_{1}$","$P_{2}$","$P_{3}$","$P_{4}$","$P_{5}$","$P_{6}$","$P_{7}$","$P_{8}$",
                                "$P_{9}$","$P_{10}$","$P_{11}$","$P_{12}$","$P_{13}$","$P_{14}$","$P_{15}$","$P_{16}$",
                                r'$P_{17}$',r'$P_{18}$',r'$P_{19}$',r'$P_{20}$',r'$P_{21}$',r'$P_{22}$',r'$P_{23}$',r'$P_{24}$'],
                        fontsize=myfontsize)
    filename = "correlations_via_corner.png"
    fullpath = os.path.join(newpath, pre_filename+filename)
    plt.savefig(fullpath)
    plt.close()

def plot_best_sky(maxlikesky,newpath, pre_filename,dv,mycomp="bulk",count=1,w=None):

    mapaxisunits="arcsec"   #...Should add as a keyword
    mapunits="Jy/Beam"      # A default unit.
    myformat="png"
    countstr = str(count); compcount = mycomp+countstr
    if mycomp == 'data' or re.match('residual*',mycomp):
        compcount = mycomp
        
    for key in maxlikesky:
        #print 'Preparing to plot the best sky model for the component ',mycomp
        #print 'As imaged for the instrument ',key
        instrument=key
        mapunits = dv[key].maps.units   # What the units should actually be.
        image = maxlikesky[key]
        wtmax = np.max(dv[key].maps.masked_wts)
        #if mapunits == 'Jy':
        #    wtmax*=1.0e6
        rmsmin= wtmax**(-0.5); mymin = -10.0*rmsmin; mymax = 10.0*rmsmin
        #print rmsmin,mymin,mymax
        cblim = False;        islog = False;        flipsign = False
        if mycomp == 'data' or re.match('residual*',mycomp):
            title = "Unsmoothed "+compcount+" model "
            cblim = True
        else:
            title = "Unfiltered "+compcount+" model "
        if mycomp == 'ptsrc':
            title = "Beamshape "+compcount+" model "
        if mycomp == 'bulk' or mycomp == 'shock':
            islog = True      # I might need to make the map positive...
            nzi = image != 0
            if np.min(image[nzi]) < 0 and np.max(image[nzi])<0:
                flipsign=True
                mymax = np.max(-1.0*image[nzi]) ### mymin = np.min(-1.0*image[nzi]);
                mymin = mymax / 1.0e4           ### I just think this is better.
            else:
                mymin = np.min(image[nzi]); mymax = np.max(image[nzi])                
            cblim=True; oldmapunits=mapunits; mapunits='-'+mapunits
        filename = pre_filename+"flux_density_"+compcount+"_skymodel"
        fullpath = os.path.join(newpath, filename+'.'+myformat)
        if type(w) == type(None):
            plot_sky_map(fullpath,image,title,mapaxisunits,mapunits)
        else:
            wti = w[key]; mapmask=0 # define some BS variable for mapmask.
            si.plot_image_wcs(image, mapmask, wti,dpi=200,myfontsize=15,zoom=True,filename=filename,
                              plotmask=False,savedir=newpath,cblim=cblim,mymin=mymin,mymax=mymax,mytitle=title,
                              target=dv[key].maps.name,format=myformat,islog=islog,cbtit=mapunits,cbar=True,
                              subtitle=key,cmap='viridis',flipsign=flipsign)
##################################################################################################

        if mycomp == 'bulk' or mycomp == 'shock':
            cblim=False; flipsign = False; mapunits = oldmapunits
        if mycomp != 'mnlvl':
            #print 'Proceeding to filter, in some manner, the map'
            if mycomp == 'ptsrc':
                beam_conv = maxlikesky[key]
            else:
                beam_conv = ip.conv_inst_beam(maxlikesky[key],dv[instrument].mapping.pixsize,instrument=instrument)
            image = beam_conv
            title = "Smoothed "+compcount+" model "
            if mycomp != 'data' and not(re.match('residual*',mycomp)):
                gc_model=ip.apply_xfer(beam_conv, dv[instrument].mapping.tab,instrument=instrument)
                image = gc_model
                title = "Filtered "+compcount+" model "

            mapaxisunits="arcsec"   #...Should add as a keyword
            filename = pre_filename+"flux_density_"+compcount+"_filtered"
            if mycomp == 'data' or re.match('residual*',mycomp):
                filename = pre_filename+"flux_density_"+compcount+"_smoothed"
                mymin = -5.0*rmsmin; mymax = 5.0*rmsmin
            
            fullpath = os.path.join(newpath, filename+'.'+myformat)
            if type(w) == type(None):
                plot_sky_map(fullpath,image,title,mapaxisunits,mapunits)
            else:
                wti = w[key]; mapmask=0 # define some BS variable for mapmask.
                si.plot_image_wcs(image, mapmask, wti,dpi=200,myfontsize=15,zoom=True,filename=filename,
                                  plotmask=False,savedir=newpath,cblim=cblim,mymin=mymin,mymax=mymax,mytitle=title,
                                  target=dv[key].maps.name,format=myformat,islog=False,cbtit=mapunits,cbar=True,
                                  subtitle=key,cmap='viridis')
        ### (end if mycomp!= mnlvl:)
##################################################################################################

        if mycomp == 'data' or re.match('residual*',mycomp):
            weightmap = dv[instrument].maps.wts
            #if mapunits == 'Jy':
            #    weightmap= dv[instrument].maps.wts *1.0e6
            wt_conv = ip.conv_inst_beam(weightmap,dv[instrument].mapping.pixsize,instrument=instrument)
            #bv = rdi.get_beamvolume(instrument) / (dv[instrument].mapping.pixsize**2)
            bv = dv[instrument].bv / (dv[instrument].mapping.pixsize**2)
            #print bv
            #import pdb;pdb.set_trace()

            ### Is this at all right?
            wt_conv*= bv.value # Weight has increased with smoothing...
            nzi = (wt_conv > 0)
            inv_rms_map = (wt_conv)**(0.5)
            snr_map = np.zeros(inv_rms_map.shape)
            snr_map[nzi] = beam_conv[nzi] * inv_rms_map[nzi]
            title = "Smoothed "+compcount+" SNR map "
            filename = pre_filename+"flux_density_"+compcount+"_smoothed_SNR"
            fullpath = os.path.join(newpath, filename+'.'+myformat)
            mapunits='SNR'
            if type(w) == type(None):
                plot_sky_map(fullpath,snr_map,title,mapaxisunits,mapunits)
            else:
                wti = w[key]; mapmask=0 # define some BS variable for mapmask.
                si.plot_image_wcs(snr_map, mapmask, wti,dpi=200,myfontsize=15,zoom=True,filename=filename,
                                  plotmask=False,savedir=newpath,cblim=False,mytitle=title,
                                  target=dv[key].maps.name,format=myformat,islog=False,cbtit=mapunits,cbar=True,
                                  subtitle=key,cmap='bwr')

        return image
    #### End of loop.

def plot_sky_map(fullpath,image,title,mapaxisunits,mapunits,format='png'):

    plt.figure(2,figsize=(20,12));        plt.clf()
    plt.imshow(image) # May want to add "extent"
    plt.title(title)
    plt.xlabel(mapaxisunits,fontsize=myfontsize); plt.ylabel(mapaxisunits,fontsize=myfontsize)
    cbar = plt.colorbar();    cbar.set_label(mapunits,fontsize=myfontsize)
    plt.savefig(fullpath,format=format)
    plt.close()
    
def overplot_input(hk,efv,cluster):

    nw = True
#    if efv.sector == 'se':
#        nw = False
    
    if nw == True:
        radarr = np.append(hk.hk_ins.rads_nw, 500.0)
        eden = hk.hk_ins.eden_nw; temp = hk.hk_ins.temp_nw
        ealp = hk.hk_ins.ealp_nw; talp = hk.hk_ins.talp_nw
        title= "Abell 2146; Northwest (upstream) shock"
    else:
        radarr = np.append(hk.hk_ins.rads_se, 500.0)
        eden = hk.hk_ins.eden_se; temp = hk.hk_ins.temp_se
        ealp = hk.hk_ins.ealp_se; talp = hk.hk_ins.talp_se
        title= "Abell 2146; Southeast (bow) shock"

    radarr = (radarr / cluster.scale).value
        
    for idx in range(len(radarr)-1):
        if idx == 0:
            rmin = 5.0
            rmax = radarr[1]
            darr = np.array([(rmin/rmax)**ealp[idx],1.0])*eden[idx]
            tarr = np.array([(rmin/rmax)**talp[idx],1.0])*temp[idx]
        else:
            rmin = radarr[idx]
            rmax = radarr[idx+1]
            darr = np.array([1.0,(rmax/rmin)**ealp[idx]])*eden[idx]
            tarr = np.array([1.0,(rmax/rmin)**talp[idx]])*temp[idx]
            
        plt.plot([rmin,rmax],darr*tarr,"g",label = "Input Pressure",fontsize=myfontsize)

    return title

def overplot_russell(efv,cluster):

    nw = True
#    if efv.sector == 'se':
#        nw = False

    Rfile = '/home/romero/Python/model_creation/Abell_2146_Russell_parresults_Fig11.dat'
    prof = np.loadtxt(Rfile, comments='#')

    myind = np.arange(9,14,1)
    if nw == True:
        myind = np.arange(9)
        
    rads = prof[myind,0]; rerr= prof[myind,1]; xerr=[rerr,rerr]
    temp = prof[myind,2]; tlow= prof[myind,3]; thig= prof[myind,4]
    dens = prof[myind,5]; dlow= prof[myind,6]; dhigh=prof[myind,7]

    pres = temp*dens
    plow= np.sqrt((tlow*dens)**2 + (temp*dlow)**2)
    phigh= np.sqrt((thig*dens)**2 + (temp*dhig)**2)
    yerr = [plow,phig]
    
    plt.errorbar(rads,pres,xerr=xerr,yerr=yerr,fmt='dr',label="Russell et al. 2012",capsize=3)
    #import pdb;pdb.set_trace()
    #cv.print_attributes(fit_params)
#    psolns = (efv.paramdict[inst][compname]/(efv.Pdl2y)).to(efv.punits)

