import matplotlib.pyplot as plt
import retrieve_data_info as rdi
import max_like_fitting as mlf
import instrument_processing as ip
import numpy as np
import os
import corner
import astropy.units as u          # Install astropy
import collect_variables as cv

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
#    plot_steps(sampler,hk.cfp,hk.hk_outs.newpath,hk.hk_outs.prefilename+tstr)
### I need to develop a way to sift through multiple bulk profiles, multiple shocks.
    plot_pres_bins(arcbins,efv,hk,hk.cluster,tstr,overlay=overlay)
    pres_samples = ((efv.samples/(efv.Pdl2y)).to(efv.punits)).value # In keV cm**-3
    plot_correlations(pres_samples,hk.hk_outs.newpath, hk.hk_outs.prefilename+tstr)
    maxlikesky=mlf.get_best_comp_maps(efv,hk,dv)
    import pdb;pdb.set_trace()
    plot_best_bulk(maxlikesky,hk.hk_outs.newpath, hk.hk_outs.prefilename+tstr,dv)

def plot_results(hk,dv,efv,sampler,hdu,overlay=None):

    tstr = 'Test_Run_'
    #if hk.cfp.testmode == False:
    #    tstr = 'Full_Run_'
    tstr = hk.cfp.testmode
    tstr = tstr+'_'+efv.tag

    ibulk=0   # I'm not sure what to do here.
    for myinst in hk.instruments:
        ### Set up the necessary variables for the model map(s) of each instrument
        #nx,ny = dv[myinst].mapping.radmap.shape
        #maps[myinst]=np.zeros((nx,ny))+ pos[posind]
        #posind+=1
        #hdu[myinst] = []

        r_bins= hk.cfp.bulkarc[ibulk] # Bulkarc is defined in radians.
        arcbins = r_bins * (u.rad).to("arcsec")
        #arcbins = (r_bins*u.kpc/hk.cluster.d_a)*(u.rad).to(efv.runits)
        #arcbins = arcbins.value #...*sigh*
        plot_steps(sampler,hk.cfp,hk.hk_outs.newpath,hk.hk_outs.prefilename+tstr)
        
        ###############################################################################
        ### Need to list keys in efv.paramdict[inst]
        mycomponent = 'Bulk'; mycount=1
        plot_pres_bins(arcbins,efv,hk,hk.cluster,tstr+'bulk_',overlay=overlay,inst=myinst,count=mycount,
                       component=mycomponent)
        ###############################################################################
        mycomponent = 'Shock'; mycount=1
        plot_pres_bins(arcbins,efv,hk,hk.cluster,tstr+'shock_',overlay=overlay,inst=myinst,count=mycount,
                       component=mycomponent)
        ###############################################################################
        ###
        pres_samples = ((efv.samples/(efv.Pdl2y)).to(efv.punits)).value # In keV cm**-3
        plot_correlations(pres_samples,hk.hk_outs.newpath, hk.hk_outs.prefilename+tstr)
        maxlikesky=mlf.get_best_comp_maps(efv,hk,dv,myinst,"bulk"+str(ibulk+1),hdu)
        plot_best_bulk(maxlikesky,hk.hk_outs.newpath, hk.hk_outs.prefilename+tstr,dv)

def plot_steps(sampler,fit_params,newpath,pre_filename):
    stepmap = plt.figure(1)
    for i in range(fit_params.ndim):
        ax = stepmap.add_subplot(fit_params.ndim+1,1,i+1)
        ax.plot(np.array([sampler.chain[:,j,i] for j in range(fit_params.nsteps)]),"k")
        #    ax.get_yaxis().set_visible(False) # Maybe works?
        ax.set_ylabel((r'$P_{0}$',r'$P_{1}$',r'$P_{2}$',r'$P_{3}$',r'$P_{4}$',
                       r'$P_{5}$',r'$P_{6}$',r'$P_{7}$',r'$P_{8}$',r'$P_{9}$',r'$P_{10}$',r'$P_{11}$')[i],)
    ax = stepmap.add_subplot(fit_params.ndim+1,1,fit_params.ndim+1)
    ax.plot(np.array([sampler._lnprob[:,j] for j in range(fit_params.nsteps)]),"b")
    myyr = [np.min(sampler._lnprob[:,fit_params.burn_in*2:])*1.1,np.max(sampler._lnprob[:,fit_params.burn_in*2:])*0.9]
    ax.set_ylim(myyr)
    ax.set_ylabel(r'$ln(\mathcal{L}$)')
  
    plt.xlabel('Steps')
    filename = "step.png"
    fullpath = os.path.join(newpath, pre_filename+filename)
#2106_MUSTANG_6_B_Real_200S_40B_ML-NO_PP-NO_POWER_20W/
#    fullpath='/home/romero/Results_Python/plots_to_show/MUSTANG_Real_step.png'
    plt.savefig(fullpath)

### Now, radarr should be able to be removed (1 Mar 2018)
def plot_pres_bins(radarr,efv,hk,cluster,tstr,center=True,overlay=None,inst=None,count=1,
                   component='Bulk'):

    fit_params=hk.cfp; plt.figure(2);    plt.clf()
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
    if center == True:
        lradarr = np.append(0,binsolns[:-1])
        cents   = (lradarr + binsolns)/2.0
        offlow  = cents - lradarr
        offhig  = binsolns - cents
        xerr    = [offlow,offhig]
        #lradarr = np.append(0,radarr[:-1])
        #cents = (lradarr + radarr)*(u.rad).to("arcsec")/2.0
        #offlow = cents - lradarr*(u.rad).to("arcsec")
        #offhig = radarr*(u.rad).to("arcsec") - cents
        #xerr = [offlow,offhig]
#################################################################
    print cents
    #import pdb; pdb.set_trace()
    plt.errorbar(cents,valstoplot,xerr=xerr,yerr=errstoplot,fmt='.',label="Deprojected profile",capsize=5)
#plt.plot(range_r, clj1226.pprofile(range_r, P0true, rptrue, 
#         atrue,ctrue,btrue), "b", label = "Nominal Profile")
#if map_type != "BOLOCAM":
#    plt.axvline(FoV.value,color="r", linestyle ="dashed")
    rin = FoV.value; rout = fwhm.value; axcol = 'r'
    if hk.hk_ins.name == "abell_2146":
        rin = (hk.hk_ins.rads_nw[1]*u.kpc/cluster.scale).value
        rout= (hk.hk_ins.rads_nw[2]*u.kpc/cluster.scale).value

    if overlay == 'input':
        overplot_input(hk,efv,cluster)

    if overlay == 'russell':
        overplot_russell(efv,cluster)
        
    plt.axvline(rin,color=axcol, linestyle ="dashed")
    plt.axvline(rout,color=axcol, linestyle ="dashed")
    plt.yscale("log")
    plt.xscale("log")
    plt.xlabel("Radius ("+runits+")")
    plt.ylabel("Pressure ("+punits+")")
    plt.title(cluster.name)
    plt.grid()
    filename = tstr+"pressure.png"
    fullpath = os.path.join(hk.hk_outs.newpath,hk.hk_outs.prefilename+filename)
    plt.savefig(fullpath)
    filename = tstr+"pressure.eps"
    fullpath = os.path.join(hk.hk_outs.newpath,hk.hk_outs.prefilename+filename)
    plt.savefig(fullpath,format='eps')

def plot_correlations(samples,newpath, pre_filename):

    plt.figure(2)
    plt.clf()
    fig = corner.corner(samples, bins = 45,quantiles=[0.16,0.50,0.84],
                        labels=["$P_{1}$","$P_{2}$","$P_{3}$","$P_{4}$","$P_{5}$","$P_{6}$","$P_{7}$","$P_{8}$",
                                "$P_{9}$","$P_{10}$","$P_{11}$","$P_{12}$"])
    filename = "correlations_via_corner.png"
    fullpath = os.path.join(newpath, pre_filename+filename)
    plt.savefig(fullpath)

def plot_best_bulk(maxlikesky,newpath, pre_filename,dv,count=1):

    mapaxisunits="arcsec"   #...Should add as a keyword
    for key in maxlikesky:
        instrument=key
        plt.figure(2)
        plt.clf()
        plt.imshow(maxlikesky[key]) # May want to add "extent"
        plt.title("Model Compton y")
        #plt.xlabel("arcsec")
        #plt.ylabel("arcsec")
        plt.xlabel(mapaxisunits)
        plt.ylabel(mapaxisunits)
        cbar = plt.colorbar() 
        cbar.set_label('mJy / beam')
        filename = "flux_density_skymodel.png"
        fullpath = os.path.join(newpath, pre_filename+filename)
        plt.savefig(fullpath)

        beam_conv = ip.conv_inst_beam(maxlikesky[key],dv[instrument].mapping.pixsize,instrument=instrument)
        gc_model=ip.apply_xfer(beam_conv, dv[instrument].mapping.tab,instrument=instrument)
    
        mapaxisunits="arcsec"   #...Should add as a keyword
        plt.figure(2)
        plt.clf()
        plt.imshow(gc_model) # May want to add "extent"
        plt.title("Model Compton y")
        #plt.xlabel("arcsec")
        #plt.ylabel("arcsec")
        plt.xlabel(mapaxisunits)
        plt.ylabel(mapaxisunits)
        cbar = plt.colorbar() 
        cbar.set_label('mJy / beam')
        filename = "flux_density_filtered.png"
        fullpath = os.path.join(newpath, pre_filename+filename)
        plt.savefig(fullpath)
        
    #### End of loop.

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
            
        plt.plot([rmin,rmax],darr*tarr,"g",label = "Input Pressure")

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

