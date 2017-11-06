import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import tSZ_spectrum as tsz
import kSZ_spectrum as ksz
import astropy.units as u
import my_astro_defaults as mad
import astropy.constants as const
#import radec_image as radec
import importlib

############################################################################

            
def inst_params(instrument):

    if instrument == "MUSTANG":
        fwhm1 = 8.7*u.arcsec  # arcseconds
        norm1 = 0.94          # normalization
        fwhm2 = 28.4*u.arcsec # arcseconds
        norm2 = 0.06          # normalization
        fwhm  = 9.0*u.arcsec
        smfw  = 10.0*u.arcsec
        freq  = 90.0*u.gigahertz # GHz
        FoV   = 42.0*u.arcsec #
        
    if instrument == "MUSTANG2":
        fwhm1 = 8.7*u.arcsec  # arcseconds
        norm1 = 0.94          # normalization
        fwhm2 = 28.4*u.arcsec # arcseconds
        norm2 = 0.06          # normalization
        fwhm  = 9.0*u.arcsec
        smfw  = 10.0*u.arcsec
        freq  = 90.0*u.gigahertz # GHz
        FoV   = 4.25*u.arcmin * (u.arcmin).to("arcsec")
        
    if instrument == "NIKA":
        fwhm1 = 8.7*2.0*u.arcsec  # arcseconds
        norm1 = 0.94     # normalization
        fwhm2 = 28.4*2.0*u.arcsec # arcseconds
        norm2 = 0.06     # normalization
        fwhm  = 18.0*u.arcsec
        smfw  = 10.0*u.arcsec
        freq  = 150.0*u.gigahertz    # GHz
        FoV   = 2.15*u.arcmin * (u.arcmin).to("arcsec")
        
    if instrument == "NIKA2":
        fwhm1 = 8.7*2.0*u.arcsec  # arcseconds
        norm1 = 0.94     # normalization
        fwhm2 = 28.4*2.0*u.arcsec # arcseconds
        norm2 = 0.06     # normalization
        fwhm  = 18.0*u.arcsec
        smfw  = 10.0*u.arcsec
        freq  = 150.0*u.gigahertz    # GHz
        FoV   = 6.5*u.arcmin * (u.arcmin).to("arcsec")
        
    if instrument == "BOLOCAM":
        fwhm1 = 8.7*7.0*u.arcsec  # arcseconds
        norm1 = 0.94     # normalization
        fwhm2 = 28.4*7.0*u.arcsec # arcseconds
        norm2 = 0.06     # normalization
        fwhm  = 58.0*u.arcsec
        smfw  = 60.0*u.arcsec
        freq  = 140.0*u.gigahertz    # GHz
        FoV   = 8.0*u.arcmin * (u.arcmin).to("arcsec")

#    else:
#        fwhm1=9.0*u.arcsec ; norm1=1.0
#        fwhm2=30.0*u.arcsec ; norm2=0.0
#        fwhm = 9.0*u.arcsec ; smfw = 10.0*u.arcsec
#        freq = 90.0*u.gigahertz 
#        FoV   = 1.0*u.arcmin * (u.arcmin).to("arcsec")
#        
#    import pdb; pdb.set_trace()

    return fwhm1,norm1,fwhm2,norm2,fwhm,smfw,freq,FoV

def inst_bp(instrument,array="2"):
    """
    Returns a frequency and bandpass array
    
    Parameters
    __________
    instrument : MUSTANG, MUSTANG2, BOLOCAM, NIKA, or NIKA2
    currently only MUSTANG2 and NIKA2 are supported

    Returns
    -------
    -> farr   - The frequency array (in GHz)
    -> band   - the bandpass, with aperture efficiency applied
    """

    if instrument == "MUSTANG2" or instrument == "MUSTANG":
        srms = (300*u.um).to("m")        # surface RMS (microns)
        ### Reference: https://science.nrao.edu/facilities/gbt/proposing/GBTpg.pdf
        EA90 = 0.36   # Aperture efficiency at 90 GHz
        ### The beam efficiencies should be taken as 1.37* Aperture Efficiency
        R90  = np.exp(-4.0*np.pi*(srms/(const.c/(9.0e10*u.s**-1))).value)    #
        Gnot = EA90/R90                   # Unphysical, but see documentation...
        if instrument == "MUSTANG2":
            flow = 75.0   # GHz
            fhig = 105.0  # GHz
        else:
            flow = 82.5   # GHz
            fhig = 97.5  # GHz
            
        farr = np.arange(flow,fhig,1.0)  # frequency array.
        tran = farr*0.0 + 1.0            # Let the transmission be unity everywhere.
        Larr = const.c.value/(farr*1.0e9) # Keep calm and carry on.        
        Ruze = Gnot * np.exp(-4.0*np.pi*(srms.value)/Larr)
        NRuz = Ruze / np.max(Ruze)        # Normalize it
        band = tran * Ruze                # Bandpass, with (unnormalized) Ruze efficiency
       
    if instrument == "NIKA2" or instrument == "NIKA":
        caldir='/home/romero/NIKA2/NIKA_SVN/Processing/Pipeline/Calibration/BP/'
        bpfile=caldir+'Transmission_2017_Jan_NIKA2_v1.fits'
        hdulist = fits.open(bpfile)

        if array == "1H":      # 1mm (260 GHz) array, Horizontal Polarization
            tbdata = hdulist[1].data # 1H
            freq = tbdata.field(0)
            tran = tbdata.field(1)
            erro = tbdata.field(2)
            atmt = tbdata.field(3)
            cfreq1h = np.sum(freq*tran)/np.sum(tran)
        
        if array == "1V":     # 1mm (260 GHz) array, Vertical Polarization
            tbdata = hdulist[2].data # 1V
            freq = tbdata.field(0)
            tran = tbdata.field(1)
            erro = tbdata.field(2)
            atmt = tbdata.field(3)
            cfreq1v = np.sum(freq*tran)/np.sum(tran)
        
        if array == "2":       # 2mm (150 GHz) array
            tbdata = hdulist[3].data # 2
            freq = tbdata.field(0)
            tran = tbdata.field(1)
            erro = tbdata.field(2)
            atmt = tbdata.field(3)
            cfreq2 = np.sum(freq*tran)/np.sum(tran)

        ### Trim the zero-frequency listing, if any.
        gi=np.where(freq > 0)
        freq = freq[gi]
        tran = tran[gi]
        erro = erro[gi]
        atmt = atmt[gi]
        
### Calculate Aperture efficiencies from information found at:
### http://www.iram.es/IRAMES/mainwiki/Iram30mEfficiencies
        Beff = 0.630         # at 210 GHz
        Aeff = Beff/1.27     # See text on webpage
        srms = (66.0*u.um).to("m")        # surface RMS (microns)
        R210 = np.exp(-4.0*np.pi*(srms/(const.c/(2.1e11*u.s**-1))).value)    #
        Gnot = Aeff/R210                   # Unphysical, but see documentation...

        Larr = const.c.value/(freq*1.0e9) # Keep calm and carry on.        
        Ruze = Gnot * np.exp(-4.0*np.pi*(srms.value)/Larr)
        NRuz = Ruze / np.max(Ruze)        # Normalize it
        band = tran * Ruze                # Bandpass, with (unnormalized) Ruze efficiency
        farr = freq
        
#########################################################################

    return band, farr

class maps:

    def __init__(self,hk):
     
        data_map, header = fits.getdata(hk.hk_ins.fitsfile, header=True)
        wt_map = fits.getdata(hk.hk_ins.wtfile,hk.hk_ins.wtext)
        if hk.hk_ins.wtisrms == True: 
            wt_map = wt_map**(-2.0)  # I'll need to account for zeros...
        self.data   = data_map
        self.header = header
        self.wts    = wt_map
        if not(hk.hk_ins.psfile == None):
            self.ptsrc,self.pshdr = fits.getdata(hk.hk_ins.psfile, header=True)
        else:
            self.ptsrc,self.pshdr = None,None
        if not(hk.hk_ins.shfile == None):
            self.shock,self.shhdr = fits.getdata(hk.hk_ins.shfile, header=True)
        else:
            self.shock,self.shhdr = None,None
        if not(hk.hk_ins.blfile == None):
            self.blob,self.blhdr = fits.getdata(hk.hk_ins.blfile, header=True)
        else:
            self.blob,self.blhdr = None,None
        if not(hk.hk_ins.miscfile1 == None):
            self.misc1,self.misc1hdr = fits.getdata(hk.hk_ins.miscfile1, header=True)
        else:
            self.misc1,self.mis1hdr = None,None
        if not(hk.hk_ins.miscfile2 == None):
            self.misc2,self.misc2hdr = fits.getdata(hk.hk_ins.miscfile2, header=True)
        else:
            self.misc2,self.mis2hdr = None,None
        if not(hk.hk_ins.miscfile3 == None):
            self.misc3,self.misc3hdr = fits.getdata(hk.hk_ins.miscfile3, header=True)
        else:
            self.misc3,self.mis3hdr = None,None


def get_sz_bp_conversions(temp,instrument,array="2",inter=False,beta=1.0/300.0,
                          betaz=1.0/300.0,rel=True):

    szcv,szcu=mad.get_sz_values()
    freq_conv = (szcv['planck'] *1.0e9)/(szcv['boltzmann']*szcv['tcmb'])
    temp_conv = 1.0/szcv['m_e_c2']
    bv = get_beamvolume(instrument)
#    fwhm1,norm1,fwhm2,norm2,fwhm,smfw,freq,FoV = inst_params(instrument)
    band, farr = inst_bp(instrument,array)
    fstep = np.median(farr - np.roll(farr,1))
    foo = np.where(band > 0)
    band=band[foo]
    farr=farr[foo]
############################################
    llt  = temp*temp_conv  # Lower limit on temperature (well, theta)
    ult  = temp*temp_conv  # Upper limit on temperature
    st   = temp            # Temperature (theta) step
    flow = np.min(farr)
 #   fhigh= np.ceil((np.max(farr)-flow)/fstep)*fstep+flow+fstep/100.0
    fhigh= np.max(farr)
    sx   = fstep*freq_conv
    llx  = flow*freq_conv
    nste = (fhigh-flow)/fstep
    ulx  = fhigh*freq_conv + sx/(2.0*nste)
    
#### The old way:
#    temparr,freqarr,conv = tsz.tSZ_conv_range(tlow,thigh,tstep,flow,fhigh,fstep)
#### The new way (for now):

    if inter == True:    
        tarr = np.arange(llt,ult,st)
        xarr = np.arange(llx,ulx,sx)
        T = tsz.itoh_2004_int(tarr,xarr)
    else:
        tarr, xarr, T = tsz.tSZ_conv(llt,llx,ult,ulx,st,sx)
        
    tarr, xarr, K = ksz.kSZ_conv(beta,betaz,llt,llx,ult,ulx,st,sx,rel=rel)
    
#    import pdb; pdb.set_trace()
### Let's check that the frequency spacing is close to what was given in the
### bandpass retreival:
    fq = np.abs(farr*freq_conv - xarr)/(farr*freq_conv)
    if np.max(fq) > 0.05:
        raise Exception("Frequency arrays differ significantly.")
## Else (implicit here):
    TY = T/tarr    # Divide by tarr (thetae) to get proper conversion units
    KY = K/tarr    # Divide by tarr (thetae) to get proper conversion units
    bT = np.sum(TY*band)/np.sum(band)  # Average over the bandpass
    bK = np.sum(KY*band)/np.sum(band)  # Average over the bandpass

    JypB = tsz.Jyperbeam_factors(bv)
    tSZ_JyBeam_per_y = JypB * bT  # Just multiply by Compton y to get Delta I (tSZ)
    kSZ_JyBeam_per_t = JypB * bK  # Just multiply by tau (of electrons) to get Delta I (kSZ)

    print 'Assuming a temperature of ',temp,' keV, we find the following.'
    print 'To go from Compton y to Jy/Beam, multiply by: ', tSZ_JyBeam_per_y
    print 'To go from tau to Jy/Beam (kSZ), multiply by: ', kSZ_JyBeam_per_t
    
#    import pdb; pdb.set_trace()
        
    return tSZ_JyBeam_per_y.value, kSZ_JyBeam_per_t.value

def get_maps_and_info(instrument,target,real=True):
    """
    Would like to deprecate this.
    """

    fitsfile, wtfile, wtext, wtisrms, tab = get_fits_files(instrument,target,real=True)
    data_map, header = fits.getdata(fitsfile, header=True)
    wt_map = fits.getdata(wtfile,wtext)
    if wtisrms == True:
        wt_map = 1.0/wt_map**2

    image_data, ras, decs, hdr, pixs = get_astro(fitsfile)
    w = WCS(fitsfile)
    
    return data_map, wt_map, header, ras, decs, pixs, w, tab

#def get_map_info(file):
#
#    image_data, ras, decs, hdr, pixs = get_astro(file)

def get_beamvolume(instrument):

    fwhm1,norm1,fwhm2,norm2,fwhm,smfw,freq,FoV = inst_params(instrument)

    sc = 2.0 * np.sqrt(2.0*np.log(2))
    sig1 = fwhm1/sc
    sig2 = fwhm2/sc
    bv1 = 2.0*np.pi * norm1*sig1**2
    bv2 = 2.0*np.pi * norm2*sig2**2
    beamvolume = bv1 + bv2  # In units of FWHM**2 (should be arcsec squared) 

    return beamvolume

def get_sz_conversion(temp,instrument,beta=0.0,betaz=0.0):
    """
    Returns the tSZ and kSZ conversions to Jy/beam for a given instrument,
    for a given electron temperature, velocity relative to the speed of
    light (beta), and beta along the line of sight (betaz).
    
    Parameters
    __________
    instrument : MUSTANG, BOLOCAM, or NIKA
    target: The target of the cluster
    real: True or False (do you want to analyze/fit real data or
             'virtual' (simulated) data?).import numpy as np

    Returns
    -------
    CtSZ, CkSZ
    """

    bv = get_beamvolume(instrument)
    fwhm1,norm1,fwhm2,norm2,fwhm,smfw,freq,FoV = inst_params(instrument)
 #   conv = tsz.tSZ_conv_single(temp,freq.value)
    bpsr = tsz.intensity_factors(freq.value,bv)
    szcv,szcu=mad.get_sz_values()
    freq_conv = (szcv['planck'] *1.0e9)/(szcv['boltzmann']*szcv['tcmb'])
    temp_conv = 1.0/szcv['m_e_c2']
    conv = tsz.itoh_2004_int(temp*temp_conv,freq*freq_conv)

    print conv,bpsr
    yJyBeam = conv.item(0)*bpsr/(temp*temp_conv) # Conversion from Compton y to Jy/beam

    return yJyBeam

###################################################################
### Under verification below.

def astro_from_hdr(hdr):
    
    xsz = hdr['naxis1']
    ysz = hdr['naxis2']
    xar = np.outer(np.arange(xsz),np.zeros(ysz)+1.0)
    yar = np.outer(np.zeros(xsz)+1.0,np.arange(ysz))
    ####################
    
    xcen = hdr['CRPIX1']
    ycen = hdr['CRPIX2']
    dxa = xar - xcen
    dya = yar - ycen
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
    
    return ras, decs, pixs

def get_astro(file):

    hdu = fits.open(file)
    hdr = hdu[0].header
    image_data = hdu[0].data

    ras, decs, pixs = astro_from_hdr(hdr)

    return image_data, ras, decs, hdr, pixs

class astrometry:

    def __init__(self,hdr):

        ras,decs,pixs = astro_from_hdr(hdr)
        self.ras = ras
        self.decs= decs
        self.pixs= pixs
        
def get_xfer(inputs):

    if inputs.tabformat == 'ascii':
        tab = np.loadtxt(inputs.tabfile, comments=inputs.tabcomments)
    if inputs.tabformat == 'fits':
        ktab = fits.getdata(inputs.tabfile)
        xtab = fits.getdata(inputs.tabfile,ext=1)
        tab = np.vstack((ktab,xtab))
        
    if inputs.tabdims == '1D':
        if inputs.instrument == "MUSTANG" or inputs.instrument == "MUSTANG2":
            tab = tab.T            # Transpose the table.
#        import pdb;pdb.set_trace()
        if inputs.tabextend==True:
            tdim = tab.shape
#            import pdb;pdb.set_trace()
            pfit = np.polyfit(tab[0,tdim[1]/2:],tab[1,tdim[1]/2:],1)
            addt = np.max(tab[0,:]) * np.array([2.0,4.0,8.0,16.0,32.0])
            extt = np.polyval(pfit,addt)
            foo = np.stack((addt,extt))
            tab = np.concatenate((tab,foo),axis=1)
#            newt = [np.append(tab[0,:],addt),np.append(tab[1,:],extt)]
            
    return tab

############################################################################

