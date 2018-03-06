import my_astro_defaults as mad
import cosmolopy.distance as cd
import astropy.units as u
from astropy.coordinates import Angle
import numpy as np
import astropy.constants as const
import retrieve_data_info as rdi   # Finds and loads data files (hard-coded)
import pressure_profiles as pp     # Creates radius + pressure arrays.
from astropy.wcs import WCS

#        import pdb; pdb.set_trace()

def clust_info(name=None):
    #######################################################
    ### CLASH CLUSTERS: (RA / DEC FROM ACCEPT)
    ### ACCEPT: http://www.pa.msu.edu/astro/MC2/accept/
    if name == 'a1835':              # From Mantz+ 2010
        z=0.253                      # Redshift
        ra = Angle('14h1m1.951s')    # Hour angle
        dec= Angle('2d52m43.18s')    # Degrees
        M_500=1.2e15                 # Solar masses
        Tx = 9.0                     # keV
    if name == 'a611':               # From Mantz+ 2010
        z=0.288                      # Redshift
        ra = Angle('8h0m56.832s')    # Hour angle
        dec= Angle('36d3m24.09s')    # Degrees
        M_500=7.4e14                 # Solar masses
        Tx = 6.8                     # keV
    if name == 'm1115':              # From Mantz+ 2010
        z=0.355                      # Redshift
        ra = Angle('11h15m52.048s')  # Hour angle
        dec= Angle('1d29m56.56s')    # Degrees
        M_500=8.6e14                 # Solar masses
        Tx = 9.2                     # keV
    if name == 'm0429':              # From Mantz+ 2010
        z=0.399                      # Redshift
        ra = Angle('4h29m36.088s')   # Hour angle
        dec= Angle('-2d53m9.02s')    # Degrees
        M_500=5.8e14                 # Solar masses
        Tx = 8.3                     # keV
    if name == 'm1206':              # From Mantz+ 2010
        z=0.439                      # Redshift
        ra = Angle('12h6m12.276s')   # Hour angle
        dec= Angle('-8d48m2.40s')    # Degrees
        M_500=1.9e15                 # Solar masses
        Tx = 10.7                    # keV
    if name == 'm0329':              # From Mantz+ 2010
        z=0.450                      # Redshift
        ra = Angle('3h29m41.681s')   # Hour angle
        dec= Angle('-2d11m47.67s')   # Degrees
        M_500=7.9e14                 # Solar masses
        Tx = 6.3                     # keV
    if name == 'rxj1347':            # From Mantz+ 2010
        z=0.451                      # Redshift
        ra = Angle('13h47m30.593s')  # Hour angle
        dec= Angle('-11d45m10.05s')  # Degrees
        M_500 = 2.2e15               # Solar masses
        Tx = 10.8                    # keV
    if name == 'm1311':              # From Mantz+ 2010
        z=0.494                      # Redshift
        ra = Angle('13h11m1.665s')   # Hour angle
        dec= Angle('-3d10m39.50s')   # Degrees
        M_500=3.9e14                 # Solar masses
        Tx = 6.0                     # keV
    if name == 'm1423':              # From Mantz+ 2010
        z=0.543                      # Redshift
        ra = Angle('14h23m47.9s')    # Hour angle
        dec= Angle('24d4m43s')       # Degrees
        M_500=6.6e14                 # Solar masses
        Tx = 6.9                     # keV
    if name == 'm1149':              # From Mantz+ 2010
        z=0.544                      # Redshift
        ra = Angle('11h49m35.856s')  # Hour angle
        dec= Angle('22d23m55.02s')   # Degrees
        M_500=1.9e15                 # Solar masses
        Tx = 8.5                     # keV
    if name == 'm0717':              # From Mantz+ 2010
        z=0.546                      # Redshift
        ra = Angle('7h17m31.654s')   # Hour angle
        dec= Angle('37d45m18.52s')   # Degrees
        M_500=2.5e15                 # Solar masses
        Tx = 11.8                    # keV
    if name == 'm0647':              # From Mantz+ 2010
        z=0.591                      # Redshift
        ra = Angle('6h47m50.029s')   # Hour angle
        dec= Angle('70d14m49.66s')   # Degrees
        M_500=1.1e15                 # Solar masses
        Tx = 11.5                    # keV
    if name == 'm0744':              # From Mantz+ 2010
        z=0.698                      # Redshift
        ra = Angle('7h44m52.802s')   # Hour angle
        dec= Angle('39d27m24.41s')   # Degrees
        M_500=1.3e15                 # Solar masses
        Tx = 8.1                     # keV
    if name == 'clj1226':            # From Mantz+ 2010
        z=0.888                      # Redshift
        ra = Angle('12h26m58.373s')  # Hour angle
        dec= Angle('33d32m47.36s')   # Degrees
        M_500=7.8e14                 # Solar masses
        Tx = 12.0                    # keV
    if name == 'ms0735':            # From Mantz+ 2010
        z=0.216                      # Redshift
        ra = Angle('07h41m50.2s')  # Hour angle
        dec= Angle('74d14m51.0s')   # Degrees
        M_500=4.8e14                 # Solar masses
        Tx = 10.0                    # keV
    #######################################################
    ### Do any conversions desired before returning values
    ra = ra.to("deg")
    
    return z, ra, dec, M_500, Tx
        
class cluster:

    def __init__(self,hk=None,name="Unknown"):
        """
        Returns a structure with important variables related to the cluster parameters 
        (which physically described the cluster), as well as parameters related to our
        viewing of the cluster (e.g. angular diameter) that depend on cosmology.

        Parameters
        __________
        name      - The name of the cluster
        M_500     - The mass enclosed within R_500
        z         - The redshift of the cluster
        ra        - The Right Ascenscion (in hours or degrees; degrees preferred)
        dec       - The Declination (in degrees)

        Returns
        -------
        The structure cluster
        """

    ### Get general cosmological parameters:
        
        cosmo,mycosmo=mad.get_cosmology()
        sz_constants=mad.get_sz_values()

        if not(hk == None):
            z,ra,dec,M_500,Tx = hk.hk_ins.z,hk.hk_ins.ra,hk.hk_ins.dec,\
                                hk.hk_ins.M_500,hk.hk_ins.Tx
            name = hk.hk_ins.name
        else:
            z, ra, dec, M_500, Tx = clust_info(name) #hk.hk_ins.name

        H0 = mycosmo['H0']
        h_70 = mycosmo['h_70'].value
        self.name = name
        self.z = z
        self.E = (mycosmo['omega_m']*(1 + self.z)**3 + mycosmo['omega_l'])**0.5
        self.H = H0 * self.E
        self.dens_crit = (3 * (self.H)**2)/(8 * np.pi * const.G)
        self.M_500 = M_500 * const.M_sun
        self.P_500 =(1.65 * 10**-3) *((self.E)**(8./3)) *((
            self.M_500 * h_70 )/((3*10**14)  * const.M_sun)
            )**(2./3.)*(h_70)**2  *u.keV /u.cm**3
        self.R_500 =(3 * self.M_500/(4 * np.pi  * 500 * self.dens_crit))**(1/3.)
        self.R_500 = self.R_500.to(u.kpc)
        self.R_max = 5 * self.R_500
        self.d_a = cd.angular_diameter_distance(self.z, **cosmo) *u.Mpc
        self.d_a = self.d_a.to("kpc") #convert to kpc (per radian)
        self.scale = self.d_a*(u.arcsec).to("rad") / u.arcsec
        self.theta_500 = (self.R_500 *u.rad/ self.d_a).to("arcmin")
        self.theta_max = (self.R_max *u.rad /self.d_a).to("arcmin")
        self.ra_deg=ra.to("deg")
        self.dec_deg=dec
        self.Tx = Tx

class mapping:

    def __init__(self,cluster,hk,fit_params,mycosmo,tSZ,kSZ,instrument,pixsize=1.0,
                 ltrmax=5.0,N=120):
        """
        Returns a structure with important variables related to gridding a map.
        In this sense, I have called it "astrometry", even though every variable
        may not classically fit within astrometry.

        Parameters
        __________
        cluster   - A structure, from the class/routing above in this file.
        pixsize   - The pixel size in arcseconds, but without units in Python
        ltrmax    - How many times "theta_max" should the profile be calculated for?
                    [The default is 5.]
        N         - The number of (logarithmic) points for which to calculate the
                    integrated profile.

        Returns
        -------
        The structure mapping
        """
        ############################################################################################
        ### First we do need to calculate a few variables
        
        rintmax=cluster.R_500.value
        theta_min=(0.2 * u.arcsec).to("radian").value
        ltm=np.log10(theta_min/2.0)          # Log of Theta Min (radians)
        ### Theta_max is already 5* R_500
        ltx=np.log10(cluster.theta_max.to("radian").value)  # Log of Theta Max (radians)
        w = WCS(hk.hk_ins.fitsfile)
        x0,y0=w.wcs_world2pix(cluster.ra_deg,cluster.dec_deg,0)
        image_data, ras, decs, hdr, pixs = rdi.get_astro(hk.hk_ins.fitsfile)
        xymap = pp.get_xymap(image_data,pixs,xcentre=x0.item(0),ycentre=y0.item(0))
        arcmap = pp.get_radial_map(image_data,pixs,xcentre=x0.item(0),ycentre=y0.item(0))  # In arcseconds
        x,y = xymap
        angmap = np.arctan2(y,x)
#        import pdb; pdb.set_trace()
        radmap = (arcmap*u.arcsec).to("rad").value            # In radians
        rmval = radmap; bi=np.where(rmval < theta_min); rmval[bi]=theta_min
        radmap = rmval
        r_bins, pressure = pp.get_default_a10_profile(instrument,cluster,mycosmo,nbins=fit_params.bins)
        ############################################################################################
        ### And now we can put variables into our structure.
        ######################################################### But first, another quick comment:
        ### ltrmax = Log(Theta_Range)_max - for LTRMAX*R_500
        ### Thus, the default is 5*R_500
        self.theta_min=theta_min
        self.theta_max= (15.0* u.arcmin).to("radian").value  # 15' in radians.
        #self.theta_range = np.logspace(ltm,ltx+np.log10(ltrmax), N)
        self.theta_range = np.logspace(ltm,ltx, N)
        self.w=w
        self.ra=cluster.ra_deg
        self.dec=cluster.dec_deg
        self.pixsize=pixs
        self.x0=x0
        self.y0=y0
        self.tSZ=tSZ
        self.kSZ=kSZ
        self.tab=rdi.get_xfer(hk)
        self.tabdims=hk.hk_ins.tabdims     # Better to include this variable, if necessary
        self.r_bins=r_bins           # Need just the value (do not want a "quantity")
        self.a10_pressure = pressure # Need just the value (do not want a "quantity")        
        self.radmap=radmap           # Need just the value (do not want a "quantity")
        self.arcmap=arcmap
        self.angmap=angmap
        self.xymap=xymap
        self.instrument=instrument


    
