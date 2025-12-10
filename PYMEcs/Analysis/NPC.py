import matplotlib.pyplot as plt
import numpy as np
from PYMEcs.pyme_warnings import warn
import logging

logger = logging.getLogger(__name__)

piover4 = np.pi/4.0

# from circle_fit import taubinSVD
# def fitcirc(x,y,sigma=None):
#     pcs = np.vstack((x,y)).T
#     xc, yc, r, sigma = taubinSVD(pcs)
#     return (xc, yc, r, sigma)

# def centreshift(x,y,xc,yc):
#     return (x-xc,y-yc)

# def plot_segments(rad,rotang=0):
#     for i in range(8):
#         ang = i * np.pi / 4.0
#         plt.plot([0,rad*np.cos(ang+rotang)],[0,rad*np.sin(ang+rotang)],'b')

# def phi_from_coords(xn,yn):
#     phis = np.arctan2(yn,xn)
#     return phis

# def rot_coords(xn,yn,ang):
#     c, s = np.cos(ang), np.sin(ang)
#     R = np.array(((c, -s), (s, c)))
#     pcs = np.vstack((xn,yn))
#     crot = R @ pcs
    
#     return (crot.T[:,0],crot.T[:,1])

from circle_fit import taubinSVD
def fitcirc(x,y,sigma=None):
    pcs = np.vstack((x,y)).T
    xc, yc, r, sigma = taubinSVD(pcs)
    return (xc, yc, r, sigma)

def centreshift(x,y,xc,yc):
    return (x-xc,y-yc)

def plot_segments(rad,rotang=0,ax=None):
    if ax is None:
        ax = plt.gca()
    for i in range(8):
        ang = i * np.pi / 4.0
        ax.plot([0,rad*np.cos(ang+rotang)],[0,rad*np.sin(ang+rotang)],'b')

def phi_from_coords(xn,yn):
    phis = np.arctan2(yn,xn)
    return phis

def r_from_coords(xn,yn):
    r = np.sqrt(yn*yn+xn*xn)
    return r

# this implementation seems broken; retire for now
def estimate_rotation(xn,yn,spacing=1.0,mode='abs',do_plot=False,secondpass=False):
    n_ang = int(45.0/spacing)
    rotrad = np.arange(n_ang)*piover4/n_ang
    phis = phi_from_coords(xn,yn)
    r = r_from_coords(xn,yn)
    frac, integ = np.modf((np.pi+phis)/piover4)
    if mode == 'abs':
        sqdiff = [ np.sum(np.abs(rd - frac*piover4)) for rd in rotrad]
    elif mode == 'square':
        sqdiff = [ np.sum((rd - frac*piover4)**2) for rd in rotrad]
    else:
        raise RuntimeError("unknown estimation mode %s" % mode)
        
    if do_plot:
        fig,ax = plt.subplots(2,1)
        ax[0].scatter(np.degrees(rotrad),sqdiff)
        ax[1].scatter(np.degrees(phis),r,alpha=0.4)

    indmin = np.argmin(sqdiff)
    radmin = rotrad[indmin]
    radrot = piover4/2.0-radmin

    if secondpass:
        xn2, yn2 = rot_coords(xn,yn,radrot)
        radrot2 = estimate_rotation(xn2,yn2,spacing=spacing,mode=mode,do_plot=do_plot,secondpass=False)
        radrot = radrot+radrot2
    
    return radrot

# FIXED up optimal rotation estimator
# we calculate a metric that looks at angles of events and is designed to have the smallest penalty (zero)
# at the center of a pi/4 (45 deg) segment and increases linearly or squarely towards the edges of the
# segments
# we do this by calculating the fractional part of the angle modulo pi/4 and subtract 0.5 so that the center
# gets a value of 0 and edges go to +- 0.5; we then take abs or squares of this "angle penalty" and sum these up
# we do this by rotating angles by a range of 0..pi/4 and then find the rotation angle that minimises the total penalty
# NOTE: important property of a well working routine: if we rotate the data the determined "optimal rotation"
# should shift linearly with this external rotation
def estimate_rotation2(xn,yn,spacing=1.0,mode='abs',do_plot=False):
    n_ang = int(45.0/spacing)
    rotrad = np.arange(n_ang)*piover4/n_ang
    phis = phi_from_coords(xn,yn)
    r = r_from_coords(xn,yn)

    sqdiff = []
    if mode == 'abs':
        for rd in rotrad:
            frac, integ = np.modf((np.pi+phis+rd)/piover4)
            sqdiff.append(np.sum(np.abs(frac-0.5)))
    elif mode == 'square':
        for rd in rotrad:
            frac, integ = np.modf((np.pi+phis+rd)/piover4)
            sqdiff.append(np.sum((frac-0.5)**2))
    else:
        raise RuntimeError("unknown estimation mode %s" % mode)
        
    if do_plot:
        fig,ax = plt.subplots(2,1)
        ax[0].scatter(np.degrees(rotrad),sqdiff)
        ax[1].scatter(np.degrees(phis),r,alpha=0.4)
        ax[1].set_ylim(0,80)

    indmin = np.argmin(sqdiff)
    radmin = rotrad[indmin]

    return radmin

def rot_coords(xn,yn,ang):
    c, s = np.cos(ang), np.sin(ang)
    R = np.array(((c, -s), (s, c)))
    pcs = np.vstack((xn,yn))
    crot = R @ pcs
    
    return (crot.T[:,0],crot.T[:,1])

def rfilt(xn,yn,r0,dr=25.0):
    
    ri = np.sqrt(xn*xn+yn*yn)
    rimask = (ri >r0-dr)*(ri < r0+dr)
    return (xn[rimask],yn[rimask])

def estimate_nlabeled(x,y,r0=None,nthresh=10,dr=30.0,rotation=None,
                      do_plot=False,secondpass=False,fitmode='abs',return_radius=False,return_bysegments=False):
    if r0 is None:
        xc, yc, r0, sigma = fitcirc(x,y)
        xn, yn = centreshift(x, y, xc, yc)
        if rotation is None:
            radrot = estimate_rotation2(xn,yn,mode=fitmode)
        else:
            radrot = rotation
        xr1,yr1 = rot_coords(xn,yn,radrot)
        xr, yr = rfilt(xr1,yr1,r0,dr=dr)
        if secondpass:
            xc2, yc2, r0, sigma = fitcirc(xr,yr)
            xn, yn = centreshift(xr, yr, xc2, yc2)
            if rotation is None:
                radrot = estimate_rotation2(xn,yn,mode=fitmode)
            else:
                radrot = rotation
            xr,yr = rot_coords(xn,yn,radrot)
    else:
        if rotation is None:
            radrot = estimate_rotation2(x,y,mode=fitmode)
        else:
            radrot = rotation
        xr1,yr1 = rot_coords(x,y,radrot)
        xr, yr = rfilt(xr1,yr1,r0,dr=dr)

    phis = phi_from_coords(xr,yr)
    phibinedges = -np.pi + piover4*np.arange(9)
    nhist,be = np.histogram(phis,bins=phibinedges)
    Nlabeled = np.sum(nhist>=nthresh)
    
    if do_plot:
        segment_radius = 80.0
        fig, axs = plt.subplots(2)
        # first subplot
        axs[0].set_aspect('equal')
        axs[0].scatter(xr,yr,s=10,alpha=0.4,edgecolors='none')
        axs[0].scatter([0],[0],marker='+')
        plot_segments(segment_radius,ax=axs[0])
        cir2 = plt.Circle((0, 0), r0, color='r',fill=False)
        axs[0].add_patch(cir2)
        from matplotlib.patches import Wedge
        phibinedges_deg = np.degrees(phibinedges)

        phibincenters = 0.5*(phibinedges[0:-1]+phibinedges[1:])
        textradius = segment_radius + 13
        for i,phic in enumerate(phibincenters):
            xt,yt = (textradius*np.cos(phic),textradius*np.sin(phic))
            axs[0].text(xt,yt,str(i+1),horizontalalignment='center', # we number segments from 1 to 8
                        verticalalignment='center',alpha=0.4)
        axs[0].set_xlim(-100,100)
        axs[0].set_ylim(-100,100)
        for i in range(nhist.size):
            if nhist[i] >= nthresh:
                axs[0].add_patch(Wedge(
                    (0, 0),                # (x,y)
                    segment_radius,        # radius
                    phibinedges_deg[i],    # theta1 (in degrees)
                    phibinedges_deg[i+1],  # theta2
                    color="r", alpha=0.1))
        axs[0].set_title('NPC Segments = %d, r0 = %.1f nm\nEvent threshold = %d, mode = %s, rot= %.2f rad' % (Nlabeled,r0,nthresh,fitmode,radrot))
        axs[0].invert_yaxis() # the y axis direction seems inverted WRT PYMEVisualise, so try to make equal
                
        # second suplot
        def radtosegno(rad):
            return (rad + np.pi + 0.5*piover4) / piover4

        def segnotorad(sec):
            return -np.pi + 0.5*piover4 + sec * piover4

        axs[1].hist(phis,bins=phibinedges)
        axs[1].plot([phibinedges[0],phibinedges[-1]],[nthresh,nthresh],'r--')
        axs[1].set_xlabel('Angle range $\phi$, $\pi/4$ per segment (radians -$\pi,\cdots,\pi$)')
        axs[1].set_ylabel('Events in segment')
        secax = axs[1].secondary_xaxis('top', functions=(radtosegno, segnotorad))
        secax.set_xlabel('segment number')
        plt.tight_layout()

    if return_radius:
        return (Nlabeled,r0)
    elif return_bysegments:
        return (Nlabeled,nhist)
    else:
        return Nlabeled

from scipy.special import binom
from scipy.optimize import curve_fit

def pn(k,n,p):
    return (binom(n,k)*(np.power(p,k)*np.power((1-p),(n-k))))

def pnpc(k,plabel):
    pbright = 1-pn(0,4,plabel)
    p_k_npc = pn(k,8,pbright)
    
    return p_k_npc

# this is the formula for the 16-spot 3D arrangement, with 2 chances to label per spot
def pnpc3d(k,plabel):
    pbright = 1-pn(0,2,plabel)
    p_k_npc = pn(k,16,pbright)
    
    return p_k_npc

def pnpc3dc(kfit,plabel):
    krange = np.arange(17,dtype='i')
    # important point: we must always first evealuate the probability expressions at the
    # canonical points (0-16), form the cumluative sum and only then
    # interpolate onto the coordinates where the fit is tested in a second step
    pc = np.cumsum(pnpc3d(krange,plabel))
    return np.interp(kfit,krange,pc)

def prangeNPC3D():
    krange = np.arange(17,dtype='i')
    prange=0.1*np.arange(1,10)

    probs = {}
    probs['krange'] = krange
    for p in prange:
        probs[p] = pnpc3d(krange,p)

    return probs


def npclabel_fit(nphist,sigma=None):
    npnormed = nphist/nphist.sum()
    ks = np.arange(9)
    popt, pcov = curve_fit(pnpc, ks, npnormed, sigma=sigma, method='lm', p0=[0.3])
    perr = np.sqrt(np.diag(pcov))
    nlabels_fit = pnpc(ks,popt[0])
    n_labels_scaled = nphist.sum()*nlabels_fit

    return (popt[0],n_labels_scaled,perr[0])

from PYMEcs.misc.utils import get_timestamp_from_filename
def plotcdf_npc3d(nlab,plot_as_points=True,timestamp=None,thresh=None,return_data=False):
    pr = prangeNPC3D()
    for p in pr.keys():
        if p != 'krange':
            plt.plot(pr['krange'],np.cumsum(pr[p]),label="p=%.1f" % p,alpha=0.5)
    if timestamp is None:
        labelexp = 'experiment'
    else:
        labelexp = 'exp %s' % timestamp

    if thresh is not None:
        labelexp = "thresh %d, %s" % (thresh,labelexp)

    # make sure our bin centers are integer spaced from 0 to 16
    histret = plt.hist(nlab,bins=np.arange(17)+0.5,density=True,
                       histtype="step", cumulative=1, label=labelexp, alpha=0.3)
    if plot_as_points:
        histn = histret[0]
        histctr = 0.5*(histret[1][1:]+histret[1][0:-1])
        plt.scatter(histctr,histn)

    # ensure we only have integer major ticks in the plot
    ax = plt.gca()
    ax.xaxis.get_major_locator().set_params(integer=True)
    
    popt,perr, pcbfx, pcbestfit = npclabel_fit3D(histctr,histn)
    plt.plot(pcbfx,pcbestfit,'--')
    plt.legend()
    
    plt.title("NPC 3D analysis using %d NPCs, LE = %.1f %% +- %.1f %%" %
              (nlab.size,100.0*popt,100.0*perr))
    plt.xlabel("N labeled")
    plt.ylabel("CDF")

    if return_data:
        return (histctr,histn)

def npclabel_fit3D(histx,histv,sigma=0.1):
    popt, pcov = curve_fit(pnpc3dc, histx, histv, sigma=sigma, method='lm', p0=[0.4])
    perr = np.sqrt(np.diag(pcov))
    krange = np.arange(17,dtype='i')
    pcumulative = pnpc3dc(krange,popt[0])

    return (popt[0],perr[0], krange, pcumulative)


#################
# NPC 3D Analysis
#################

def to3vecs(x,y,z):
    return np.stack((x,y,z),axis=1)

def xyzfrom3vec(v):
    return (v[:,0],v[:,1],v[:,2])

from scipy.spatial.transform import Rotation as R
from scipy.interpolate import RegularGridInterpolator
from scipy.signal import fftconvolve

def fpinterpolate(fp3d,x,y,z,method='linear', bounds_error=True, fill_value=np.nan):
    # V[i,j,k] = 100*x[i] + 10*y[j] + z[k]
    fpinterp = RegularGridInterpolator((x,y,z), fp3d, method=method, bounds_error=bounds_error, fill_value=fill_value)
    return fpinterp

# variation on makeNPC function in SimuFLUX by Marin & Ries (https://github.com/ries-lab/SimuFLUX)
def makeNPC(center=[0,0,0], R=50, copynumber=32, dz=50, rotation=0, twistangle=np.pi/32, shiftangle=np.pi/16, dR=3, dzpair=3.2):
    if not isinstance(center, np.ndarray):
        center = np.array(center, dtype=float)

    dphi = np.pi/4
    maxphi = 2 * np.pi - dphi - rotation
    nphi = int(maxphi / dphi)
    phi = np.arange(nphi+1)*dphi
    v0 = np.zeros_like(phi)

    if copynumber == 8:
        phiall = phi
        zall = v0
        Rall = v0 + R
    elif copynumber == 16:
        phiall = np.hstack([phi, -phi+twistangle])
        zall = np.hstack([v0+dz/2, v0-dz/2])
        Rall = np.hstack([v0,v0]) + R
    elif copynumber == 32:
        phiall = np.hstack([phi, phi+shiftangle, -phi+twistangle, -phi + twistangle - shiftangle])
        zall = np.hstack([v0+dz/2+dzpair/2, v0+dz/2-dzpair/2, v0-dz/2-dzpair/2, v0-dz/2+dzpair/2])
        Rall = np.hstack([v0+R+dR, v0+R, v0+R+dR, v0+R])
    else:
        raise ValueError("NPC copy number: must be 8, 16 or 32, is: %d" % copynumber)

    posnpc = np.vstack([Rall*np.cos(phiall), Rall*np.sin(phiall), zall]).T
    posnpc += center

    return posnpc

def npctemplate_detailed(x3d,y3d,z3d,npcgeometry,sigma=5.0):
    vol = np.zeros_like(x3d)
    d0, h0 = npcgeometry
    npcpts = makeNPC(R=d0/2, dz=h0)
    x=npcpts[:,0]
    y=npcpts[:,1]
    z=npcpts[:,2]
    for i in range(x.size):
        vol += np.exp(-((x3d-x[i])**2+(y3d-y[i])**2+(z3d-z[i])**2)/(2*sigma*sigma))

    return vol

# we may rewrite this for our purpose if bounds violations become a problem
# code from https://stackoverflow.com/questions/21670080/how-to-find-global-minimum-in-python-optimization-with-bounds
class RandomDisplacementBounds(object):
    """random displacement with bounds"""
    def __init__(self, xmin, xmax, stepsize=0.5):
        self.xmin = xmin
        self.xmax = xmax
        self.stepsize = stepsize

    def __call__(self, x):
        """take a random step but ensure the new position is within the bounds"""
        while True:
            # this could be done in a much more clever way, but it will work for example purposes
            xnew = x + np.random.uniform(-self.stepsize, self.stepsize, np.shape(x))
            if np.all(xnew < self.xmax) and np.all(xnew > self.xmin):
                break
        return xnew

# # define the new step taking routine and pass it to basinhopping
# take_step = RandomDisplacementBounds(xmin, xmax)
# result = basinhopping(f, x0, niter=100, minimizer_kwargs=minimizer_kwargs,
#                       take_step=take_step)

maxshift = 50.0
class LLmaximizerNPC3D(object):
    # the bgprop value needs a little more thought, it could be specific for this set of parameters
    def __init__(self, npcgeometry, extent_nm=150.0, voxelsize_nm=2.0, sigma=5.0, bgprob=1e-9,volcallback=None):
        self.x = np.arange(-extent_nm/2.0, extent_nm/2.0+1.0, voxelsize_nm, dtype='f')
        self.y = self.x.copy()
        self.z = self.x.copy()
        self.npcgeometry = npcgeometry
        self.sigma = sigma
        self.bgprob = bgprob
        self.volcallback = volcallback
        
        x2d,y2d = np.meshgrid(self.x,self.y)
        x3d,y3d, z3d = np.meshgrid(self.x,self.y,self.z)

        if volcallback is not None:
            if not callable(volcallback):
                raise RuntimeError("volcallback option is not callable")
            self.fpg3d = volcallback(x3d,y3d,z3d,npcgeometry,sigma=sigma)
        else:
            eps=15.0 # we have removed the eps variable as input, set it here
            d0, h0 = npcgeometry # diameter of ring and ring spacing
            self.circ2d = (x2d**2 + y2d**2 -0.25*d0**2 <= eps**2) & (x2d**2 + y2d**2 -0.25*d0**2 >= -eps**2)
            self.g3d = np.exp(-(x3d**2+y3d**2+z3d**2)/2.0/(sigma**2))
            self.fp3d = np.zeros_like(x3d)
            idz = np.argmin(np.abs(self.z-h0/2))
            self.fp3d[:,:,idz] = self.circ2d
            idz = np.argmin(np.abs(self.z-(-h0/2)))
            self.fp3d[:,:,idz] = self.circ2d
            self.fpg3d = np.clip(fftconvolve(self.fp3d,self.g3d,mode='same'),0,None)

        self.fpg3d += bgprob # add background probability
        self.fpg3d /= self.fpg3d.sum()
        self.nllfp = -np.log10(self.fpg3d)
        self.nllfpi = None
        self.interpolator()

        self.points = None
        self.pars0 = (0.,0.,0., 0., 0., 100.0, 100.0) # shift_x, shift_y, shift_z, angle_around_z, angle_around_y, scale_xy, scale_z
        self.bounds = (
            (-maxshift,maxshift), # p[0]
            (-maxshift,maxshift), # p[1]
            (-maxshift,maxshift), # p[2]
            (-90.0,90.0), # p[3]
            (-35.0,35.0), # p[4]
            (80.0,120.0), # p[5] - limit to 20% variation to avoid misfits
            (50.0,150.0) # p[6]  - limit to 50% variation to avoid misfits
        )

    def registerPoints(self,pts): # register candidate points for fitting
        self.points = pts
        self.transform_coords(self.pars0)

    def fit(self,method='XXX'): # run the maximumLL fit
        pass

    def fitpars(self): # get the best fit parameters
        pass

    def LLcalc(self,params=None): # return log-likelihood for given parameter set; if None use best fit params
        pass

    def interpolator(self):
        if self.nllfpi is None:
            self.nllfpi = fpinterpolate(self.nllfp, self.x, self.y, self.z,bounds_error=False,fill_value=15.0) # need to think about the fill_value, it could need tweaking

        return self.nllfpi

    def lleval(self,pars):
        c3d = self.transform_coords(pars)
        llvals = self.nllfpi((c3d[:,0],c3d[:,1],c3d[:,2]))
        return llvals.sum()

    def transform_coords(self,pars):
        if self.points is None:
            raise RuntimeError("no valid points, please register points first")
        c3d = self.points + [pars[0],pars[1],pars[2]] # pars[0:3] should be vector offset
        self.c3dr = R.from_euler('zy', [pars[3],pars[4]], degrees=True).apply(c3d)
        self.c3dr[:,0:2] *= 0.01*pars[5]
        self.c3dr[:,2]   *= 0.01*pars[6]
        self._lastpars = pars
        return self.c3dr

    def transform_coords_inv(self,pars): # this is apparently currently not used and needs to be debugged before use
        if 'c3dr' not in dir(self) or self.c3dr is None:
            raise RuntimeError("need transformed points to start with")
        c3dr = self.c3dr.copy()
        c3dr[:,0:2] /= 0.01*pars[5]
        c3dr[:,2]   /= 0.01*pars[6]
        c3di = R.from_euler('zy', [pars[3],pars[4]], degrees=True).inv().apply(c3dr)
        c3di -= [pars[0],pars[1],pars[2]]
        self.c3di = c3di
        return self.c3di

    def plot_points(self,mode='transformed',external_pts=None,axes=None,p0=None): # supported modes should be 'original', 'transformed', 'both', external
        if mode == 'transformed':
            x,y,z = xyzfrom3vec(self.c3dr)
        elif mode == 'original':
            x,y,z = xyzfrom3vec(self.points)
        elif mode == 'both':
            if p0 is not None: # in this mode we consider a p0 as initial transform if provided
                plast = self._lastpars # perhaps better to use opt_result.x?
                self.transform_coords(p0)
                x,y,z = xyzfrom3vec(self.c3dr)
                self.transform_coords(plast) # restore transformed coords to what they were at the start
            else:
                x,y,z = xyzfrom3vec(self.points)
            x1,y1,z1 = xyzfrom3vec(self.c3dr)
        elif mode == 'external':
            if external_pts is None:
                raise RuntimeError("with mode='external' need external points but none supplied")
            x,y,z = xyzfrom3vec(external_pts)
        else:
            raise RuntimeError("unknown mode %s" % mode)
        
        if mode == 'both':
            if axes is None:
                fig, (axt,axb) = plt.subplots(2,3)
            else:
                (axt,axb) = axes
        else:
            if axes is None:
                fig, axt = plt.subplots(1,3,figsize=(6.4,2.4))
            else:
                axt = axes

        axt[0].cla()
        axt[0].imshow(self.fpg3d.sum(axis=2).T,extent=[self.x.min(), self.x.max(), self.y.min(), self.y.max()])
        axt[0].scatter(x,y,c='orange',s=10)
        axt[0].set_aspect('equal')
        axt[0].set_title('x-y')
        axt[1].cla()
        axt[1].imshow(self.fpg3d.sum(axis=1).T,extent=[self.x.min(), self.x.max(), self.z.min(), self.z.max()])
        axt[1].scatter(x,z,c='orange',s=10)
        axt[1].set_aspect('equal')
        axt[1].set_title('x-z')
        axt[2].cla()
        axt[2].imshow(self.fpg3d.sum(axis=0).T,extent=[self.y.min(), self.y.max(), self.z.min(), self.z.max()])
        axt[2].scatter(y,z,c='orange',s=10)
        axt[2].set_aspect('equal')
        axt[2].set_title('y-z')

        if mode == 'both':
            axb[0].cla()
            axb[0].imshow(self.fpg3d.sum(axis=2).T,extent=[self.x.min(), self.x.max(), self.y.min(), self.y.max()])
            axb[0].scatter(x1,y1,c='orange',s=10)
            axb[0].set_aspect('equal')
            axb[0].set_title('x-y')
            axb[1].cla()
            axb[1].imshow(self.fpg3d.sum(axis=1).T,extent=[self.x.min(), self.x.max(), self.z.min(), self.z.max()])
            axb[1].scatter(x1,z1,c='orange',s=10)
            axb[1].set_aspect('equal')
            axb[1].set_title('x-z')
            axb[2].cla()
            axb[2].imshow(self.fpg3d.sum(axis=0).T,extent=[self.y.min(), self.y.max(), self.z.min(), self.z.max()])
            axb[2].scatter(y1,z1,c='orange',s=10)
            axb[2].set_aspect('equal')
            axb[2].set_title('y-z')


    def function_to_minimize(self):
        def minfunc(p):
            return self.lleval(p)

        return minfunc
    
    # minimize the negative log likelihood
    def nllminimize(self,p0=(0,0,0,0,0,100.0,100.0),method='L-BFGS-B'):
        from scipy.optimize import minimize
        self.p0 = p0
        self.minmethod = method
        self.opt_result = minimize(self.function_to_minimize(),p0,method=method,bounds=self.bounds)

    def nll_basin_hopping(self,p0,method='L-BFGS-B',bounds=None):
        from scipy.optimize import basinhopping
        self.p0 = p0 # we record as p0 since pars0 is used at the time of "registerPoints"; will cause issues as nllm object is reused
        self.minmethod = "basinhopping with %s" % method
        if bounds is None:
            bounds=self.bounds # default bounds
        minimizer_kwargs = dict(method=method, bounds=bounds)
        self.opt_result = basinhopping(self.function_to_minimize(), p0, minimizer_kwargs=minimizer_kwargs)

    def pprint_lastpars(self):
        print("Origin: %s" % self._lastpars[0:3])
        print("Angles: %d rot-z, %d rot-y" % tuple(np.round(self._lastpars[3:5])))
        print("Ring diam: %d, ring spacing: %d" % tuple(np.round(np.array(self.npcgeometry)*100.0/np.array(self._lastpars[5:]))))

class NPC3D(object):
    def __init__(self, points=None, pipeline=None, objectID=None, zclip=None, offset_mode='mean'):
        self.points = points
        if pipeline is not None:
            if objectID is None:
                raise RuntimeError("need an objectID to set points from pipeline, None was given")
            npcidx = pipeline['objectID'] == objectID
            self.points = to3vecs(pipeline['x'][npcidx],pipeline['y'][npcidx],pipeline['z'][npcidx])
            self.t = pipeline['t'][npcidx]
            self.objectID = objectID
        self.npts = None
        if self.points is not None:
            self.normalize_points(zclip=zclip, mode=offset_mode)
        self.transformed_pts = None
        self.opt_result = None
        self.filtered_pts = None
        self.fitted = False

    def normalize_points(self,zclip=None,mode='mean'):
        if mode == 'mean':
            self.offset = self.points.mean(axis=0)[None,:]
        elif mode == 'median':
            self.offset = np.median(self.points,axis=0)[None,:]
        else:
            raise RuntimeError("unknown mode '%s', should be mean or median" % mode)
        npts = self.points - self.offset
        nt = self.t
        if not zclip is None:
            zgood = (npts[:,2] > -zclip)*(npts[:,2] < zclip)
            npts = npts[zgood,:]
            nt = nt[zgood]
        self.npts = npts
        self.nt = nt

    def fitbymll(self,nllminimizer,plot=True,printpars=True,axes=None,preminimizer=None,axespre=None):
        nllm = nllminimizer
        self.nllminimizer = nllm
        self.preminimizer = preminimizer

        if preminimizer is not None:
            preminimizer.registerPoints(self.npts)
            preminimizer.nll_basin_hopping(p0=(0,0,0,0,0,100.0,100.0))
            self.opt_result_pre = preminimizer.opt_result
            self.transformed_pts_pre = preminimizer.c3dr
            self.bounds_pre = preminimizer.bounds
            # in the second stage llm minimizing stage start with best fit from previous fit as p0
            # and allow mainly variation in rotation angles 
            # for other parameters only allow deviation from robust fitting in quite narrow range
            p0 = self.opt_result_pre.x
            dc = 3.0 # max deviation in coordinates (in nm)
            dperc = 5.0 # max deviation in scaling percentage
            bounds = (
                (p0[0]-dc,p0[0]+dc), # p[0]
                (p0[1]-dc,p0[1]+dc), # p[1]
                (p0[2]-dc,p0[2]+dc), # p[2]
                (-90.0,90.0), # p[3]
                (-35.0,35.0), # p[4]
                (p0[5]-dperc,p0[5]+dperc), # p[5] - limit to 20% variation to avoid misfits
                (p0[6]-dperc,p0[6]+dperc) # p[6]  - limit to 20% variation to avoid misfits
            )
            logger.info("p0 %s" % p0)
            logger.info("bounds %s" % repr(bounds))

            nllm.registerPoints(self.npts)
            nllm.nll_basin_hopping(p0=p0,bounds=bounds)
            self.opt_result = nllm.opt_result
            self.transformed_pts = nllm.c3dr
            self.bounds = bounds
        else:
            nllm.registerPoints(self.npts)
            nllm.nll_basin_hopping(p0=(0,0,0,0,0,100.0,100.0)) # no bound keyword implies default bounds
            self.opt_result = nllm.opt_result
            self.transformed_pts = nllm.c3dr
            self.bounds = nllm.bounds
        self.fitted = True
        if printpars:
            nllm.pprint_lastpars()
        if plot:
            if preminimizer is not None:
                preminimizer.plot_points(mode='both',axes=axespre)
                p0 = nllm.p0
            else:
                p0 = None
            nllm.plot_points(mode='both',axes=axes,p0=p0) # if a prefit was done we use its p0

    def filter(self,axis='z',minval=0, maxval=100):
        if axis == 'x':
            coords = self.transformed_pts[:,0]
        elif axis == 'y':
            coords = self.transformed_pts[:,1]
        elif axis == 'z':
            coords = self.transformed_pts[:,2]
        else:
            raise RuntimeError("unknow axis %s requested (must be x, y or z)" % axis)

        goodidx = (coords >= minval)*(coords <= maxval)
        self.filtered_pts = self.transformed_pts[goodidx,:]
        try:
            self.filtered_t = self.nt[goodidx]
        except AttributeError: # ignore if we do not have the 'nt' attribute
            pass
        
    def plot_points(self,mode='transformed'):
        if mode == 'normalized':
            pts = self.npts
        elif mode == 'transformed':
            pts = self.transformed_pts
        elif mode == 'filtered':
            pts = self.filtered_pts
        else:
            raise RuntimeError("unknown mode %s" % mode)

        self.nllminimizer.plot_points(mode='external',external_pts=pts)
        
            
    def plot_points3D(self,mode='transformed',ax=None,with_offset=False,s=10):
        if mode == 'normalized':
            pts = self.npts
            if with_offset:
                pts = pts + self.offset
        elif mode == 'transformed':
            pts = self.transformed_pts
        elif mode == 'filtered':
            pts = self.filtered_pts
        else:
            raise RuntimeError("unknown mode %s" % mode)

        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(projection='3d')
        ax.scatter(pts[:,0], pts[:,1], pts[:,2], 'o', s=s)
        return ax

    # the nominal glyph diam and height, set by the LLM fitter npcgeometry
    # note access will only work after llm fit has taklen place!
    def get_glyph_diam(self):
        return self.nllminimizer.npcgeometry[0]

    def get_glyph_height(self):
        return self.nllminimizer.npcgeometry[1]
    
    def get_glyph(self, with_offset=True):

        def get_circ_coords(radius=50.0, npoints=25):
            angle = np.linspace( 0 , 2 * np.pi , npoints) 
            xc = radius * np.cos( angle ) 
            yc = radius * np.sin( angle )
            return (xc,yc)

        def transform_coords_invs(x,y,z,pars,offset):
            xr = x.copy()
            yr = y.copy()
            zr = z.copy()
            zr /= 0.01*pars[6]
            xr /= 0.01*pars[5]
            yr /= 0.01*pars[5]
            c3d = np.stack([xr,yr,zr],axis=1)
            c3di = R.from_euler('zy', [pars[3],pars[4]], degrees=True).inv().apply(c3d)
            c3di -= [pars[0],pars[1],pars[2]]
            c3di += offset

            return (c3di[:,0],c3di[:,1],c3di[:,2])
        
        xc, yc = get_circ_coords(radius=0.5*self.get_glyph_diam(), npoints=25)

        if with_offset:
            offset = self.offset
        else:
            offset = 0
        pars = self.opt_result.x
        glyph = {}
        x1,y1,z1 = transform_coords_invs(xc,yc,np.zeros_like(xc)-0.5*self.get_glyph_height(),pars,offset)
        glyph['circ_bot'] = to3vecs(x1,y1,z1)
        x2,y2,z2 = transform_coords_invs(xc,yc,np.zeros_like(xc)+0.5*self.get_glyph_height(),pars,offset)
        glyph['circ_top'] = to3vecs(x2,y2,z2)
        xa,ya,za = transform_coords_invs([0,0],[0,0],[-75.0,+75.0],pars,offset)
        glyph['axis'] = to3vecs(xa,ya,za)

        return glyph
        
    def plot_points3D_with_glyph(self, ax=None, with_offset=False, s=10):

        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(projection='3d')

        self.plot_points3D(mode='normalized',ax=ax,with_offset=with_offset,s=s)
        glyph = self.get_glyph(with_offset=with_offset)
        x1,y1,z1 = xyzfrom3vec(glyph['circ_bot'])
        ax.plot(x1,y1,z1,'r')
        x2,y2,z2 = xyzfrom3vec(glyph['circ_top'])
        ax.plot(x2,y2,z2,'r')
        xa,ya,za = xyzfrom3vec(glyph['axis'])
        ax.plot(xa,ya,za,'r')
        
        return ax
        

    def nlabeled(self,nthresh=1,r0=50.0,dr=25.0,do_plot=False,rotlocked=True,zrange=150.0,analysis2d=False,rotation=None):
        zrangeabs = abs(zrange)
        if rotlocked and rotation is None:
            self.filter('z',-zrangeabs,zrangeabs)
            if self.filtered_pts.size > 0:
                rotation = estimate_rotation2(self.filtered_pts[:,0],self.filtered_pts[:,1])
            else:
                rotation=None
        self.rotation = rotation # remember rotation

        if analysis2d:
            self.filter('z',-zrangeabs,zrangeabs)
            x=self.filtered_pts[:,0]
            y=self.filtered_pts[:,1]
            if self.filtered_pts.size > 0:
                self.n = estimate_nlabeled(x,y,r0=r0,dr=dr,nthresh=nthresh,do_plot=do_plot,rotation=rotation)
            else:
                self.n = 0
            return self.n
        else:
            self.filter('z',0,zrangeabs)
            if self.filtered_pts.size > 0:
                # self.plot_points('filtered')
                x=self.filtered_pts[:,0]
                y=self.filtered_pts[:,1]
                (self.n_top,self.n_top_bysegments) = estimate_nlabeled(x,y,r0=r0,dr=dr,nthresh=nthresh,do_plot=do_plot,
                                                                       rotation=rotation,return_bysegments=True)
            else:
                (self.n_top,self.n_top_bysegments) = (0,np.zeros((8)))
            self.filter('z',-zrangeabs,0)
            if self.filtered_pts.size > 0:
                # self.plot_points('filtered')
                x=self.filtered_pts[:,0]
                y=self.filtered_pts[:,1]
                (self.n_bot,self.n_bot_bysegments) = estimate_nlabeled(x,y,r0=r0,dr=dr,nthresh=nthresh,do_plot=do_plot,
                                                                       rotation=rotation,return_bysegments=True)
            else:
                (self.n_bot,self.n_bot_bysegments) = (0,np.zeros((8)))
            return (self.n_top,self.n_bot)

class NPC3DSet(object):
    def __init__(self,filename=None,zclip=75.0,offset_mode='median',NPCdiam=100.0,NPCheight=70.0,
                 foreshortening=1.0,known_number=-1,templatemode='standard',sigma=7.0):
        self.filename=filename
        self.zclip = zclip
        self.offset_mode = offset_mode
        self.npcdiam = NPCdiam
        self.npcheight = NPCheight
        self.npcs = []
        self.foreshortening=foreshortening # we just record this for reference
        # TODO: expose llm parameters to this init method as needed in practice!
        self.templatemode = templatemode
        if templatemode == 'standard':
            volcallback=None
            # sigma = 7.0 # we set this via sigma keyword now
        elif templatemode == 'detailed' or templatemode == 'twostage':
            volcallback=npctemplate_detailed
            # sigma = 7.0 # we set this via sigma keyword now
        self.llm = LLmaximizerNPC3D([self.npcdiam,self.npcheight],sigma=sigma,bgprob=1e-9,extent_nm=300.0,volcallback=volcallback)
        if templatemode == 'twostage':
            self.llmpre = LLmaximizerNPC3D([self.npcdiam,self.npcheight],sigma=sigma,bgprob=1e-9,extent_nm=300.0,volcallback=None)
        else:
            self.llmpre = None
        self.measurements = []
        self.known_number = known_number # only considered if > 0

        # v1.1
        # - has templatemode added
        # - has more complete recording of llm initialization parameters
        # v1.2
        # - twostage mode introduced with prefitting with robust template followed by detailed template
        # v1.3
        # - add time from original points
        # v1.4
        # - add bounds info to npc object when fitting
        self._version='1.4' # REMEMBER to increment version when changing this object or the underlying npc object definitions

    def registerNPC(self,npc):
        self.npcs.append(npc)

    def addNPCfromPipeline(self,pipeline,oid):
        self.registerNPC(NPC3D(pipeline=pipeline,objectID=oid,zclip=self.zclip,offset_mode=self.offset_mode))
        
    def measure_labeleff(self,nthresh=1,do_plot=False,printpars=False,refit=False):
        self.measurements = []
        if do_plot:
            fig, axes = plt.subplots(2,3)
        else:
            axes = None
        for npc in self.npcs:
            if not npc.fitted or refit:
                npc.fitbymll(self.llm,plot=do_plot,axes=axes,printpars=printpars)
            nt,nb = npc.nlabeled(nthresh=nthresh,dr=20.0)
            self.measurements.append([nt,nb])
    
    def plot_labeleff(self,thresh=None):
        from PYMEcs.misc.utils import get_timestamp_from_filename
        if len(self.measurements) < 10:
            raise RuntimeError("not enough measurements, need at least 10, got %d" %
                               len(self.measurements))
        meas = np.array(self.measurements)
        nlab = meas.sum(axis=1)
        # fill with trailing zeros if we have a known number of NPCs but have fewer measurements
        # the "missing NPCs" typically represent NPCs with no events
        if int(self.known_number) > 0 and nlab.shape[0] < int(self.known_number):
            nlab = np.pad(nlab,((0,int(self.known_number)-nlab.shape[0])))

        plt.figure()
        plotcdf_npc3d(nlab,timestamp=get_timestamp_from_filename(self.filename),thresh=thresh)

    def diam(self):
        diams = []
        for npc in self.npcs:
            diams.append(npc.get_glyph_diam()/(0.01*npc.opt_result.x[5]))
        return diams

    def height(self):
        heights = []
        for npc in self.npcs:
            heights.append(npc.get_glyph_height()/(0.01*npc.opt_result.x[6]))
        return heights

    def n_bysegments(self):
        nbs_top = []
        nbs_bot = []
        for npc in self.npcs:
            if npc.fitted:
                try:
                   nbs_top.append(npc.n_top_bysegments)
                except AttributeError:
                    pass
                try:
                   nbs_bot.append(npc.n_bot_bysegments)
                except AttributeError:
                    pass
        if len(nbs_top) == 0 and len(nbs_bot) == 0:
            return None
        else:
            return dict(top=np.array(nbs_top),bottom=np.array(nbs_bot))

    def version(self):
        return self._version
    
class NPCSetContainer(object):
    def __init__(self,npcs):
        self.npcs=npcs

    def get_npcset(self):
        if '_unpickled' in dir(self):
            warn(self._unpickled)
            return None
        return self.npcs
    
    # we add custom pickling/unpickling methods so that an npcset instance in the
    # PYME metadata won't greatly inflate the image saving footprint
    # the pickled version of the npcset should be on disk already, no need to "properly" pickle/unpickle the NPCSetContainer
    #  which was created for this very purpose, i.e. the NPCSetContainer object needs to only persist during the lifetime of the pymevis viewer
    def __getstate__(self):
        warn("NPCset is being pickled - just a dummy mostly for PYME metadata - won't be usable after unpickling")
        return 'not a valid npcset object after pickling/unpickling'
    
    def __setstate__(self, d):
        warn("NPCset is being unpickled - this is just a dummy unpickle, won't be usable after unpickling")
        self._unpickled = d

def mk_NPC_gallery(npcs,mode,zclip3d,NPCRotationAngle,xoffs=0,yoffs=0,enforce_8foldsym=False):
    x = np.empty((0))
    y = np.empty((0))
    z = np.empty((0))
    t = np.empty((0),int)
    objectID = np.empty((0),int)
    is_top = np.empty((0),int)
    segmentID = np.empty((0),int)
    phi = np.empty((0))

    if mode == 'TopOverBottom':
        gspx = 180
        gspy = 180
        gsdx = 0
        rowlength = 10
    elif mode == 'TopBesideBottom':
        gspx = 400
        gspy = 180
        gsdx = 180
        rowlength = 5

    elif mode == 'SingleAverageSBS':
        gspx = 0
        gspy = 0
        gsdx = 180
        rowlength = 5

    else:
           gspx = 0
           gspy = 0
           gsdx = 0
           rowlength = 5

    xtr = np.empty((0))
    ytr = np.empty((0))
    ztr = np.empty((0))
    objectIDtr = np.empty((0),int)
    polyidtr = np.empty((0),int)
        
    xg = []
    yg = []
    dphi = np.pi/4
    radius = 65.0
    pi_base = []
        
    for i in range(8):
        xg.extend([0,radius*np.sin(i*dphi)])
        yg.extend([0,radius*np.cos(i*dphi)])
        pi_base.extend([i+1,i+1])

    zt = np.full((16),25.0)
    zb = -np.full((16),25.0)

    polyidx = np.array(pi_base,dtype='i')
    xga = np.array(xg)
    yga = np.array(yg)

    def filtered_t(npc):
        if 'filtered_t' in dir(npc):
            return npc.filtered_t
        else:
            return range(npc.filtered_pts.shape[0])

    if enforce_8foldsym:
        angled_repeats = 8
    else:
        angled_repeats = 1

    for i,npc in enumerate(npcs.npcs):
        if not npc.fitted:
            warn("NPC not yet fitted, please call only after fitting")
            return

        gxt = (i % rowlength) * gspx
        gxb = gxt + gsdx
        gy = (i // rowlength) * gspy
            
        npc.filter('z',0,zclip3d)
        ptst = npc.filtered_pts
        tt = filtered_t(npc)
        
        npc.filter('z',-zclip3d,0)
        ptsb = npc.filtered_pts
        tb = filtered_t(npc)

        from scipy.spatial.transform import Rotation as R
        for i in range(angled_repeats):
            if npc.rotation is not None:
                if NPCRotationAngle == 'negative':
                    factor = -1.0
                elif NPCRotationAngle == 'positive':
                    factor = 1.0
                else:
                    factor = 0.0
                ptst_t = R.from_euler('z', factor*npc.rotation + i*piover4, degrees=False).apply(ptst)
                ptsb_t = R.from_euler('z', factor*npc.rotation + i*piover4, degrees=False).apply(ptsb)

            phit = phi_from_coords(ptst_t[:,0],ptst_t[:,1])
            phib = phi_from_coords(ptsb_t[:,0],ptsb_t[:,1])
        
            x = np.append(x,ptst_t[:,0] + gxt)
            y = np.append(y,ptst_t[:,1] + gy)
            z = np.append(z,ptst_t[:,2])
            phi = np.append(phi, phit)
            segmentID = np.append(segmentID, ((phit+np.pi)/piover4).astype(int))
            t = np.append(t,tt)

            x = np.append(x,ptsb_t[:,0] + gxb)
            y = np.append(y,ptsb_t[:,1] + gy)
            z = np.append(z,ptsb_t[:,2])
            phi = np.append(phi, phib)
            segmentID = np.append(segmentID, ((phib+np.pi)/piover4).astype(int))
            t = np.append(t,tb)
        
            objectID = np.append(objectID,np.full_like(ptst[:,0],npc.objectID,dtype=int))
            objectID = np.append(objectID,np.full_like(ptsb[:,0],npc.objectID,dtype=int))

            is_top = np.append(is_top,np.ones_like(phit,dtype=int))
            is_top = np.append(is_top,np.zeros_like(phib,dtype=int))

        # remaining stuff for trace dict which shows segment boundaries
        xtr = np.append(xtr,xga + gxt)
        ytr = np.append(ytr,yga + gy)
        ztr = np.append(ztr,zt)

        objectIDtr = np.append(objectIDtr, np.full_like(xga,npc.objectID,dtype=int))
        polyidtr = np.append(polyidtr, polyidx)
        polyidx += 8

        xtr = np.append(xtr,xga + gxb)
        ytr = np.append(ytr,yga + gy)
        ztr = np.append(ztr,zb)

        objectIDtr = np.append(objectIDtr, np.full_like(xga,npc.objectID,dtype=int))
        polyidtr = np.append(polyidtr, polyidx)
        polyidx += 8            

    # t = np.arange(x.size)
    A = np.full_like(x,10.0,dtype='f')
    error_x = np.full_like(x,1.0,dtype='f')
    error_y = np.full_like(x,1.0,dtype='f')
    error_z = np.full_like(x,1.0,dtype='f')

    dsdict = dict(x=x+xoffs,y=y+yoffs,z=z,
                  objectID=objectID,t=t,A=A,
                  error_x=error_x,error_y=error_y,error_z=error_z,
                  is_top=is_top,segmentID=segmentID,phi=phi)

    trdict = dict(x=xtr+xoffs,y=ytr+yoffs,z=ztr,
                  objectID=objectIDtr,polyIndex=polyidtr)
        
    from PYME.IO.tabular import DictSource
    gallery = DictSource(dsdict)
    segments = DictSource(trdict)

    return gallery,segments

def mk_npctemplates(npcs):
    x = np.empty((0))
    y = np.empty((0))
    z = np.empty((0))
    polyIndex = np.empty((0),int)
    polySize = np.empty((0),int)
    objectID = np.empty((0),int)
    NtopLabelled = np.empty((0),int)
    NbotLabelled = np.empty((0),int)
    NLabelled = np.empty((0),int)
    diams = np.empty((0),float)
    heights = np.empty((0),float)
    fitquals = np.empty((0),float)
    cx = np.empty((0),float)
    cy = np.empty((0),float)
    cz = np.empty((0),float)
    ci = 1
    for npc in npcs.npcs:
        nt, nb = (npc.n_top,npc.n_bot)
        glyph = npc.get_glyph()
        pars = npc.opt_result.x
        diam = npc.get_glyph_diam() / (0.01*pars[5])
        height = npc.get_glyph_height() / (0.01*pars[6])
        fitqual = npc.opt_result.fun/npc.npts.shape[0]
        offset = npc.offset[0]
        for poly in ['circ_bot','circ_top','axis']:
            c3 = glyph[poly]
            xg = c3[:,0]
            yg = c3[:,1]
            zg = c3[:,2]
            x = np.append(x,xg)
            y = np.append(y,yg)
            z = np.append(z,zg)
            polyIndex = np.append(polyIndex,np.full_like(xg,ci,dtype=int))
            polySize = np.append(polySize,np.full_like(xg,xg.size,dtype=int))
            ci += 1
            objectID = np.append(objectID,np.full_like(xg,npc.objectID,dtype=int))
            NtopLabelled = np.append(NtopLabelled,np.full_like(xg,nt,dtype=int))
            NbotLabelled = np.append(NbotLabelled,np.full_like(xg,nb,dtype=int))
            NLabelled = np.append(NLabelled,np.full_like(xg,nt+nb,dtype=int))               
            diams = np.append(diams,np.full_like(xg,diam,dtype=float))
            heights = np.append(heights,np.full_like(xg,height,dtype=float))
            fitquals = np.append(fitquals,np.full_like(xg,fitqual,dtype=float))
            cx = np.append(cx,np.full_like(xg,offset[0]-pars[0],dtype=float))
            cy = np.append(cy,np.full_like(xg,offset[1]-pars[1],dtype=float))
            cz = np.append(cz,np.full_like(xg,offset[2]-pars[2],dtype=float))
    t = np.arange(x.size)
    A = np.full_like(x,10.0,dtype='f')
    error_x = np.full_like(x,1.0,dtype='f')
    error_y = np.full_like(x,1.0,dtype='f')
    error_z = np.full_like(x,1.0,dtype='f')
        
    dsdict = dict(x=x,y=y,z=z,polyIndex=polyIndex,polySize=polySize,
                  NtopLabelled=NtopLabelled,NbotLabelled=NbotLabelled,NLabelled=NLabelled,
                  objectID=objectID,t=t,A=A,
                  error_x=error_x,error_y=error_y,error_z=error_z,
                  npc_height=heights,npc_diam=diams,npc_fitqual=fitquals,
                  npc_ctrx=cx,npc_ctry=cy,npc_ctrz=cz)

    from PYME.IO.tabular import DictSource
    return DictSource(dsdict)
