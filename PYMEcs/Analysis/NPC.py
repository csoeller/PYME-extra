import matplotlib.pyplot as plt
import numpy as np

piover4 = np.pi/4.0

from circle_fit import taubinSVD
def fitcirc(x,y,sigma=None):
    pcs = np.vstack((x,y)).T
    xc, yc, r, sigma = taubinSVD(pcs)
    return (xc, yc, r, sigma)

def centreshift(x,y,xc,yc):
    return (x-xc,y-yc)

def plot_segments(rad,rotang=0):
    for i in range(8):
        ang = i * np.pi / 4.0
        plt.plot([0,rad*np.cos(ang+rotang)],[0,rad*np.sin(ang+rotang)],'b')

def phi_from_coords(xn,yn):
    phis = np.arctan2(yn,xn)
    return phis

def estimate_rotation(xn,yn,spacing=1.0,mode='abs',do_plot=False):
    n_ang = int(45.0/spacing)
    rotrad = np.arange(n_ang)*piover4/n_ang
    phis = phi_from_coords(xn,yn)
    frac, integ = np.modf((np.pi+phis)/piover4)
    if mode == 'abs':
        sqdiff = [ np.sum(np.abs(rd - frac*piover4)) for rd in rotrad]
    elif mode == 'square':
        sqdiff = [ np.sum((rd - frac*piover4)**2) for rd in rotrad]
    else:
        raise RuntimeError("unknown estimation mode %s" % mode)
        
    if do_plot:
        plt.figure()
        plt.scatter(np.degrees(rotrad),sqdiff)

    indmin = np.argmin(sqdiff)
    radmin = rotrad[indmin]
    
    return piover4/2.0-radmin

def rot_coords(xn,yn,ang):
    c, s = np.cos(ang), np.sin(ang)
    R = np.array(((c, -s), (s, c)))
    pcs = np.vstack((xn,yn))
    crot = R @ pcs
    
    return (crot.T[:,0],crot.T[:,1])

piover4 = np.pi/4.0

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

def estimate_rotation(xn,yn,spacing=1.0,mode='abs',do_plot=False):
    n_ang = int(45.0/spacing)
    rotrad = np.arange(n_ang)*piover4/n_ang
    phis = phi_from_coords(xn,yn)
    frac, integ = np.modf((np.pi+phis)/piover4)
    if mode == 'abs':
        sqdiff = [ np.sum(np.abs(rd - frac*piover4)) for rd in rotrad]
    elif mode == 'square':
        sqdiff = [ np.sum((rd - frac*piover4)**2) for rd in rotrad]
    else:
        raise RuntimeError("unknown estimation mode %s" % mode)
        
    if do_plot:
        plt.figure()
        plt.scatter(np.degrees(rotrad),sqdiff)

    indmin = np.argmin(sqdiff)
    radmin = rotrad[indmin]
    
    return piover4/2.0-radmin

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

def estimate_nlabeled(x,y,r0=None,nthresh=10,dr=30.0,
                      do_plot=False,secondpass=False,fitmode='abs',return_radius=False):
    if r0 is None:
        xc, yc, r0, sigma = fitcirc(x,y)
        xn, yn = centreshift(x, y, xc, yc)
        radrot = estimate_rotation(xn,yn,mode=fitmode)
        xr1,yr1 = rot_coords(xn,yn,radrot)
        xr, yr = rfilt(xr1,yr1,r0,dr=dr)
        if secondpass:
            xc2, yc2, r0, sigma = fitcirc(xr,yr)
            xn, yn = centreshift(xr, yr, xc2, yc2)
            radrot = estimate_rotation(xn,yn,mode=fitmode)
            xr,yr = rot_coords(xn,yn,radrot)
    else:
        radrot = estimate_rotation(x,y,mode=fitmode)
        xr1,yr1 = rot_coords(x,y,radrot)
        xr, yr = rfilt(xr1,yr1,r0,dr=dr)

    phis = phi_from_coords(xr,yr)
    phibinedges = -np.pi + piover4*np.arange(9)
    nhist,be = np.histogram(phis,bins=phibinedges)
    Nlabeled = np.sum(nhist>nthresh)
    
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
            if nhist[i] > nthresh:
                axs[0].add_patch(Wedge(
                    (0, 0),                # (x,y)
                    segment_radius,        # radius
                    phibinedges_deg[i],    # theta1 (in degrees)
                    phibinedges_deg[i+1],  # theta2
                    color="r", alpha=0.1))
        axs[0].set_title('NPC Segments = %d, r0 = %.1f nm\nEvent threshold = %d, mode = %s' % (Nlabeled,r0,nthresh,fitmode))
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
def plotcdf_npc3d(nlab,plot_as_points=True,timestamp=None):
    pr = prangeNPC3D()
    for p in pr.keys():
        if p != 'krange':
            plt.plot(pr['krange'],np.cumsum(pr[p]),label="p=%.1f" % p)
    if timestamp is None:
        labelexp = 'experiment'
    else:
        labelexp = 'exp %s' % timestamp

    histret = plt.hist(nlab,range=(0,16),bins='auto',density=True,
                       histtype="step", cumulative=1, label=labelexp, alpha=0.3)
    if plot_as_points:
        histn = histret[0]
        histctr = 0.5*(histret[1][1:]+histret[1][0:-1])
        plt.scatter(histctr,histn)

    popt,perr, pcbfx, pcbestfit = npclabel_fit3D(histctr,histn)
    plt.plot(pcbfx,pcbestfit,'--')
    plt.legend()
    
    plt.title("NPC 3D analysis using %d NPCs, LE = %d %% +- %.1f %%" %
              (nlab.size,np.round(100.0*popt),100.0*perr))
    plt.xlabel("N labeled")
    plt.ylabel("CDF")

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
    def __init__(self, p0, extent_nm=150.0, voxelsize_nm=2.0, eps=15.0, sigma=5.0, bgprob=1e-9):
        self.x = np.arange(-extent_nm/2.0, extent_nm/2.0+1.0, voxelsize_nm, dtype='f')
        self.y = self.x.copy()
        self.z = self.x.copy()
        self.p0 = p0
        d0, h0 = p0 # diameter of ring and ring spacing
        x2d,y2d = np.meshgrid(self.x,self.y)
        x3d,y3d, z3d = np.meshgrid(self.x,self.y,self.z)
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
            (85.0,115.0), # p[5]
            (85.0,115.0) # p[6]
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

    def transform_coords_inv(self,pars):
        if 'c3dr' not in dir(self) or self.c3dr is None:
            raise RuntimeError("need transformed points to start with")
        c3dr = self.c3dr.copy()
        c3dr[:,0:2] /= 0.01*pars[5]
        c3dr[:,2]   /= 0.01*pars[6]
        c3di = R.from_euler('zy', [pars[3],pars[4]], degrees=True).inv.apply(c3d)
        c3di -= [pars[0],pars[1],pars[2]]
        self.c3di = c3di
        return self.c3di

    def plot_points(self,mode='transformed',external_pts=None,axes=None): # supported modes should be 'original', 'transformed', 'both', external
        if mode == 'transformed':
            x,y,z = xyzfrom3vec(self.c3dr)
        elif mode == 'original':
            x,y,z = xyzfrom3vec(self.points)
        elif mode == 'both':
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
        self.opt_result = minimize(self.function_to_minimize(),p0,method=method,bounds=self.bounds)

    def nll_basin_hopping(self,p0,method='L-BFGS-B'):
        from scipy.optimize import basinhopping
        minimizer_kwargs = dict(method=method, bounds=self.bounds)
        self.opt_result = basinhopping(self.function_to_minimize(), p0, minimizer_kwargs=minimizer_kwargs)

    def pprint_lastpars(self):
        print("Origin: %s" % self._lastpars[0:3])
        print("Angles: %d rot-z, %d rot-y" % tuple(np.round(self._lastpars[3:5])))
        print("Ring diam: %d, ring spacing: %d" % tuple(np.round(np.array(self.p0)*100.0/np.array(self._lastpars[5:]))))

class NPC3D(object):
    def __init__(self, points=None, pipeline=None, objectID=None, zclip=None):
        self.points = points
        if pipeline is not None:
            if objectID is None:
                raise RuntimeError("need an objectID to set points from pipeline, None was given")
            npcidx = pipeline['objectID'] == objectID
            self.points = to3vecs(pipeline['x'][npcidx],pipeline['y'][npcidx],pipeline['z'][npcidx])
            self.objectID = objectID
        self.npts = None
        if self.points is not None:
            self.normalize_points(zclip=zclip)
        self.transformed_pts = None
        self.opt_result = None
        self.filtered_pts = None
        self.fitted = False

    def normalize_points(self,zclip=None):
        npts = self.points - self.points.mean(axis=0)[None,:]
        if not zclip is None:
            zgood = (npts[:,2] > -zclip)*(npts[:,2] < zclip)
            npts = npts[zgood,:]
        self.npts = npts

    def fitbymll(self,nllminimizer,plot=True,printpars=True,axes=None):
        nllm = nllminimizer
        self.nllminimizer = nllm
        
        nllm.registerPoints(self.npts)
        nllm.nll_basin_hopping(p0=(0,0,0,0,0,100.0,100.0))
        self.opt_result = nllm.opt_result
        self.transformed_pts = nllm.c3dr
        self.fitted = True
        if printpars:
            nllm.pprint_lastpars()
        if plot:
            nllm.plot_points(mode='both',axes=axes)

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
        
            
    def plot_points3D(self,mode='transformed'):
        if mode == 'normalized':
            pts = self.npts
        elif mode == 'transformed':
            pts = self.transformed_pts
        elif mode == 'filtered':
            pts = self.filtered_pts
        else:
            raise RuntimeError("unknown mode %s" % mode)

        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        ax.scatter(pts[:,0], pts[:,1], pts[:,2], 'o')
        return ax
        
        
    def nlabeled(self,nthresh=1,r0=50.0,dr=25.0,do_plot=False):
        self.filter('z',0,150)
        if self.filtered_pts.size > 0:
            # self.plot_points('filtered')
            x=self.filtered_pts[:,0]
            y=self.filtered_pts[:,1]
            self.n_top = estimate_nlabeled(x,y,r0=r0,dr=dr,nthresh=nthresh-1,do_plot=do_plot)
        else:
            self.n_top = 0
        self.filter('z',-150,0)
        if self.filtered_pts.size > 0:
            # self.plot_points('filtered')
            x=self.filtered_pts[:,0]
            y=self.filtered_pts[:,1]
            self.n_bot = estimate_nlabeled(x,y,r0=r0,dr=dr,nthresh=nthresh-1,do_plot=do_plot)
        else:
            self.n_bot = 0
        return (self.n_top,self.n_bot)

class NPC3DSet(object):
    def __init__(self,filename=None):
        self.filename=filename
        self.npcs = []
        self.llm = LLmaximizerNPC3D([100.0,70.0],eps=15.0,sigma=7.0,bgprob=1e-9,extent_nm=300.0)
        self.measurements = []

    def registerNPC(self,npc):
        self.npcs.append(npc)

    def addNPCfromPipeline(self,pipeline,oid):
        self.registerNPC(NPC3D(pipeline=pipeline,objectID=oid,zclip=75))
        
    def measure_labeleff(self,nthresh=1,do_plot=False,printpars=False,refit=False):
        self.measurements = []
        if do_plot:
            fig, axes = plt.subplots(2,3)
        else:
            axes = None
        for npc in self.npcs:
            if not npc.fitted and not refit:
                npc.fitbymll(self.llm,plot=do_plot,axes=axes,printpars=printpars)
            nt,nb = npc.nlabeled(nthresh=nthresh,dr=20.0)
            self.measurements.append([nt,nb])
    
    def plot_labeleff(self):
        from PYMEcs.misc.utils import get_timestamp_from_filename
        if len(self.measurements) < 10:
            raise RuntimeError("not enough measurements, need at least 10, got %d" %
                               len(self.measurements))
        meas = np.array(self.measurements)
        nlab = meas.sum(axis=1)

        plt.figure()
        plotcdf_npc3d(nlab,timestamp=get_timestamp_from_filename(self.filename))

