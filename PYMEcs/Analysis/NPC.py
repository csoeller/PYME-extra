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

def estimate_nlabeled(x,y,nthresh=10,do_plot=False,secondpass=False,fitmode='abs'):
    xc, yc, r0, sigma = fitcirc(x,y)
    xn, yn = centreshift(x, y, xc, yc)
    radrot = estimate_rotation(xn,yn,mode=fitmode)
    xr1,yr1 = rot_coords(xn,yn,radrot)
    xr, yr = rfilt(xr1,yr1,r0,dr=30.0)
    if secondpass:
        xc2, yc2, r0, sigma = fitcirc(xr,yr)
        xn, yn = centreshift(xr, yr, xc2, yc2)
        radrot = estimate_rotation(xn,yn,mode=fitmode)
        xr,yr = rot_coords(xn,yn,radrot)
    
    phis = phi_from_coords(xr,yr)
    phibinedges = -np.pi + piover4*np.arange(9)
    nhist,be = np.histogram(phis,bins=phibinedges)
    Nlabeled = np.sum(nhist>nthresh)
    
    if do_plot:
        fig, axs = plt.subplots(2)
        axs[0].set_aspect('equal')
        axs[0].scatter(xr,yr,s=5)
        axs[0].scatter([0],[0],marker='+')
        plot_segments(80.0,ax=axs[0])
        cir2 = plt.Circle((0, 0), r0, color='r',fill=False)
        axs[0].add_patch(cir2)
        axs[0].set_title('N = %d, r0 = %.1f' % (Nlabeled,r0))
        axs[1].hist(phis,bins=phibinedges)
        axs[1].plot([phibinedges[0],phibinedges[-1]],[nthresh,nthresh],'r--')
        plt.tight_layout()

    return Nlabeled

from scipy.special import binom
from scipy.optimize import curve_fit

def pn(k,n,p):
    return (binom(n,k)*(np.power(p,k)*np.power((1-p),(n-k))))

def pnpc(k,plabel):
    pbright = 1-pn(0,4,plabel)
    p_k_npc = pn(k,8,pbright)
    
    return p_k_npc

def npclabel_fit(nphist,sigma=None):
    npnormed = nphist/nphist.sum()
    ks = np.arange(9)
    popt, pcov = curve_fit(pnpc, ks, npnormed, sigma=sigma, method='lm', p0=[0.3])
    perr = np.sqrt(np.diag(pcov))
    nlabels_fit = pnpc(ks,popt[0])
    n_labels_scaled = nphist.sum()*nlabels_fit

    return (popt[0],n_labels_scaled,perr[0])
