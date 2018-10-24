import numpy as np
from scipy import signal
from scipy import stats

def sigmaline(L):
    thresh = (0.2071 * np.sqrt(L) + 1.9102) / (1.2071 * np.sqrt(L) + 0.9102)
    return np.minimum(thresh,1)

def tukey2d(shape,fraction=0.5):
    tx = signal.tukey(shape[0],fraction)
    ty = signal.tukey(shape[1],fraction)

    #tY, tX = np.meshgrid(ty,tx)
    #return tY*tX

    return np.outer(tx,ty)

def spectrum_mean_over_R(pspec,vszx,vszy,binwidth=None, debug = False):
    
    # get size N and M of x and y dims
    N, M = pspec.shape[0:2]
    
    # set up spatial frequency coordinates in x and y [make sure we get the sequence of the two axes right]
    # put on meshgrid and get muR from this
    mux = (np.arange(N) - N/2)/float(N*vszx)
    muy = (np.arange(M) - M/2)/float(M*vszy)
    
    muY, muX = np.meshgrid(muy,mux)
    
    if debug:
        print(N,M)
        print(muX.shape)
    
    muR = np.sqrt(muX*muX+muY*muY)
    
    # calculate bins using the specified binwidth (with suitable default)
    # note that it makes sense to specify binwidth in multiples of dmu
    #      with minimum of 1*dmu (i.e. min of allowed binwidth is 1)
    K = min(N,M)
    mumult = max(int(binwidth),1)
    bins = mumult*np.arange(K/2/mumult)/float(K*vszx) # there would be an issue if vszx ne vszy
    
    # calculate binned_statistics over spectrum
    # the spectrum must be realvalued, either powerspec or real/imaginary part
    means, binedges, bnum = stats.binned_statistic(muR.ravel(),pspec.ravel(),statistic='sum', bins=bins)
    # return bincenters and means for each bin
    binctrs = 0.5*(binedges[1:]+binedges[:-1])
    return (binctrs,means)

def frc(i0,i1,vszx,vszy,muwidth = 2):
    
    t2d = tukey2d(i0.shape,0.25)
    I0 = np.fft.fftshift(np.fft.fftn(i0*t2d))
    I1 = np.fft.fftshift(np.fft.fftn(i1*t2d))
    
    CC = np.real(I0 * np.conj(I1))
    PS0 = np.abs(I0)**2
    PS1 = np.abs(I1)**2
    
    bcc, mcc = spectrum_mean_over_R(CC,vszx,vszy,binwidth=muwidth)
    b0, mi0 = spectrum_mean_over_R(PS0,vszx,vszy,binwidth=muwidth)
    b1, mi1 = spectrum_mean_over_R(PS1,vszx,vszy,binwidth=muwidth)
    # count the number of pixels contributing to each ring
    b2, L = spectrum_mean_over_R(np.ones(PS1.shape),vszx,vszy,binwidth=muwidth)
    
    # in principle should check that bcc, b0, b1 have the same bin locations
    frcv = mcc/np.sqrt(mi0*mi1)
    
    return (bcc,frcv,L)

class FRCplotter:
    def __init__(self, dsviewer):
        self.dsviewer = dsviewer
        dsviewer.AddMenuItem('Experimental>Analysis', 'FRC of image pair', self.OnFRC)

    def OnFRC(self, event=None):
        image = self.dsviewer.image
        im0 = image.data[:,:,:,0].squeeze()
        im1 = image.data[:,:,:,1].squeeze()

        mdh = self.dsviewer.image.mdh
        vx = 1e3*mdh['voxelsize.x']
        vy = 1e3*mdh['voxelsize.y']

        freqs,frc1,L = frc(im0,im1,vx,vy,muwidth = 2)
        import matplotlib.pyplot as plt
        
        plt.figure()
        plt.plot(freqs,frc1)
        plt.plot(freqs,sigmaline(L))
        plt.plot(freqs,1/7.0*np.ones(freqs.shape))
        plt.plot(freqs,np.zeros(freqs.shape),'--')
        plt.xlim(0,freqs[-1])
        plt.xlabel('spatial frequency (nm^-1)')
        plt.ylabel('FRC values')
        plt.show()

def Plug(dsviewer):
    """Plugs this module into the gui"""
    dsviewer.frcplt = FRCplotter(dsviewer)
