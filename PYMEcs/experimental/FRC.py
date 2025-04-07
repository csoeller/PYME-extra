import numpy as np
from scipy import signal
from scipy import stats
import wx
from PYMEcs.pyme_warnings import warn

##################
# FRC.py
#
# Copyright Christian Soeller, 2018
# c.soeller@gmail.com
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
##################

# in this implementation we follow the principles as coded in the Fiji FRC plugin coded by Alex Herbert
# the source code is at https://c4science.ch/source/ijp-frc/browse/master/src/main/java/ch/epfl/biop/frc/FRC.java
# it was crucial to have some implementation to look at since not all details of the FRC approach are made
# clear in the FRC papers, for example that the real part of F(i0)*conj(F(i1)) is being used (although it turns out
# that F(i0)*conj(F(i1)) should be real given the symmetry of the Fourier Transform of real valued function)

# also the use of a tukey window is reproduced (described already in the
# Nieuwenhuizen et al 2013 paper)

# the half-bit line formula was copied although
# it would be nice to figure out how the constants where arrived at
# Note there is some formula in the Nat Meth paper supplementary by
# Nieuwenhuizen et al 2013, but that is more general and we need to
# figure out if this one is a version of the general formula with
# specific parameter settings

def sigmaline(L):
    thresh = (0.2071 * np.sqrt(L) + 1.9102) / (1.2071 * np.sqrt(L) + 0.9102)
    return np.minimum(thresh,1)

def tukey2d(shape,fraction=0.5):
    try:
        from scipy.signal import tukey
    except ImportError:
        from scipy.signal.windows import tukey # from some version onwards only available from the scipy.signal.windows namespace
        
    tx = tukey(shape[0],fraction)
    ty = tukey(shape[1],fraction)

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

import statsmodels.api as sm
lowess = sm.nonparametric.lowess

def padsquare(image,newsize=None):
    N, M = image.shape[0:2]
    if newsize is None:
        newsize = max(N,M)
    K = newsize # less typing
    
    if N != newsize or M != newsize:
        if N>newsize or M>newsize:
            raise RuntimeError('trying to embed image (%d,%d) into smaller container (%d,%d)' % (N,M,newsize,newsize)) 
        
        newim = np.zeros((K,K))
        # note that these divisions need to be integer divisons as we do index arithmetic
        startn = (K-N) // 2
        startm = (K-M) // 2
        newim[startn:startn+N,startm:startm+M] = image
        return newim
    else:
        return image

def frc(i0,i1,vszx,vszy,muwidth = 2, zeropad = False, lowessFraction = 1/20.0):
    
    t2d = tukey2d(i0.shape,0.25)

    if zeropad:
        im0 = padsquare(i0*t2d)
        im1 = padsquare(i1*t2d)
    else:
        im0 = i0*t2d
        im1 = i1*t2d

    I0 = np.fft.fftshift(np.fft.fftn(im0))
    I1 = np.fft.fftshift(np.fft.fftn(im1))
    
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

    smoothed = lowess(frcv, bcc, frac=lowessFraction, it=3, delta=0.0,
                      is_sorted=False, missing='drop', return_sorted=False)
    
    return (bcc,frcv,smoothed,L)

def frc_from_image(image,channels,muwidth=2,zeropad=True,lowessFraction=5.0*0.01):
    im0 = image.data_xytc[:,:,:,channels[0]].squeeze()
    im1 = image.data_xytc[:,:,:,channels[1]].squeeze()

    mdh = image.mdh
    vx = 1e3*mdh['voxelsize.x']
    vy = 1e3*mdh['voxelsize.y']
    freqs,frc1,smoothed,L = frc(im0,im1,vx,vy,muwidth=2,zeropad=True,
                                lowessFraction=lowessFraction)
    halfbit = sigmaline(L)
    fhb = zc.zerocross1d(freqs,smoothed-halfbit)
    f7= zc.zerocross1d(freqs,smoothed-1.0/7.0)
    
    return (freqs,frc1,smoothed,fhb,f7,halfbit)

def frc_plot(freqs,frc1,smoothed,fhb,f7,halfbit,chanNames=['block0','block1'],
             showGrid=True,showHalfbitThreshold=False):
    import matplotlib.pyplot as plt
        
    from mpl_toolkits.axes_grid1 import host_subplot
    import mpl_toolkits.axisartist as AA

    plt.figure()
    ax = host_subplot(111, axes_class=AA.Axes)
    ax.plot(freqs,frc1)
    ax.plot(freqs,smoothed)
    if showHalfbitThreshold:
        ax.plot(freqs,halfbit)
    ax.plot(freqs,1/7.0*np.ones(freqs.shape))

    if len(f7) > 0:
        ax.plot([f7[0],f7[0]],[0,1],'--')
        ax.text(0.55, 0.8,'res-1/7     = %3.1f nm' % (1.0/f7[0]),horizontalalignment='left',
                     verticalalignment='center', transform=plt.gca().transAxes)
    else:
        ax.text(0.55, 0.8,'res-1/7 - no intercept',horizontalalignment='left',
                    verticalalignment='center', transform=plt.gca().transAxes)

    if showHalfbitThreshold:
        if len(fhb) > 0:
            ax.plot([fhb[0],fhb[0]],[0,1],'--')
            ax.text(0.55, 0.9,'res-halfbit = %3.1f nm' % (1.0/fhb[0]),horizontalalignment='left',
                    verticalalignment='center', transform=plt.gca().transAxes)
        else:
            ax.text(0.55, 0.9,'res-halfbit - no intercept',horizontalalignment='left',
                    verticalalignment='center', transform=plt.gca().transAxes)

    ax.plot(freqs,np.zeros(freqs.shape),'--')
    ax.set_xlim(0,freqs[-1])
    ax.set_xlabel('spatial frequency (nm^-1)')
    ax.set_ylabel('FRC values')
    xt = np.array([10., 15, 20, 30, 40, 50, 75, 150])
    ft = 1.0 / xt

    ax2 = ax.twin()  # ax2 is responsible for "top" axis and "right" axis
    ax2.set_xticks(ft[::-1])
    ax2.set_xticklabels(['%d' % xi for xi in xt[::-1]],rotation=90)
    ax2.set_xlabel('resolution (nm)')
    if showGrid:
        ax2.grid(True)
        
    ax2.axis["right"].major_ticklabels.set_visible(False)
    ax2.axis["top"].major_ticklabels.set_visible(True)

    plt.title("FRC for channels %s and %s" % (chanNames[0],chanNames[1]), y=1.08)
    plt.show()


def save_vol_mrc(data, grid_spacing, outfilename, origin=None, overwrite=False):
    import mrcfile # will bomb unless this is installed
    
    data = data.astype('float32') 
    with mrcfile.new(outfilename, overwrite=overwrite) as mrc:
        mrc.set_data(data)
        mrc.voxel_size = grid_spacing
        if origin is not None:
            mrc.header.origin = origin
        mrc.update_header_from_data()
        mrc.update_header_stats()

import PYMEcs.Analysis.zerocross as zc
from traits.api import HasTraits, Str, Int, CStr, List, Enum, Float, Bool

class FRCsettings(HasTraits):
    ZeroPadding = Bool(True)
    ShowHalfbitThreshold = Bool(False)
    ShowGrid = Bool(True)
    LowessSmoothingPercentage = Float(5.0)
    
class FRCplotter:
    def __init__(self, dsviewer):
        self.dsviewer = dsviewer
        dsviewer.AddMenuItem('Experimental>Analysis', 'FRC of image pair', self.OnFRC)
        dsviewer.AddMenuItem('Experimental>Analysis', 'save last FRC curves', self.OnFRCSave)
        dsviewer.AddMenuItem('Experimental>Analysis', 'adjust FRC settings', self.OnFRCSettings)
        dsviewer.AddMenuItem('Experimental>Analysis', 'save MRC volumes for FSC', self.OnFSCsave_as_MRC)
        
        self.frcSettings = FRCsettings()
        self.lastFRC = None
        
    def OnFRCSettings(self, event=None):
        if self.frcSettings.configure_traits(kind='modal'):
            pass

    def OnFRC(self, event=None):
        from PYME.DSView.modules.coloc import ColocSettingsDialog
        image = self.dsviewer.image
        voxelsize = image.voxelsize

        try:
            names = image.mdh.getEntry('ChannelNames')
        except:
            names = ['Channel %d' % n for n in range(image.data_xytc.shape[3])]

        with ColocSettingsDialog(self.dsviewer, voxelsize[0], names, show_bins=False) as dlg:
            if dlg.ShowModal() != wx.ID_OK:
                return
            chans = dlg.GetChans()

        chanNames = [names[chans[0]],names[chans[1]]]

        freqs,frc1,smoothed,fhb,f7,halfbit = frc_from_image(image,chans,muwidth=2,
                                                    zeropad=self.frcSettings.ZeroPadding,
                                                    lowessFraction=self.frcSettings.LowessSmoothingPercentage*0.01)
        frc_plot(freqs,frc1,smoothed,fhb,f7,halfbit,chanNames,
                 showHalfbitThreshold=self.frcSettings.ShowHalfbitThreshold,
                 showGrid=self.frcSettings.ShowGrid)
        
        self.lastFRC = dict(frequencies=freqs,frc_curve=frc1,frc_smoothed=smoothed,frchbval=fhb,frcval=f7,frc_halfbit=halfbit)

    def OnFRCSave(self, event=None):
        if self.lastFRC is None:
            warn("could not find data from prior FRC Plot; cannot scave, please plot first")
            return
        with wx.FileDialog(self.dsviewer, 'Save FRC data as ...',
                           wildcard='CSV (*.csv)|*.csv',
                           style=wx.FD_SAVE) as fdialog:
            if fdialog.ShowModal() != wx.ID_OK:
                return
            else:
                fpath = fdialog.GetPath()

        lfrc = self.lastFRC
        import pandas as pd
        df = pd.DataFrame.from_dict(dict(freqs=lfrc['frequencies'],frc=lfrc['frc_curve'],frcsmooth=lfrc['frc_smoothed']))
        df.to_csv(fpath,index=False)
        
    def OnFSCsave_as_MRC(self, event=None):
        from PYME.DSView.modules.coloc import ColocSettingsDialog
        im = self.dsviewer.image
        try:
            names = im.mdh.getEntry('ChannelNames')
        except:
            names = ['Channel %d' % n for n in range(im.data_xyztc.shape[4])]

        with ColocSettingsDialog(self.dsviewer, im.mdh.voxelsize.x, names, show_bins=False) as dlg:
            if dlg.ShowModal() != wx.ID_OK:
                return
            chans = dlg.GetChans()

        chanNames = [names[chans[0]],names[chans[1]]]

        vol1 = im.data_xyztc[:,:,:,0,chans[0]].squeeze()
        vol2 = im.data_xyztc[:,:,:,0,chans[1]].squeeze()

        MINFLUXts = im.mdh.get('Source.MINFLUX.TimeStamp')
        if MINFLUXts is not None:
            defaultFiles = ["%s-%s.mrc" % (MINFLUXts,names[chans[0]]),
                            "%s-%s.mrc" % (MINFLUXts,names[chans[1]])]
        else:
            defaultFiles = ['','']

        for i,vol in enumerate([vol1,vol2]):
            with wx.FileDialog(self.dsviewer, 'Save channel %d as ...' % i,
                                    wildcard='MRC (*.mrc)|*.mrc',
                                    defaultFile=defaultFiles[i],
                                    style=wx.FD_SAVE) as fdialog:
                if fdialog.ShowModal() != wx.ID_OK:
                    return
                fpath = fdialog.GetPath()
            save_vol_mrc(vol.T,im.mdh.voxelsize_nm.x,fpath)

        # we also write a JSON file of the metadata so that we can be sure about voxelsizes etc
        # currently just the stem of the second MRC file name with json extension
        from pathlib import Path
        p = Path(fpath)
        p.with_suffix('.json').write_text(im.mdh.to_JSON())
        
def Plug(dsviewer):
    """Plugs this module into the gui"""
    dsviewer.frcplt = FRCplotter(dsviewer)
