from PYME.IO.MetaDataHandler import NestedClassMDHandler
from PYME.IO.image import ImageStack
from PYME.DSView import ViewIm3D
import numpy as np
import wx
from PYMEcs.misc.guiMsgBoxes import Warn
from scipy.ndimage import gaussian_filter

from traits.api import HasTraits, Str, Int, CStr, List, Enum, Float
#from traitsui.api import View, Item, Group
#from traitsui.menu import OKButton, CancelButton, OKCancelButtons

class procPtsChoice(HasTraits):
    minimal_distance_in_nm = Float(30)
    sigma_for_Gaussian_blur = Float(2)


class procPts:
    """
GUI class to process point positions produced by blobfind 
    """
    def __init__(self, dsviewer):
        self.dsviewer = dsviewer
        self.procPtsSel = procPtsChoice()

        dsviewer.AddMenuItem('Experimental>Object Processing',
                             'Set all pixels of detected object centers in new image',
                             self.OnProcPoints)

    def OnProcPoints(self, event=None):
        try:
            points = self.dsviewer.view.points
        except AttributeError:
            Warn(None, 'no object locations found')
            return

        if len(self.dsviewer.view.points) < 1:
            Warn(None, 'object location list empty')
            return
        
        if not self.procPtsSel.configure_traits(kind='modal'):
            return
        
        xp = np.rint(self.dsviewer.view.points[:,0])
        yp = np.rint(self.dsviewer.view.points[:,1])

        mindistnm = self.procPtsSel.minimal_distance_in_nm # make this an option
        mindistpix = mindistnm / (1e3*self.dsviewer.image.mdh.voxelsize.x)
        # here we need some code to remove points with NND < mindist
        xd = np.subtract.outer(xp,xp)
        yd = np.subtract.outer(yp,yp)
        d = np.sqrt(xd**2+yd**2)
        np.fill_diagonal(d,1e6)
        dmin = d.min(0)
        xp2 = xp[dmin >= mindistpix]
        yp2 = yp[dmin >= mindistpix]
        
        imd = np.zeros(self.dsviewer.image.data.shape[0:2])
        imd[xp2.astype('i'),yp2.astype('i')] = 1
        sigma = self.procPtsSel.sigma_for_Gaussian_blur # in future may make this a config parameter
        imf = gaussian_filter(imd,sigma)
        
        im = ImageStack(imf, titleStub = 'object positions')
        mdh2 = NestedClassMDHandler(self.dsviewer.image.mdh)
        im.mdh.copyEntriesFrom(mdh2)

        if self.dsviewer.mode == 'visGUI':
            mode = 'visGUI'
        else:
            mode = 'lite'

        dv = ViewIm3D(im, mode=mode, glCanvas=self.dsviewer.glCanvas,
                      parent=wx.GetTopLevelParent(self.dsviewer))


def Plug(dsviewer):
    procPts(dsviewer)
