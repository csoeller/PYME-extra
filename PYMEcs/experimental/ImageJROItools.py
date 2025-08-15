import logging
logger = logging.getLogger(__file__)
import wx
import numpy as np
import roifile as rf # new dependency on roifile (available from pypi)
import skimage as ski

from PYME.DSView import ViewIm3D
from PYME.IO.image import ImageStack

class ImageJROItools:

    def __init__(self, dsviewer):
        self.dsviewer = dsviewer
        dsviewer.AddMenuItem('Experimental>ROIs',
                          'Generate ROI mask from Fiji ROISet',
                          self.OnROISet,
                          helpText='Load at FIJI ROISet (e.g. ROISet.zip) and construct a ROI mask from it')

    def OnROISet(self, event=None):
        roi_filename = wx.FileSelector('Load ROI set...',
                                   wildcard="ROISet files (*.zip)|*.zip|ROI files (*.roi)|*.roi", 
                                   flags = wx.FD_OPEN | wx.FD_FILE_MUST_EXIST)
        rois = rf.ImagejRoi.fromfile(roi_filename)
        roimask = np.zeros(self.dsviewer.image.data_xytc.shape[0:2],dtype='int')
        roiszx,roiszy = roimask.shape
        counter = 1
        for roi in rois:
            if roi.roitype in [rf.ROI_TYPE.RECT,
                               rf.ROI_TYPE.OVAL,
                               rf.ROI_TYPE.POLYGON,
                               rf.ROI_TYPE.FREEHAND]: # this needs to be relaxed pretty soonish!
                coords = np.round(roi.coordinates()).astype('int')
                r = coords[:,0]
                c = coords[:,1]
                rr, cc = ski.draw.polygon(r, c)
                roimask[np.clip(cc,0,roiszx-1),np.clip(rr,0,roiszy-1)] = counter
                counter += 1
        ROImaskImg = ImageStack(roimask, titleStub = 'ROI region mask')
        ROImaskImg.mdh.copyEntriesFrom(self.dsviewer.image.mdh)
        # add entries for the mask
        ROImaskImg.mdh['ROISet'] = roi_filename

        ViewIm3D(ROImaskImg, mode='visGUI', title='ROI mask',
                 glCanvas=self.dsviewer.glCanvas, parent=self.dsviewer)

def Plug(dsviewer):
    """Plugs this module into the gui"""
    ImageJROItools(dsviewer)
