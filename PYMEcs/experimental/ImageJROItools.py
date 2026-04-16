import logging
logger = logging.getLogger(__file__)
import wx
import numpy as np
import roifile as rf # new dependency on roifile (available from pypi)
import skimage as ski

from PYME.DSView import ViewIm3D
from PYME.IO.image import ImageStack
from PYMEcs.pyme_warnings import warn

class ImageJROItools:

    def __init__(self, dsviewer):
        self.dsviewer = dsviewer
        self.image = self.dsviewer.image
        
        dsviewer.AddMenuItem('Experimental>ROIs',
                          'Generate ROI mask from Fiji ROISet',
                          self.OnROISet,
                          helpText='Load at FIJI ROISet (e.g. ROISet.zip) and construct a ROI mask from it')
        dsviewer.AddMenuItem('Experimental>Overlays',
                          'Obtain Origin Info from bioformats metadata',
                          self.OnBFOrigin,
                          helpText='works only for data loaded from MSR files; get origin info from metadata; used for image overlays for MINFLUX data')

    def OnBFOrigin(self, event=None):
        mdh = self.image.mdh
        # now modified to pure obf version
        # i.e. no dependency on bffile
        
        # try:
        #     import bffile  # noqa: F401 - only test availability here
        # except ImportError:
        #     warn('bffile not available, required to retrieve bioformats metadata - install bffile to use this functionality')
        #     return

        ds = self.image.data
        try:
            obf = ds.obf
        except AttributeError:
            warn('not OBF data, giving up')
            return
        try:
            series_no = ds.stack_number
        except AttributeError:
            warn('stack number should be in datasource, not found - giving up')
            return
   
        #fn = self.image.filename.split(':')[0]
        #from bffile import BioFile
        #bf = BioFile(fn)
        #bf.open()
        n_images = len(obf.stacks)
        if series_no >= n_images:
            warn('stack number should be less than number of images %d, but is %d - giving up' % (n_images,series_no))
            return
        #img_meta = bf.core_metadata(series=series_no)
        #ds.img_meta = img_meta
        #offsets = img_meta.series_metadata['Offsets']
        offsets = obf.stacks[series_no].offsets
        mdh['Origin.x'] = 1e9*offsets[0] # convert to nm
        mdh['Origin.y'] = 1e9*offsets[1] # convert to nm

        zsize = self.image.data_xyztc.shape[2]
        if zsize > 1:
            origin_z = - self.image.voxelsize_nm.z * zsize/2 # so that center of stack will appear at ~0 in VisGUI
        else:
            origin_z = 0
        mdh['Origin.z'] = origin_z

        warn("Informational message only: \n\nadded image origin at\n\n    %.1f, %.1f\n\nfrom MSR metdata" % (mdh['Origin.x'],mdh['Origin.y']))
        
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
