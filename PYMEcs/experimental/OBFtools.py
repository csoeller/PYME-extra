import logging
logger = logging.getLogger(__file__)
import wx
import numpy as np

from PYMEcs.pyme_warnings import warn

class OBFtools:

    def __init__(self, dsviewer):
        self.dsviewer = dsviewer
        self.image = self.dsviewer.image

        dsviewer.AddMenuItem('Experimental>OBF-MSR',
                             'Obtain Origin Info from OBF metadata (for .obf and .msr sourced data)',
                             self.OnOBFOrigin,
                             helpText='works only for data loaded from MSR or OBF files; get origin info from metadata; used for image overlays for MINFLUX data')

    def OnOBFOrigin(self, event=None):
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

def Plug(dsviewer):
    """Plugs this module into the gui"""
    OBFtools(dsviewer)
