from traits.api import HasTraits, Str, Int, CStr, List, Enum, Float, Bool
from traitsui.api import View, Item, Group
from traitsui.menu import OKButton, CancelButton, OKCancelButtons
from PYME.DSView.dsviewer import View3D
from PYME.IO.MetaDataHandler import NestedClassMDHandler

import numpy as np
from PYMEcs.Analysis.Sofi import calcCorrelates, calcCorrelatesI
import sys

# plugin interface to David's SOFI code
# current version allows integer zooming but no interface yet to calcCorrelatesZ
# this may better be addressed once the SOFI code data model modernisation has been considered

class SOFIconfig(HasTraits):
        numberOfOrders = Int(3) # so far I have not seen to much reason to include higher orders
        startAtFrame = Int(50)
        stopAtFrame = Int(sys.maxsize)
        filterHalfWidth = Int(25)
        useZooming = Bool(False)
        zoomFactor = Int(1)


from PYME.DSView.modules._base import Plugin
class SOFIengine(Plugin):
    def __init__(self, dsviewer):
        Plugin.__init__(self, dsviewer)
        self.sconf = SOFIconfig(stopAtFrame=self.dsviewer.image.data.shape[2])
        
        dsviewer.AddMenuItem('Experimental>Sofi', "Calculate SOFI moments", self.OnCalcSofi)

    def OnCalcSofi(self, event):
        if not self.sconf.configure_traits(kind='modal'):
            return
        sc = self.sconf
        data = self.dsviewer.image.data # need to check how to best switch to new data model
        if sc.useZooming:
            # note: we need to modify result voxel sizes if zooming!
            corrs, means = calcCorrelatesI(data,
                                           nOrders=sc.numberOfOrders,
                                           zoom=sc.zoomFactor,
                                           startAt=sc.startAtFrame,
                                           stopAt=sc.stopAtFrame,
                                           filtHalfWidth=sc.filterHalfWidth)
        else:
            corrs, means = calcCorrelates(data,
                                          nOrders=sc.numberOfOrders,
                                          startAt=sc.startAtFrame,
                                          stopAt=sc.stopAtFrame,
                                          filtHalfWidth=sc.filterHalfWidth)


        cmdh = NestedClassMDHandler(self.dsviewer.image.mdh)
        cmdh['SofiProcessing.StartFrame'] = sc.startAtFrame
        cmdh['SofiProcessing.StopFrame'] = sc.stopAtFrame
        cmdh['SofiProcessing.FilterHalfWidth'] = sc.filterHalfWidth
        cmdh['SofiProcessing.NumberOfOrders'] = sc.numberOfOrders
        cmdh['SofiProcessing.Mode'] = 'moments'
        if sc.useZooming: # is this the right way to handle voxel sizes? more suitable API interface available?
            cmdh['voxelsize.x'] = self.dsviewer.image.mdh['voxelsize.x'] / float(sc.zoomFactor)
            cmdh['voxelsize.y'] = self.dsviewer.image.mdh['voxelsize.y'] / float(sc.zoomFactor)
            cmdh['SofiProcessing.ZoomFactor'] = sc.zoomFactor
            
        mmdh = NestedClassMDHandler(self.dsviewer.image.mdh)
        mmdh['SofiProcessing.StartFrame'] = sc.startAtFrame
        mmdh['SofiProcessing.StopFrame'] = sc.stopAtFrame
        mmdh['SofiProcessing.Mode'] = 'mean'
        if sc.useZooming:
            mmdh['voxelsize.x'] = self.dsviewer.image.mdh['voxelsize.x'] / float(sc.zoomFactor)
            mmdh['voxelsize.y'] = self.dsviewer.image.mdh['voxelsize.y'] / float(sc.zoomFactor)
            mmdh['SofiProcessing.ZoomFactor'] = sc.zoomFactor

        View3D(corrs, titleStub='SOFI correlates', mdh=cmdh)
        View3D(means, titleStub='Data mean', mdh=mmdh)


def Plug(dsviewer):
    return SOFIengine(dsviewer)
