from traits.api import HasTraits, Str, Int, CStr, List, Enum, Float, Bool
from traitsui.api import View, Item, Group
from traitsui.menu import OKButton, CancelButton, OKCancelButtons
from PYME.DSView.dsviewer import View3D

import numpy as np
from PYMEcs.Analysis.Sofi import calcCorrelates, calcCorrelatesI
import sys

class SOFIconfig(HasTraits):
        numberOfOrders = Int(5)
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
        data = self.dsviewer.image.data # need to check how to best switch to new data handling
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

        View3D(corrs, titleStub='SOFI correlates', mdh=self.dsviewer.image.mdh)
        View3D(means, titleStub='Data mean', mdh=self.dsviewer.image.mdh)


def Plug(dsviewer):
    return SOFIengine(dsviewer)
