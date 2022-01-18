import PYMEcs.Analysis.stackTracker as st
from PYMEcs.misc.guiMsgBoxes import Warn

from traits.api import HasTraits, Str, Int, CStr, List, Enum, Float, Bool

class CalcZfConfig(HasTraits):
    offsetFromCenter = Int(0)
    planesToUse = Int(11)
    
from PYME.DSView.modules._base import Plugin
class CalcZf(Plugin):
    def __init__(self, dsviewer):
        Plugin.__init__(self, dsviewer)
        self.czfconf = CalcZfConfig()
        
        dsviewer.AddMenuItem('Experimental', "Calculate z-factor", self.OnCalcZf)
        dsviewer.AddMenuItem('Experimental', "z-factor calculation settings", self.OnCalcZfSettings)

    def OnCalcZfSettings(self, event):
        self.czfconf.configure_traits(kind='modal')
        
    def OnCalcZf(self, event):
        # first a few checks that the data set is suitable
        # do we need further checks?
        from math import isclose
        if not isclose(self.dsviewer.image.voxelsize_nm.z,50.0,abs_tol=1.0):
            Warn(self.dsviewer,'z voxelsize needs to be 50 nm, is %.1f nm' % self.dsviewer.image.voxelsize_nm.z)
            return

        use_zszh = self.czfconf.planesToUse // 2
        offset = self.czfconf.offsetFromCenter
        need_planes = 2*use_zszh+1 + abs(offset)
        
        zsz = self.dsviewer.image.data_xyztc.shape[2]
        zszh = zsz // 2
        if 2*zszh+1 != zsz:
            Warn(self.dsviewer,'z dimension size should be odd number, is %d' % zsz)
            return
        if zsz < need_planes:
            Warn(self.dsviewer,'need at least %d planes, stack has %d planes' % (need_planes,zsz))
            return

        dataset = self.dsviewer.image.data_xyztc[:,:,zszh+offset-use_zszh:zszh+offset+use_zszh+1,0,0].squeeze()
        
        tstack = st.initialise_data(dataset,4,self.dsviewer.image.voxelsize_nm)
        dxnm,dynm,dznm = st.get_shifts_from_stackobject(tstack)

        st.fit_and_plot_zf(dxnm,dynm,dznm,tstack)

def Plug(dsviewer):
    return CalcZf(dsviewer)
