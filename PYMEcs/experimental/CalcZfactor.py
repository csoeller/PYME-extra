import PYMEcs.Analysis.stackTracker as st
from PYMEcs.misc.guiMsgBoxes import Warn

from PYME.DSView.modules._base import Plugin
class CalcZf(Plugin):
    def __init__(self, dsviewer):
        Plugin.__init__(self, dsviewer)
        
        dsviewer.AddMenuItem('Experimental', "Calculate z-factor", self.OnCalcZf)


    def OnCalcZf(self, event):
        # TODO: first need to do add a few checks that the data set is suitable!
        from math import isclose
        if not isclose(self.dsviewer.image.voxelsize_nm.z,50.0,abs_tol=1.0):
            Warn(self.dsviewer,'z voxelsize needs to be 50 nm, is %.1f nm' % self.dsviewer.image.voxelsize_nm.z)
            return
        
        dataset = self.dsviewer.image.data_xyztc[:,:,:,0,0].squeeze()
        tstack = st.initialise_data(dataset,4,self.dsviewer.image.voxelsize_nm)
        dxnm,dynm,dznm = st.get_shifts_from_stackobject(tstack)

        st.fit_and_plot_zf(dxnm,dynm,dznm,tstack)

def Plug(dsviewer):
    return CalcZf(dsviewer)
