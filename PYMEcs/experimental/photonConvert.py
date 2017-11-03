import wx
from  wx.lib.dialogs import ScrolledMessageDialog

from PYME.localization.remFitBuf import CameraInfoManager
import PYME.Analysis.gen_sCMOS_maps as gmaps
from PYME.IO.MetaDataHandler import NestedClassMDHandler
from PYME.IO.image import ImageStack
from PYME.DSView import ViewIm3D

# import PYMEcs.Analysis.pyNCS.denoisetools as ncs

def Warn(parent, message, caption = 'Warning!'):
    dlg = wx.MessageDialog(parent, message, caption, wx.OK | wx.ICON_WARNING)
    dlg.ShowModal()
    dlg.Destroy()

class photonConverter:
    """
GUI class to convert a single PYME frame to photoelectron units
    """
    def __init__(self, dsviewer):
        self.dsviewer = dsviewer
        self.do = dsviewer.do
        self.image = dsviewer.image
        self.ci = CameraInfoManager()
       
        dsviewer.AddMenuItem('Experimental>Maps',
                             'Convert current frame to photo-electron counts',
                             self.OnPhotonConvert)
        dsviewer.AddMenuItem('Experimental>Maps',
                             'Show dark map',
                             self.OnShowDark)

    def check_mapexists(self, mdh, type = 'dark'):
        import os
        if type == 'dark':
            id = 'Camera.DarkMapID'
        else:
            id = 'Camera.VarianceMapID'
            
        mapPath = gmaps.mkDefaultPath(type,mdh,create=False)
        if os.path.exists(mapPath):
            mdh[id] = mapPath
            return mapPath
        else:
            return None

    def OnShowDark(self, event=None):
        mdh2 = NestedClassMDHandler(self.image.mdh)
        # overwrite the map location with default maps if exist
        if self.check_mapexists(mdh2,type='dark') is None:
            Warn(None,'no suitable map in default location')
            return

        darkf = self.ci.getDarkMap(mdh2)
        im = ImageStack(darkf, titleStub = 'Dark Map')
        im.mdh.copyEntriesFrom(mdh2)
        #im.mdh['Parent'] = self.image.filename
        #im.mdh['Units'] = 'PhotoElectrons'

        if self.dsviewer.mode == 'visGUI':
            mode = 'visGUI'
        else:
            mode = 'lite'

        dv = ViewIm3D(im, mode=mode, glCanvas=self.dsviewer.glCanvas, parent=wx.GetTopLevelParent(self.dsviewer))
        

        
    def OnPhotonConvert(self, event=None):

        # we try color channel 0; should only be done on monochrome anyway
        curFrame = self.image.data[:,:, self.do.zp, 0].squeeze()

        # this makes a new metadata structure that copies all entries from the argument
        mdh2 = NestedClassMDHandler(self.image.mdh)
        # overwrite the map location with default maps if exist
        self.check_mapexists(mdh2,type='dark')

        darkf = self.ci.getDarkMap(mdh2)
        corrFrame = float(mdh2['Camera.ElectronsPerCount'])*self.ci.correctImage(mdh2, curFrame)

        im = ImageStack(corrFrame, titleStub = 'Frame %d in photoelectron units' % self.do.zp)
        im.mdh.copyEntriesFrom(mdh2)
        im.mdh['Parent'] = self.image.filename
        im.mdh['Units'] = 'PhotoElectrons'

        if self.dsviewer.mode == 'visGUI':
            mode = 'visGUI'
        else:
            mode = 'lite'

        dv = ViewIm3D(im, mode=mode, glCanvas=self.dsviewer.glCanvas, parent=wx.GetTopLevelParent(self.dsviewer))

    # def OnListMaps(self, event=None):
      
    #   with open(self.getLogFileName(),"r") as f:
    #      txt = "\n".join(f.readlines())
    #   dlg = ScrolledMessageDialog(self.visFr, txt, "VisGUI Error Output", size=(900,400),
    #                               style=wx.RESIZE_BORDER | wx.DEFAULT_DIALOG_STYLE )
    #   dlg.ShowModal()
    #   dlg.Destroy()

def Plug(dsviewer):
    dsviewer.photonConv = photonConverter(dsviewer)
