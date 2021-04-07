import wx
from  wx.lib.dialogs import ScrolledMessageDialog

from PYME.localization.remFitBuf import CameraInfoManager
import PYME.Analysis.gen_sCMOS_maps as gmaps
from PYME.IO.MetaDataHandler import NestedClassMDHandler
from PYME.IO.image import ImageStack
from PYME.DSView import ViewIm3D
import numpy as np

from PYMEcs.misc.mapUtils import defaultCalibrationDir,defaultMapName,installedCams,install_map,\
    checkAndInstallMap,installMapsFrom,getInstalledMapList,check_mapexists,mk_compositeMap,\
    addMap2composite,export_mapFromComposite,_getDefaultMap,get_dark_default,get_variance_default,\
    get_flatfield_default

import logging
logger = logging.getLogger(__name__)

# FUNCTIONAILITY that would be good to add
# - check all installed maps for sanity etc
# - install a single map file (rather than whole directory)
# - force option for install maps commands
# generally: better API for map functions than current ones in
#        PYME.Analysis.gen_sCMOS_maps
#        CameraInfoManager in PYME.localization.remFitBuf



from traits.api import HasTraits, Str, Int, CStr, List, Enum, Float
from traitsui.api import View, Item, Group
from traitsui.menu import OKButton, CancelButton, OKCancelButtons

from PYMEcs.misc.guiMsgBoxes import Warn

class cameraChoice(HasTraits):
    _clist = List([])
    Camera = Enum(values='_clist')

    traits_view = View(Group(Item(name = 'Camera'),
                              label = 'Select Camera',
                             show_border = True),
                       buttons = OKCancelButtons)

    def add_cams(self,camlist):
        for cam in camlist:
            if cam not in self._clist:
                self._clist.append(cam)

class meanVarianceCalc(HasTraits):
    Start = Int(0)
    End = Int(-1)


class mapTools:
    """
GUI class to supply various map tools
    """
    def __init__(self, dsviewer):
        self.dsviewer = dsviewer
        self.do = dsviewer.do
        self.image = dsviewer.image
        self.ci = CameraInfoManager()
        self.loadedMaps = {}
       
        dsviewer.AddMenuItem('Experimental>Map Tools',
                             'Convert current frame to photo-electron counts',
                             self.OnPhotonConvert)
        dsviewer.AddMenuItem('Experimental>Map Tools',
                             'Calculate Mean and Variance of frame sequence',
                             self.OnMeanVariance)
        dsviewer.AddMenuItem('Experimental>Map Tools',
                             'Show dark map',
                             self.OnShowDark)
        dsviewer.AddMenuItem('Experimental>Map Tools',
                             'Show variance map',
                             self.OnShowVariance)
        dsviewer.AddMenuItem('Experimental>Map Tools',
                             'Show flatfield map',
                             self.OnShowFlatField)
        dsviewer.AddMenuItem('Experimental>Map Tools',
                             'List installed maps',
                             self.OnListMaps)
        dsviewer.AddMenuItem('Experimental>Map Tools',
                             'Install maps to system calibration directory',
                             self.OnInstallMapsToSystem)
        dsviewer.AddMenuItem('Experimental>Map Tools',
                             'Copy maps from system to user directory',
                             self.OnCopyMapsToUserDir)

            
    def getMapSafely(self, type='dark'):
        try:
            mdh = self.dsviewer.LMAnalyser.analysisController.analysisMDH
        except AttributeError:
            mdh = self.image.mdh
        try:
            theMap = self.ci.getDarkMap(mdh)
            self.loadedMaps[type] = mdh
        except IOError:
            logger.exception('Dark map specified but not found, falling back on defaults')
            (theMap,mdh2) = _getDefaultMap(self.ci,mdh,type,return_loadedmdh=True)
            self.loadedMaps[type] = mdh2
        return theMap

        
    def showMap(self, type='dark'):
        mdh2 = NestedClassMDHandler(self.image.mdh)
        # # overwrite the map location with default maps if exist
        # if check_mapexists(mdh2,type=type) is None:
        #     Warn(None,'no suitable map in default location')
        #     return

        check_mapexists(mdh2,type=type)
        
        theMap = self.getMapSafely(type)
        twoDMap = np.ones(self.image.data.shape[0:2])*theMap # promote to 2D if necessary
        
        im = ImageStack(twoDMap, titleStub = '%s Map' % type.capitalize())
        im.mdh.copyEntriesFrom(mdh2)

        if self.dsviewer.mode == 'visGUI':
            mode = 'visGUI'
        else:
            mode = 'lite'

        dv = ViewIm3D(im, mode=mode, glCanvas=self.dsviewer.glCanvas, parent=wx.GetTopLevelParent(self.dsviewer))
        

    def OnShowDark(self, event=None):
        self.showMap(type='dark')

    def OnShowVariance(self, event=None):
        self.showMap(type='variance')

    def OnShowFlatField(self, event=None):
        self.showMap(type='flatfield')

    def OnMeanVariance(self, event=None):
        mvChoice = meanVarianceCalc(End=self.image.data.shape[2]-1)
        if not mvChoice.configure_traits(kind='modal'):
            return
        
        m, v = gmaps._meanvards(self.image.data,start=mvChoice.Start,end=mvChoice.End)
        mmdh = NestedClassMDHandler(self.image.mdh)
        mmdh.setEntry('Analysis.name', 'mean-variance')
        mmdh.setEntry('Analysis.start', mvChoice.Start)
        mmdh.setEntry('Analysis.end', mvChoice.End)
        mmdh.setEntry('Analysis.resultname', 'mean')
        mmdh.setEntry('Analysis.units', 'ADU')
        
        imm = ImageStack(m, mdh=mmdh, titleStub = 'Mean')

        vmdh = NestedClassMDHandler(mmdh)
        vmdh.setEntry('Analysis.resultname', 'variance')
        vmdh.setEntry('Analysis.units', 'ADU^2')

        imv = ImageStack(v, mdh=vmdh, titleStub = 'Variance')
        
        if self.dsviewer.mode == 'visGUI':
            mode = 'visGUI'
        else:
            mode = 'lite'

        dv = ViewIm3D(imm, mode=mode, glCanvas=self.dsviewer.glCanvas,
                      parent=wx.GetTopLevelParent(self.dsviewer))
        dv = ViewIm3D(imv, mode=mode, glCanvas=self.dsviewer.glCanvas,
                      parent=wx.GetTopLevelParent(self.dsviewer))
        
    def OnPhotonConvert(self, event=None):

        # we try color channel 0; should only be done on monochrome anyway
        curFrame = self.image.data[:,:, self.do.zp, 0].squeeze()

        # this makes a new metadata structure that copies all entries from the argument
        mdh2 = NestedClassMDHandler(self.image.mdh)
        # some old files do not have a camera serialname
        # fake one, which ensures no map is found and we get uniform maps
        try:
            t = mdh2['Camera.SerialNumber']
        except AttributeError:
             mdh2['Camera.SerialNumber'] = 'XXXXX'

        # overwrite the map location with default maps if exist
        check_mapexists(mdh2,type='dark')
        check_mapexists(mdh2,type='flatfield')

        #darkf = self.ci.getDarkMap(mdh2)
        darkf = self.getMapSafely(type='dark')
        corrFrame = float(mdh2['Camera.ElectronsPerCount'])*self.ci.correctImage(mdh2, curFrame)/mdh2.getEntry('Camera.TrueEMGain')

        im = ImageStack(corrFrame, titleStub = 'Frame %d in photoelectron units' % self.do.zp)
        im.mdh.copyEntriesFrom(mdh2)
        im.mdh['Parent'] = self.image.filename
        im.mdh['Units'] = 'PhotoElectrons'
        im.mdh['Camera.ElectronsPerCount'] = 1.0
        im.mdh['Camera.TrueEMGain'] = 1.0
        im.mdh['Camera.ADOffset'] = 0

        if self.dsviewer.mode == 'visGUI':
            mode = 'visGUI'
        else:
            mode = 'lite'

        dv = ViewIm3D(im, mode=mode, glCanvas=self.dsviewer.glCanvas, parent=wx.GetTopLevelParent(self.dsviewer))

        
    def OnListMaps(self, event=None):
      
      maps = getInstalledMapList()
      if len(maps) > 0:
          dlg = ScrolledMessageDialog(self.dsviewer, "\n".join(maps), "Installed maps", size=(900,400),
                                      style=wx.RESIZE_BORDER | wx.DEFAULT_DIALOG_STYLE )
          dlg.ShowModal()
          dlg.Destroy()
      else:
          Warn(None,'no suitable maps found')

        
    def OnInstallMapsToSystem(self, event=None):
        instmsg = 'Install maps from user directory...'
        fdialog = wx.DirDialog(None, instmsg,
                               style=wx.DD_DEFAULT_STYLE|wx.DD_DIR_MUST_EXIST|wx.DD_CHANGE_DIR)

        if fdialog.ShowModal() == wx.ID_OK:
            dirSelection = fdialog.GetPath().encode()
            fdialog.Destroy()
        else:
            fdialog.Destroy()
            return        

        inst, msg = installMapsFrom(dirSelection, calibrationDir=defaultCalibrationDir)
        
        dlg = ScrolledMessageDialog(self.dsviewer, msg, instmsg, size=(900,400),
                                    style=wx.RESIZE_BORDER | wx.DEFAULT_DIALOG_STYLE )
        dlg.ShowModal()
        dlg.Destroy()

        
    def OnCopyMapsToUserDir(self, event=None):
        import os
        instmsg = 'Copy maps to user directory...'
        cdict = {os.path.basename(camdir) : camdir for camdir in installedCams()}
        cdict['All Cameras'] = defaultCalibrationDir
        
        cChoice = cameraChoice()
        cChoice.add_cams(sorted(cdict.keys()))
        if not cChoice.configure_traits(kind='modal'):
            return
        camdir = cdict[cChoice.Camera]

        fdialog = wx.DirDialog(None, instmsg,
                               style=wx.DD_DEFAULT_STYLE|wx.DD_DIR_MUST_EXIST|wx.DD_CHANGE_DIR)

        if fdialog.ShowModal() == wx.ID_OK:
            dirSelection = fdialog.GetPath().encode()
            fdialog.Destroy()
        else:
            fdialog.Destroy()
            return

        inst, msg = installMapsFrom(camdir, calibrationDir=dirSelection)
        
        dlg = ScrolledMessageDialog(self.dsviewer, msg, instmsg, size=(900,400),
                                    style=wx.RESIZE_BORDER | wx.DEFAULT_DIALOG_STYLE )
        dlg.ShowModal()
        dlg.Destroy()

  

def Plug(dsviewer):
    dsviewer.mapTool = mapTools(dsviewer)
