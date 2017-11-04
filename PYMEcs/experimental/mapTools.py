import wx
from  wx.lib.dialogs import ScrolledMessageDialog

from PYME.localization.remFitBuf import CameraInfoManager
import PYME.Analysis.gen_sCMOS_maps as gmaps
from PYME.IO.MetaDataHandler import NestedClassMDHandler
from PYME.IO.image import ImageStack
from PYME.DSView import ViewIm3D
from PYME.IO.FileUtils import nameUtils

# import PYMEcs.Analysis.pyNCS.denoisetools as ncs

from traits.api import HasTraits, Str, Int, CStr, List, Enum, Float
from traitsui.api import View, Item, Group
from traitsui.menu import OKButton, CancelButton, OKCancelButtons

class cameraChoice(HasTraits):
    clist = List([])
    Camera = Enum(values='clist')

    traits_view = View(Group(Item(name = 'Camera'),
                             label = 'Select Camera',
                             show_border = True),
                       buttons = OKCancelButtons)

    def add_cams(self,camlist):
        for cam in camlist:
            if cam not in self.clist:
                self.clist.append(cam)


def Warn(parent, message, caption = 'Warning!'):
    dlg = wx.MessageDialog(parent, message, caption, wx.OK | wx.ICON_WARNING)
    dlg.ShowModal()
    dlg.Destroy()

defaultCalibrationDir = nameUtils.getCalibrationDir('',create=False)

# FIXME: needs to learn to deal with flatfield maps
def defaultMapName(source, createPath=False, calibrationDir=defaultCalibrationDir):
    resname = source.mdh.getOrDefault('Analysis.resultname',None)
    if resname is None:
        return None

    if resname == 'mean':
        maptype = 'dark'
    else:
        maptype = 'variance'

    mapname = gmaps.mkDefaultPath(maptype, source.mdh, create=createPath, calibrationDir=calibrationDir)
    return mapname

# return a list of existing camera directories that have tif files inside (assuming these are camera maps)
def installedCams(calibrationDir=defaultCalibrationDir):
    from glob import glob
    import os

    camdirs = []
    for (dirpath, dirnames, filenames) in os.walk(calibrationDir):
        camdirs.extend(dirnames)
        break
    fulldirs = [os.path.join(calibrationDir,cdir) for cdir in camdirs]
    validdirs = [cdir for cdir in fulldirs if glob(os.path.join(cdir, '*.tif'))]

    return validdirs


# source passed as PYME ImageStack
def install_map(source, calibrationDir=defaultCalibrationDir):
    """Installs a map file to a calibration directory. By default uses the system claibration directory."""

    import os
    if source.mdh.getOrDefault('Analysis.name', '') != 'mean-variance':
        msg = 'map %s, Analysis.name is not equal to "mean-variance" - probably not a map' % source.filename
        return msg

    validROIHeight = source.mdh.getOrDefault('Analysis.valid.ROIHeight',
                                             source.mdh['Camera.ROIHeight'])
    validROIWidth = source.mdh.getOrDefault('Analysis.valid.ROIWidth',
                                             source.mdh['Camera.ROIWidth'])
    if not (validROIHeight == source.mdh['Camera.ROIHeight']
            and validROIWidth == source.mdh['Camera.ROIWidth']):
        msg = 'Partial (ROI based) maps cannot be installed to a calibration directory'
        return msg

    if source.mdh.getOrDefault('Analysis.isuniform', False):
        msg = 'Uniform maps cannot be installed to ba calibration directory'
        return msg

    if source.mdh['Analysis.resultname'] == 'mean':
        maptype = 'dark'
    else:
        maptype = 'variance'

    mapname = gmaps.mkDefaultPath(maptype, source.mdh, create=True, calibrationDir=calibrationDir)
    if os.path.isfile(mapname):
        msg = 'map %s exists, not overwriting' % mapname
        return msg
    
    source.Save(filename=mapname)
    return None
    
# attempt to install a map in a calibration dir but check a few things
# 1) does it have the .tif extension?
# 2) can it be opened as a tiffstack
# 3) if all ok so far pass on to install_map which does additional checks; note that
#    this function does not overwrite existing maps at the destination
def checkAndInstallMap(mapf, calibrationDir=defaultCalibrationDir):
    import os
    inst = 0
    
    ext = os.path.splitext(mapf)[-1].lower()
    if ext != ".tif":
        msg = 'asked to install %s, not a tif extension' % mapf
        return (inst, msg)
    try:
        source = ImageStack(filename=mapf)
    except:
        msg = 'asked to install %s, could not open as PYME ImageStack, not a map?' % mapf
        return (inst, msg)

    msg = install_map(source, calibrationDir=calibrationDir)
    if msg is None:
        msg = 'installed map %s in default location' % mapf
        msg += "\n\t-> %s" % defaultMapName(source, calibrationDir=calibrationDir)
        inst = 1

    return (inst, msg)


# install maps, potentially several
# if fromfile is a directory attempt to install all maps below that directory
# if fromfile is a file, attempt to install that single file
def installMapsFrom(fromfile, calibrationDir=defaultCalibrationDir):
    from glob import glob
    import os
    
    msgs = []
    ntotal = 0
    if os.path.isdir(fromfile):
        # this is a directory, walk it
        msgs.append("Installing maps from directory %s -> System Folder\n" % fromfile)
        mapfiles = [y for x in os.walk(fromfile) for y in glob(os.path.join(x[0], '*.tif'))]
        for mapf in mapfiles:
            ninst, msg = checkAndInstallMap(mapf, calibrationDir=calibrationDir)
            ntotal += ninst
            msgs.append(msg)
    else:
        ninst, msg = checkAndInstallMap(fromfile, calibrationDir=calibrationDir)
        ntotal = 1
        msgs.append(msg)
    msgs.append("\ninstalled %d maps" % ntotal)
    return ntotal, "\n".join(msgs)

    
def getInstalledMapList():
    from PYME.IO.FileUtils import nameUtils
    from glob import glob
    import os
    rootdir = defaultCalibrationDir
    result = [y for x in os.walk(rootdir) for y in glob(os.path.join(x[0], '*.tif'))]

    return result

class mapTools:
    """
GUI class to supply various map tools
    """
    def __init__(self, dsviewer):
        self.dsviewer = dsviewer
        self.do = dsviewer.do
        self.image = dsviewer.image
        self.ci = CameraInfoManager()
       
        dsviewer.AddMenuItem('Experimental>Map Tools',
                             'Convert current frame to photo-electron counts',
                             self.OnPhotonConvert)
        dsviewer.AddMenuItem('Experimental>Map Tools',
                             'Show dark map',
                             self.OnShowDark)
        dsviewer.AddMenuItem('Experimental>Map Tools',
                             'Show variance map',
                             self.OnShowVariance)
        dsviewer.AddMenuItem('Experimental>Map Tools',
                             'List installed maps',
                             self.OnListMaps)
        dsviewer.AddMenuItem('Experimental>Map Tools',
                             'Install maps to system calibration directory',
                             self.OnInstallMapsToSystem)
        dsviewer.AddMenuItem('Experimental>Map Tools',
                             'Copy maps from system to user directory',
                             self.OnCopyMapsToUserDir)

        
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


    def showMap(self, type='dark'):
        mdh2 = NestedClassMDHandler(self.image.mdh)
        # overwrite the map location with default maps if exist
        if self.check_mapexists(mdh2,type=type) is None:
            Warn(None,'no suitable map in default location')
            return

        if type == 'dark':
            darkf = self.ci.getDarkMap(mdh2)
        else:
            darkf = self.ci.getVarianceMap(mdh2)

        im = ImageStack(darkf, titleStub = '%s Map' % type.capitalize())
        im.mdh.copyEntriesFrom(mdh2)
        #im.mdh['Parent'] = self.image.filename
        #im.mdh['Units'] = 'PhotoElectrons'

        if self.dsviewer.mode == 'visGUI':
            mode = 'visGUI'
        else:
            mode = 'lite'

        dv = ViewIm3D(im, mode=mode, glCanvas=self.dsviewer.glCanvas, parent=wx.GetTopLevelParent(self.dsviewer))
        

    def OnShowDark(self, event=None):
        self.showMap(type='dark')

    def OnShowVariance(self, event=None):
        self.showMap(type='variance')

        
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
