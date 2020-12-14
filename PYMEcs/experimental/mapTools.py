import wx
from  wx.lib.dialogs import ScrolledMessageDialog

from PYME.localization.remFitBuf import CameraInfoManager
import PYME.Analysis.gen_sCMOS_maps as gmaps
from PYME.IO.MetaDataHandler import NestedClassMDHandler
from PYME.IO.image import ImageStack
from PYME.DSView import ViewIm3D
from PYME.IO.FileUtils import nameUtils
import numpy as np

# FUNCTIONAILITY that would be good to add
# - check all installed maps for sanity etc
# - install a single map file (rather than whole directory)
# - force option for install maps commands
# generally: better API for map functions than current ones in
#        PYME.Analysis.gen_sCMOS_maps
#        CameraInfoManager in PYME.localization.remFitBuf


# the NCS functionality needs the pyNCS in the path, either linked
# into the PYMEcs.Analysis directory (e.g. via symlink) or as pyNCS in the PYTHONPATH
# the code can be obtained from https://github.com/HuanglabPurdue/NCS,
#   in the python3-6 directory
# my testing shows that the implementation runs fine under python-2.7
try:
    import PYMEcs.Analysis.pyNCS.denoisetools as ncs
except ImportError:
    try:
        import pyNCS.denoisetools as ncs
    except ImportError:
        ncs = None

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


class ncsSelect(HasTraits):
    winSize = Enum(64,128,256)
    Rs = Int(8) # IMPORTANT: looks like the pyncs code assumes Rs is a divider of the imgsz (see imagsz2 below)!!
    Lambda_nm = Float(690) # emission wavelength
    NA = Float(1.49)     # objective NA
    iterations = Int(15)
    alpha = Float(0.2)

defaultCalibrationDir = nameUtils.getCalibrationDir('',create=False)

def defaultMapName(source, createPath=False, calibrationDir=defaultCalibrationDir):
    resname = source.mdh.getOrDefault('Analysis.resultname',None)
    if resname is None:
        return None

    if resname == 'mean':
        if source.mdh.getOrDefault('Analysis.FlatField',False):
            maptype = 'flatfield'
        else:
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
        if source.mdh.getOrDefault('Analysis.FlatField',False):
            maptype = 'flatfield'
        else:
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
        msgs.append("Installing maps from directory %s -> %s Folder\n" % (fromfile,calibrationDir))
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


def check_mapexists(mdh, type = 'dark'):
    import os
    if type == 'dark':
        id = 'Camera.DarkMapID'
    elif type == 'variance':
        id = 'Camera.VarianceMapID'
    elif type == 'flatfield':
        id = 'Camera.FlatfieldMapID'
    else:
        raise RuntimeError('unknown map type %s' % type)
        
    mapPath = gmaps.mkDefaultPath(type,mdh,create=False)
    if os.path.exists(mapPath):
        mdh[id] = mapPath
        return mapPath
    else:
        return None

def mk_compositeMap(sourcemap):
    # make new composite map using sourcemap to populate
    # check source is a map
    # composite map has first channel with map data, second channel is mask where map is valid
    # returns composite map (identified by its metadata)
    pass

def addMap2composite(map,compMap):
    # insert map according to valid ROI into composite map
    # update the valid channel accordingly
    # possibly perform a number of checks
    pass

def export_mapFromComposite(compMap):
    # export a normal map from the composite map
    # in the exported map we set the valid ROI to the whole chip area
    # return as image?
    pass


class mapTools:
    """
GUI class to supply various map tools
    """
    def __init__(self, dsviewer):
        self.dsviewer = dsviewer
        self.do = dsviewer.do
        self.image = dsviewer.image
        self.ci = CameraInfoManager()
        self.ncsSel = ncsSelect() # by making it part of the object we retain parameters across invocations
       
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
        dsviewer.AddMenuItem('Experimental>Map Tools',
                             'NCS denoising of small square ROI',
                             self.OnNCSDenoise)

    def OnNCSDenoise(self, event=None):
        ncsSel = self.ncsSel
        if not ncsSel.configure_traits(kind='modal'):
            return
        
        mdh = self.image.mdh
        img = self.image.data[:,:, self.do.zp, 0].squeeze() # current frame

        # code below from pyncs.addnoise
        # I is presumably average photon number; bg is background photons/pixel
        # idealimg = np.abs(normimg)*I+bg
        # poissonimg = np.random.poisson(idealimg)

        # UNITS: scmosimg (ADUs) = electrons * [gain-as-ADU/electrons / flatfield] + readnoise(ADUs) + offset
        # from which follows gain_for_pyNCS = [gain-as-ADU/electrons / flatfield]
        # scmosimg = poissonimg*gainmap + np.sqrt(varmap)*np.random.randn(R,R)
        # scmosimg += offset

        if check_mapexists(mdh,'dark') is None:
            Warn(None, 'ncsimage: no dark map found')
            return
        else:
            dark = self.ci.getDarkMap(mdh)

        if check_mapexists(mdh,'variance') is None:
            Warn(None, 'ncsimage: no variance map found')
            return
        else:
            # we need to convert to units of ADU^2
            variance = self.ci.getVarianceMap(mdh)/(mdh['Camera.ElectronsPerCount']**2)

        if check_mapexists(mdh,'flatfield') is None:
            gain = 1.0/mdh['Camera.ElectronsPerCount']*np.ones_like(variance)
        else:
            # code above argues we need to divide by flatfield
            gain = 1.0/(mdh['Camera.ElectronsPerCount']*self.ci.getFlatfieldMap(mdh))
        
        # the slice code needs a little bit of further checking
        isz = ncsSel.winSize
        imshape = img.shape
        xstart = self.do.xp - isz/2
        if xstart < 0:
            xstart = imshape[0]/2 - isz/2
        ystart = self.do.yp - isz/2
        if ystart < 0:
            ystart = imshape[1]/2 - isz/2

        sliceSquare = np.s_[xstart:xstart+isz,ystart:ystart+isz] # either on crosshairs or in the centre
        roi = [[xstart,xstart+isz],[ystart,ystart+isz],[self.do.zp,self.do.zp]]
        
        var_sl = variance[sliceSquare]
        dark_sl = dark[sliceSquare]
        gain_sl = gain[sliceSquare]
        
        # next few lines from NCSdemo_experiment which show how raw cmos data is pre-corrected
        #  apply gain and offset correction
        #  N = subims.shape[0]
        #  imsd = (subims-np.tile(suboffset,(N,1,1)))/np.tile(subgain,(N,1,1))
        #  imsd[imsd<=0] = 1e-6
        
        # therefore this needs to be in photoelectrons
        imgcorr = mdh['Camera.ElectronsPerCount']*self.ci.correctImage(mdh, img)
        imgc_sl = imgcorr[sliceSquare].squeeze()

        imgc_sl_T = imgc_sl[:,:,None].transpose((2,0,1))
        ret = np.clip(imgc_sl_T,1e-6,None,out=imgc_sl_T) # clip inplace at 1e-6

        Rs = ncsSel.Rs # IMPORTANT: looks like the pyncs code assumes Rs is a divider of the imgsz (see imagsz2 below)!!
        Pixelsize = mdh['voxelsize.x']
        Lambda = ncsSel.Lambda_nm / 1e3 # emission wavelength
        NA = ncsSel.NA     # objective NA
        iterationN = ncsSel.iterations
        alpha = ncsSel.alpha

        if ncs is not None:
            out = ncs.reducenoise(Rs,imgc_sl_T,var_sl,gain_sl,isz,Pixelsize,NA,Lambda,alpha,iterationN)
        else:
            out = imgc_sl_T # no reduction performed

        # now we need code to show this image and make it possible to save that
        disp_img = np.dstack([imgc_sl_T.squeeze(), out.squeeze()])
        im = ImageStack(disp_img, titleStub = '%d pixel crop of Frame %d denoised' % (isz,self.do.zp))
        im.mdh.copyEntriesFrom(mdh)

        # NCS parameters and crop info
        im.mdh['Parent'] = self.image.filename
        im.mdh['Processing.Units'] = 'PhotoElectrons'
        im.mdh['Processing.Type'] = 'NCS Denoising'
        im.mdh['Processing.NCS.alpha'] = ncsSel.alpha
        im.mdh['Processing.NCS.iterations'] = ncsSel.iterations
        im.mdh['Processing.NCS.Rs'] = ncsSel.Rs
        im.mdh['Processing.NCS.LambdaNM'] = ncsSel.Lambda_nm
        im.mdh['Processing.NCS.NA'] = ncsSel.NA
        im.mdh['Processing.CropROI'] = roi
        im.mdh['Processing.Comment'] = 'First frame: original (photoelectrons), second frame: denoised'
        
        vx, vy, vz = self.image.voxelsize
        ox, oy, oz = self.image.origin
        im.mdh['Origin.x'] = ox + roi[0][0]*vx
        im.mdh['Origin.y'] = oy + roi[1][0]*vy
        im.mdh['Origin.z'] = oz

        if self.dsviewer.mode == 'visGUI':
            mode = 'visGUI'
        else:
            mode = 'lite'

        dv = ViewIm3D(im, mode=mode, glCanvas=self.dsviewer.glCanvas, parent=wx.GetTopLevelParent(self.dsviewer))        
                            
            
    def showMap(self, type='dark'):
        mdh2 = NestedClassMDHandler(self.image.mdh)
        # overwrite the map location with default maps if exist
        if check_mapexists(mdh2,type=type) is None:
            Warn(None,'no suitable map in default location')
            return

        if type == 'dark':
            darkf = self.ci.getDarkMap(mdh2)
        elif type == 'variance':
            darkf = self.ci.getVarianceMap(mdh2)
        elif type == 'flatfield':
            darkf = self.ci.getFlatfieldMap(mdh2)
        else:
            raise RuntimeError('unknown type %s' % type)
            
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

        darkf = self.ci.getDarkMap(mdh2)
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
