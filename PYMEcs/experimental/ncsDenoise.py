from PYME.DSView import ViewIm3D
from PYME.localization.remFitBuf import CameraInfoManager

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

from PYMEcs.misc.mapUtils import check_mapexists

class ncsSelect(HasTraits):
    winSize = Enum(64,128,256)
    Rs = Int(8) # IMPORTANT: looks like the pyncs code assumes Rs is a divider of the imgsz (see imagsz2 below)!!
    Lambda_nm = Float(690) # emission wavelength
    NA = Float(1.49)     # objective NA
    iterations = Int(15)
    alpha = Float(0.2)

class ncsDenoiser:
    """
GUI class to supply various map tools
    """
    def __init__(self, dsviewer):
        self.dsviewer = dsviewer
        self.do = dsviewer.do
        self.image = dsviewer.image
        self.ncsSel = ncsSelect() # by making it part of the object we retain parameters across invocations
        self.ci = CameraInfoManager()

        dsviewer.AddMenuItem('Experimental',
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
                            
def Plug(dsviewer):
    dsviewer.ncsDenoiser = ncsDenoiser(dsviewer)
