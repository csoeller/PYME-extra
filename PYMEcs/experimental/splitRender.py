from PYME.LMVis.renderers import TriangleRenderer, renderMetadataProviders
from PYME.IO.image import GeneratedImage, ImageBounds
# from PYME.LMVis import visHelpers
from PYME.LMVis import statusLog
# from PYME.IO import tabular

from PYME.IO import MetaDataHandler
import numpy as np

from traits.api import HasTraits, Str, Int, CStr, List, Enum, Float, Bool
from traitsui.api import View, Item, Group
from traitsui.menu import OKButton, CancelButton, OKCancelButtons

class TimeBlock(HasTraits):
    BlockSize = Int(100)   
    Colour = Enum(values='clist')
    clist = List([])

    traits_view = View(Group(Item(name = 'BlockSize'),
                             Item(name = 'Colour')),
                       title = 'Select',
                       buttons = OKCancelButtons
    )

    def add_keys(self,chans):
        for chan in chans:
            if chan not in self.clist:
                self.clist.append(chan)


class SplitRenderer(TriangleRenderer):
    """A renderer that specialises the TriangleRenderer to generate two images
       by splitting the events into two groups. It does not iterate over colour channels but rather will do 
       a single colour which must be selected prior to calling this renderer.
    """
    
    # def __init__(self, visFr, pipeline, mainWind = None, splitSettings = None):
    #     super(SplitRenderer, self).__init__(visFr, pipeline, mainWind)
    #     self.splitSettings = splitSettings
        
    def Generate(self, settings, splitSettings):
        mdh = MetaDataHandler.NestedClassMDHandler()
        mdh['Rendering.Method'] = self.name
        if 'imageID' in self.pipeline.mdh.getEntryNames():
            mdh['Rendering.SourceImageID'] = self.pipeline.mdh['imageID']
        mdh['Rendering.SourceFilename'] = getattr(self.pipeline, 'filename', '')
        mdh['Rendering.NEventsRendered'] = len(self.pipeline[self.pipeline.keys()[0]]) # in future good to use colourfilter for per channel info?
        mdh.Source = MetaDataHandler.NestedClassMDHandler(self.pipeline.mdh)

        if splitSettings is not None:
            for key in splitSettings:
                mdh['Rendering.SplitSettings.%s' % key] = splitSettings[key]
            
        for cb in renderMetadataProviders:
            cb(mdh)

        pixelSize = settings['pixelSize']

        status = statusLog.StatusLogger('Generating %s Image ...' % self.name)
        
        imb = self._getImBounds(*settings.get('zBounds', [None, None]))

        #record the pixel origin in nm from the corner of the camera for futrue overlays
        if 'Source.Camera.ROIPosX' in mdh.getEntryNames():
            #a rendered image with information about the source ROI
            voxx, voxy = 1e3 * mdh['Source.voxelsize.x'], 1e3 * mdh['Source.voxelsize.y']

            ox = (mdh['Source.Camera.ROIPosX'] - 1) * voxx + imb.x0
            oy = (mdh['Source.Camera.ROIPosY'] - 1) * voxy + imb.y0
            if 'Source.Positioning.PIFoc' in mdh.getEntryNames():
                oz = mdh['Source.Positioning.PIFoc'] * 1e3
            else:
                oz = imb.z0
        else:
            ox = imb.x0
            oy = imb.y0
            oz = imb.z0

        mdh['Origin.x'] = ox
        mdh['Origin.y'] = oy
        mdh['Origin.z'] = oz

        colours = settings['colours']
        curCol = self.colourFilter.currentColour
        if curCol is None:
            colnames = ['Everything_even','Everything_odd']
        else:
            colnames = ['%s_even' % curCol,'%s_odd' % curCol]

        try:
            oldBlock = self.pipeline.filterKeys['timeBlock']
        except KeyError, AttributeError:
            oldBlock = None
        
        ims = []

        for filter in [(-0.1,0.5),(0.5,1.5)]:
            self.pipeline.filterKeys['timeBlock'] = filter
            self.pipeline.Rebuild()
            ims.append(np.atleast_3d(self.genIm(settings, imb, mdh)))
        
        # needs to reset timeBlock in some fashion
        if oldBlock is None:
            del self.pipeline.filterKeys['timeBlock']
        else:
            self.pipeline.filterKeys['timeBlock'] = oldBlock
        self.pipeline.Rebuild()

        return GeneratedImage(ims, imb, pixelSize, settings['zSliceThickness'], colnames, mdh=mdh)

    def GenerateGUI(self, event=None, splitSettings=None):
        import wx
        from PYME.LMVis import genImageDialog
        from PYME.DSView import ViewIm3D
        
        jitVars = ['1.0']
        jitVars += self.colourFilter.keys()

        self.genMeas = self.pipeline.GeneratedMeasures.keys()
        if not 'neighbourDistances' in self.genMeas:
            self.genMeas.append('neighbourDistances')
            
        if not 'neighbourErrorMin' in self.genMeas:
            self.genMeas.append('neighbourErrorMin')
            
        jitVars += self.genMeas
        
        
        if 'z' in self.pipeline.keys():
            zvals = self.pipeline['z']
        else:
            zvals = None

        dlg = genImageDialog.GenImageDialog(self.mainWind, mode=self.mode, defaultPixelSize=self._defaultPixelSize, colours=self.colourFilter.getColourChans(), zvals = zvals, jitterVariables = jitVars, jitterVarDefault=self._getDefaultJitVar(jitVars), jitterVarDefaultZ=self._getDefaultZJitVar(jitVars))
        ret = dlg.ShowModal()

        #bCurr = wx.BusyCursor()

        if ret == wx.ID_OK:
            im = self.Generate(dlg.get_settings(),splitSettings=splitSettings)
            imfc = ViewIm3D(im, mode='visGUI', title='Generated %s - %3.1fnm bins' % (self.name, im.pixelSize),
                            glCanvas=self.visFr.glCanvas, parent=self.mainWind)
        else:
            imfc = None

        dlg.Destroy()
        return imfc

class splitRenderPlugin:
    def __init__(self, visFr):
        self.visFr = visFr
        self.pipeline = visFr.pipeline
        self.renderer = SplitRenderer(visFr, visFr.pipeline)
        self.blockSize = 100

        visFr.AddMenuItem('Experimental>Rendering',
                          'Split Render by Time Blocks',
                          self.OnSplitRender,
                          helpText='this renders a 2 channel image with events split by time blocks that can be used to evaluate the FRC; for multi-colour images the correct colour needs to be selected')

    def addTimeBlock(self):
        tb = TimeBlock(BlockSize=self.blockSize)
        tb.add_keys(['Everything']+self.pipeline.colourFilter.getColourChans())
        if tb.configure_traits(kind='modal'):
            psd = self.pipeline.selectedDataSource
            blockSel = np.mod((psd['t']/tb.BlockSize).astype('int'),2)
            self.blockSize = tb.BlockSize # remember choice
            self.colourChoice = tb.Colour
            self.pipeline.selectedDataSource.addColumn('timeBlock',blockSel)
            self.pipeline.Rebuild()
            return True
        else:
            return False

    def OnSplitRender(self,event=None):
        if not self.addTimeBlock():
            return
        oldCol = self.pipeline.colourFilter.currentColour
        if self.colourChoice == 'Everything':
            cChoice = None
        else:
            cChoice = self.colourChoice

        if oldCol != cChoice:
            self.pipeline.colourFilter.setColour(cChoice)
        self.renderer.GenerateGUI(splitSettings={'Mode':'TimeBlock','ChunkSize':self.blockSize})
        if oldCol != cChoice:
            self.pipeline.colourFilter.setColour(oldCol)
            self.pipeline.Rebuild()

def Plug(visFr):
    '''Plugs this module into the gui'''
    visFr.splitRender = splitRenderPlugin(visFr)
