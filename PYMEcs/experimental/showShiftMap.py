import numpy as np
import logging
logger = logging.getLogger(__file__)

# note: to turn this into a standalone shiftmap display tool, add code to open a shiftfield:
# import PYME.Acquire.Hardware.splitter as sp
# dx2,dy2 = sp.LoadShiftField(h5file_or_sf_file)
#
# if it is purely a shiftfield file then we have to supply some additional info, e.g. voxelsizes
# and ROI size for the splitter

class ShowMap:
    """

    """
    def __init__(self, visFr):
        self.visFr = visFr
        self.pipeline = visFr.pipeline

        visFr.AddMenuItem('Experimental>ShiftMap', 'Show Shiftmap', self.OnShowShiftMap,
                          helpText='Show a shiftmap from metadata info')
        visFr.AddMenuItem('Experimental>ShiftMap', 'Shiftmap as image', self.OnShiftMapAsImage,
                          helpText='Show a shiftmap from metadata info and display as image')

    def OnShowShiftMap(self, event=None):
        from PYME.Analysis.points import twoColour, twoColourPlot
        mdh = self.pipeline.mdh
        vs = [mdh['voxelsize.x']*1e3, mdh['voxelsize.y']*1e3, mdh['voxelsize.z']*1e3]
        dx = mdh.getEntry('chroma.dx')
        dy = mdh.getEntry('chroma.dy')

        shape = mdh['Splitter.Channel0ROI'][2:]

        twoColourPlot.PlotShiftField2(dx,dy,shape,vs)

    def OnShiftMapAsImage(self, event=None):
        mdh = self.pipeline.mdh
        voxelsize = [mdh['voxelsize.x']*1e3, mdh['voxelsize.y']*1e3, mdh['voxelsize.z']*1e3]
        spx = mdh.getEntry('chroma.dx')
        spy = mdh.getEntry('chroma.dy')

        shape = mdh['Splitter.Channel0ROI'][2:]
        xi, yi = np.meshgrid(np.arange(0, shape[0]*voxelsize[0], 100), np.arange(0, shape[1]*voxelsize[1], 100))
        xin = xi.ravel()
        yin = yi.ravel()
        dx = spx.ev(xin[:], yin[:]).reshape(xi.shape)
        dy = spy.ev(xin[:], yin[:]).reshape(xi.shape)
        
        from PYME.DSView.dsviewer import View3D
        image = np.stack([dx[::-1, :],dy[::-1, :]],axis=-1)
        View3D(image,  parent=self.visFr, titleStub='ShiftMap')

def Plug(visFr):
    """Plugs this module into the gui"""
    visFr.showShiftMap = ShowMap(visFr)

