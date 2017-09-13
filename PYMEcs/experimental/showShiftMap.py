import numpy as np
import logging
logger = logging.getLogger(__file__)

class ShowMap:
    """

    """
    def __init__(self, visFr):
        self.visFr = visFr
        self.pipeline = visFr.pipeline

        visFr.AddMenuItem('Experimental>ShiftMap', 'Show Shiftmap', self.OnShowShiftMap,
                          helpText='Show a shiftmap from metadata info')

    def OnShowShiftMap(self, event=None):
        from PYME.Analysis.points import twoColour, twoColourPlot
        mdh = self.pipeline.mdh
        vs = [mdh['voxelsize.x']*1e3, mdh['voxelsize.y']*1e3, mdh['voxelsize.z']*1e3]
        dx = mdh.getEntry('chroma.dx')
        dy = mdh.getEntry('chroma.dy')

        shape = mdh['Splitter.Channel0ROI'][2:]

        twoColourPlot.PlotShiftField2(dx,dy,shape,vs)

def Plug(visFr):
    """Plugs this module into the gui"""
    visFr.showShiftMap = ShowMap(visFr)

