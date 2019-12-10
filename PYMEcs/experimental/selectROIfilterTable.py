import numpy as np
from PYME.recipes.tablefilters import FilterTable

class SelectROIFT:
    def __init__(self, visFr):
        self.visFr = visFr

        visFr.AddMenuItem('Experimental>View', "Add ROI FilterTable module from selection", self.OnSelROIFT)

    def OnSelROIFT(self, event):
        try:
            #old glcanvas
            x0, y0 = self.visFr.glCanvas.selectionStart
            x1, y1 = self.visFr.glCanvas.selectionFinish
        except AttributeError:
            #new glcanvas
            x0, y0 = self.visFr.glCanvas.selectionSettings.start
            x1, y1 = self.visFr.glCanvas.selectionSettings.finish

        filters = {}
        filters['x'] = [min(x0, x1), max(x0, x1)]
        filters['y'] = [min(y0, y1), max(y0,y1)]

        recipe = self.visFr.pipeline.recipe
        ftable = FilterTable(recipe, inputName=self.visFr.pipeline.selectedDataSourceKey,
                                  outputName='selectedROI', filters=filters)
        if not ftable.configure_traits(kind='modal'):
            return

        recipe.add_module(ftable)
        recipe.execute()

def Plug(visFr):
    """Plugs this module into the gui"""
    SelectROIFT(visFr)
