from traitsui.api import View, Item, Group
from traitsui.menu import OKButton, CancelButton, OKCancelButtons
#from PYME.DSView.dsviewer import View3D
#from PYME.IO.MetaDataHandler import NestedClassMDHandler

#import numpy as np

from PYME.DSView.modules._base import Plugin
class SIMPLER():
    def __init__(self, visFr):
        self.visFr = visFr
             
        visFr.AddMenuItem('Experimental>Corrections', "Filter events for SIMPLER", self.OnFilterSIMPLER)
        visFr.AddMenuItem('Experimental>Corrections', "Cluster modes for SIMPLER N0 distribution", self.OnModesN0SIMPLER)


    def OnFilterSIMPLER(self, event):
        from PYME.recipes.tablefilters import FilterTable
        from PYME.recipes.tracking import FindClumps
        recipe = self.visFr.pipeline.recipe
        pipeline = self.visFr.pipeline
        recipe.add_modules_and_execute([FindClumps(recipe, inputName=pipeline.selectedDataSourceKey,
                                                   outputName='with_clumps',
                                                   timeWindow=0,
                                                   clumpRadiusVariable='error_x',
                                                   clumpRadiusScale=2.0),
                                        FilterTable(recipe, inputName='with_clumps',
                                                    outputName='simpler_filtered',
                                                    filters={'clumpEdge' : [-0.1,0.1]})])
        self.visFr.pipeline.selectDataSource('simpler_filtered')

    def OnModesN0SIMPLER(self, event):
        

def Plug(visFr):
    return SIMPLER(visFr)
