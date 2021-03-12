import numpy as np
from PYMEcs.recipes.localisations import ClusterTimeRange, ValidClumps
from PYME.recipes.tablefilters import FilterTable
from PYME.recipes.localisations import DBSCANClustering

class SpecLabeling:
    def __init__(self, visFr):
        self.visFr = visFr

        visFr.AddMenuItem('Experimental>Corrections', "Filter specific labeling (DNA-PAINT)", self.OnSpecFilter)
        visFr.AddMenuItem('Experimental>Corrections', "Apply specific labeling to non-clumped data (DNA-PAINT)", self.OnSpecFilterNC)

    def OnSpecFilter(self, event):
        from PYMEcs.recipes import localisations
        recipe = self.visFr.pipeline.recipe
        bigLimit = 1e6 # size for filters big enough to capture everything
        
        # check 'coalesced' is available
        recipe.add_module(DBSCANClustering(recipe,inputName='coalesced', outputName='dbscanClustered', columns=['x', 'y'],
                                            searchRadius=50, minClumpSize=3, clumpColumnName='dbscanClumpID'))
        filters={'dbscanClumpID' : [0.5,1e6]}
        recipe.add_module(FilterTable(recipe, inputName='dbscanClustered',
                                  outputName='validCluster', filters={'dbscanClumpID' : [0.5,bigLimit]}))
        
        recipe.add_module(ClusterTimeRange(recipe, inputName='validCluster',
                                           outputName='withTrange', IDkey='dbscanClumpID'))
        recipe.add_module(FilterTable(recipe, inputName='withTrange',
                                  outputName='specLabeling', filters={'trange' : [2500,bigLimit]}))
        
        recipe.execute()
        self.visFr.pipeline.selectDataSource('specLabeling')

        
    def OnSpecFilterNC(self, event):
        from PYMEcs.recipes import localisations
        recipe = self.visFr.pipeline.recipe

        recipe.add_module(ValidClumps(inputName='with_clumps',
                                      inputValid='specLabeling',
                                      outputName='with_validClumps'))
        recipe.add_module(FilterTable(recipe, inputName='with_validClumps',
                                  outputName='specLabelingNC', filters={'validID' : [0.5,1.5]}))
        recipe.add_module(FilterTable(recipe, inputName='with_validClumps',
                                  outputName='nonSpecNC', filters={'validID' : [-0.5,0.5]}))
        recipe.execute()
        self.visFr.pipeline.selectDataSource('specLabelingNC')


def Plug(visFr):
    """Plugs this module into the gui"""
    SpecLabeling(visFr)
