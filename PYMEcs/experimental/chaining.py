import numpy as np

class ParticleTracker2:
    def __init__(self, visFr):
        self.visFr = visFr

        visFr.AddMenuItem('Experimental>Deprecated>Chaining', "Clump consecutive appearances", self.OnCoalesce)

    def OnCoalesce(self, event):
        from PYMEcs.recipes import localisations
        recipe = self.visFr.pipeline.recipe

        from PYMEcs.misc.utils import unique_name
        outputName = unique_name('coalesced',self.visFr.pipeline.dataSources.keys())
        recipe.add_module(localisations.MergeClumpsTperiod(recipe, inputName='with_clumps', outputName=outputName))
    
        recipe.execute()
        self.visFr.pipeline.selectDataSource(outputName)


def Plug(visFr):
    """Plugs this module into the gui"""
    ParticleTracker2(visFr)
