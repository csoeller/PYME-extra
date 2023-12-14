import numpy as np

class ParticleTracker2:
    def __init__(self, visFr):
        self.visFr = visFr

        visFr.AddMenuItem('Experimental>Deprecated>Chaining', "Clump consecutive appearances", self.OnCoalesce)

    def OnCoalesce(self, event):
        from PYMEcs.recipes import localisations
        recipe = self.visFr.pipeline.recipe
    
        recipe.add_module(localisations.MergeClumpsTperiod(recipe, inputName='with_clumps', outputName='coalesced'))
    
        recipe.execute()
        self.visFr.pipeline.selectDataSource('coalesced')
        #self.visFr.CreateFoldPanel() #TODO: can we capture this some other way?


def Plug(visFr):
    """Plugs this module into the gui"""
    ParticleTracker2(visFr)
