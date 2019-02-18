class TimedSpecies:
    def __init__(self, visFr):
        self.visFr = visFr
        self.pipeline = visFr.pipeline
        
        visFr.AddMenuItem('Experimental', "&Sequential Imaging - Timed species assignment",self.OnTimedSpecies)

    def OnTimedSpecies(self,event):
        from PYMEcs.recipes import localisations
        recipe = self.pipeline.recipe
        tsm = localisations.TimedSpecies(recipe,
                                         inputName=self.pipeline.selectedDataSourceKey,
                                         outputName='timedSpecies')
        if not tsm.configure_traits(kind='modal'):
            return
        recipe.add_module(tsm)
        recipe.execute()

        self.pipeline.selectDataSource('timedSpecies')
        self.pipeline.Rebuild()
 
def Plug(visFr):
    '''Plugs this module into the gui'''
    visFr.timedSpecies = TimedSpecies(visFr)
