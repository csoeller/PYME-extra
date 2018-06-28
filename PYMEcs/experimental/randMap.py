import numpy as np

class RandMap:
    """

    """
    def __init__(self, visFr):
        self.visFr = visFr
        self.pipeline = visFr.pipeline

        visFr.AddMenuItem('Experimental',
                          'Add random value column',
                          self.OnRandMap,
                          helpText='the added random value column can be used to select a fraction of all events')

    def OnRandMap(self, event=None):
        self.pipeline.selectedDataSource.setMapping('randVal','0*x+np.random.rand(x.size)')
        self.pipeline.Rebuild()

def Plug(visFr):
    '''Plugs this module into the gui'''
    visFr.randMap = RandMap(visFr)
