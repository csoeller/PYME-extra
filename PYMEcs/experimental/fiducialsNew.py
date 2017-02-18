import numpy as np

class FiducialTracker:
    """

    """
    def __init__(self, visFr):
        self.visFr = visFr
        self.pipeline = visFr.pipeline
    
        visFr.AddMenuItem('Experimental>Fiducials', 'Add mean fiducial track', self.OnFiducialTrack,
                          helpText='Add mean fiducial track')
        visFr.AddMenuItem('Experimental>Fiducials', "Plot Fiducial Track", self.OnPlotFiducial,
                          helpText='Plot mean fiducial tracks for all available dims')

    def OnFiducialTrack(self, event=None):
        """

        """
        from PYMEcs.recipes import localisations

        fTracker = localisations.FiducialTrack(filterMethod = 'Gaussian')
        if fTracker.configure_traits(kind='modal'):
            namespace = {fTracker.inputName: self.pipeline}
            fTracker.execute(namespace)

            for fiducial in ['fiducial_%s' % dim for dim in ['x','y','z']]:
                try:
                    self.pipeline.addColumn(fiducial, namespace[fTracker.outputName][fiducial])
                except:
                    pass

    def OnPlotFiducial(self, event):
        pipeline = self.visFr.pipeline
        t = pipeline['t']
        x = pipeline['fiducial_x']
        y = pipeline['fiducial_y']
        tu,idx = np.unique(t.astype('int'), return_index=True)
        xu = x[idx]
        yu = y[idx]
        import matplotlib.pyplot as plt
        plt.figure()
        plt.plot(tu, xu)
        plt.plot(tu, yu)

def Plug(visFr):
    """Plugs this module into the gui"""
    FiducialTracker(visFr)
