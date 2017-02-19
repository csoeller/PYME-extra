import numpy as np

# interpolate the key from the source to the selected datasource of the pipeline
def finterpDS(pipeline,sourcep,key):
    tsource, idx = np.unique(sourcep['t'], return_index=True)
    fsource = sourcep[key][idx]
    fDS = np.interp(pipeline.selectedDataSource['t'], tsource, fsource)
    return fDS

class FiducialTracker:
    """

    """
    def __init__(self, visFr):
        self.visFr = visFr
        self.pipeline = visFr.pipeline
    
        visFr.AddMenuItem('Experimental>Fiducials', 'Add mean fiducial track', self.OnFiducialTrack,
                          helpText='Add mean fiducial track')
        visFr.AddMenuItem('Experimental>Fiducials', 'New DS with mean fiducial track applied',
                          self.OnFiducialCorrectDS,
                          helpText='Apply mean fiducial track')
        visFr.AddMenuItem('Experimental>Fiducials', "Plot Fiducial Track", self.OnPlotFiducial,
                          helpText='Plot mean fiducial tracks for all available dims')

    def OnFiducialTrack(self, event=None):
        """

        """
        from PYMEcs.recipes import localisations

        fTracker = localisations.FiducialTrack(filterMethod = 'Gaussian')
        if fTracker.configure_traits(kind='modal'):
            # we call this with the pipeline to allow filtering etc
            namespace = {fTracker.inputName: self.pipeline}
            fTracker.execute(namespace)

            # the fiducial needs to be entered for the whole data source
            # otherwise we have an issue that fiducial data is not available
            # when filters are changed; this makes the code a bit ugly
            ds = namespace[fTracker.outputName]
            for fiducial in ['fiducial_%s' % dim for dim in ['x','y','z']]:
                if fiducial in ds.keys():
                    self.pipeline.selectedDataSource.addColumn(fiducial,
                                                               finterpDS(self.pipeline,
                                                                         ds,
                                                                         fiducial))

    def OnFiducialCorrectDS(self, event=None):
        """

        """
        from PYMEcs.recipes.localisations import FiducialApply
        recipe = self.pipeline.recipe
        recipe.add_module(FiducialApply(recipe, inputName=self.pipeline.selectedDataSourceKey,
                                        outputName='with_fiducial'))
        recipe.execute()
        self.pipeline.selectDataSource('with_fiducial')

        self.visFr.RefreshView()
        self.visFr.CreateFoldPanel()


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
