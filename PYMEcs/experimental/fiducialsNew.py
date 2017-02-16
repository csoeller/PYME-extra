class FiducialTracker:
    """

    """
    def __init__(self, visFr):
        self.visFr = visFr
        self.pipeline = visFr.pipeline
    
        visFr.AddMenuItem('Experimental>Fiducials', 'Add fiducial track', self.OnFiducialTrack,
                          helpText='Add fiducial track')

    def OnFiducialTrack(self, event=None):
        """

        """
        from PYMEcs.recipes import localisations

        fTracker = localisations.FiducialTrack()
        # set trait defaults for our specific application
        fTracker.filterMethod = 'Gaussian'
        if fTracker.configure_traits(kind='modal'):
            namespace = {fTracker.inputName: self.pipeline}
            fTracker.execute(namespace)

            for fiducial in ['fiducial_%s' % dim for dim in ['x','y','z']]:
                try:
                    self.pipeline.addColumn(fiducial, namespace[fTracker.outputName][fiducial])
                except:
                    pass


def Plug(visFr):
    """Plugs this module into the gui"""
    FiducialTracker(visFr)
