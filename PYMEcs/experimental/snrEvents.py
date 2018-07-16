import numpy as np
# you may need a lot more imports depending what functionality your require in your plugin
class SNRcalculator:
    """
    A plugin, very simple to demonstrate the concept. Also providing a simple
    measure of some kind of SNR, the formula used is probably debatable.
    For example, low background estimates cause very high SNRs which may or may not
    be reasonable given the uncertainty in determining the background etc
    """
    def __init__(self, visFr):
        self.visFr = visFr
        self.pipeline = visFr.pipeline

        visFr.AddMenuItem('Experimental>ExtraColumns', 'Add SNR property', self.OnAddSNR,
                          helpText='Add an event property that provides some measure of SNR for events (from background and amplitude)')
        

    def OnAddSNR(self, event=None):
        """
        this function adds an 'SNR' property to events - there could be some discussion how that is actually best calculated
        """

        # the formula below is very adhoc
        # I am not even sure this is remotely correct nor the best way
        # so use only as a basis for experimentation and/or better plugins
        A = self.pipeline['A']
        bg = self.pipeline['fitResults_background']
        snr = np.sqrt(np.clip(A,0,None)/np.abs(bg))
        self.pipeline.addColumn('SNR', snr)

        self.pipeline.Rebuild()
        self.visFr.CreateFoldPanel() # to make, for example, new columns show up in filter column selections

            
def Plug(visFr):
    """Plugs this module into the gui"""
    visFr.snrCalc = SNRcalculator(visFr)
