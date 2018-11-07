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

        # the original formula was very adhoc
        # we have now changed to a different way
        # here we use an approach derived from a formula from Tang et al,
        # Scientific Reports | 5:11073 | DOi: 10.1038/srep11073
        mdh = self.pipeline.mdh
        # there is an issue if we don't have the nPhotons property FIXME!
        nph = self.pipeline['nPhotons']
        bgraw = self.pipeline['fitResults_background']
        bgph = np.clip((bgraw-mdh['Camera.ADOffset'])*mdh['Camera.ElectronsPerCount'],1,None)
        
        npixroi = (2*mdh.getOrDefault('Analysis.ROISize',5) + 1)**2
        snr = 1.0/npixroi * np.clip(nph,0,None)/np.sqrt(bgph)
        self.pipeline.addColumn('SNR', snr)

        self.pipeline.Rebuild()
        self.visFr.CreateFoldPanel() # to make, for example, new columns show up in filter column selections

            
def Plug(visFr):
    """Plugs this module into the gui"""
    visFr.snrCalc = SNRcalculator(visFr)
