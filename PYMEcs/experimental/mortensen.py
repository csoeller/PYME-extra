import numpy as np
# you may need a lot more imports depending what functionality your require in your plugin
class MortensenFormula:
    """
    A plugin, very simple to demonstrate the concept. Also providing a simple
    measure of some kind of SNR, the formula used is probably debatable.
    For example, low background estimates cause very high SNRs which may or may not
    be reasonable given the uncertainty in determining the background etc
    """
    def __init__(self, visFr):
        self.visFr = visFr
        self.pipeline = visFr.pipeline

        visFr.AddMenuItem('Experimental>ExtraColumns', 'Add Mortensen Formula', self.OnAddMort,
                          helpText='Add an event property that provides an estimate by the Mortensen Formula (from background and amplitude)')
        visFr.AddMenuItem('Experimental>ExtraColumns', 'Plot Mortensen Error', self.OnPlotMort,
                          helpText='Scatterplot estimate by the Mortensen Formula')
        

    def OnAddMort(self, event=None):
        """
        this function adds a 'mortensenError' property to events - there could be some discussion how that is actually best calculated
        """
        import math
        mdh = self.pipeline.mdh
        # the formula below is very adhoc
        # I am not even sure this is remotely correct nor the best way
        # so use only as a basis for experimentation and/or better plugins
        N = self.pipeline['nPhotons']
        Nb = mdh['Camera.ElectronsPerCount'] * np.maximum(0,self.pipeline['fitResults_background']- mdh['Camera.ADOffset'])
        a = 1e3*mdh['voxelsize.x']
        siga = np.sqrt(self.pipeline['sig']**2+a*a/12.0)

        emort = siga /np.sqrt(N) * np.sqrt(16.0/9.0 + 8*math.pi*siga*siga*Nb/(N*a*a))
        emort_nobg = siga /np.sqrt(N) * np.sqrt(16.0/9.0)
        
        cb = 1.0*Nb/N
        self.pipeline.addColumn('cb_estimate', cb)
        self.pipeline.addColumn('mortensenError',emort)
        self.pipeline.addColumn('mortensenErrorNoBG',emort_nobg)
        self.pipeline.addColumn('backgroundPhotons',Nb)
        
        self.pipeline.Rebuild()
        self.visFr.CreateFoldPanel() # to make, for example, new columns show up in filter column selections


    def OnPlotMort(self, event=None):
        import matplotlib.pyplot as plt
        pipeline = self.pipeline
        plt.figure()
        ebg = plt.scatter(pipeline['error_x'],pipeline['mortensenError'],
                    c='g',alpha=0.5)
        enobg = plt.scatter(pipeline['error_x'],pipeline['mortensenErrorNoBG'],
                    c='r',alpha=0.5)
        plt.legend((ebg,enobg),('error with bg','error assuming zero bg'))
        plt.plot([0,20],[0,20])
        plt.xlabel('Fit error x')
        plt.ylabel('Error from Mortensen Formula')
        
def Plug(visFr):
    """Plugs this module into the gui"""
    visFr.mortForm = MortensenFormula(visFr)
