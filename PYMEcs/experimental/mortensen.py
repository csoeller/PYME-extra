import numpy as np
# you may need a lot more imports depending what functionality your require in your plugin

from traits.api import HasTraits, Str, Int, CStr, List, Enum, Float
#from traitsui.api import View, Item, Group
#from traitsui.menu import OKButton, CancelButton, OKCancelButtons

import PYMEcs.misc.shellutils as su

class PlotOptions(HasTraits):
    plotMode = Enum(('Compare with and without background',
                     'Colour errors by photon number',
                     'Scatter Density Plot'))

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

        visFr.AddMenuItem('Experimental>ExtraColumns>Errors', 'Add Mortensen Formula', self.OnAddMort,
                          helpText='Add an event property that provides an estimate by the Mortensen Formula (from background and amplitude)')
        visFr.AddMenuItem('Experimental>ExtraColumns>Errors', 'Plot Mortensen Error', self.OnPlotMort,
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

        err = pipeline['mortensenError']
        err1 = np.percentile(err,1)
        err99 = np.percentile(err,99)

        errnbg = pipeline['mortensenErrorNoBG']
        errnbg1 = np.percentile(errnbg,1)
        errnbg99 = np.percentile(errnbg,99)

        popt = PlotOptions()
        if popt.configure_traits(kind='modal'):
            if popt.plotMode == 'Compare with and without background':
                plt.figure()
                ebg = plt.scatter(pipeline['error_x'],pipeline['mortensenError'],
                                  c='g',alpha=0.5)
                enobg = plt.scatter(pipeline['error_x'],pipeline['mortensenErrorNoBG'],
                                    c='r',alpha=0.5)
                plt.legend((ebg,enobg),('error with bg','error assuming zero bg'))
                plt.plot([errnbg1,err99],[errnbg1,err99])
                plt.xlabel('Fit error x')
                plt.ylabel('Error from Mortensen Formula')
            elif popt.plotMode == 'Colour errors by photon number':
                nph = pipeline['nPhotons']
                nph5 = np.percentile(nph,5)
                nph95 = np.percentile(nph,95)
                plt.figure()
                plt.scatter(pipeline['error_x'],pipeline['mortensenError'],
                            c=nph,vmin=nph5,vmax=nph95,cmap=plt.cm.jet)
                plt.plot([err1,err99],[err1,err99])
                plt.xlabel('Fit error x')
                plt.ylabel('Error from Mortensen Formula')
                plt.title('error coloured with nPhotons')
                plt.colorbar()
            elif popt.plotMode == 'Scatter Density Plot':
                plt.figure()
                su.scatterdens(pipeline['error_x'],pipeline['mortensenError'],
                               subsample=0.2, xlabel='Fit error x',
                               ylabel='Error from Mortensen Formula',s=20)
                plt.plot([err1,err99],[err1,err99])
                
def Plug(visFr):
    """Plugs this module into the gui"""
    visFr.mortForm = MortensenFormula(visFr)
