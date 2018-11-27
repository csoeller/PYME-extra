import numpy as np
# you may need a lot more imports depending what functionality your require in your plugin

from traits.api import HasTraits, Str, Int, CStr, List, Enum, Float
#from traitsui.api import View, Item, Group
#from traitsui.menu import OKButton, CancelButton, OKCancelButtons

from PYMEcs.misc.guiMsgBoxes import Warn
import PYMEcs.misc.shellutils as su

class PlotOptions(HasTraits):
    plotMode = Enum(('Compare with and without background',
                     'Colour errors by photon number',
                     'Scatter Density Plot',
                     'Plot z errors (colour by photon number)'))

# We use the formula (S30) from Mortensen et al, 2010 [1] which provides a nice closed-form expression for the localisation error;
# note that this is likely a lower bound on actual error due to simplifying model assumptions (no read noise etc) 

# Reference

# 1. Optimized localization analysis for single-molecule tracking and super-resolution microscopy.
#    Kim I Mortensen, L Stirling Churchman, James A Spudich, and Henrik Flyvbjerg.
#    Nat Meth, 2017 vol. 18 (5) pp. 377-381.
#    http://www.nature.com/doifinder/10.1038/nmeth.1447

class MortensenFormula:
    """
    A plugin that calculates errors according to the Mortensen formula using photon number and background
    estimates that event analysis has provided.
    """
    def __init__(self, visFr):
        self.visFr = visFr
        self.pipeline = visFr.pipeline
        self.plotOptions = PlotOptions()

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
        # we think we do not need to subtract the camera offset
        Nb = mdh['Camera.ElectronsPerCount'] * np.maximum(0,self.pipeline['fitResults_background'] / mdh['Camera.TrueEMGain'])
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

        popt = self.plotOptions
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
            elif popt.plotMode == 'Plot z errors (colour by photon number)':
                if 'fitError_z0' not in pipeline.keys():
                    Warn('No z error -  works only with fitError_z0 property')
                    return
                nph = pipeline['nPhotons']
                nph5 = np.percentile(nph,5)
                nph95 = np.percentile(nph,95)
                plt.figure()
                plt.scatter(pipeline['fitError_z0'],pipeline['mortensenError'],
                            c=nph,vmin=nph5,vmax=nph95,cmap=plt.cm.jet)
                plt.plot([err1,err99],[err1,err99])
                plt.xlabel('Fit error z')
                plt.ylabel('Error from Mortensen Formula')
                plt.title('error coloured with nPhotons')
                plt.colorbar()
                
def Plug(visFr):
    """Plugs this module into the gui"""
    visFr.mortForm = MortensenFormula(visFr)
