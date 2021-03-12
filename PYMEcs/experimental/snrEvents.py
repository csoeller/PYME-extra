import numpy as np
from PYMEcs.misc.guiMsgBoxes import Warn
import wx
import os

def populate_fresults_wholep(fitMod,pipeline,bgzero=True):
    r = np.zeros(pipeline['x'].size, fitMod.fresultdtype)
    for k in pipeline.keys():
        f = k.split('_')
        if len(f) == 1:
            try:
                r[f[0]] = pipeline[k]
            except ValueError:
                pass
        elif len(f) == 2:
            try:
                r[f[0]][f[1]] = pipeline[k]
            except ValueError:
                pass
        elif len(f) == 3:
            try:
                r[f[0]][f[1]][f[2]] = pipeline[k]
            except ValueError:
                pass
        else:
            raise RuntimeError('more fields than expected: %d' % len(f))

    if bgzero:
        r['fitResults']['bg'] = 0
        r['fitResults']['br'] = 0
    
    return r


from PYME.IO.MetaDataHandler import NestedClassMDHandler
def genFitImage(fitMod,fr,mdh,psfname=None):
    mdh2 = NestedClassMDHandler(mdh)
    if psfname is not None:
        mdh2['PSFFile'] = psfname
    fitim = fitMod.genFitImage(fr,mdh2)

    return fitim

def get_photons(fitim,mdh):
    nph = fitim.sum()*mdh.getEntry('Camera.ElectronsPerCount')/mdh.getEntry('Camera.TrueEMGain')
    return nph


def nPhotons(fitMod,fr,mdh,psfname=None,nmax=100,progressBar=None,updateStep=100):
    mdh2 = NestedClassMDHandler(mdh)
    if psfname is not None:
        mdh2['PSFFile'] = psfname
    npoints = min(fr.shape[0],nmax)
    nph = np.zeros((npoints))
    us = int(updateStep)
    for i in range(npoints):
        nph[i] = get_photons(genFitImage(fitMod,fr[i],mdh2,psfname=None), mdh2)
        if (progressBar is not None) and ((i % us) == 0):
            progressBar.Update(100.0*i/float(npoints))
            wx.Yield()
    return nph


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
        visFr.AddMenuItem('Experimental>Biplane', 'Add Biplane nPhotons', self.OnAddnphBPlane,
                          helpText='Add nPhotons for Biplane fitting')
        visFr.AddMenuItem('Experimental>Biplane', 'Plot Biplane zErrors', self.OnPlotErrBPlane,
                          helpText='Plot Biplane fitting z-errors')


    def OnAddSNR(self, event=None):
        """
        this function adds an 'SNR' property to events - there could be some discussion how that is actually best calculated
        """
        from PYMEcs.recipes import localisations

        if False: # this way does not make a special recipe module but uses the generic Mapping module and adds a few columns
            from PYME.IO import tabular
            from PYME.recipes import tablefilters
            # the original formula was very adhoc
            # we have now changed to a different way
            # here we use an approach derived from a formula from Tang et al,
            # Scientific Reports | 5:11073 | DOi: 10.1038/srep11073
            mdh = self.pipeline.mdh
            pipeline = self.pipeline
            # there is an issue if we don't have the nPhotons property FIXME!
            nph = self.pipeline['nPhotons']
            bgraw = self.pipeline['fitResults_background']
            bgph = np.clip((bgraw)*mdh['Camera.ElectronsPerCount']/mdh.getEntry('Camera.TrueEMGain'),1,None)
        
            npixroi = (2*mdh.getOrDefault('Analysis.ROISize',5) + 1)**2
            snr = 1.0/npixroi * np.clip(nph,0,None)/np.sqrt(bgph)
            #dirty copy of the pipeline output

            recipe = pipeline.recipe
            mapp = tablefilters.Mapping(recipe,inputName=pipeline.selectedDataSourceKey,
                                  outputName='snr')
            recipe.add_module(mapp)
            recipe.execute()

            snrmap = recipe.namespace['snr']
            snrmap.addColumn('SNR', snr)
            snrmap.addColumn('backgroundPhotons',bgph)

        else:
            recipe = self.pipeline.recipe
            snr = localisations.SnrCalculation(recipe,inputName=self.pipeline.selectedDataSourceKey,
                                               outputName='snr')
            recipe.add_module(snr)
            recipe.execute()
            
        self.pipeline.selectDataSource('snr')

        self.pipeline.Rebuild()
        #self.visFr.CreateFoldPanel() # we should not need this anymore


    def OnAddnphBPlane(self, event=None):
        """
        this function adds nPhoton property to events for biplane
        """
        fdialog = wx.FileDialog(None, 'Please select PSF to use ...',
                                #defaultDir=os.path.split(self.image.filename)[0],
                                wildcard='PSF Files|*.psf|TIFF files|*.tif', style=wx.FD_OPEN)
        succ = fdialog.ShowModal()
        if (succ == wx.ID_OK):
            psfn = filename = fdialog.GetPath()
            mdh = self.pipeline.mdh
            pipeline = self.pipeline
            if mdh.getEntry('Analysis.FitModule') not in ['SplitterFitInterpBNR']:
                Warn('Plugin works only for Biplane analysis')
                return
            fitMod = __import__('PYME.localization.FitFactories.' +
                                self.pipeline.mdh.getEntry('Analysis.FitModule'),
                                fromlist=['PYME', 'localization', 'FitFactories'])
            fr = populate_fresults_wholep(fitMod, pipeline)
            progress = wx.ProgressDialog("calculating photon numbers",
                                         "calculating...", maximum=100, parent=None,
                                         style=wx.PD_SMOOTH|wx.PD_AUTO_HIDE)
            nph = nPhotons(fitMod, fr, mdh, psfname=psfn, nmax=1e6,
                           progressBar=progress, updateStep = 100)
            progress.Destroy()
            self.pipeline.addColumn('nPhotons', nph)
            self.pipeline.addColumn('fitResults_background', pipeline['fitResults_bg']+pipeline['fitResults_br'])
            self.pipeline.addColumn('sig',float(137.0)+np.zeros_like(pipeline['x'])) # this one is a straight kludge for mortensenError
            self.pipeline.Rebuild()
            self.visFr.CreateFoldPanel()


    def OnPlotErrBPlane(self, event=None):
        mdh = self.pipeline.mdh
        pipeline = self.pipeline
        if mdh.getEntry('Analysis.FitModule') not in ['SplitterFitInterpBNR']:
            Warn('Plugin works only for Biplane analysis')
            return
        bgraw = pipeline['fitResults_bg']+pipeline['fitResults_br']
        bgph = bgraw * pipeline.mdh.getEntry('Camera.ElectronsPerCount')/pipeline.mdh.getEntry('Camera.TrueEMGain')

        import matplotlib.pyplot as plt
        plt.figure()
        plt.scatter(pipeline['nPhotons'],pipeline['fitError_z0'],s=10,c=bgph)
        plt.colorbar()
        plt.xlim(0,None)
        plt.xlabel('Photons')
        plt.ylabel('z error (nm)')
        
# example use after adding the nPhotons entry:
# def bg(pipeline):
#     return (pipeline['fitResults_bg']+pipeline['fitResults_br'])*pipeline.mdh.getEntry('Camera.ElectronsPerCount')/pipeline.mdh.getEntry('Camera.TrueEMGain')
# scatter(pipeline['nPhotons'],pipeline['fitError_z0'],s=10,c=bg(pipeline))
# colorbar()            
            
def Plug(visFr):
    """Plugs this module into the gui"""
    visFr.snrCalc = SNRcalculator(visFr)
