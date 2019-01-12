import numpy as np

import logging
logger = logging.getLogger(__file__)

# interpolate the key from the source to the selected datasource of the pipeline
def finterpDS(pipeline,sourcep,key):
    tsource, idx = np.unique(sourcep['t'], return_index=True)
    fsource = sourcep[key][idx]
    fDS = np.interp(pipeline.selectedDataSource['t'], tsource, fsource)
    return fDS

def zshift(t,data,navg=50):
    ti,idx = np.unique(t.astype('int'),return_index=True)
    di = data[idx]
    nm = min(navg,di.shape[0])
    offset = di[0:nm].mean()
    return data - offset

from traits.api import HasTraits, Str, Int, CStr, List, Enum, Float
#from traitsui.api import View, Item, Group
#from traitsui.menu import OKButton, CancelButton, OKCancelButtons

class SetZPars(HasTraits):
    scaleFactor = Float(-1e3)
    shiftFrames = Int(0)

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
        visFr.AddMenuItem('Experimental>Fiducials', "Set Z Parameters", self.OnSetZPars,
                          helpText='Set shift and scale parameters for driftz track')
        visFr.AddMenuItem('Experimental>Fiducials', "Set Z drift (from aligned driftz)", self.OnSetZDrift,
                          helpText='Set Z drift compensation from scaled and aligned driftz track')
        visFr.AddMenuItem('Experimental>Fiducials', "Clear Z driftz mapping", self.clearDriftZ,
                          helpText='Remove Z drift mapping by popping any mapping for z')
        
        self.scaleFactor = -1e3
        self.shiftFrames = 0
        self.zeroAlignFrames = 200
        self.zDrift = None

    def OnFiducialTrack(self, event=None):
        """

        """
        from PYMEcs.recipes import localisations

        if False:
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
                pds = self.pipeline.selectedDataSource
                isfid = np.zeros(len(pds['x']), dtype='i')
                isfid[self.pipeline.filter.Index] = ds['isFiducial']
                pds.addColumn('isFiducial',isfid)
        else:
            recipe = self.pipeline.recipe
            recipe.add_module(localisations.FiducialTrack(recipe, inputName=self.pipeline.selectedDataSourceKey,
                                                          outputName='with_fiducial'))
            recipe.execute()
            self.pipeline.selectDataSource('with_fiducial')

    def OnFiducialCorrectDS(self, event=None):
        """

        """
        from PYMEcs.recipes.localisations import FiducialApply
        recipe = self.pipeline.recipe
        recipe.add_module(FiducialApply(recipe, inputName='with_fiducial',
                                        outputName='corrected_from_fiducial'))
        recipe.execute()
        self.pipeline.selectDataSource('corrected_from_fiducial')

        #self.visFr.RefreshView()
        #self.visFr.CreateFoldPanel()


    def OnPlotFiducial(self, event):
        import PYMEnf.DriftCorrection.compactFit as cf
        
        pipeline = self.visFr.pipeline
        t = pipeline['t']
        x = pipeline['fiducial_x']
        y = pipeline['fiducial_y']
        z = pipeline['fiducial_z']
        
        tu,idx = np.unique(t.astype('int'), return_index=True)
        xu = x[idx]
        yu = y[idx]
        zu = z[idx]

        hasdp = True
        try:
            driftPane = self.visFr.driftPane
        except:
            hasdp = False

        import matplotlib.pyplot as plt
        plt.figure()
        plt.plot(tu, xu, label='x')
        plt.plot(tu, yu, label='y')
        plt.plot(tu, zu, label='z')
        if hasdp:
            if 'driftx' in driftPane.dp.driftExprX:
                indepVars = { 't': t, 'driftx': pipeline['driftx'], 'drifty': pipeline['drifty'] }
                dx,dy,tt = cf.xyDriftCurves(driftPane.dp.driftCorrFcn,driftPane.dp.driftCorrParams,indepVars,t)
                plt.plot(tt,-zshift(tt,dx), '--', label='x-drift')
                plt.plot(tt,-zshift(tt,dy), '--', label='y-drift')
                ti = np.arange(tt.min(),tt.max(),dtype=t.dtype)
                tu,iu = np.unique(t,return_index=True)
                dzi = zshift(ti,np.interp(ti,tu,pipeline['driftz'][iu]),navg=self.zeroAlignFrames)
                dzir = self.scaleFactor*np.roll(dzi,self.shiftFrames)
                self.zDrift = [ti, dzir]
                plt.plot(ti, dzir, '--', label='z-drift')
        plt.legend()


    def OnSetZPars(self, event=None):
        setPar = SetZPars(scaleFactor=self.scaleFactor,shiftFrames=self.shiftFrames,zeroAlignFrames = self.zeroAlignFrames)
        if setPar.configure_traits(kind='modal'):
            self.scaleFactor = setPar.scaleFactor
            self.shiftFrames = setPar.shiftFrames
            self.zeroAlignFrames = setPar.zeroAlignFrames


    def OnSetZDrift(self, event=None):
        if self.zDrift is None:
            logger.error('No zDrift found - cannot correct drift')
            return
        t, dz = self.zDrift
        self.visFr.pipeline.mapping.dz = np.interp(self.visFr.pipeline.mapping['t'], t, dz)
        self.visFr.pipeline.mapping.setMapping('z', 'z - dz')

        self.visFr.pipeline.ClearGenerated()

    def clearDriftZ(self, event=None):
        try:
            self.visFr.pipeline.mapping.mappings.pop('z')
        except KeyError:
            pass

        self.visFr.pipeline.ClearGenerated()
        
def Plug(visFr):
    """Plugs this module into the gui"""
    FiducialTracker(visFr)
