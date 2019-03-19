import numpy as np
import matplotlib.pyplot as plt

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
        visFr.AddMenuItem('Experimental>Fiducials', "Compare fiducial and drift", self.fiducialvsdrift,
                          helpText='Compare fiducial and drift information')
        
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


    def fiducialvsdrift(self, event=None):
        from PYMEcs.misc.guiMsgBoxes import Warn
        from scipy.optimize import leastsq
        import PYMEcs.misc.shellutils as su
        
        dfunc = lambda p, v: -100.0*p[0]*v[0]-100.0*p[1]*v[1]
        efunc = lambda p, fx, fy, dx, dy: np.append(fx-dfunc(p[:2],[dx,dy]),
                                                fy-dfunc(p[2:],[dx,dy]))
        zfunc = lambda p, dz: p[0]*dz + p[1]
        ezfunc = lambda p, fz, dz: fz-zfunc(p,dz)

        pipeline = self.pipeline
        if 'corrected_fiducials' not in pipeline.dataSources:
            Warn(self.visFr,"no 'corrected_fiducials' data source")
            return

        if 'driftx' not in pipeline.keys():
            Warn(self.visFr,"no 'driftx' property")
            return
        
        fids = pipeline.dataSources['corrected_fiducials']

        tuq, idx = np.unique(fids['t'], return_index=True)
        fidz = fids['fiducial_z'][idx]
        fidy = fids['fiducial_y'][idx]
        fidx = fids['fiducial_x'][idx]

        # what do we do when these do not exist?
        # answer: we may have to interpolate onto the times from the normal pipeline -> check that approach

        tup, idxp = np.unique(pipeline['t'], return_index=True)
        dxp = pipeline['driftx'][idxp]
        dyp = pipeline['drifty'][idxp]
        dzp = pipeline['driftz'][idxp]
    
        dx = np.interp(tuq, tup, dxp)
        dy = np.interp(tuq, tup, dyp)
        dz = np.interp(tuq, tup, dzp)
    
        #dy = fids['drifty'][idx]
        #dx = fids['driftx'][idx]

        fx = su.zs(fidx)
        fy = su.zs(fidy)
        fz = su.zs(fidz)
    
        dxx = su.zs(dx)
        dyy = su.zs(dy)
        p,suc = leastsq(efunc,np.zeros(4),args=(fx,fy,dxx,dyy))

        pz,sucz = leastsq(ezfunc,[-1e3,0],args=(fz,dz))

        plt.figure()
        plt.plot(tuq,fx,label='fiducial x')
        plt.plot(tuq,dfunc(p[:2],[dxx,dyy]),label='best fit x drift')

        plt.plot(tuq,fy,label='fiducial y')
        plt.plot(tuq,dfunc(p[2:],[dxx,dyy]),label='best fit y drift')
        plt.legend()
        plt.xlabel('Time (frames)')
        plt.ylabel('Drift (nm)')
        plt.title("Best fit params (a11 %.2f,a12 %.2f,a21 %.2f,a22 %.2f): " % tuple(p.tolist()))
        
        plt.figure()
        plt.plot(tuq,fz,label='fiducial z')
        plt.plot(tuq,zfunc(pz,dz),label='best fit z drift')
        plt.legend()
        plt.xlabel('Time (frames)')
        plt.ylabel('Drift (nm)')
        plt.title("Best fit parameters (zfactor %.2f, zoffs %.2f): " % tuple(pz.tolist()))

        plt.figure()
        plt.plot(tuq,fx-dfunc(p[:2],[dxx,dyy]),label='x difference')
        plt.plot(tuq,fy-dfunc(p[2:],[dxx,dyy]),label='y difference')
        plt.plot(tuq,fz-zfunc(pz,dz),label='z difference')
        plt.legend()
        plt.xlabel('Time (frames)')
        plt.ylabel('Drift (nm)')

def Plug(visFr):
    """Plugs this module into the gui"""
    FiducialTracker(visFr)
