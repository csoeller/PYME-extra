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
    
        visFr.AddMenuItem('Experimental>Deprecated>Fiducials', 'Add mean fiducial track', self.OnFiducialTrack,
                          helpText='Add mean fiducial track')
        visFr.AddMenuItem('Experimental>Deprecated>Fiducials', 'New DS with mean fiducial track applied',
                          self.OnFiducialCorrectDS,
                          helpText='Apply mean fiducial track')
        visFr.AddMenuItem('Experimental>Deprecated>Fiducials', "Plot Fiducial Track", self.OnPlotFiducial,
                          helpText='Plot mean fiducial tracks for all available dims')
        visFr.AddMenuItem('Experimental>Deprecated>Fiducials', "Set Z Parameters", self.OnSetZPars,
                          helpText='Set shift and scale parameters for driftz track')
        visFr.AddMenuItem('Experimental>Deprecated>Fiducials', "Set Z drift (from aligned driftz)", self.OnSetZDrift,
                          helpText='Set Z drift compensation from scaled and aligned driftz track')
        visFr.AddMenuItem('Experimental>Deprecated>Fiducials', "Clear Z driftz mapping", self.clearDriftZ,
                          helpText='Remove Z drift mapping by popping any mapping for z')
        visFr.AddMenuItem('Experimental>Deprecated>Fiducials', "Diagnose Fiducials", lambda e: self.fiducial_diagnosis(),
                          helpText='Diagnose quality of fiducial correction')
        visFr.AddMenuItem('Experimental>Deprecated>Fiducials', "Compare fiducial and drift", self.fiducialvsdrift,
                          helpText='Compare fiducial and drift information')
        visFr.AddMenuItem('Experimental>Corrections',"Fiducial - extract fiducial track and correct datasource",
                          self.OnFiducialCorrectNew,
                          helpText='faster tracking by inserting FiducialTrack and FiducialApplyFromFiducials modules')
        
        self.scaleFactor = -1e3
        self.shiftFrames = 0
        self.zeroAlignFrames = 200
        self.zDrift = None


    def OnFiducialCorrectNew(self, event=None):
        from PYMEcs.recipes import localisations
        recipe = self.pipeline.recipe
        # change defaults back to median filter
        ftrack = localisations.FiducialTrack(recipe, inputName='Fiducials',
                                             filterMethod='Median',
                                             outputName='fiducialAdded')
        if not ftrack.configure_traits(kind='modal'):
            return
        recipe.add_module(ftrack)
        recipe.add_module(localisations.FiducialApplyFromFiducials(recipe, inputData=self.pipeline.selectedDataSourceKey,
                                                                   inputFiducials='fiducialAdded',
                                                                   outputName='fiducialApplied',
                                                                   outputFiducials='corrected_fiducials'))
        recipe.execute()
        self.pipeline.selectDataSource('fiducialApplied')
            
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


    def fiducial_diagnosis(self):
        import numpy as np
        import matplotlib.pyplot as plt
        from PYME.Analysis.points.fiducials import FILTER_FUNCS
        pipeline = self.pipeline
        
        fids = pipeline.dataSources['corrected_fiducials']
        if 'clumpIndex' in fids.keys():
            ci = fids['clumpIndex']

            #cis = np.arange(1, ci.max())

            #clump_lengths = [(ci == c).sum() for c in cis]

            #largest_clump = ci == cis[np.argmax(clump_lengths)]

            #x_c = fids['x'][largest_clump]
            #y_c = fids['y'][largest_clump]

            f1 = plt.figure()
            a1 = plt.axes()
            plt.title('Y residuals')
            plt.grid()
            f2 = plt.figure()
            plt.title('X residuuals')
            a2 = plt.axes()
            plt.grid()
            f3 = plt.figure()
            plt.title('Z residuuals')
            a3 = plt.axes()
            plt.grid()

            sel_ids = np.unique(fids['fiducialID'])
            for pos in range(sel_ids.shape[0]):
                a1.text(-250, pos * 50, "%d" % sel_ids[pos])
                a2.text(-250, pos * 50, "%d" % sel_ids[pos])
                a3.text(-250, pos * 150, "%d" % sel_ids[pos])
            
            for i in range(1, ci.max()+1):
                mask = fids['clumpIndex'] == i
                if mask.sum() > 0:
                    f_id = fids['fiducialID'][mask][0]
                    pos = np.where(sel_ids==f_id)[0][0] # position in set of selected ids

                    fid_m = fids['fiducialID'] == f_id

                    ym = fids['y'][fid_m].mean()
                    xm = fids['x'][fid_m].mean()
                    zm = fids['z'][fid_m].mean()

                    # also plot a filtered version to see the trend in the noisy trace
                    yfilt = FILTER_FUNCS['Median'](fids['t'][mask],{'y':fids['y'][mask]},13)
                    xfilt = FILTER_FUNCS['Median'](fids['t'][mask],{'x':fids['x'][mask]},13)
                    zfilt = FILTER_FUNCS['Median'](fids['t'][mask],{'z':fids['z'][mask]},13)

                    a1.plot(fids['t'][mask], fids['y'][mask] - ym + pos * 50,
                            color=plt.cm.hsv( (i % 20.0)/20.))
                    a1.plot(fids['t'][mask], yfilt['y'] - ym + pos * 50, '--',
                            color='#b0b0b0', alpha=0.7)

                    a2.plot(fids['t'][mask], fids['x'][mask] - xm + pos * 50,
                            color=plt.cm.hsv((i % 20.0) / 20.))
                    a2.plot(fids['t'][mask], xfilt['x'] - xm + pos * 50, '--',
                            color='#b0b0b0', alpha=0.7)

                    a3.plot(fids['t'][mask], fids['z'][mask] - zm + pos * 150,
                            color=plt.cm.hsv((i % 20.0) / 20.))
                    a3.plot(fids['t'][mask], zfilt['z'] - zm + pos * 150, '--',
                            color='#b0b0b0', alpha=0.7)

        else:
            # stuff to do when no clumpIndex
            f1 = plt.figure()
            a1 = plt.axes()
            plt.title('Y residuals')
            plt.grid()
            f2 = plt.figure()
            plt.title('X residuuals')
            a2 = plt.axes()
            plt.grid()
            f3 = plt.figure()
            plt.title('Z residuuals')
            a3 = plt.axes()
            plt.grid()

            ym = fids['y'].mean()
            xm = fids['x'].mean()
            zm = fids['z'].mean()

            # also plot a filtered version to see the trend in the noisy trace
            yfilt = FILTER_FUNCS['Median'](fids['t'],{'y':fids['y']},13)
            xfilt = FILTER_FUNCS['Median'](fids['t'],{'x':fids['x']},13)
            zfilt = FILTER_FUNCS['Median'](fids['t'],{'z':fids['z']},13)

            a1.plot(fids['t'], fids['y'] - ym)
            a1.plot(fids['t'], yfilt['y'] - ym, '--',
                    color='#b0b0b0', alpha=0.7)
            
            a2.plot(fids['t'], fids['x'] - xm)
            a2.plot(fids['t'], xfilt['x'] - xm, '--',
                    color='#b0b0b0', alpha=0.7)
            
            a3.plot(fids['t'], fids['z'] - zm)
            a3.plot(fids['t'], zfilt['z'] - zm, '--',
                    color='#b0b0b0', alpha=0.7)

        # plot the trace derived from the fiducials
        tuq, idx = np.unique(fids['t'], return_index=True)
        fidz = fids['fiducial_z'][idx]
        fidy = fids['fiducial_y'][idx]
        fidx = fids['fiducial_x'][idx]
        plt.figure()
        plt.subplot(311)
        plt.plot(tuq, -fidz, label = 'fiducial z')
        plt.title('Fiducial z')
        plt.subplot(312)
        plt.plot(tuq, -fidx, label = 'fiducial x')
        plt.title('Fiducial x')
        plt.subplot(313)
        plt.plot(tuq, -fidy, label = 'fiducial y')
        plt.title('Fiducial y')


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
        plt.plot(tuq,fy-dfunc(p[2:],[dxx,dyy])+50,label='y difference (+50 nm offset)')
        plt.plot(tuq,fz-zfunc(pz,dz)+100,label='z difference (+100 nm offset)')
        plt.legend()
        plt.xlabel('Time (frames)')
        plt.ylabel('Drift (nm)')
        plt.grid()

def Plug(visFr):
    """Plugs this module into the gui"""
    FiducialTracker(visFr)
