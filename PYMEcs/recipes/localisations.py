from PYME.recipes.base import register_module, ModuleBase, Filter
from PYME.recipes.traits import Input, Output, Float, Enum, CStr, Bool, Int, List, DictStrStr, DictStrList, ListFloat, ListStr

import numpy as np
import pandas as pd
from PYME.IO import tabular

import logging
logger = logging.getLogger(__file__)

@register_module('FiducialTrack')
class FiducialTrack(ModuleBase):
    """
    Extract average fiducial track from input pipeline

    Parameters
    ----------

        radiusMultiplier: this number is multiplied with error_x to obtain search radius for clustering
        timeWindow: the window along the time dimension used for clustering
        filterScale: the size of the filter kernel used to smooth the resulting average fiducial track
        filterMethod: enumrated choice of filter methods for smoothing operation (Gaussian, Median or Uniform kernel)

    Notes
    -----

    Output is a new pipeline with added fiducial_x, fiducial_y columns

    """
    import PYMEcs.Analysis.trackFiducials as tfs
    inputName = Input('filtered')
    
    radiusMultiplier = Float(5.0)
    timeWindow = Int(25)
    filterScale = Float(11)
    filterMethod = Enum(tfs.FILTER_FUNCS.keys())
    clumpMinSize = Int(50)
    singleFiducial = Bool(True)

    outputName = Output('fiducialAdded')

    def execute(self, namespace):
        import PYMEcs.Analysis.trackFiducials as tfs

        inp = namespace[self.inputName]
        mapped = tabular.mappingFilter(inp)

        if self.singleFiducial:
            # if all data is from a single fiducial we do not need to align
            # we then avoid problems with incomplete tracks giving rise to offsets between
            # fiducial track fragments
            align = False
        else:
            align = True

        t, x, y, z, isFiducial = tfs.extractTrajectoriesClump(inp, clumpRadiusVar = 'error_x',
                                                              clumpRadiusMultiplier=self.radiusMultiplier,
                                                              timeWindow=self.timeWindow, clumpMinSize=self.clumpMinSize,
                                                              align=align)
        rawtracks = (t, x, y, z)
        tracks = tfs.AverageTrack(inp, rawtracks, filter=self.filterMethod,
                                         filterScale=self.filterScale,align=align)
        
        # add tracks for all calculated dims to output
        for dim in tracks.keys():
            mapped.addColumn('fiducial_%s' % dim, tracks[dim])
        mapped.addColumn('isFiducial', isFiducial)
        
        # propogate metadata, if present
        try:
            mapped.mdh = inp.mdh
        except AttributeError:
            pass

        namespace[self.outputName] = mapped

    @property
    def hide_in_overview(self):
        return ['columns']

@register_module('FiducialApply')
class FiducialApply(ModuleBase):
    inputName = Input('filtered')
    outputName = Output('fiducialApplied')

    def execute(self, namespace):
        inp = namespace[self.inputName]
        mapped = tabular.mappingFilter(inp)

        for dim in ['x','y','z']:
            try:
                mapped.addColumn(dim, inp[dim]-inp['fiducial_%s' % dim])
            except:
                logger.warn('Could not set dim %s' % dim)

        # propogate metadata, if present
        try:
            mapped.mdh = inp.mdh
        except AttributeError:
            pass
        
        namespace[self.outputName] = mapped


@register_module('MergeClumpsTperiod')
class MergeClumpsTperiod(ModuleBase):
    """
    Create a new mapping object which derives mapped keys from original ones.
    Also adds the time period of bursts by adding tmin and tmax columns for each clump.
    """
    inputName = Input('clumped')
    outputName = Output('merged')
    labelKey = CStr('clumpIndex')

    def execute(self, namespace):
        from PYME.Analysis.points.DeClump import pyDeClump
        from PYME.Analysis.points.DeClump import deClump as deClumpC
        
        inp = namespace[self.inputName]

        grouped = pyDeClump.mergeClumps(inp, labelKey=self.labelKey)

        # work out tmin and tmax
        I = np.argsort(inp[self.labelKey])
        sorted_src = {k: inp[k][I] for k in [self.labelKey,'t']}
        # tmin and tmax - tentative addition
        NClumps = int(np.max(sorted_src[self.labelKey]) + 1)
        tmin = deClumpC.aggregateMin(NClumps, sorted_src[self.labelKey].astype('i'), sorted_src['t'].astype('f'))
        tmax = -deClumpC.aggregateMin(NClumps, sorted_src[self.labelKey].astype('i'), -1.0*sorted_src['t'].astype('f'))
        grouped.addColumn('tmin',tmin)
        grouped.addColumn('tmax',tmax)
        
        try:
            grouped.mdh = inp.mdh
        except AttributeError:
            pass

        namespace[self.outputName] = grouped

# interpolate the key from the source to the selected target of the pipeline
def finterpDS(target,source,key):
    tsource, idx = np.unique(source['t'], return_index=True)
    fsource = source[key][idx]
    fDS = np.interp(target['t'], tsource, fsource)
    return fDS

@register_module('FiducialApplyFromFiducials')
class FiducialApplyFromFiducials(ModuleBase):
    inputData = Input('filtered')
    inputFiducials = Input('Fiducials')
    outputName = Output('fiducialApplied')
    outputFiducials = Output('corrected_fiducials')
    
    def execute(self, namespace):
        inp = namespace[self.inputData]
        fiducial = namespace[self.inputFiducials]
        
        mapped = tabular.mappingFilter(inp)
        out_f = tabular.mappingFilter(fiducial)
        
        for dim in ['x','y','z']:
            fiducial_dim = finterpDS(inp,fiducial,'fiducial_%s' % dim)
            
            mapped.addColumn('fiducial_%s' % dim,fiducial_dim)
            mapped.addColumn(dim, inp[dim]-fiducial_dim)

            out_f.setMapping(dim, '{0} - fiducial_{0}'.format(dim))

        # propogate metadata, if present
        try:
            mapped.mdh = inp.mdh
            out_f.mdh = fiducial.mdh
        except AttributeError:
            pass

        namespace[self.outputName] = mapped
        namespace[self.outputFiducials] = out_f


@register_module('ClusterTimeRange')
class ClusterTimeRange(ModuleBase):
    
    inputName = Input('dbscanClustered')
    IDkey = CStr('dbscanClumpID')
    outputName = Output('withTrange')

    def execute(self, namespace):
        from scipy.stats import binned_statistic
        
        inp = namespace[self.inputName]
        mapped = tabular.mappingFilter(inp)

        ids = inp[self.IDkey]
        t = inp['t']
        maxid = int(ids.max())
        edges = -0.5+np.arange(maxid+2)
        resmin = binned_statistic(ids, t, statistic='min', bins=edges)
        resmax = binned_statistic(ids, t, statistic='max', bins=edges)
        trange = resmax[0][ids] - resmin[0][ids] + 1

        mapped.addColumn('trange', trange)
        
        # propogate metadata, if present
        try:
            mapped.mdh = inp.mdh
        except AttributeError:
            pass
        
        namespace[self.outputName] = mapped

@register_module('ClusterStats')
class ClusterStats(ModuleBase):
    
    inputName = Input('with_clumps')
    IDkey = CStr('clumpIndex')
    StatMethod = Enum(['std','min','max', 'mean', 'median', 'count', 'sum'])
    StatKey = CStr('x')
    outputName = Output('withClumpStats')

    def execute(self, namespace):
        from scipy.stats import binned_statistic
        
        inp = namespace[self.inputName]
        mapped = tabular.mappingFilter(inp)

        ids = inp[self.IDkey] # I imagine this needs to be an int type key
        prop = inp[self.StatKey]
        maxid = int(ids.max())
        edges = -0.5+np.arange(maxid+2)
        resstat = binned_statistic(ids, prop, statistic=self.StatMethod, bins=edges)

        mapped.addColumn(self.StatKey+"_"+self.StatMethod, resstat[0][ids])
        
        # propogate metadata, if present
        try:
            mapped.mdh = inp.mdh
        except AttributeError:
            pass
        
        namespace[self.outputName] = mapped

    @property
    def _key_choices(self):
        #try and find the available column names
        try:
            return sorted(self._parent.namespace[self.inputName].keys())
        except:
            return []

    @property
    def default_view(self):
        from traitsui.api import View, Group, Item
        from PYME.ui.custom_traits_editors import CBEditor

        return View(Item('inputName', editor=CBEditor(choices=self._namespace_keys)),
                    Item('_'),
                    Item('IDkey', editor=CBEditor(choices=self._key_choices)),
                    Item('StatKey', editor=CBEditor(choices=self._key_choices)),
                    Item('StatMethod'),
                    Item('_'),
                    Item('outputName'), buttons=['OK'])

@register_module('ValidClumps')
class ValidClumps(ModuleBase):
    
    inputName = Input('with_clumps')
    inputValid = Input('valid_clumps')
    IDkey = CStr('clumpIndex')
    outputName = Output('with_validClumps')

    def execute(self, namespace):
        
        inp = namespace[self.inputName]
        valid = namespace[self.inputValid]
        mapped = tabular.mappingFilter(inp)

        # note: in coalesced data the clumpIndices are float!
        # this creates issues in comparisons unless these are converted to int before comparisons are made!!
        # that is the reason for the rint and astype conversions below
        ids = np.rint(inp[self.IDkey]).astype('i')
        validIDs = np.in1d(ids,np.unique(np.rint(valid[self.IDkey]).astype('i')))
        
        mapped.addColumn('validID', validIDs.astype('f')) # should be float or int?
        
        # propogate metadata, if present
        try:
            mapped.mdh = inp.mdh
        except AttributeError:
            pass
        
        namespace[self.outputName] = mapped

@register_module('CopyMapped')
class CopyMapped(ModuleBase):
    inputName = Input('filtered')
    outputName = Output('filtered-copy')

    def execute(self, namespace):
        inp = namespace[self.inputName]
        mapped = tabular.mappingFilter(inp)
        namespace[self.outputName] = mapped

@register_module('QindexScale')
class QindexScale(ModuleBase):
    inputName = Input('qindex')
    outputName = Output('qindex-calibrated')
    qIndexkey = CStr('qIndex')
    qindexValue = Float(1.0)
    NEquivalent = Float(1.0)

    def execute(self, namespace):
        inp = namespace[self.inputName]
        mapped = tabular.mappingFilter(inp)
        qkey = self.qIndexkey
        scaled = inp[qkey]
        qigood = inp[qkey] > 0
        scaled[qigood] = inp[qkey][qigood] * self.NEquivalent / self.qindexValue

        self.newKey = '%sCal' % qkey
        mapped.addColumn(self.newKey, scaled)
        namespace[self.outputName] = mapped

    @property
    def _key_choices(self):
        #try and find the available column names
        try:
            return sorted(self._parent.namespace[self.inputName].keys())
        except:
            return []

    @property
    def default_view(self):
        from traitsui.api import View, Group, Item
        from PYME.ui.custom_traits_editors import CBEditor

        return View(Item('inputName', editor=CBEditor(choices=self._namespace_keys)),
                    Item('_'),
                    Item('qIndexkey', editor=CBEditor(choices=self._key_choices)),
                    Item('qindexValue'),
                    Item('NEquivalent'),
                    Item('_'),
                    Item('outputName'), buttons=['OK'])

@register_module('QindexRatio')
class QindexRatio(ModuleBase):
    inputName = Input('qindex')
    outputName = Output('qindex-calibrated')
    qIndexDenom = CStr('qIndex1')
    qIndexNumer = CStr('qIndex2')
    qIndexRatio = CStr('qRatio')

    def execute(self, namespace):
        inp = namespace[self.inputName]
        mapped = tabular.mappingFilter(inp)
        qkey1 = self.qIndexDenom
        qkey2 = self.qIndexNumer

        v2 = inp[qkey2]
        ratio = np.zeros_like(v2,dtype='float64')
        qigood = v2 > 0
        ratio[qigood] = inp[qkey1][qigood] / v2[qigood]

        mapped.addColumn(self.qIndexRatio, ratio)
        namespace[self.outputName] = mapped

    @property
    def _key_choices(self):
        #try and find the available column names
        try:
            return sorted(self._parent.namespace[self.inputName].keys())
        except:
            return []

    @property
    def default_view(self):
        from traitsui.api import View, Group, Item
        from PYME.ui.custom_traits_editors import CBEditor

        return View(Item('inputName', editor=CBEditor(choices=self._namespace_keys)),
                    Item('_'),
                    Item('qIndexDenom', editor=CBEditor(choices=self._key_choices)),
                    Item('qIndexNumer', editor=CBEditor(choices=self._key_choices)),
                    Item('qIndexRatio'),
                    Item('_'),
                    Item('outputName'), buttons=['OK'])


@register_module('ObjectVolume')
class ObjectVolume(ModuleBase):
    inputName = Input('objectID')
    outputName = Output('volumes')

    def execute(self, namespace):
        from PYMEcs.Analysis.objectVolumes import objectVolumes
        inp = namespace[self.inputName]
        mapped = tabular.mappingFilter(inp)

        volumes = objectVolumes(np.vstack([inp[k] for k in ('x','y')]).T,inp['objectID'])

        mapped.addColumn('volumes', volumes)
        namespace[self.outputName] = mapped

def uniqueByID(ids,column):
    uids, idx = np.unique(ids.astype('int'), return_index=True)
    ucol = column[idx]
    valid = uids > 0
    return uids[valid], ucol[valid]

@register_module('ScatterbyID')         
class ScatterbyID(ModuleBase):
    """Take just certain columns of a variable"""
    inputName = Input('measurements')
    IDkey = CStr('objectID')
    xkey = CStr('qIndex')
    ykey = CStr('objArea')
    outputName = Output('outGraph')
    
    def execute(self, namespace):
        meas = namespace[self.inputName]
        ids = meas[self.IDkey]
        uid, x = uniqueByID(ids,meas[self.xkey])
        uid, y = uniqueByID(ids,meas[self.ykey])
        
        import pylab
        pylab.figure()
        pylab.scatter(x,y)
        
        pylab.grid()
        pylab.xlabel(self.xkey)
        pylab.ylabel(self.ykey)
        #namespace[self.outputName] = out

    @property
    def _key_choices(self):
        #try and find the available column names
        try:
            return sorted(self._parent.namespace[self.inputName].keys())
        except:
            return []

    @property
    def default_view(self):
        from traitsui.api import View, Group, Item
        from PYME.ui.custom_traits_editors import CBEditor

        return View(Item('inputName', editor=CBEditor(choices=self._namespace_keys)),
                    Item('_'),
                    Item('IDkey', editor=CBEditor(choices=self._key_choices)),
                    Item('xkey', editor=CBEditor(choices=self._key_choices)),
                    Item('ykey', editor=CBEditor(choices=self._key_choices)),
                    Item('_'),
                    Item('outputName'), buttons=['OK'])

@register_module('HistByID')         
class HistByID(ModuleBase):
    """Plot histogram of a column by ID"""
    inputName = Input('measurements')
    IDkey = CStr('objectID')
    histkey = CStr('qIndex')
    outputName = Output('outGraph')
    nbins = Int(50)
    minval = Float(float('nan'))
    maxval = Float(float('nan'))
    
    def execute(self, namespace):
        import math
        meas = namespace[self.inputName]
        ids = meas[self.IDkey]
            
        uid, valsu = uniqueByID(ids,meas[self.histkey])
        if math.isnan(self.minval):
            minv = valsu.min()
        else:
            minv = self.minval
        if math.isnan(self.maxval):
            maxv = valsu.max()
        else:
            maxv = self.maxval

        import matplotlib.pyplot as plt
        plt.figure()
        plt.hist(valsu,self.nbins,range=(minv,maxv))
        plt.xlabel(self.histkey)

    @property
    def _key_choices(self):
        #try and find the available column names
        try:
            return sorted(self._parent.namespace[self.inputName].keys())
        except:
            return []

    @property
    def default_view(self):
        from traitsui.api import View, Group, Item
        from PYME.ui.custom_traits_editors import CBEditor

        return View(Item('inputName', editor=CBEditor(choices=self._namespace_keys)),
                    Item('_'),
                    Item('IDkey', editor=CBEditor(choices=self._key_choices)),
                    Item('histkey', editor=CBEditor(choices=self._key_choices)),
                    Item('nbins'),
                    Item('minval'),
                    Item('maxval'), 
                    Item('_'),
                    Item('outputName'), buttons=['OK'])

    # a version of David's module which we include here so that we can test/hack a few things
@register_module('DBSCANClustering2')
class DBSCANClustering2(ModuleBase):
    """
    Performs DBSCAN clustering on input dictionary

    Parameters
    ----------

        searchRadius: search radius for clustering
        minPtsForCore: number of points within SearchRadius required for a given point to be considered a core point

    Notes
    -----

    See `sklearn.cluster.dbscan` for more details about the underlying algorithm and parameter meanings.

    """
    import multiprocessing
    inputName = Input('filtered')

    columns = ListStr(['x', 'y', 'z'])
    searchRadius = Float()
    minClumpSize = Int()
    numberOfJobs = Int(max(multiprocessing.cpu_count()-1,1)) # this is a feature of the latest dbscan in scipy
    
    outputName = Output('dbscanClustered')

    def execute(self, namespace):
        from sklearn.cluster import dbscan

        inp = namespace[self.inputName]
        mapped = tabular.mappingFilter(inp)

        # Note that sklearn gives unclustered points label of -1, and first value starts at 0.
        try:
            core_samp, dbLabels = dbscan(np.vstack([inp[k] for k in self.columns]).T,
                                         self.searchRadius, self.minClumpSize, n_jobs=self.numberOfJobs)
            multiproc = True
        except:
            core_samp, dbLabels = dbscan(np.vstack([inp[k] for k in self.columns]).T,
                                         self.searchRadius, self.minClumpSize)
            multiproc = False

        if multiproc:
            logger.info('using dbscan multiproc version')
        else:
            logger.info('falling back to dbscan single-threaded version')

            # shift dbscan labels up by one to match existing convention that a clumpID of 0 corresponds to unclumped
        mapped.addColumn('dbscanClumpID', dbLabels + 1)

        # propogate metadata, if present
        try:
            mapped.mdh = inp.mdh
        except AttributeError:
            pass

        namespace[self.outputName] = mapped

    @property
    def hide_in_overview(self):
        return ['columns']

@register_module('SnrCalculation')
class SnrCalculation(ModuleBase):
    inputName = Input('filtered')
    outputName = Output('snr')
   
    def execute(self, namespace):
        inp = namespace[self.inputName]
        mapped = tabular.mappingFilter(inp)

        if 'mdh' not in dir(inp):
            raise RuntimeError('SnrCalculation needs metadata')
        else:
            mdh = inp.mdh

        nph = inp['nPhotons']
        bgraw = inp['fitResults_background']
        bgph = np.clip((bgraw)*mdh['Camera.ElectronsPerCount']/mdh.getEntry('Camera.TrueEMGain'),1,None)
        
        npixroi = (2*mdh.getOrDefault('Analysis.ROISize',5) + 1)**2
        snr = 1.0/npixroi * np.clip(nph,0,None)/np.sqrt(bgph)

        mapped.addColumn('SNR', snr)
        mapped.addColumn('backgroundPhotons',bgph)
        
        mapped.mdh = inp.mdh

        namespace[self.outputName] = mapped

@register_module('TimedSpecies')
class TimedSpecies(ModuleBase):
    inputName = Input('filtered')
    outputName = Output('timedSpecies')
    Species_1_Name = CStr('Species1')
    Species_1_Start = Float(0)
    Species_1_Stop = Float(1e6)
    
    Species_2_Name = CStr('')
    Species_2_Start = Float(0)
    Species_2_Stop = Float(0)
    
    Species_3_Name = CStr('')
    Species_3_Start = Float(0)
    Species_3_Stop = Float(0)
    
    def execute(self, namespace):
        inp = namespace[self.inputName]
        mapped = tabular.mappingFilter(inp)
        timedSpecies = self.populateTimedSpecies()
        
        mapped.addColumn('ColourNorm', np.ones_like(mapped['t'],'float'))
        for species in timedSpecies:
            mapped.addColumn('p_%s' % species['name'],
                             (mapped['t'] >= species['t_start'])*
                             (mapped['t'] < species['t_end']))

        if 'mdh' in dir(inp):
            mapped.mdh = inp.mdh
            mapped.mdh['TimedSpecies'] = timedSpecies

        namespace[self.outputName] = mapped

    def populateTimedSpecies(self):
        ts = []
        if self.Species_1_Name:
            ts.append({'name'   : self.Species_1_Name,
                       't_start': self.Species_1_Start,
                       't_end'  : self.Species_1_Stop})

        if self.Species_2_Name:
            ts.append({'name'   : self.Species_2_Name,
                       't_start': self.Species_2_Start,
                       't_end'  : self.Species_2_Stop})

        if self.Species_3_Name:
            ts.append({'name'   : self.Species_3_Name,
                       't_start': self.Species_3_Start,
                       't_end'  : self.Species_3_Stop})

        return ts

import wx
def populate_fresults(fitMod,inp,bgzero=True):
    r = np.zeros(inp['x'].size, fitMod.fresultdtype)
    for k in inp.keys():
        f = k.split('_')
        if len(f) == 1:
            try:
                r[f[0]] = inp[k]
            except ValueError:
                pass
        elif len(f) == 2:
            try:
                r[f[0]][f[1]] = inp[k]
            except ValueError:
                pass
        elif len(f) == 3:
            try:
                r[f[0]][f[1]][f[2]] = inp[k]
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


@register_module('BiplanePhotons')
class BiplanePhotons(ModuleBase):

    inputName = Input('filtered')
    outputName = Output('withPhotons')

    def execute(self, namespace):
        inp = namespace[self.inputName]
        mapped = tabular.mappingFilter(inp)

        fdialog = wx.FileDialog(None, 'Please select PSF to use ...',
                                #defaultDir=os.path.split(self.image.filename)[0],
                                wildcard='PSF Files|*.psf|TIFF files|*.tif', style=wx.FD_OPEN)
        succ = fdialog.ShowModal()
        if (succ == wx.ID_OK):
            psfn = filename = fdialog.GetPath()
           
            mdh = inp.mdh
            if mdh.getEntry('Analysis.FitModule') not in ['SplitterFitInterpBNR']:
                Warn('Plugin works only for Biplane analysis')
                return
            fitMod = __import__('PYME.localization.FitFactories.' +
                                mdh.getEntry('Analysis.FitModule'),
                                fromlist=['PYME', 'localization', 'FitFactories'])

            fr = populate_fresults(fitMod, inp)
            progress = wx.ProgressDialog("calculating photon numbers",
                                         "calculating...", maximum=100, parent=None,
                                         style=wx.PD_SMOOTH|wx.PD_AUTO_HIDE)
            nph = nPhotons(fitMod, fr, mdh, psfname=psfn, nmax=1e6,
                           progressBar=progress, updateStep = 100)
            progress.Destroy()
            mapped.addColumn('nPhotons', nph)
            mapped.addColumn('fitResults_background', inp['fitResults_bg']+inp['fitResults_br'])
            #mapped.addColumn('sig',float(137.0)+np.zeros_like(inp['x'])) # this one is a straight kludge for mortensenError


        # propogate metadata, if present
        try:
            mapped.mdh = inp.mdh
        except AttributeError:
            pass
        
        namespace[self.outputName] = mapped

@register_module('ClusterModes')
class ClusterModes(ModuleBase):
    
    inputName = Input('dbscanClustered')
    IDkey = CStr('dbscanClumpID')
    outputName = Output('with_clusterModes')
    PropertyKey = CStr('nPhotons')

    def execute(self, namespace):
        from PYMEcs.Analysis.Simpler import clusterModes
        
        inp = namespace[self.inputName]
        cmodes = tabular.mappingFilter(inp)

        ids = inp[self.IDkey] # I imagine this needs to be an int type key
        props = inp[self.PropertyKey]

        cm, ce, ccx, ccy = clusterModes(inp['x'],inp['y'],ids,props)
        cmodes.addColumn('clusterMode',cm)
        cmodes.addColumn('clusterModeError',ce)
        cmodes.addColumn('clusterCentroid_x',ccx)
        cmodes.addColumn('clusterCentroid_y',ccy)
        
        # propogate metadata, if present
        try:
            cmodes.mdh = inp.mdh
        except AttributeError:
            pass
        
        namespace[self.outputName] = cmodes

    @property
    def _key_choices(self):
        #try and find the available column names
        try:
            return sorted(self._parent.namespace[self.inputName].keys())
        except:
            return []

    @property
    def default_view(self):
        from traitsui.api import View, Group, Item
        from PYME.ui.custom_traits_editors import CBEditor

        return View(Item('inputName', editor=CBEditor(choices=self._namespace_keys)),
                    Item('_'),
                    Item('IDkey', editor=CBEditor(choices=self._key_choices)),
                    Item('PropertyKey', editor=CBEditor(choices=self._key_choices)),
                    Item('_'),
                    Item('outputName'), buttons=['OK'])


@register_module('SIMPLERzgenerator')
class SIMPLERzgenerator(ModuleBase):

    inputName = Input('filtered')
    outputName = Output('with_simpler_z')
    df_in_nm = Float(88.0)
    alphaf = Float(0.9)
    N0_scale_factor = Float(1.0)
    
    def execute(self, namespace):
        inp = namespace[self.inputName]
        mapped = tabular.mappingFilter(inp)

        from six.moves import cPickle
        fdialog = wx.FileDialog(None, 'Please select N0 map to use ...',
                    wildcard='N0 map file (*.n0m)|*.n0m', style=wx.FD_OPEN)
        succ = fdialog.ShowModal()
        if (succ == wx.ID_OK):
            nmFilename = fdialog.GetPath()
            with open(nmFilename, 'rb') as fid:
                n0m,bb = cPickle.load(fid)
        
            try:
                N0 = n0m(inp['x'],inp['y'],grid=False) # this should ensure N) is floating point type
            except TypeError:
                N0 = n0m(inp['x'],inp['y'])
            N0 *= self.N0_scale_factor
            N = inp['nPhotons']
            NoverN0 = N/N0
            simpler_z = self.df_in_nm*np.log(self.alphaf/(NoverN0 - (1 - self.alphaf)))
            simpler_z[np.isnan(simpler_z)] = -100.0

            mapped.addColumn('N0', N0)
            mapped.addColumn('NoverN0', NoverN0)
            mapped.addColumn('z', simpler_z)

        # propogate metadata, if present
        try:
            mapped.mdh = inp.mdh
        except AttributeError:
            pass
        
        namespace[self.outputName] = mapped
