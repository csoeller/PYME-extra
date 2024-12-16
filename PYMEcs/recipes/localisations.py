from PYME.recipes.base import register_module, ModuleBase, Filter
from PYME.recipes.traits import Input, Output, Float, Enum, CStr, Bool, Int, List, DictStrStr, DictStrList, ListFloat, ListStr, FileOrURI

import numpy as np
import pandas as pd
from PYME.IO import tabular
from PYME.warnings import warn

import logging
logger = logging.getLogger(__file__)

@register_module('NNdist')
class NNdist(ModuleBase):

    inputName = Input('coalesced')
    outputName = Output('withNNdist')

    def execute(self, namespace):
        inp = namespace[self.inputName]
        mapped = tabular.MappingFilter(inp)

        from scipy.spatial import KDTree
        coords = np.vstack([inp[k] for k in ['x','y','z']]).T
        tree = KDTree(coords)
        dd, ii = tree.query(coords,k=3)
        mapped.addColumn('NNdist', dd[:,1])
        mapped.addColumn('NNdist2', dd[:,2])
        
        try:
            mapped.mdh = inp.mdh
        except AttributeError:
            pass
        
        namespace[self.outputName] = mapped

@register_module('ClumpFromTID')
class ClumpFromTID(ModuleBase):
    """
    Generate a clump index from the tid field (Trace ID from MINFLUX data)
    generates contiguous indices, preferable if existing IDs not contiguous
    (as otherwise will create empty clumps with zero clumpsize after merging)

    Parameters
    ----------
    inputName: string - name of the data source containing a field named 'tid'
    
    Returns
    -------
    outputLocalizations : tabular.MappingFilter that contains new clumpIndex and clumpSize fields (will overwrite existing in inputName if present)
    
    """
    inputName = Input('with_TID')
    outputName = Output('with_clumps')

    def execute(self, namespace):
        inp = namespace[self.inputName]
        mapped = tabular.MappingFilter(inp)

        # we replace the non-sequential trace ids from MINFLUX data field TID with a set of sequential ids
        # this works better for clumpIndex assumptions in the end (we think)
        uids,revids,idcounts = np.unique(inp['tid'],return_inverse=True,return_counts=True)
        ids = np.arange(1,uids.size+1,dtype='int32')[revids]
        counts = idcounts[revids]

        mapped.addColumn('clumpIndex', ids)
        mapped.addColumn('clumpSize', counts)
        
        try:
            mapped.mdh = inp.mdh
        except AttributeError:
            pass
        
        namespace[self.outputName] = mapped

@register_module('NNdistMutual')
class NNdistMutual(ModuleBase):

    inputChan1 = Input('channel1')
    inputChan2 = Input('channel2')
    outputName = Output('c1withNNdist')

    def execute(self, namespace):
        inpc1 = namespace[self.inputChan1]
        inpc2 = namespace[self.inputChan2]
        mapped = tabular.MappingFilter(inpc1)

        from scipy.spatial import KDTree
        coords1 = np.vstack([inpc1[k] for k in ['x','y','z']]).T
        coords2 = np.vstack([inpc2[k] for k in ['x','y','z']]).T
        tree = KDTree(coords2)
        dd, ii = tree.query(coords1,k=2)
        mapped.addColumn('NNdistMutual', dd[:,0])
        mapped.addColumn('NNdistMutual2', dd[:,1])
        
        try:
            mapped.mdh = inpc1.mdh
        except AttributeError:
            pass
        
        namespace[self.outputName] = mapped

@register_module('NNfilter')
class NNfilter(ModuleBase):

    inputName = Input('withNNdist')
    outputName = Output('NNfiltered')
    nnMin = Float(0)
    nnMax = Float(1e5)
    

    def execute(self, namespace):
        inp = namespace[self.inputName]
        mapped = tabular.MappingFilter(inp)

        nn_good = (inp['NNdist'] > self.nnMin) * (inp['NNdist'] < self.nnMax)
        nn2_good = (inp['NNdist2'] > self.nnMin) * (inp['NNdist2'] < self.nnMax)
        nn_class = 1.0*nn_good + 2.0*(np.logical_not(nn_good)*nn2_good)
        mapped.addColumn('NNclass', nn_class)

        filterKeys = {'NNclass': (0.5,2.5)}
        filtered = tabular.ResultsFilter(mapped, **filterKeys)
        
        try:
            filtered.mdh = inp.mdh
        except AttributeError:
            pass
        
        namespace[self.outputName] = filtered

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
        mapped = tabular.MappingFilter(inp)

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
        mapped = tabular.MappingFilter(inp)

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
    Also adds the time period of bursts by adding clumpTmin, clumpTmax amd clumpLength columns.
    """
    inputName = Input('clumped')
    outputName = Output('merged')
    labelKey = CStr('clumpIndex')

    def execute(self, namespace):
        from PYME.Analysis.points.DeClump import pyDeClump
        from PYME.Analysis.points.DeClump import deClump as deClumpC
        
        inp = namespace[self.inputName]

        grouped = pyDeClump.mergeClumps(inp, labelKey=self.labelKey)

        # we need this because currently the addColumn method of DictSrc is broken
        def addColumn(dictsrc,name, values):
            if not isinstance(values, np.ndarray):
                raise TypeError('New column "%s" is not a numpy array' % name)

            if not len(values) == len(dictsrc):
                raise ValueError('Columns are different lengths')

            dictsrc._source[name] = values # this was missing I think

        # work out tmin and tmax
        I = np.argsort(inp[self.labelKey])
        sorted_src = {k: inp[k][I] for k in [self.labelKey,'t']}
        # tmin and tmax - tentative addition
        NClumps = int(np.max(sorted_src[self.labelKey]) + 1)
        tmin = deClumpC.aggregateMin(NClumps, sorted_src[self.labelKey].astype('i'), sorted_src['t'].astype('f'))
        tmax = -deClumpC.aggregateMin(NClumps, sorted_src[self.labelKey].astype('i'), -1.0*sorted_src['t'].astype('f'))
        if '_source' in dir(grouped): # appears to be a DictSource which currently has a broken addColumn method
            addColumn(grouped,'clumpTmin',tmin) # use our fixed function to add a column
            addColumn(grouped,'clumpTmax',tmax)
            addColumn(grouped,'clumpLength',tmax-tmin+1) # this only works if time is in frame units (otherwise a value different from +1 is needed)
        else:
            grouped.addColumn('clumpTmin',tmin)
            grouped.addColumn('clumpTmax',tmax)
            grouped.addColumn('clumpLength',tmax-tmin+1)
            
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
        
        mapped = tabular.MappingFilter(inp)
        out_f = tabular.MappingFilter(fiducial)
        
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
        mapped = tabular.MappingFilter(inp)

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
        mapped = tabular.MappingFilter(inp)

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
        mapped = tabular.MappingFilter(inp)

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
        mapped = tabular.MappingFilter(inp)
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
        mapped = tabular.MappingFilter(inp)
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
        mapped = tabular.MappingFilter(inp)
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
        mapped = tabular.MappingFilter(inp)

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


@register_module('subClump')
class subClump(ModuleBase):
    """
    Groups clumps into smaller sub clumps of given max sub clump size
    
    Parameters
    ----------

        labelKey: datasource key serving as input clump index

        subclumpMaxSize: target size of new subclumps - see also comments in the code

        clumpColumnName: column name of generated sub clump index

        sizeColumnName: column name of sub clump sizes
    
    """
    inputName = Input('with_clumps')
    labelKey = CStr('clumpIndex')
    subclumpMaxSize = Int(4)
    
    clumpColumnName = CStr('subClumpID')
    sizeColumnName = CStr('subClumpSize')
    outputName = Output('with_subClumps')

    def execute(self, namespace):

        inp = namespace[self.inputName]
        mapped = tabular.MappingFilter(inp)

        # subid calculations
        subids = np.zeros_like(inp['x'],dtype='int')
        clumpids = inp[self.labelKey]
        uids,revids,counts = np.unique(clumpids,return_inverse=True, return_counts=True)
        # the clumpsz per se should not be needed, we just use the counts below
        # clumpsz = counts[revids] # NOTE: could the case ID=0 be an issue? These are often points not part of clusters etc

        # the code below has the task to split "long" clumps into subclumps
        # we use the following strategy:
        #    - don't split clumps with less than 2*scMaxSize events
        #    - if clumps are larger than that break the clump into subclumps
        #       - each subclump has at least scMaxSize events
        #       - if there are "extra" events at the end that would not make a full scMaxSize clump
        #            then add these to the otherwise previous subclump
        #       - as a result generated subclumps have therefore from scMaxSize to 2*scMaxSize-1 events
        # NOTE 1: this strategy is not the only possible one, reassess as needed
        # NOTE 2: we do NOT explicitly sort by time of the events in the clump when grouping into sub clumps
        #         (should not matter but think about again)
        # NOTE 3: this has not necessarily been tuned for efficiency - tackle if needed
        
        scMaxSize = self.subclumpMaxSize
        curbaseid = 1
        for i,id in enumerate(uids):
            if id > 0:
                if counts[i] > 2*scMaxSize-1:
                    cts = counts[i]
                    ctrange = np.arange(cts,dtype='int') # a range we can use in splitting the events up into subclumps
                    subs = cts // scMaxSize # the number of subclumps we are going to make
                    # the last bit of the expression below ensures that events that would not make a full size subClump
                    # get added to the last of the subclumps; a subclump has there from scMaxSize to 2*scMaxSize-1 events
                    subids[revids == i] = curbaseid + ctrange // scMaxSize - (ctrange >= (subs * scMaxSize))
                    curbaseid += subs
                else:
                    subids[revids == i] = curbaseid
                    curbaseid += 1

        suids,srevids,scounts = np.unique(subids,return_inverse=True, return_counts=True) # generate counts for new subclump IDs
        subclumpsz = scounts[srevids] # NOTE: could the case ID==0 be an issue? These are often points not part of clusters etc
        # end subid calculations
        
        mapped.addColumn(str(self.clumpColumnName), subids)
        mapped.addColumn(str(self.sizeColumnName), subclumpsz)
        
        # propogate metadata, if present
        try:
            mapped.mdh = inp.mdh
        except AttributeError:
            pass

        namespace[self.outputName] = mapped

        
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
    searchRadius = Float(10)
    minClumpSize = Int(1)
    
    #exposes sklearn parallelism. Recipe modules are generally assumed
    #to be single-threaded. Enable at your own risk
    multithreaded = Bool(False)
    numberOfJobs = Int(max(multiprocessing.cpu_count()-1,1))
    
    clumpColumnName = CStr('dbscanClumpID')
    sizeColumnName = CStr('dbscanClumpSize')
    outputName = Output('dbscanClustered')

    def execute(self, namespace):
        from sklearn.cluster import dbscan
        from scipy.stats import binned_statistic

        inp = namespace[self.inputName]
        mapped = tabular.MappingFilter(inp)

        # Note that sklearn gives unclustered points label of -1, and first value starts at 0.
        if self.multithreaded:
            core_samp, dbLabels = dbscan(np.vstack([inp[k] for k in self.columns]).T,
                                         eps=self.searchRadius, min_samples=self.minClumpSize, n_jobs=self.numberOfJobs)
        else:
            #NB try-catch from Christians multithreaded example removed as I think we should see failure here
            core_samp, dbLabels = dbscan(np.vstack([inp[k] for k in self.columns]).T,
                                     eps=self.searchRadius, min_samples=self.minClumpSize)

        # shift dbscan labels up by one to match existing convention that a clumpID of 0 corresponds to unclumped
        dbids = dbLabels + 1
        maxid = int(dbids.max())
        edges = -0.5+np.arange(maxid+2)
        resstat = binned_statistic(dbids, np.ones_like(dbids), statistic='sum', bins=edges)

        mapped.addColumn(str(self.clumpColumnName), dbids)
        mapped.addColumn(str(self.sizeColumnName),resstat[0][dbids])
        
        # propogate metadata, if present
        try:
            mapped.mdh = inp.mdh
        except AttributeError:
            pass

        namespace[self.outputName] = mapped

    
    @property
    def hide_in_overview(self):
        return ['columns']

    def _view_items(self, params=None):
        from traitsui.api import Item, TextEditor
        return [Item('columns', editor=TextEditor(auto_set=False, enter_set=True, evaluate=ListStr)),
                    Item('searchRadius'),
                    Item('minClumpSize'),
                    Item('multithreaded'),
                    Item('numberOfJobs'),
                    Item('clumpColumnName'),
                    Item('sizeColumnName'),]


@register_module('SnrCalculation')
class SnrCalculation(ModuleBase):
    inputName = Input('filtered')
    outputName = Output('snr')
   
    def execute(self, namespace):
        inp = namespace[self.inputName]
        mapped = tabular.MappingFilter(inp)

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
        mapped = tabular.MappingFilter(inp)
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
        mapped = tabular.MappingFilter(inp)

        fdialog = wx.FileDialog(None, 'Please select PSF to use ...',
                                #defaultDir=os.path.split(self.image.filename)[0],
                                wildcard='PSF Files|*.psf|TIFF files|*.tif', style=wx.FD_OPEN)
        succ = fdialog.ShowModal()
        if (succ == wx.ID_OK):
            psfn = filename = fdialog.GetPath()
           
            mdh = inp.mdh
            if mdh.getEntry('Analysis.FitModule') not in ['SplitterFitInterpBNR']:
                warn('Plugin works only for Biplane analysis')
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

@register_module('NPCAnalysisByID')
class NPCAnalysisByID(ModuleBase):
    
    inputName = Input('processed')
    outputName = Output('withNPCanalysis')

    IDkey = CStr('objectID')
    SegmentThreshold = Int(10)
    SecondPass = Bool(False)
    FitMode = Enum(['abs','square'])

    def execute(self, namespace):
        from PYMEcs.Analysis.NPC import estimate_nlabeled
        npcs = namespace[self.inputName]
        ids = npcs[self.IDkey]
        x = npcs['x']
        y = npcs['y']

        uids,revids,idcounts = np.unique(ids,return_inverse=True,return_counts=True)

        npcnlabeled = np.zeros_like(uids,dtype = 'int')
        npcradius = np.zeros_like(uids,dtype = 'float')
        for i, id in enumerate(uids):
            if (id > 0) and (idcounts[i] > 0):
                idx_thisid = ids == id
                xid = x[idx_thisid]
                yid = y[idx_thisid]
                npcnlabeled[i],npcradius[i] = estimate_nlabeled(xid,yid,nthresh=self.SegmentThreshold,
                                                   do_plot=False,secondpass=self.SecondPass,
                                                   fitmode=self.FitMode, return_radius=True)

        npcnall = npcnlabeled[revids] # map back on all events
        npcradall = npcradius[revids] # map back on all events
        mapped = tabular.MappingFilter(npcs)
        mapped.addColumn('NPCnlabeled', npcnall)
        mapped.addColumn('NPCradius', npcradall)

        try:
            mapped.mdh = NestedClassMDHandler(npcs.mdh)
        except AttributeError:
            mapped.mdh = NestedClassMDHandler() # make empty mdh

        mapped.mdh['NPCAnalysis.EventThreshold'] = self.SegmentThreshold
        mapped.mdh['NPCAnalysis.RotationAlgorithm'] = self.FitMode
        mapped.mdh['NPCAnalysis.SecondPass'] = self.SecondPass
        
        namespace[self.outputName] = mapped

        
from scipy.stats import binned_statistic
from scipy.signal import savgol_filter
from scipy.interpolate import CubicSpline

# a function that is supposed to return a smoothed site trajectory as a function
# as input we need site coordinates
# the choice of savgol_filter for smoothing and CubicSpline for interpolation is a little arbitrary for now
# and can be in future adjusted as needed
# the bins need to be chosen in a robust manner - FIX
def smoothed_site_func(t,coord_site,statistic='mean',bins=75,sgwindow_length=10,sgpolyorder=6):
    csitem, tbins, binnum = binned_statistic(t,coord_site,statistic=statistic,bins=bins)
    tmid = 0.5*(tbins[0:-1]+tbins[1:])
    return CubicSpline(tmid,savgol_filter(csitem,10,6))

from scipy.optimize import curve_fit
def gfit_func(x, a0, a1, a2):
    z = (x - a1) / a2
    y = a0 * np.exp(-z**2 / 2)
    return y

def gaussfit(x,y,params=None):
    try:
        parameters, covar = curve_fit(gfit_func, x, y, params)
    except RuntimeError:
        return None, None
    else:
        return parameters, covar

def fit_binedges(tms,delta_s=100,width_s=500):
    tmin = tms.min()
    tmax = tms.max()
    delta_ms = 1e3*delta_s
    width_ms = 1e3*width_s
    numsamp = (tmax-tmin)/delta_ms
    t_samp = np.linspace(tmin,tmax,int(numsamp+0.5))
    edgeleft = np.maximum(t_samp - 0.5*width_ms, tmin)
    edgeright = np.minimum(t_samp + 0.5*width_ms, tmax)
    return (edgeleft, edgeright,t_samp)

def fit_sitedrift(tms,xsite,delta_s=100,width_s=500,bins=20):
    bel,ber,t_samp=fit_binedges(tms,delta_s=delta_s,width_s=width_s)
    pars = np.zeros((3,bel.size))
    errs = np.zeros_like(pars)
    for i in range(bel.size):
        tgood = (tms >= bel[i]) & (tms<=ber[i])
        hv, be = np.histogram(xsite[tgood],bins=bins)
        bctrs = 0.5*(be[:-1]+be[1:])
        par,cov = gaussfit(bctrs,hv,params=(10.,1.,3.)) # reasonable starting values
        pars[:,i] = par
        errs[:,i] = np.sqrt(np.diag(cov))

    return t_samp, pars, errs

def site_fit_gaussian(tms,xsite,delta_s=100,width_s=500):
    t_samp, pars, errs = fit_sitedrift(tms,xsite,delta_s=delta_s,width_s=width_s)
    ifunc = CubicSpline(t_samp, pars[1,:])
    return ifunc

class ModuleBaseMDHmod(ModuleBase):
    """
    Exactly like ModuleBase but with a modified execute method
    that allows passing an 'mdh' key in the dict to modify mdh from run method
    """
    # NOTE: we override the 'execute' method here because we want to meddle with the metadata from within the run method;
    #       a bit of a hack but mainly used
    #       for experimentation in this case; we may solve this differently in future if really needed
    def execute(self, namespace):
        """
        takes a namespace (a dictionary like object) from which it reads its inputs and
        into which it writes outputs

        NOTE: This was previously the function to define / override to make a module work. To support automatic metadata propagation
        and reduce the ammount of boiler plate, new modules should override the `run()` method instead.
        """
        from PYME.IO import MetaDataHandler
        inputs = {k: namespace[v] for k, v in self._input_traits.items()}

        ret = self.run(**inputs)
        mdhret = None # the default unless 'mdh' is in returned dictironary

        # convert output to a dictionary if needed
        if isinstance(ret, dict):
            out = {k : ret[v] for v, k in self._output_traits.items()}
            if 'mdh' in ret:
                mdhret = ret['mdh']
        elif isinstance(ret, List):
            out = {k : v  for k, v in zip(self.outputs, ret)} #TODO - is this safe (is ordering consistent)
        else:
            # single output
            if len(self.outputs) > 1:
                raise RuntimeError('Module has multiple outputs, but .run() returns a single value')

            out = {list(self.outputs)[0] : ret}

        # complete metadata (injecting as appropriate)
        mdhin = MetaDataHandler.DictMDHandler(getattr(list(inputs.values())[0], 'mdh', None))
        
        mdh = MetaDataHandler.DictMDHandler()
        self._params_to_metadata(mdh)

        for v in out.values():
            if v is None:
                continue
            if getattr(v, 'mdh', None) is None:
                v.mdh = MetaDataHandler.DictMDHandler()

            v.mdh.mergeEntriesFrom(mdhin) #merge, to allow e.g. voxel size overrides due to downsampling
            #print(v.mdh, mdh)
            v.mdh.copyEntriesFrom(mdh) # copy / overwrite with module processing parameters
            if mdhret is not None:
                v.mdh.copyEntriesFrom(mdhret)

        namespace.update(out)


@register_module('OrigamiSiteTrack')
class OrigamiSiteTrack(ModuleBaseMDHmod):
    """
    Recipe module aimed at analysing oregami MINFLUX datasets. More docs to be added.

    Inputs
    ------
    inputClusters: name of tabular input containing already
        "labeled" clusters used for lining up
    inputSites: name of tabular input containing already
        coalesced clusters that should match the inputClusters in terms of label etc
    inputAllPoints: [optional] tabular input that contains typically
        "all" data to which the calculated drift track should be applied

    Outputs
    -------
    outputName: name of tabular output that will contain
        coordinate corrected
        version of inputClusters. Also contains a wealth of additional info, such
        as new error estimates based on site "scatter" and preserves the
        original and also precorrected versions of those properties.
    outputAllPoints: [optional] tabular output that contains the
        coordinate corrected (drift-corrected) version
        of inputAllPoints. Will only be generated if optional
        input inputAllPoints was supplied
    
    """
    inputClusters = Input('siteclusters')
    inputSites = Input('sites')
    inputAllPoints = Input('')  # optional input to apply the correction to
    outputName = Output('corrected_siteclusters')
    outputAllPoints = Output('') # optional output if our inputAllPoints was actually used
    
    labelKey = CStr('dyeID',
                    desc="property name of the ID indentifying site clusters; needs to match the ID name generated by DBSCAN clustering module") # should really be siteID
    smoothingBinWidthsSeconds = Float(200,label='temporal binning (s)',
                                      desc="parameter that sets the temporal binning (in units of s) when" +
                                      " estimating a smoothed drift trajectory from the data")
    savgolWindowLength = Int(10,label='savgol filter window length',
                             desc="the window_length argument of the scipy.signal.savgol_filter function used in smoothing to obtain the drift trajectory")
    savgolPolyorder = Int(6,label='savgol filter polynomial order',
                          desc="the polyorder argument of the scipy.signal.savgol_filter function used in smoothing to obtain the drift trajectory")
    binnedStatistic = Enum(['mean','median','Gaussian'],
                           desc="statistic for smoothing when using binned_statistic function on data to obtain the drift trajectory by time windowing of the 'siteclouds'")
    gaussianBinSizeSeconds = Float(500)
    
    def run(self, inputClusters, inputSites,inputAllPoints=None):
        site_id = self.labelKey
        idsunique = inputSites[site_id].astype('i')
        # centroid coordinates
        xc = inputSites['x']
        yc = inputSites['y']
        zc = inputSites['z']

        # actual localisation coordinates
        x = inputClusters['x']
        y = inputClusters['y']
        z = inputClusters['z']
        t = inputClusters['t']

        ids = inputClusters[site_id]

        xsite = np.zeros_like(x)
        ysite = np.zeros_like(y)
        zsite = np.zeros_like(z)
        xerr = np.zeros_like(x)
        xerrnc = np.zeros_like(x)
        yerr = np.zeros_like(x)
        yerrnc = np.zeros_like(x)
        zerr = np.zeros_like(x)
        zerrnc = np.zeros_like(x)
        
        for j,id in enumerate(idsunique):
            idx = ids == id
            xsite[idx] = x[idx]-xc[j]
            ysite[idx] = y[idx]-yc[j]
            zsite[idx] = z[idx]-zc[j]
            xerrnc[idx] = np.std(x[idx])
            yerrnc[idx] = np.std(y[idx])
            zerrnc[idx] = np.std(z[idx])

        trange = t.max() - t.min()
        delta_ms = self.smoothingBinWidthsSeconds * 1e3 # 200 s 
        nbins = int(trange/delta_ms)

        if self.binnedStatistic in ['mean','median']:
            # we should check that with this choice of nbins we get no issues! (i.e. too few counts in some bins)
            c_xsite = smoothed_site_func(t,xsite,bins=nbins,statistic=self.binnedStatistic,
                                         sgwindow_length=self.savgolWindowLength,sgpolyorder=self.savgolPolyorder)
            c_ysite = smoothed_site_func(t,ysite,bins=nbins,statistic=self.binnedStatistic,
                                         sgwindow_length=self.savgolWindowLength,sgpolyorder=self.savgolPolyorder)
            c_zsite = smoothed_site_func(t,zsite,bins=nbins,statistic=self.binnedStatistic,
                                         sgwindow_length=self.savgolWindowLength,sgpolyorder=self.savgolPolyorder)
        else:
            c_xsite = site_fit_gaussian(t,xsite,delta_s=self.smoothingBinWidthsSeconds,width_s=self.gaussianBinSizeSeconds)
            c_ysite = site_fit_gaussian(t,ysite,delta_s=self.smoothingBinWidthsSeconds,width_s=self.gaussianBinSizeSeconds)
            c_zsite = site_fit_gaussian(t,zsite,delta_s=self.smoothingBinWidthsSeconds,width_s=self.gaussianBinSizeSeconds)

        for j,id in enumerate(idsunique):
            idx = ids == id
            xerr[idx] = np.std(x[idx]-c_xsite(t[idx]))
            yerr[idx] = np.std(y[idx]-c_ysite(t[idx]))
            zerr[idx] = np.std(z[idx]-c_zsite(t[idx]))
            
        # note to self: we shall preserve the site coordinates in the new data source
        # new properties to create: [xyz]_site, [xyz]_ori, new [xyz]

        mapped_ds = tabular.MappingFilter(inputClusters)
        mapped_ds.setMapping('x_ori', 'x')
        mapped_ds.setMapping('y_ori', 'y')
        mapped_ds.setMapping('z_ori', 'z')

        mapped_ds.setMapping('error_x_ori', 'error_x')
        mapped_ds.setMapping('error_y_ori', 'error_y')
        mapped_ds.setMapping('error_z_ori', 'error_z')
        
        # mapped_ds.addColumn('x_ori', x)
        # mapped_ds.addColumn('y_ori', y)
        # mapped_ds.addColumn('z_ori', z)
        
        mapped_ds.addColumn('x_site_nc', xsite)
        mapped_ds.addColumn('y_site_nc', ysite)
        mapped_ds.addColumn('z_site_nc', zsite)

        mapped_ds.addColumn('x_site', xsite-c_xsite(t))
        mapped_ds.addColumn('y_site', ysite-c_ysite(t))
        mapped_ds.addColumn('z_site', zsite-c_zsite(t))

        mapped_ds.addColumn('error_x_nc', xerrnc)
        mapped_ds.addColumn('error_y_nc', yerrnc)
        mapped_ds.addColumn('error_z_nc', zerrnc)

        mapped_ds.addColumn('error_x', xerr)
        mapped_ds.addColumn('error_y', yerr)
        mapped_ds.addColumn('error_z', zerr)
        
        mapped_ds.addColumn('x', x-c_xsite(t))
        mapped_ds.addColumn('y', y-c_ysite(t))
        mapped_ds.addColumn('z', z-c_zsite(t))

        if 'driftx' in inputClusters.keys():
            mapped_ds.setMapping('driftx_ori', 'driftx')
            mapped_ds.setMapping('drifty_ori', 'drifty')
            mapped_ds.setMapping('driftz_ori', 'driftz')

        mapped_ds.addColumn('driftx', c_xsite(t))
        mapped_ds.addColumn('drifty', c_ysite(t))
        mapped_ds.addColumn('driftz', c_zsite(t))

        if inputAllPoints is not None:
            mapped_ap = tabular.MappingFilter(inputAllPoints)
            # actual localisation coordinates
            x = inputAllPoints['x']
            y = inputAllPoints['y']
            z = inputAllPoints['z']
            t = inputAllPoints['t']

            mapped_ap.addColumn('x', x-c_xsite(t))
            mapped_ap.addColumn('y', y-c_ysite(t))
            mapped_ap.addColumn('z', z-c_zsite(t))

            mapped_ap.addColumn('driftx', c_xsite(t))
            mapped_ap.addColumn('drifty', c_ysite(t))
            mapped_ap.addColumn('driftz', c_zsite(t))
        else:
            # how to deal with an "optional" output
            # this would be a dummy assignment in the absence of inputAllPoints
            # mapped_ap = tabular.MappingFilter(inputClusters)
            mapped_ap = None # returning none in this case seems better and appears to work

        # for now we disable this way of passing the correction function via metadata
        if False:
            from PYME.IO import MetaDataHandler

            mdh = MetaDataHandler.DictMDHandler()
            mdhin = inputClusters.mdh
        
            def checkmdh(item,mdh):
                if item in mdh:
                    return "%s_1" % item
                else:
                    return item
            
            mdh[checkmdh('Processing.DriftFuncs.X',mdhin)] = c_xsite
            mdh[checkmdh('Processing.DriftFuncs.Y',mdhin)] = c_ysite
            mdh[checkmdh('Processing.DriftFuncs.Z',mdhin)] = c_zsite
            mdh[checkmdh('Processing.DriftFuncs.TRange',mdhin)] = [t.min(),t.max()]
        
        return {'outputName': mapped_ds, 'outputAllPoints' : mapped_ap, 'mdh' : None } # pass proper mdh instead of None if metadata output needed


import numpy as np
import scipy.special

@register_module('MINFLUXcolours')
class MINFLUXcolours(ModuleBase):

    inputLocalizations = Input('localizations')
    output = Output('localizations_mcolour') # localisations with MINFLUX colour info

    dcrisgfrac = Bool(False)
    
    def run(self,inputLocalizations):
    
        mapped_ds = tabular.MappingFilter(inputLocalizations)
        mapped_ds.setMapping('A','nPhotons')

        if self.dcrisgfrac:
            mapped_ds.setMapping('fitResults_Ag','dcr*nPhotons')
            mapped_ds.setMapping('fitResults_Ar','(1-dcr)*nPhotons')
        else:
            mapped_ds.setMapping('fitResults_Ag','nPhotons*dcr/(1+dcr)')
            mapped_ds.setMapping('fitResults_Ar','nPhotons/(1+dcr)')


        mapped_ds.setMapping('fitError_Ag','1*sqrt(fitResults_Ag/1)')
        mapped_ds.setMapping('fitError_Ar','1*sqrt(fitResults_Ar/1)')

        mapped_ds.setMapping('gFrac', 'fitResults_Ag/(fitResults_Ag + fitResults_Ar)')
        mapped_ds.setMapping('error_gFrac',
                      'sqrt((fitError_Ag/fitResults_Ag)**2 + (fitError_Ag**2 + fitError_Ar**2)/(fitResults_Ag + fitResults_Ar)**2)' +
                      '*fitResults_Ag/(fitResults_Ag + fitResults_Ar)')

        sg = mapped_ds['fitError_Ag']
        sr = mapped_ds['fitError_Ar']
        g = mapped_ds['fitResults_Ag']
        r = mapped_ds['fitResults_Ar']
        I = mapped_ds['A']

        colNorm = np.sqrt(2 * np.pi) * sg * sr / (2 * np.sqrt(sg ** 2 + sr ** 2) * I) * (
            scipy.special.erf((sg ** 2 * r + sr ** 2 * (I - g)) / (np.sqrt(2) * sg * sr * np.sqrt(sg ** 2 + sr ** 2)))
            - scipy.special.erf((sg ** 2 * (r - I) - sr ** 2 * g) / (np.sqrt(2) * sg * sr * np.sqrt(sg ** 2 + sr ** 2))))

        mapped_ds.addColumn('ColourNorm', colNorm)

        return mapped_ds

from pathlib import Path
def check_mbm_name(mbmfilename,timestamp,endswith='__MBM-beads'):
    if timestamp is None:
        return True
    mbmp = Path(mbmfilename)
    
    # should return False if warning is necessary
    return mbmp.stem.startswith(timestamp) and (mbmp.stem.endswith(endswith) or mbmp.suffix.endswith('.zip'))

def get_bead_dict_from_mbm(mbm):
    beads = {}
    raw_beads = mbm._raw_beads
    for bead in raw_beads:
        beads[bead] = {}
        beads[bead]['x'] = 1e9*raw_beads[bead]['pos'][:,0]
        beads[bead]['y'] = 1e9*raw_beads[bead]['pos'][:,1]
        beads[bead]['z'] = 1e9*raw_beads[bead]['pos'][:,2]
        beads[bead]['A'] = raw_beads[bead]['str']
        beads[bead]['t'] = np.asarray(1e3*raw_beads[bead]['tim'],dtype=int)
        beads[bead]['tim'] = raw_beads[bead]['tim']
        if 'tid' in raw_beads[bead].dtype.fields:
            beads[bead]['tid'] = raw_beads[bead]['tid']
        if 'gri' in raw_beads[bead].dtype.fields:
            beads[bead]['tid'] = raw_beads[bead]['gri']
            

    x = np.empty((0))
    y = np.empty((0))
    z = np.empty((0))
    t = np.empty((0),int)
    tid = np.empty((0),int)
    beadID = np.empty((0),int)
    objectID = np.empty((0),int)
    A = np.empty((0))
    tim = np.empty((0))
    good = np.empty((0),int)

    for bead in beads:
        beadid = int(bead[1:])
        beadisgood = mbm.beadisgood[bead]
        # print('beadid %d' % beadid)
        x = np.append(x,beads[bead]['x'])
        y = np.append(y,beads[bead]['y'])
        z = np.append(z,beads[bead]['z'])
        t = np.append(t,beads[bead]['t'])
        tim = np.append(tim,beads[bead]['tim'])
        tid = np.append(tid,beads[bead]['tid'])
        A = np.append(A,beads[bead]['A'])
        beadID = np.append(beadID,np.full_like(beads[bead]['x'],beadid,dtype=int))
        # objectIDs start at 1, beadIDs can start at 0
        objectID = np.append(objectID,np.full_like(beads[bead]['x'],beadid+1,dtype=int))
        good = np.append(good,np.full_like(beads[bead]['x'],beadisgood,dtype=int))

    return dict(x=x,y=y,z=z,t=t,tim=tim,beadID=beadID,tid=tid,
                A=A,good=good,objectID=objectID)

# remove this once sessionpath PR is accepted!!
try:
    from PYME.LMVis.sessionpaths import register_path_modulechecks
    register_path_modulechecks('PYMEcs.MBMcorrection','mbmfile','mbmsettings')
except ImportError:
    pass

# code on suggestion from https://blog.finxter.com/5-best-ways-to-compute-the-hash-of-a-python-tuple/
import hashlib
def tuple_hash(tuple_obj):
  hasher = hashlib.sha256()
  hasher.update(repr(tuple_obj).encode())
  return hasher.hexdigest()

@register_module('MBMcorrection')
class MBMcorrection(ModuleBaseMDHmod):
    inputLocalizations = Input('localizations')
    output = Output('mbm_corrected')
    outputTracks = Output('mbm_tracks')
    outputTracksCorr = Output('mbm_tracks_corrected')
    
    mbmfile = FileOrURI('')
    mbmsettings = FileOrURI('')
    mbmfilename_checks = Bool(True)
    
    Median_window = Int(5)
    MBM_lowess_fraction = Float(0.1,label='lowess fraction for MBM smoothing',
                                desc='lowess fraction used for smoothing of mean MBM trajectories (default 0.1); 0 = no smoothing')
    MBM_beads = List() # this is a dummy to make sure older style PVS files are read ok - TODO: get rid off at some stage!!!
    _MBM_beads = List() # this one does the real work and with the leading "_" is NOT treated as a parameter that "fires" the module! 

    _mbm_allbeads = List()

    _mbm_cache = {}
    _lowess_cache = {}

    def lowess_cachetuple(self):
        mbm = self.getmbm()
        return (Path(self.mbmfile).name,self.MBM_lowess_fraction,self.Median_window,str(mbm.beadisgood))
    
    def lowess_cachekey(self):
        return tuple_hash(self.lowess_cachetuple())

    def lowess_cachehit(self):
        cachekey = self.lowess_cachekey()
        if cachekey in self._lowess_cache:
            logger.debug("CACHEHIT from in-memory-cache")
            return self._lowess_cache[cachekey]
        elif self.lowess_chacheread():
            logger.debug("CACHEHIT from file-cache")
            return self._lowess_cache[cachekey]
        else:
            logger.debug("NO CACHEHIT!!! cache tuple: %s" % str(self.lowess_cachetuple()))
            return None
        
    def lowess_cachestore(self,axis,value):
        cachekey = self.lowess_cachekey()
        if cachekey not in self._lowess_cache:
            self._lowess_cache[cachekey] = {}
        self._lowess_cache[cachekey][axis] = value

    def lowess_cachesave(self):
        cachehit = self.lowess_cachehit()
        if cachehit is not None:
            fpath = self.lowess_cachefilepath()
            np.savez(str(fpath),**cachehit)

    def lowess_chacheread(self):
        cachekey = self.lowess_cachekey()
        fpath = self.lowess_cachefilepath()
        if fpath.exists():
            self._lowess_cache[cachekey] = np.load(str(fpath))
            return True
        else:
            return False

    def lowess_cachefilepath(self):
        cachekey = self.lowess_cachekey()
        from pathlib import Path
        fdir = Path(self.mbmfile).parent
        fpath = fdir / (".mbm_lowess_%s.npz" % cachekey)
        return fpath
    
    def getmbm(self):
        mbmkey = self.mbmfile
        if mbmkey in self._mbm_cache:
            return self._mbm_cache[mbmkey]
        else:
            return None

    # cached lowess calc with time coordinate retrieval
    def lowess_calc(self,axis):
        mbm = self.getmbm()
        axismean = mbm.mean(axis)
        axismean_g = axismean[~np.isnan(axismean)]
        t_g = mbm.t[~np.isnan(axismean)]
        cachehit = self.lowess_cachehit() # LOWESS CACHE OP
        if self.MBM_lowess_fraction > 1e-5:
            from statsmodels.nonparametric.smoothers_lowess import lowess
            if cachehit is not None and axis in cachehit: # not all axes may have been calculated yet, so check for axis!
                axismean_sm = cachehit[axis] # LOWESS CACHE OP
            else:
                axismean_sm = lowess(axismean_g, t_g, frac=self.MBM_lowess_fraction,
                                     return_sorted=False)
                self.lowess_cachestore(axis,axismean_sm) # LOWESS CACHE OP
            return (t_g,axismean_sm)
        else:
            return (t_g,axismean_g)
    
    def run(self,inputLocalizations):
        import json
        from PYME.IO import unifiedIO
        from pathlib import Path
        from PYMEcs.Analysis.MBMcollection import MBMCollectionDF

        mapped_ds = tabular.MappingFilter(inputLocalizations)

        if self.mbmfile != '':
            mbmkey = self.mbmfile
            if mbmkey not in self._mbm_cache.keys():
                from PYMEcs.IO.MINFLUX import foreshortening
                mbm = MBMCollectionDF(name=Path(self.mbmfile).stem,filename=self.mbmfile,
                                      foreshortening=foreshortening)
                self._mbm_allbeads = [bead for bead in mbm.beadisgood if np.sum(np.logical_not(np.isnan(mbm.beads['x'][bead]))) > 1] # exclude "empty" trajectories
                self._MBM_beads = self._mbm_allbeads
                self._mbm_cache[mbmkey] = mbm
            else:
                mbm = self._mbm_cache[mbmkey]

            mbmsettingskey = self.mbmsettings
            if mbmsettingskey != '' and mbmsettingskey not in self._mbm_cache.keys():
                s = unifiedIO.read(self.mbmsettings)
                mbmconf = json.loads(s)

                for bead in mbmconf['beads']:
                    mbm.beadisgood[bead] =  mbmconf['beads'][bead]
                self._MBM_beads = [bead for bead in mbmconf['beads'].keys() if mbmconf['beads'][bead] and bead in self._mbm_allbeads]
                self._mbm_cache[mbmsettingskey] = mbmconf
            # in a second pass (if just loaded from mbmconf) set beadisgood to the useful self._MBM_beads subset
            for bead in mbm.beadisgood:
                mbm.beadisgood[bead] = bead in self._MBM_beads

            mbm.median_window = self.Median_window

            bead_ds_dict = get_bead_dict_from_mbm(mbm)
            tnew = 1e-3*inputLocalizations['t']
            mbmcorr = {}
            mbmtrack_corr = {}
            for axis in ['x','y','z']:
                tsm, axismean_sm = self.lowess_calc(axis) # cache based lowess calculation
                axis_interp = np.interp(tnew,tsm,axismean_sm)
                mbmcorr[axis] = axis_interp
                mbmtrack_corr[axis] = np.interp(bead_ds_dict['tim'],tsm,axismean_sm)
                axis_nc = "%s_nc" % axis
                mapped_ds.addColumn('mbm%s' % axis, mbmcorr[axis])
                if axis_nc in inputLocalizations.keys():
                    mapped_ds.addColumn(axis,inputLocalizations[axis_nc] - mbmcorr[axis])

            if self.mbmfilename_checks:
                if not check_mbm_name(self.mbmfile,inputLocalizations.mdh.get('MINFLUX.TimeStamp')):
                    warn("check MBM filename (%s) vs Series timestamp (%s)" %
                         (Path(self.mbmfile).name,inputLocalizations.mdh.get('MINFLUX.TimeStamp')))
                if mbmsettingskey != '' and not check_mbm_name(self.mbmsettings,inputLocalizations.mdh.get('MINFLUX.TimeStamp'),endswith='npz-settings'):
                    warn("check MBM settings filename (%s) vs Series timestamp (%s)" %
                         (Path(self.mbmsettings).name,inputLocalizations.mdh.get('MINFLUX.TimeStamp')))

            from PYME.IO import MetaDataHandler
            MBMmdh = MetaDataHandler.DictMDHandler()
            MBMmdh['Processing.MBMcorrection.mbm'] = mbm

            from PYME.IO.tabular import DictSource
            mbmtrack_ds = DictSource(bead_ds_dict)
            mbmtrack_dscorr = tabular.MappingFilter(mbmtrack_ds)
            for axis in ['x','y','z']:
                mbmtrack_dscorr.addColumn(axis,mbmtrack_ds[axis] - mbmtrack_corr[axis])
            
            return dict(output=mapped_ds, outputTracks=mbmtrack_ds,
                        outputTracksCorr=mbmtrack_dscorr, mdh=MBMmdh)
            # mapped_ds.mbm = mbm # attach mbm object to the output

        return dict(output=mapped_ds, outputTracks=None,
                        outputTracksCorr=None, mdh=None)  # return empty slots for outputs etc in this case
    
    @property
    def default_view(self):
        from traitsui.api import View, Group, Item, CheckListEditor
        from PYME.ui.custom_traits_editors import CBEditor

        return View(Item('inputLocalizations', editor=CBEditor(choices=self._namespace_keys)),
                    Item('_'),
                    Item('mbmfile'),
                    Item('mbmsettings'),
                    Item('mbmfilename_checks'),
                    Item('_MBM_beads', editor=CheckListEditor(values=self._mbm_allbeads,cols=4),
                         style='custom',
                         ),
                    Item('Median_window'),
                    Item('MBM_lowess_fraction'),
                    Item('_'),
                    Item('output'),
                    Item('outputTracks'),
                    Item('outputTracksCorr'),
                    buttons=['OK'])


try:
    import alphashape
except ImportError:
    has_ashape = False
else:
    has_ashape = True

from scipy.spatial import ConvexHull

def shape_measure_alpha(points,alpha=0.01):
    alpha_shape = alphashape.alphashape(points[:,0:2], alpha)
    # alpha_vol = alphashape.alphashape(points, 1e-3) # errors with 3D shapes occasionally; try to fix this to convex_hull_like

    area = alpha_shape.area
    # vol = alpha_vol.volume
    vol = ConvexHull(points[:,0:3]).volume # we fal back to convex hull to avoid issues
    nlabel = points.shape[0]

    if alpha_shape.geom_type == 'MultiPolygon':
        polys = []
        for geom in alpha_shape.geoms:
            pc = []
            for point in geom.exterior.coords:
                pc.append(point)
            polys.append(np.array(pc))        
    elif alpha_shape.geom_type == 'Polygon':
        polys = [np.array(alpha_shape.boundary.coords)]
    else:
        raise RuntimeError("Got alpha shape that is bounded by unknown geometry, got geom type %s with alpha = %.2f" %
                           (alpha_shape.geom_type,alpha))

    return (area,vol,nlabel,polys)

def shape_measure_convex_hull(points):
    shape2d = ConvexHull(points[:,0:2])
    area = shape2d.volume
    poly0 = points[shape2d.vertices,0:2]
    polys = [np.append(poly0,poly0[0,None,:],axis=0)] # close the polygon
    shape3d = ConvexHull(points[:,0:3])
    vol = shape3d.volume
    nlabel = points.shape[0]

    return (area,vol,nlabel,polys)

import PYME.config
def shape_measure(points,alpha=0.01,useAlphaShapes=True):
    if has_ashape and useAlphaShapes:
        return shape_measure_alpha(points,alpha=alpha)
    else:
        return shape_measure_convex_hull(points)

@register_module('SiteDensity')
class SiteDensity(ModuleBase):
    """Documentation to be added."""
    inputLocalisations = Input('clustered')
    outputName = Output('cluster_density')
    outputShapes = Output('cluster_shapes')
    
    IDName = CStr('dbscanClumpID')
    clusterSizeName = CStr('dbscanClumpSize')
    alpha = Float(0.01) # trying to ensure wqe stay with single polygon boundaries by default
    use_alpha_shapes = Bool(True,
                            desc="use alpha shapes IF available, use convex hull if false (on systems without alphashape module falls back to convex hull in any case)")

    def run(self, inputLocalisations):

        ids = inputLocalisations[self.IDName]
        uids = np.unique(ids)
        sizes = inputLocalisations[self.clusterSizeName] # really needed?

        # stddevs
        sx = np.zeros_like(ids,dtype=float)
        sy = np.zeros_like(ids,dtype=float)
        sz = np.zeros_like(ids,dtype=float)

        # centroids
        cx = np.zeros_like(ids,dtype=float)
        cy = np.zeros_like(ids,dtype=float)
        cz = np.zeros_like(ids,dtype=float)

        area = np.zeros_like(ids,dtype=float)
        vol = np.zeros_like(ids,dtype=float)
        dens = np.zeros_like(ids,dtype=float)
        nlabel = np.zeros_like(ids,dtype=float)

        polx = np.empty((0))
        poly = np.empty((0))
        polz = np.empty((0))
        polarea = np.empty((0))
        clusterid = np.empty((0),int)
        polid = np.empty((0),int)
        polstdz = np.empty((0))

        polidcur = 1 # we keep our own list of poly ids since we can get multiple polygons per cluster
        for id in uids:
            roi = ids==id

            roix = inputLocalisations['x'][roi]
            roiy = inputLocalisations['y'][roi]
            roiz = inputLocalisations['z'][roi]
            roi3d=np.stack((roix,roiy,roiz),axis=1)

            arear,volr,nlabelr,polys = shape_measure(roi3d, self.alpha, useAlphaShapes=self.use_alpha_shapes)
            
            #alpha_shape = ashp.alphashape(roi3d[:,0:2], self.alpha)
            #alpha_vol = ashp.alphashape(roi3d, self.alpha)

            sxr = np.std(roix)
            syr = np.std(roiy)
            szr = np.std(roiz)
            sx[roi] = sxr
            sy[roi] = syr
            sz[roi] = szr

            cxr = np.mean(roix)
            cyr = np.mean(roiy)
            czr = np.mean(roiz)
            #stdy.append(syr)
            #stdz.append(szr)
            cx[roi] = cxr
            cy[roi] = cyr
            cz[roi] = czr

            area[roi] = arear
            vol[roi] = volr
            nlabel[roi] = nlabelr
            dens[roi] = nlabelr/(arear/1e6)

            # some code to process the polygons
            for pol in polys:
                polxr = pol[:,0]
                polx = np.append(polx,polxr)
                poly = np.append(poly,pol[:,1])
                polz = np.append(polz,np.full_like(polxr,czr))
                polid = np.append(polid,np.full_like(polxr,polidcur,dtype=int))
                clusterid = np.append(clusterid,np.full_like(polxr,id,dtype=int))
                polarea = np.append(polarea,np.full_like(polxr,arear))
                polstdz = np.append(polstdz,np.full_like(polxr,szr))
                polidcur += 1

        mapped_ds = tabular.MappingFilter(inputLocalisations)
        mapped_ds.addColumn('clst_cx', cx)
        mapped_ds.addColumn('clst_cy', cy)
        mapped_ds.addColumn('clst_cz', cz)

        mapped_ds.addColumn('clst_stdx', sx)
        mapped_ds.addColumn('clst_stdy', sy)
        mapped_ds.addColumn('clst_stdz', sz)

        mapped_ds.addColumn('clst_area', area)
        mapped_ds.addColumn('clst_vol', vol)
        mapped_ds.addColumn('clst_size', nlabel)
        mapped_ds.addColumn('clst_density', dens)

        from PYME.IO.tabular import DictSource
        dspoly = DictSource(dict(x=polx,
                                 y=poly,
                                 z=polz,
                                 polyIndex=polid,
                                 polyArea=polarea,
                                 clusterID=clusterid,
                                 stdz=polstdz))
        
        return dict(outputName=mapped_ds,outputShapes=dspoly)
