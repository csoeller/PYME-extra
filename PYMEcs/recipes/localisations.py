from PYME.recipes.base import register_module, ModuleBase, Filter
from PYME.recipes.traits import Input, Output, Float, Enum, CStr, Bool, Int, List, DictStrStr, DictStrList, ListFloat, ListStr, FileOrURI

import numpy as np
import pandas as pd
from PYME.IO import tabular

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
                    Item('clumpColumnName'),]


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
    
@register_module('OrigamiSiteTrack')
class OrigamiSiteTrack(ModuleBase):
    """Recipe module aimed at analysing oregami MINFLUX datasets. More docs to be added.""" ### fix this!!!!
    inputClusters = Input('siteclusters')
    inputSites = Input('sites')
    outputName = Output('corrected_siteclusters')
    labelKey = CStr('dyeID') # should really be siteID
    smoothingBinWidthsSeconds = Float(200,label='temporal binning (s)',
                                      desc="parameter that sets the temporal binning (in units of s) when" +
                                      " estimating a smoothed drift trajectory from the data")
    savgolWindowLength = Int(10,label='savgol filter window length')
    savgolPolyorder = Int(6,label='savgol filter polynomial order')
    binnedStatistic = Enum(['mean','median'],desc="statistic for smoothing when using binned_statistic function on data")
    
    def run(self, inputClusters, inputSites):
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

        # we should check that with this choice of nbins we get no issues! (i.e. too few counts in some bins)
        c_xsite = smoothed_site_func(t,xsite,bins=nbins,statistic=self.binnedStatistic,
                                     sgwindow_length=self.savgolWindowLength,sgpolyorder=self.savgolPolyorder)
        c_ysite = smoothed_site_func(t,ysite,bins=nbins,statistic=self.binnedStatistic,
                                     sgwindow_length=self.savgolWindowLength,sgpolyorder=self.savgolPolyorder)
        c_zsite = smoothed_site_func(t,zsite,bins=nbins,statistic=self.binnedStatistic,
                                     sgwindow_length=self.savgolWindowLength,sgpolyorder=self.savgolPolyorder)

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
        
        return mapped_ds
