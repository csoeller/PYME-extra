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

    outputName = Output('fiducialAdded')

    def execute(self, namespace):
        import PYMEcs.Analysis.trackFiducials as tfs

        inp = namespace[self.inputName]
        mapped = tabular.mappingFilter(inp)

        tracks = tfs.extractAverageTrajectory(inp, clumpRadiusVar = 'error_x',
                                              clumpRadiusMultiplier=self.radiusMultiplier,
                                              timeWindow=self.timeWindow, filter=self.filterMethod,
                                              filterScale=self.filterScale)

        # add tracks for all calculated dims to output
        for dim in tracks.keys():
            mapped.addColumn('fiducial_%s' % dim, tracks[dim])

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
    qindexValue = Float(1.0)
    NEquivalent = Float(1.0)

    def execute(self, namespace):
        inp = namespace[self.inputName]
        mapped = tabular.mappingFilter(inp)

        scaled = inp['qIndex']
        qigood = inp['qIndex'] > 0
        scaled[qigood] = inp['qIndex'][qigood] * self.NEquivalent / self.qindexValue

        mapped.addColumn('qIndexCalibrated', scaled)
        namespace[self.outputName] = mapped


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


# a version of David's module which we include here so that we can test/hack a few things
@register_module('DBSCANClustering')
class DBSCANClustering(ModuleBase):
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
    numberOfJobs = Int(max(multiprocessing.cpu_count()-1,1))
    
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
    def pipeline_view(self):
        from traitsui.api import View, Group, Item
        from PYME.ui.custom_traits_editors import CBEditor

        modname = 'mylabel'

        return View(Group(Item('parameter', editor=CBEditor(choices=self._key_choices)), label=modname))

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
    def pipeline_view(self):
        from traitsui.api import View, Group, Item
        from PYME.ui.custom_traits_editors import CBEditor

        modname = 'mylabel'

        return View(Group(Item('parameter', editor=CBEditor(choices=self._key_choices)), label=modname))

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
