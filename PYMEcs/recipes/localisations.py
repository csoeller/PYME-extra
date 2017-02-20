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
