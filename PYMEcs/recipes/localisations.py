from PYME.recipes.base import register_module, ModuleBase, Filter
from PYME.recipes.traits import Input, Output, Float, Enum, CStr, Bool, Int, List, DictStrStr, DictStrList, ListFloat, ListStr

import numpy as np
import pandas as pd
from PYME.IO import tabular

@register_module('FiducialTrack')
class FiducialTrack(ModuleBase):
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

        # shift dbscan labels up by one to match existing convention that a clumpID of 0 corresponds to unclumped
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
