import numpy as np
import sys
from scipy import ndimage

import logging
logger = logging.getLogger(__file__)

def zshift(data,navg=100):
    nm = min(navg,data.shape[0])
    offset = data[0:nm].mean()
    return data - offset

def myfilter(d,width=11):
    return ndimage.gaussian_filter(d,width)

class ClusterTracker:
    """

    """
    def __init__(self, visFr):
        self.visFr = visFr
        self.pipeline = visFr.pipeline
        self.clusterTracks = []

        visFr.AddMenuItem('Experimental>Clusters', 'DBSCAN Clump', self.OnClumpDBSCAN,
                          helpText='Calculate ClumpIndex using DBSCAN algorithm')
        visFr.AddMenuItem('Experimental>Clusters', 'Track Clumps', self.OnTrackClumps,
                          helpText='extract the tracks for all clusters (clumps) that we found')
        visFr.AddMenuItem('Experimental>Clusters', 'Plot Tracks', self.OnShowTracks,
                          helpText='plot tracks of clusters (clumps) that we found')
        visFr.AddMenuItem('Experimental>Clusters', 'Plot Tracks Filtered', self.OnShowTracksFiltered,
                          helpText='plot filtered tracks of clusters (clumps) that we found')
        visFr.AddMenuItem('Experimental>Clusters', 'Clear Tracks', self.OnClearTracks,
                          helpText='clear tracks from memory')
        

    def OnClumpDBSCAN(self, event=None):
        """
        Runs sklearn DBSCAN clustering algorithm on pipeline filtered results using the GUI defined in the DBSCAN
        recipe module.

        Args are user defined through GUI
            searchRadius: search radius for clustering
            minClumpSize: number of points within eps required for a given point to be considered a core point

        This version is generally used to identify clumps identifying fiduciaries and therefore the
        default searchRadius is set fairly generous.
        """
        from PYME.recipes import localisations

        clumper = localisations.DBSCANClustering(minClumpSize = 50, searchRadius = 20.0)
        if clumper.configure_traits(kind='modal'):
            namespace = {clumper.inputName: self.pipeline}
            clumper.execute(namespace)

            self.pipeline.addColumn(clumper.outputName, namespace[clumper.outputName]['dbscanClumpID'])

    def OnTrackClumps(self, event=None):
        pipeline = self.pipeline
        from PYME.recipes import localisations
        clumper = localisations.DBSCANClustering()
        clusterID = clumper.outputName

        if not clusterID in pipeline.keys():
            logger.error('Cannot find column %s in pipeline' % clusterID)
            return

        clusterIDs = pipeline[clusterID].astype('int')
        idmax = max(clusterIDs)
       
        for id in range(1,idmax+1):
            thiscluster = clusterIDs == id 
            t_id = pipeline['t'][thiscluster]
            x_id = pipeline['x'][thiscluster]
            y_id = pipeline['y'][thiscluster]
            I = np.argsort(t_id)
            self.clusterTracks.append([t_id[I],x_id[I],y_id[I]])

    def OnShowTracks(self, event=None):
        import matplotlib.pyplot as plt
        if len(self.clusterTracks) > 0:
            plt.figure()
            for entry in self.clusterTracks:
                t,x,y = entry
                plt.plot(t,zshift(x))
            plt.figure()
            for entry in self.clusterTracks:
                t,x,y = entry
                plt.plot(t,zshift(y))

    def OnShowTracksFiltered(self, event=None):
        import matplotlib.pyplot as plt
        if len(self.clusterTracks) > 0:
            plt.figure()
            for entry in self.clusterTracks:
                t,x,y = entry
                plt.plot(t,myfilter(zshift(x)))
            plt.figure()
            for entry in self.clusterTracks:
                t,x,y = entry
                plt.plot(t,myfilter(zshift(y)))

    def OnClearTracks(self, event=None):
        self.clusterTracks = []

            
def Plug(visFr):
    """Plugs this module into the gui"""
    visFr.clusterTracker = ClusterTracker(visFr)

