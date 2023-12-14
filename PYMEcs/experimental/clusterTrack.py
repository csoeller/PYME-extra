import numpy as np
import sys
from scipy import ndimage
import matplotlib.pyplot as plt

import logging
logger = logging.getLogger(__file__)

from traits.api import HasTraits, Str, Int, CStr, List, Enum, Float
from traitsui.api import View, Item, Group
from traitsui.menu import OKButton, CancelButton, OKCancelButtons

class FilterChoice(HasTraits):
    windowSize = Int(11)
    filterType = Enum(['Gaussian','Median'])
    funcmap = {
        'Gaussian' : ndimage.gaussian_filter1d,
        'Median'   : ndimage.median_filter}
    
    def get_filter(self):
        
        def filterfunc(data):
            return self.funcmap[self.filterType](data,self.windowSize)

        return filterfunc
    
class GetTime(HasTraits):
    alignmentTime = Int(0)
    averagePeriod = Int(50)

def zeroshift(t,data,navg=50, alignmentTime=0):
    ti,idx = np.unique(t.astype('int'),return_index=True)
    di = data[idx]
    if alignmentTime<0:
        alignmentTime = 0
    nmin = min(alignmentTime,di.shape[0])
    nmax = min(alignmentTime+navg,di.shape[0])
    offset = di[nmin:nmax].mean()
    return data - offset

def atLeastRange(rmin=-20,rmax=-20):
    ymin,ymax = plt.ylim()
    if ymin > rmin:
        plt.ylim(ymin=rmin)
    if ymax < rmax:
        plt.ylim(ymax=rmax)

class ClusterTracker:
    """

    """
    def __init__(self, visFr):
        self.visFr = visFr
        self.pipeline = visFr.pipeline
        self.clusterTracks = []
        self.alignmentTime = 0
        self.averagePeriod = 50

        visFr.AddMenuItem('Experimental>Deprecated>Clusters', 'DBSCAN Clump', self.OnClumpDBSCAN,
                          helpText='Calculate ClumpIndex using DBSCAN algorithm')
        visFr.AddMenuItem('Experimental>Deprecated>Clusters', 'Track Clumps', self.OnTrackClumps,
                          helpText='extract the tracks for all clusters (clumps) that we found')
        visFr.AddMenuItem('Experimental>Deprecated>Clusters', 'Plot Tracks', self.OnShowTracks,
                          helpText='plot tracks of clusters (clumps) that we found')
        visFr.AddMenuItem('Experimental>Deprecated>Clusters', 'Plot Tracks Filtered', self.OnShowTracksFiltered,
                          helpText='plot filtered tracks of clusters (clumps) that we found')
        visFr.AddMenuItem('Experimental>Deprecated>Clusters', 'Clear Tracks', self.OnClearTracks,
                          helpText='clear tracks from memory')
        visFr.AddMenuItem('Experimental>Deprecated>Clusters', 'Set Alignment Time', self.OnSetAlignmentTime,
                          helpText='set alignment time')
        

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
        from PYMEcs.recipes import localisations

        clumper = localisations.DBSCANClustering2(minClumpSize = 50, searchRadius = 20.0)
        if clumper.configure_traits(kind='modal'):
            namespace = {clumper.inputName: self.pipeline}
            clumper.execute(namespace)

            self.pipeline.addColumn(clumper.outputName, namespace[clumper.outputName]['dbscanClumpID'])

    def OnTrackClumps(self, event=None):
        pipeline = self.pipeline
        from PYMEcs.recipes import localisations
        clumper = localisations.DBSCANClustering2()
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
            z_id = pipeline['z'][thiscluster]
            I = np.argsort(t_id)
            self.clusterTracks.append([t_id[I],x_id[I],y_id[I],z_id[I]])


    def OnShowTracks(self, event=None):
        import matplotlib.pyplot as plt
        if len(self.clusterTracks) > 0:
            navg = self.averagePeriod
            atime = self.alignmentTime
            plt.figure()
            for entry in self.clusterTracks:
                t,x,y,z = entry
                plt.plot(t,zeroshift(t,x, navg, atime))
            plt.title('x tracks')
            atLeastRange(-20,20)
            plt.figure()
            for entry in self.clusterTracks:
                t,x,y,z = entry
                plt.plot(t,zeroshift(t,y, navg, atime))
            plt.title('y tracks')
            atLeastRange(-20,20)
            plt.figure()
            for entry in self.clusterTracks:
                t,x,y,z = entry
                plt.plot(t,zeroshift(t,z, navg, atime))
            plt.title('z tracks')
            atLeastRange(-20,20)

    def OnShowTracksFiltered(self, event=None):
        import matplotlib.pyplot as plt
        fc = FilterChoice()
        if len(self.clusterTracks) > 0 and fc.configure_traits(kind='modal'):
            navg = self.averagePeriod
            atime = self.alignmentTime
            filterfunc = fc.get_filter()
            plt.figure()
            for entry in self.clusterTracks:
                t,x,y,z = entry
                plt.plot(t,filterfunc(zeroshift(t,x, navg, atime)))
            plt.title('x tracks')
            atLeastRange(-20,20)
            plt.figure()
            for entry in self.clusterTracks:
                t,x,y,z = entry
                plt.plot(t,filterfunc(zeroshift(t,y, navg, atime)))
            plt.title('y tracks')
            atLeastRange(-20,20)
            plt.figure()
            for entry in self.clusterTracks:
                t,x,y,z = entry
                plt.plot(t,filterfunc(zeroshift(t,z, navg, atime)))
            plt.title('z tracks')
            atLeastRange(-20,20)

    def OnClearTracks(self, event=None):
        self.clusterTracks = []

    def OnSetAlignmentTime(self, event=None):
        gtime = GetTime()
        if gtime.configure_traits(kind='modal'):
            self.alignmentTime = gtime.alignmentTime
            self.averagePeriod = gtime.averagePeriod
            
def Plug(visFr):
    """Plugs this module into the gui"""
    visFr.clusterTracker = ClusterTracker(visFr)

