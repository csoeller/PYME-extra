import matplotlib.pyplot as plt
import numpy as np

def plot_errors(pipeline):
    curds = pipeline.selectedDataSourceKey
    pipeline.selectDataSource('coalesced_nz')
    p = pipeline
    plt.figure()
    plt.subplot(221)
    if 'error_z' in pipeline.keys():
        plt.boxplot([p['error_x'],p['error_y'],p['error_z']],labels=['error_x','error_y','error_z'])
    else:
        plt.boxplot([p['error_x'],p['error_y']],labels=['error_x','error_y'])
    plt.ylabel('loc error - coalesced (nm)')
    pipeline.selectDataSource('with_clumps')
    plt.subplot(222)
    plt.boxplot([p['fbg']],labels=['background'])
    uids, idx = np.unique(p['clumpIndex'],return_index=True)
    plt.subplot(223)
    if 'error_z' in pipeline.keys():
        plt.boxplot([p['error_x'][idx],p['error_y'][idx],p['error_z'][idx]],
                    labels=['error_x','error_y','error_z'])
    else:
        plt.boxplot([p['error_x'][idx],p['error_y'][idx]],labels=['error_x','error_y'])
    plt.ylabel('loc error - raw (nm)')
    plt.subplot(224)
    bp_dict = plt.boxplot([p['nPhotons']],labels=['photons'])
    for line in bp_dict['medians']:
        # get position data for median line
        x, y = line.get_xydata()[0] # top of median line
        # overlay median value
        plt.text(x, y, '%.0f' % y,
                 horizontalalignment='right') # draw above, centered
    plt.tight_layout()
    pipeline.selectDataSource(curds)
    

def plot_cluster_analysis(pipeline, ds='dbscanClustered'):
    curds = pipeline.selectedDataSourceKey
    pipeline.selectDataSource(ds)
    p = pipeline
    uids, cts = np.unique(p['dbscanClumpID'], return_counts=True)
    ctsgt1 = cts[cts > 1.1]
    pipeline.selectDataSource(curds)
    plt.figure()
    plt.subplot(221)
    h = plt.hist(cts,bins=15)
    plt.xlabel('Cluster Size')
    plt.plot([np.mean(cts),np.mean(cts)],[0,h[0].max()])
    plt.plot([np.median(cts),np.median(cts)],[0,h[0].max()],'--')
    plt.subplot(222)
    h = plt.hist(ctsgt1,bins=15)
    plt.xlabel('Cluster Size (clusters > 1)')
    plt.plot([np.mean(ctsgt1),np.mean(ctsgt1)],[0,h[0].max()])
    plt.plot([np.median(ctsgt1),np.median(ctsgt1)],[0,h[0].max()],'--')
    plt.subplot(223)
    bp_dict = plt.boxplot([cts,ctsgt1],labels=['cluster size','clusters > 1'], showmeans=True)
    for line in bp_dict['means']:
        # get position data for median line
        x, y = line.get_xydata()[0] # top of median line
        # overlay median value
        plt.text(x, y, '%.1f' % y,
                 horizontalalignment='center') # draw above, centered

    plt.tight_layout()
    pipeline.selectDataSource(curds)


class MINFLUXanalyser():
    def __init__(self, visFr):
        self.visFr = visFr

        visFr.AddMenuItem('Experimental>MINFLUX', "Localisation Error analysis", self.OnErrorAnalysis)
        visFr.AddMenuItem('Experimental>MINFLUX', "Cluster sizes - 3D", self.OnCluster3D)
        visFr.AddMenuItem('Experimental>MINFLUX', "Cluster sizes - 2D", self.OnCluster2D)

    def OnErrorAnalysis(self, event):
        plot_errors(self.visFr.pipeline)

    def OnCluster3D(self, event):
        plot_cluster_analysis(self.visFr.pipeline, ds='dbscanClustered')

    def OnCluster2D(self, event):
        plot_cluster_analysis(self.visFr.pipeline, ds='dbscan2D')

        
def Plug(visFr):
    return MINFLUXanalyser(visFr)
