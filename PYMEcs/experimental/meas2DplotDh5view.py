class Meas2DPlotter:
    """

    """
    def __init__(self, dsviewer):
        self.dsviewer = dsviewer
        dsviewer.AddMenuItem('Experimental>Meas2D', 'Plot cluster measurements', self.OnPlotClust)
        dsviewer.AddMenuItem('Experimental>Meas2D', 'Plot per cluster colocalisation', self.OnClusterColoc)


    def OnPlotClust(self, event=None):
        nodata = False
        try:
            pipeline = self.dsviewer.pipeline
        except AttributeError:
            nodata = True
        else:
            if 'area' not in pipeline.keys():
                nodata = True
        if nodata:
            print('no area column found - returning')
            return
        
        mdh = self.dsviewer.image.mdh
        vx = 1e3*mdh['voxelsize.x']
        vy = 1e3*mdh['voxelsize.y']
        RyRsz = 900.0 # RyR footprint in nm
        
        ryrs = pipeline['area']*vx*vy/RyRsz
        ryrmean = ryrs.mean()
        
        import matplotlib.pyplot as plt
        # plot data and fitted curves
        plt.figure()
        plt.hist(ryrs, bins=20)
        plt.xlabel('RyRs')
        plt.ylabel('Frequency')
        plt.title('RyR size distribution (mean %.1f RyRs, %d clusters)' % (ryrmean,ryrs.shape[0]))
        plt.show()

    def OnClusterColoc(self, event=None):
        nodata = False
        try:
            pipeline = self.dsviewer.pipeline
        except AttributeError:
            nodata = True
        else:
            if 'area' not in pipeline.keys():
                nodata = True
            if 'mean_intensity' not in pipeline.keys():
                nodata = True
        if nodata:
            print('no area column found - returning')
            return
        
        mdh = self.dsviewer.image.mdh
        vx = 1e3*mdh['voxelsize.x']
        vy = 1e3*mdh['voxelsize.y']
        RyRsz = 900.0 # RyR footprint in nm
        ryrmin = 5 # minimal size of cluster to check
        
        ryrs = pipeline['area']*vx*vy/RyRsz
        ryrgtmin = ryrs > ryrmin

        import matplotlib.pyplot as plt
        import numpy as np
        plt.figure()
        plt.scatter(ryrs[ryrgtmin],pipeline['mean_intensity'][ryrgtmin],
                    facecolors='lightgray', edgecolors='black')
        plt.xlabel('RyRs')
        plt.ylabel('fractional coloc')
        plt.title('area fraction of IP3R in RyR clusters')
        plt.ylim(0,1)

        fracs1 = pipeline['mean_intensity'][ryrgtmin]
        meanclc = fracs1.mean()
        #print('mean fraction colocalising per cluster (RyRs>5): %.2f' % meanclc)

        positive1 = (fracs1 > 0.05).sum()/float(fracs1.shape[0])
        #print('fraction with positive staining: %.2f' % positive1)

        plt.figure()
        #mbins = np.arange(0,10,1.0)/10.0
        h0 = plt.hist(pipeline['mean_intensity'][ryrgtmin],bins=10)
        plt.title('mean coloc frac %.2f, cluster frac positive %.2f' % (meanclc,positive1))
        plt.xlabel('coloc fraction per cluster')
        plt.ylabel('frequency')
        plt.show()

        
def Plug(dsviewer):
    """Plugs this module into the gui"""
    dsviewer.meas2Dplt = Meas2DPlotter(dsviewer)
