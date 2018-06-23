class Meas2DPlotter:
    """

    """
    def __init__(self, dsviewer):
        self.dsviewer = dsviewer
        dsviewer.AddMenuItem('Experimental>Meas2D', 'Plot cluster measurements', self.OnPlotClust)

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

def Plug(dsviewer):
    """Plugs this module into the gui"""
    dsviewer.meas2Dplt = Meas2DPlotter(dsviewer)
