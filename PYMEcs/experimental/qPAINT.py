import numpy as np
import sys
from scipy import ndimage

import logging
logger = logging.getLogger(__file__)

class QPCalc:
    """

    """
    def __init__(self, visFr):
        self.visFr = visFr
        self.pipeline = visFr.pipeline
        self.qpMeasurements = None

        visFr.AddMenuItem('Experimental>qPAINT', 'objID by DBSCAN clumping', self.OnClumpObjects,
                          helpText='Calculate objID using DBSCAN algorithm')
        visFr.AddMenuItem('Experimental>qPAINT', "Measure nonzero object ID dark times",self.OnMeasureTau)
        visFr.AddMenuItem('Experimental>qPAINT', "Measure nonzero object ID volumes (area for 2D data)",self.OnMeasureVol)
        visFr.AddMenuItem('Experimental>qPAINT', "Calibrate qIndex",self.OnQindexCalibrate)
        visFr.AddMenuItem('Experimental>qPAINT', "Plot qIndex histogram",self.OnQindexHist)
        visFr.AddMenuItem('Experimental>qPAINT', "Plot qIndex vs Area/Volume",self.OnPlotQIvsVol)

    def OnClumpObjects(self, event=None):
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

        clumper = localisations.DBSCANClustering(minClumpSize = 20, searchRadius = 20.0)
        if clumper.configure_traits(kind='modal'):
            namespace = {clumper.inputName: self.pipeline}
            clumper.execute(namespace)

            self.pipeline.addColumn('objectID', namespace[clumper.outputName]['dbscanClumpID'])

    def OnMeasureTau(self, event):
        from PYMEcs.Analysis import fitDarkTimes

        # chans = self.pipeline.colourFilter.getColourChans()

        ids = np.unique(self.pipeline['objectID'].astype('int'))

        pipeline = self.pipeline
        
        self.qpMeasurements, tau1, qidx, ndt = fitDarkTimes.measureObjectsByID(self.pipeline, set(ids))
        pipeline.addColumn('taudark',tau1)
        pipeline.addColumn('NDarktimes',ndt)
        pipeline.addColumn('qIndex',qidx)

        # pipeline.Rebuild()

    def OnMeasureVol(self, event):
        from PYMEcs.recipes import localisations

        VolMeasurer = localisations.ObjectVolume()
        if VolMeasurer.configure_traits(kind='modal'):
            # we call this with the pipeline to allow filtering etc
            namespace = {VolMeasurer.inputName: self.pipeline}
            VolMeasurer.execute(namespace)

            self.pipeline.addColumn('volumes', namespace[VolMeasurer.outputName]['volumes'])
        
    def OnQindexHist(self, event):
        ids, idx = np.unique(self.pipeline['objectID'].astype('int'), return_index=True)
        try:
            qidx = self.pipeline['qIndexCalibrated'][idx]
            isCalibrated = True
        except:
            qidx = self.pipeline['qIndex'][idx]
            isCalibrated = False

        if qidx.shape[0] > 200:
            nbins = 50
        else:
            nbins = 20
            
        import matplotlib.pyplot as plt
        plt.figure()
        plt.hist(qidx[qidx > 0],nbins)
        if isCalibrated:
            plt.title('Calibrated QIndex from %d objects' % qidx.shape[0])
        else:
            plt.title('Uncalibrated QIndex from %d objects' % qidx.shape[0])
        plt.xlabel('qIndex')

    def OnQindexCalibrate(self, event):
        from PYMEcs.recipes import localisations

        QIScaler = localisations.QindexScale()
        if QIScaler.configure_traits(kind='modal'):
            # we call this with the pipeline to allow filtering etc
            namespace = {QIScaler.inputName: self.pipeline}
            QIScaler.execute(namespace)

            self.pipeline.addColumn('qIndexCalibrated', namespace[QIScaler.outputName]['qIndexCalibrated'])
        

    def OnPlotQIvsVol(self, event):
        if 'qIndexCalibrated' in self.pipeline.keys():
            qi = self.pipeline['qIndexCalibrated']
        else:
            qi = self.pipeline['qIndex']

        ids, idx = np.unique(self.pipeline['objectID'], return_index=True)
        qiu = qi[idx]
        vols = self.pipeline['volumes'][idx]
        
        import matplotlib.pyplot as plt
        plt.figure()        
        plt.scatter(qiu, vols)
        if 'qIndexCalibrated' in self.pipeline.keys():
            plt.xlabel('Calibrated qIndex')
        else:
            plt.xlabel('qIndex (a.u.)')
        plt.ylabel('Area (nm^2)')

def Plug(visFr):
    """Plugs this module into the gui"""
    visFr.qPAINTCalc = QPCalc(visFr)

