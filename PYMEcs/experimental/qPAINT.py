import wx
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
        self.measurements = None

        visFr.AddMenuItem('qPAINT', 'Set objIDs by DBSCAN clumping', self.OnClumpObjects,
                          helpText='Calculate objID using DBSCAN algorithm')
        visFr.AddMenuItem('qPAINT', "Set driftpars from image",self.OnSetDriftPars)
        visFr.AddMenuItem('qPAINT', "Set objIDs from image",self.OnGetIDsfromImage)
        visFr.AddMenuItem('qPAINT', "Get Area from Label Image",self.OnAreaFromLabels)
        visFr.AddMenuItem('qPAINT', "Measure object ID dark times",self.OnMeasureTau)
        visFr.AddMenuItem('qPAINT', "Measure object ID volumes (area for 2D) from convex hull",self.OnMeasureVol)
        visFr.AddMenuItem('qPAINT', "Calibrate qIndex",self.OnQindexCalibrate)
        visFr.AddMenuItem('qPAINT', "Plot qIndex histogram",self.OnQindexHist)
        visFr.AddMenuItem('qPAINT', "Plot qIndex vs Area/Volume",self.OnPlotQIvsVol)
        visFr.AddMenuItem('qPAINT', "qDensity - qIndex to area ratio",self.OnQDensityCalc)

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

    def OnGetIDsfromImage(self, event):
        from PYME.DSView import dsviewer

        visFr = self.visFr
        pipeline = visFr.pipeline

        dlg = wx.SingleChoiceDialog(
                None, 'choose the image which contains labels', 'Use Segmentation',
                dsviewer.openViewers.keys(),
                wx.CHOICEDLG_STYLE
                )

        if dlg.ShowModal() == wx.ID_OK:
            img = dsviewer.openViewers[dlg.GetStringSelection()].image
            
            pixX = np.round((pipeline.filter['x'] - img.imgBounds.x0 )/img.pixelSize).astype('i')
            pixY = np.round((pipeline.filter['y'] - img.imgBounds.y0 )/img.pixelSize).astype('i')

            ind = (pixX < img.data.shape[0])*(pixY < img.data.shape[1])*(pixX >= 0)*(pixY >= 0)

            ids = np.zeros_like(pixX)
            #assume there is only one channel
            ids[ind] = img.data[:,:,:,0].squeeze()[pixX[ind], pixY[ind]].astype('i')

            numPerObject, b = np.histogram(ids, np.arange(ids.max() + 1.5) + .5)

            pipeline.addColumn('objectID', ids)
            pipeline.addColumn('NEvents', numPerObject[ids-1])

            pipeline.Rebuild()

        dlg.Destroy()

    def OnSetDriftPars(self, event):
        from PYME.DSView import dsviewer

        dlg = wx.SingleChoiceDialog(
                None, 'choose the image which contains drift info', 'Use Segmentation',
                dsviewer.openViewers.keys(),
                wx.CHOICEDLG_STYLE
                )
        if dlg.ShowModal() == wx.ID_OK:
            dpn = self.visFr.driftPane
            img = dsviewer.openViewers[dlg.GetStringSelection()].image
            dpn.tXExpr.SetValue(img.mdh['DriftCorrection.ExprX'])
            dpn.tYExpr.SetValue(img.mdh['DriftCorrection.ExprY'])
            dpn.tZExpr.SetValue(img.mdh['DriftCorrection.ExprZ'])
            dpn.OnDriftExprChange(None)
            destp = dpn.dp.driftCorrParams
            srcp = img.mdh['DriftCorrection.Parameters']
            for key in destp.keys():
                if key.startswith(('a','b')):
                    destp[key] = srcp[key]
            dpn.OnDriftApply(None)

    def OnAreaFromLabels(self, event):
        from PYME.DSView import dsviewer

        dlg = wx.SingleChoiceDialog(
                None, 'choose the image which contains drift info', 'Use Segmentation',
                dsviewer.openViewers.keys(),
                wx.CHOICEDLG_STYLE
                )
        if dlg.ShowModal() == wx.ID_OK:
            from PYME.recipes.measurement import Measure2D

            img = dsviewer.openViewers[dlg.GetStringSelection()].image
            # execute the measure2D module as a mini-recipe
            MeasureIt = Measure2D(measureContour=False)
            namespace = {MeasureIt.inputLabels: img, MeasureIt.inputIntensity: img}
            MeasureIt.execute(namespace)
            # save the measurements for further use
            self.measurements = namespace[MeasureIt.outputName]
            # make a new objAreas column for the pipeline
            meas = self.measurements
            areas = meas['area'] * float(img.mdh['voxelsize.x'])*float(img.mdh['voxelsize.y'])
            labels = meas['label'].astype('int')
            abyl = np.zeros(labels.max()+1)
            abyl[labels] = areas
            objIDs = self.pipeline['objectID'].astype('int')
            objAreas = np.zeros_like(objIDs,dtype='float')
            objAreas[objIDs>0] = abyl[objIDs[objIDs>0]]
            # enter the new column into the pipeline
            self.pipeline.addColumn('objArea',objAreas)

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
        
    def OnQDensityCalc(self, event):
        if 'qIndexCalibrated' in self.pipeline.keys():
            qi = self.pipeline['qIndexCalibrated']
        else:
            qi = self.pipeline['qIndex']
        
        objAreas = self.pipeline['objArea']
        valid = (qi > 0)*(objAreas > 0)
        objDensity = np.zeros_like(objAreas)
        objDensity[valid] = qi[valid]/objAreas[valid]
        self.pipeline.addColumn('objDensity',objDensity)

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

