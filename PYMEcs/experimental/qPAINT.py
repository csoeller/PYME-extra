import wx
import numpy as np
import sys
from scipy import ndimage

import logging
logger = logging.getLogger(__file__)

def uniqueByID(ids,column):
    uids, idx = np.unique(ids, return_index=True)
    ucol = column[idx]
    valid = uids > 0
    return uids[valid], ucol[valid]

def selectWithDialog(choices, message='select image from list', caption='Selection'):
    dlg = wx.SingleChoiceDialog(None, message, caption, choices, wx.CHOICEDLG_STYLE)
    if dlg.ShowModal() == wx.ID_OK:
        item = dlg.GetStringSelection()
    else:
        item = None
    dlg.Destroy()
    return item

class QPCalc:
    """

    """
    def __init__(self, visFr):
        self.visFr = visFr
        self.pipeline = visFr.pipeline
        self.qpMeasurements = None
        self.measurements = None

        visFr.AddMenuItem('qPAINT', "From Image - Set driftpars",self.OnSetDriftPars)
        visFr.AddMenuItem('qPAINT', "From Image - Set objIDs",self.OnGetIDsfromImage)
        visFr.AddMenuItem('qPAINT', "From Image - Get Areas by ID",self.OnAreaFromLabels)
        visFr.AddMenuItem('qPAINT', "Measure object ID dark times",self.OnMeasureTau)
        visFr.AddMenuItem('qPAINT', "All in 1 go: select Image, set drift, IDs, measure qindex, areas",
                          self.OnSelectImgAndProcess)

        visFr.AddMenuItem('qPAINT', itemType='separator') #--------------------------
        visFr.AddMenuItem('qPAINT', "Multicolour - set timed species from image",self.OnTimedSpeciesFromImage)
        visFr.AddMenuItem('qPAINT', "Multicolour - qindex by channel",self.OnChannelMeasureTau)

        visFr.AddMenuItem('qPAINT', itemType='separator') #--------------------------
        visFr.AddMenuItem('qPAINT', 'Points based: Set objIDs by DBSCAN clumping', self.OnClumpObjects,
                          helpText='Calculate objID using DBSCAN algorithm')
        visFr.AddMenuItem('qPAINT', "Points based: Measure object ID volumes (area if 2D) by convex hull",self.OnMeasureVol)

        visFr.AddMenuItem('qPAINT', itemType='separator') #--------------------------
        visFr.AddMenuItem('qPAINT', "Calibrate qIndex",self.OnQindexCalibrate)
        visFr.AddMenuItem('qPAINT', "qDensity - calculate qIndex to area ratio",self.OnQDensityCalc)

        visFr.AddMenuItem('qPAINT', itemType='separator') #--------------------------
        visFr.AddMenuItem('qPAINT', "Plot qIndex histogram",self.OnQindexHist)
        visFr.AddMenuItem('qPAINT', "Scatter plot by ID",self.OnScatterByID)
        visFr.AddMenuItem('qPAINT', "Save Measurements",self.OnSaveMeasurements)

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

    def OnGetIDsfromImage(self, event, img=None):
        from PYME.DSView import dsviewer

        visFr = self.visFr
        pipeline = visFr.pipeline

        if img is None:
            selection = selectWithDialog(dsviewer.openViewers.keys())
            if selection is not None:
                img = dsviewer.openViewers[selection].image

        if img is not None:            
            pixX = np.round((pipeline['x'] - img.imgBounds.x0 )/img.pixelSize).astype('i')
            pixY = np.round((pipeline['y'] - img.imgBounds.y0 )/img.pixelSize).astype('i')

            ind = (pixX < img.data.shape[0])*(pixY < img.data.shape[1])*(pixX >= 0)*(pixY >= 0)

            ids = np.zeros_like(pixX)
            #assume there is only one channel
            ids[ind] = img.data[:,:,:,0].squeeze()[pixX[ind], pixY[ind]].astype('i')

            numPerObject, b = np.histogram(ids, np.arange(ids.max() + 1.5) + .5)

            pipeline.addColumn('objectID', ids)
            pipeline.addColumn('NEvents', numPerObject[ids-1])

            #pipeline.Rebuild()

    def OnSetDriftPars(self, event, img=None):
        from PYME.DSView import dsviewer

        if img is None:
            selection = selectWithDialog(dsviewer.openViewers.keys())
            if selection is not None:
                img = dsviewer.openViewers[selection].image

        if img is not None:
            dpn = self.visFr.driftPane
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
            dpn.OnDriftExprChange(None)

    def OnTimedSpeciesFromImage(self, event, img=None):
        from PYME.DSView import dsviewer
        
        if img is None:
            selection = selectWithDialog(dsviewer.openViewers.keys())
            if selection is not None:
                img = dsviewer.openViewers[selection].image

        if img is not None:
            if 'TimedSpecies' in img.mdh.keys():
                pipeline = self.visFr.pipeline
                timedSpecies = img.mdh['TimedSpecies']
                if pipeline.selectedDataSource is not None:
                    pipeline.selectedDataSource.setMapping('ColourNorm', '1.0 + 0*t')
                    for species in timedSpecies.keys():
                        pipeline.selectedDataSource.setMapping('p_%s' % species,
                                                               '(t>= %d)*(t<%d)' % timedSpecies[species])
                        # FIXME: (1) append if exists, (2) assigning to mdh ok?
                        pipeline.selectedDataSource.mdh['TimedSpecies'] = timedSpecies
                    self.visFr.pipeline.Rebuild()
                    self.visFr.CreateFoldPanel()

                
    def OnAreaFromLabels(self, event, img=None):
        from PYME.DSView import dsviewer

        if img is None:
            selection = selectWithDialog(dsviewer.openViewers.keys())
            if selection is not None:
                img = dsviewer.openViewers[selection].image

        if img is not None:
            from PYME.recipes.measurement import Measure2D

            # execute the measure2D module as a mini-recipe
            MeasureIt = Measure2D(measureContour=False)
            namespace = {MeasureIt.inputLabels: img, MeasureIt.inputIntensity: img}
            MeasureIt.execute(namespace)
            # save the measurements for further use
            self.measurements = namespace[MeasureIt.outputName]
            # make a new objAreas column for the pipeline
            meas = self.measurements
            # currently the scale of areas is in 1000 nm^2
            areas = meas['area'] * float(img.mdh['voxelsize.x'])*float(img.mdh['voxelsize.y'])*1e3
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

    def OnChannelMeasureTau(self, event):
        from PYMEcs.Analysis import fitDarkTimes

        chan = selectWithDialog(self.pipeline.colourFilter.getColourChans())
        if chan is None:
            return

        pipeline = self.pipeline
        dispColor = self.pipeline.colourFilter.currentColour
        pipeline.colourFilter.setColour(chan)

        ids = np.unique(pipeline['objectID'].astype('int'))
        
        self.qpMeasurements, tau1, qidx, ndt = fitDarkTimes.measureObjectsByID(pipeline, set(ids))
        pipeline.addColumn('taudark_%s' % chan,tau1)
        pipeline.addColumn('NDarktimes_%s' % chan,ndt)
        pipeline.addColumn('qIndex_%s' % chan,qidx)
        # restore original display settings
        pipeline.colourFilter.setColour(dispColor)

    def OnSelectImgAndProcess(self, event):
        from PYME.DSView import dsviewer
        selection = selectWithDialog(dsviewer.openViewers.keys())
        if selection is not None:
            img = dsviewer.openViewers[selection].image
            self.OnSetDriftPars(None,img=img)
            self.OnGetIDsfromImage(None,img=img)
            self.OnAreaFromLabels(None,img=img)
            self.OnMeasureTau(None)


    def OnMeasureVol(self, event):
        from PYMEcs.recipes import localisations
        # fixme: this is currently 2D only but allows straightforward extension to 3D
        VolMeasurer = localisations.ObjectVolume()
        if VolMeasurer.configure_traits(kind='modal'):
            # we call this with the pipeline to allow filtering etc
            namespace = {VolMeasurer.inputName: self.pipeline}
            VolMeasurer.execute(namespace)
            # FIXME: scaling factor is correct for 2D only
            self.pipeline.addColumn('volume', namespace[VolMeasurer.outputName]['volume']/1e3) # area in 1000 nm^2
        
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
        

    # superseded by OnScatterByID
    def OnPlotQIvsVol(self, event):
        if 'qIndexCalibrated' in self.pipeline.keys():
            qi = self.pipeline['qIndexCalibrated']
        else:
            qi = self.pipeline['qIndex']
        if 'objArea' in self.pipeline.keys():
            vols = self.pipeline['objArea']
        else:
            vols = self.pipeline['volume']

        ids = self.pipeline['objectID']
        uids, uqi = uniqueByID(ids, qi)
        uids, uvols =  uniqueByID(ids, vols)

        import matplotlib.pyplot as plt
        plt.figure()        
        plt.scatter(uqi, uvols)
        if 'qIndexCalibrated' in self.pipeline.keys():
            plt.xlabel('Calibrated qIndex')
        else:
            plt.xlabel('qIndex (a.u.)')
        plt.ylabel('Area (1000 nm^2)')

    def OnScatterByID(self, event):
        from PYMEcs.recipes import localisations

        ScatterbyID = localisations.ScatterbyID(self)
        self.namespace = {ScatterbyID.inputName: self.pipeline}
        if ScatterbyID.configure_traits(kind='modal'):
            # we call this with the pipeline to allow filtering etc
            ScatterbyID.execute(self.namespace)

    def OnSaveMeasurements(self,event):
        from PYME.recipes import runRecipe

        if self.measurements is not None:
            filename = wx.FileSelector('Save Area measurements as ...', 
                                       wildcard="CSV files (*.csv)|*.csv|Excell files (*.xlsx)|*.xlsx|HDF5 files (*.hdf)|*.hdf", 
                                       flags = wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT)
                                   
            if not filename == '':
                runRecipe.saveOutput(self.measurements, filename)
        
        if self.qpMeasurements is not None:
            fdialog = wx.FileDialog(None, 'Save qPaint measurements ...',
                                    wildcard='Numpy array|*.npy|Tab formatted text|*.txt', style=wx.SAVE)
            succ = fdialog.ShowModal()
            if (succ == wx.ID_OK):
                outFilename = fdialog.GetPath().encode()

                if outFilename.endswith('.txt'):
                    of = open(outFilename, 'w')
                    of.write('\t'.join(self.qpMeasurements.dtype.names) + '\n')

                    for obj in self.qpMeasurements:
                        of.write('\t'.join([repr(v) for v in obj]) + '\n')
                    of.close()

                else:
                    np.save(outFilename, self.qpMeasurements) # we are assuming single channel here!

            
def Plug(visFr):
    """Plugs this module into the gui"""
    visFr.qPAINTCalc = QPCalc(visFr)

