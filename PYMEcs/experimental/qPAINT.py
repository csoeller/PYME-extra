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

def Warn(parent, message, caption = 'Warning!'):
    dlg = wx.MessageDialog(parent, message, caption, wx.OK | wx.ICON_WARNING)
    dlg.ShowModal()
    dlg.Destroy()

class QPCalc:
    """

    """
    def __init__(self, visFr):
        self.visFr = visFr
        self.pipeline = visFr.pipeline
        self.qpMeasurements = {}

        visFr.AddMenuItem('qPAINT', "Darktime distribution of selected events",self.OnDarkT)
        visFr.AddMenuItem('qPAINT', "From Image - Set driftpars",self.OnSetDriftPars)
        visFr.AddMenuItem('qPAINT', "From Image - Set objectIDs",self.OnGetIDsfromImage)
        visFr.AddMenuItem('qPAINT', "From Image - Get Areas by ID",self.OnAreaFromLabels)
        visFr.AddMenuItem('qPAINT', "Qindex - Measure object ID dark times",self.OnMeasureTau)
        visFr.AddMenuItem('qPAINT', "All in 1 go: select Image, set drift, IDs, measure qindex, areas",
                          self.OnSelectImgAndProcess)

        visFr.AddMenuItem('qPAINT', itemType='separator') #--------------------------
        visFr.AddMenuItem('qPAINT', "Multicolour - set timed species from image",self.OnTimedSpeciesFromImage)
        visFr.AddMenuItem('qPAINT', "Multicolour - qIndex by channel",self.OnChannelMeasureTau)
        visFr.AddMenuItem('qPAINT', "Multicolour - copy data source",self.OnCopyDS)
        visFr.AddMenuItem('qPAINT', "Multicolour - in 1 go: select Image, set drift, IDs, areas, species, qIndex by channel",
                          self.OnSelectImgAndProcessMulticol)

        visFr.AddMenuItem('qPAINT', itemType='separator') #--------------------------
        visFr.AddMenuItem('qPAINT', 'Points based: Set objectIDs by DBSCAN clumping', self.OnClumpObjects,
                          helpText='Calculate objectID using DBSCAN algorithm')
        visFr.AddMenuItem('qPAINT', "Points based: Measure object ID volumes (area if 2D) by convex hull",self.OnMeasureVol)

        visFr.AddMenuItem('qPAINT', itemType='separator') #--------------------------
        visFr.AddMenuItem('qPAINT', "Calibrate a qIndex (or any column)",self.OnQindexCalibrate)
        visFr.AddMenuItem('qPAINT', "Ratio two qIndices (or any two columns)",self.OnQindexRatio)

        visFr.AddMenuItem('qPAINT', itemType='separator') #--------------------------
        visFr.AddMenuItem('qPAINT', "Plot histogram of one data column",self.OnGeneralHist)
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

            pipeline.Rebuild()
            self.visFr.CreateFoldPanel() # to make, for example, new columns show up in filter columns
            

    def OnCopyDS(self, event=None):
        """

        """
        from PYMEcs.recipes.localisations import CopyMapped
        recipe = self.pipeline.recipe
        CM = CopyMapped(recipe, inputName=self.pipeline.selectedDataSourceKey,
                                        outputName='%s-copy' % self.pipeline.selectedDataSourceKey)
        if CM.configure_traits(kind='modal'):
            recipe.add_module(CM)
            recipe.execute()
            self.pipeline.selectDataSource(CM.outputName)

            self.visFr.RefreshView()
            self.visFr.CreateFoldPanel()


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
            # somewhat weired condition as sometimes TimeSpecies becomes None in the root namespce
            if img.mdh.getOrDefault('TimedSpecies',None) is not None:
                logger.debug('setting timed species from root level')
                timedSpecies = img.mdh['TimedSpecies']
            elif 'Source.TimedSpecies' in img.mdh.keys():
                logger.debug('setting timed species from source level')
                timedSpecies = img.mdh['Source.TimedSpecies']
                logger.debug('time species is %s' % repr(timedSpecies))
            else:
                return
            pipeline = self.visFr.pipeline
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
            self.qpMeasurements['area'] = namespace[MeasureIt.outputName]

            meas = self.qpMeasurements['area']
            # currently the scale of areas is in 1000 nm^2
            areas = meas['area'] * float(img.mdh['voxelsize.x'])*float(img.mdh['voxelsize.y'])*1e3

            from PYMEcs.Analysis import fitDarkTimes
            # add the area info to the qPAINT measures already available
            for chan in (self.pipeline.colourFilter.getColourChans() + ['Everything']):
                if chan in self.qpMeasurements.keys():
                    # here we copy the area over, other 2D measures could be similarly copied
                    self.qpMeasurements[chan].addNewColumnByID(meas['label'], 'area', areas)

            # make a new objAreas column for the pipeline
            #      1st we are making a lookup table
            #      for areas by label index in abyl
            labels = meas['label'].astype('int')
            abyl = np.zeros(labels.max()+1) # make it long enough to hold all labels present
            abyl[labels] = areas
            objIDs = self.pipeline['objectID'].astype('int')
            objAreas = np.zeros_like(objIDs,dtype='float')
            # now look up the area values by addressing into our lookup table abyl
            objAreas[objIDs>0] = abyl[objIDs[objIDs>0]]
            # enter the new column into the pipeline
            self.pipeline.addColumn('objArea',objAreas)

    def OnMeasureTau(self, event):
        from PYMEcs.Analysis import fitDarkTimes

        # chans = self.pipeline.colourFilter.getColourChans()

        ids = np.unique(self.pipeline['objectID'].astype('int'))
        idsnz = ids[ids > 0]
        pipeline = self.pipeline
        
        measure = fitDarkTimes.measureObjectsByID(self.pipeline, set(idsnz))
        self.qpMeasurements['Everything'] = measure

        # mapping from measurements column name to new pipeline column name to add
        colmapNames = {
            'tau1':'taudark',
            'NDarktimes':'NDarktimes',
            'qindex1':'qIndex',
            'tau1err':'taudark_error',
            'chisqr1':'taudark_chisq2'}
        
        newcols =  fitDarkTimes.retrieveMeasuresForIDs(measure,pipeline['objectID'],
                                                       columns=colmapNames.keys())       
        for sourceCol in colmapNames.keys():
            pipeline.addColumn(colmapNames[sourceCol],newcols[sourceCol])

        # pipeline.Rebuild()

    def OnChannelMeasureTau(self, event, chan=None):
        from PYMEcs.Analysis import fitDarkTimes

        if chan is None:
            chan = selectWithDialog(self.pipeline.colourFilter.getColourChans(), message='select channel')
        if chan is None:
            return

        pipeline = self.pipeline
        dispColor = self.pipeline.colourFilter.currentColour
        pipeline.colourFilter.setColour(chan)

        ids = np.unique(pipeline['objectID'].astype('int'))
        idsnz = ids[ids > 0]
        
        self.qpMeasurements[chan] = fitDarkTimes.measureObjectsByID(pipeline, set(idsnz))

        # switch back to all channels
        pipeline.colourFilter.setColour('Everything')
        colmapNames = {
            'tau1':'taudark_%s' % chan,
            'NDarktimes':'NDarktimes_%s' % chan,
            'qindex1':'qIndex_%s' % chan,
            'tau1err':'taudark_error_%s' % chan,
            'chisqr1':'taudark_chisq2_%s' % chan}
        
        newcols =  fitDarkTimes.retrieveMeasuresForIDs(self.qpMeasurements[chan],pipeline['objectID'],
                                                       columns=colmapNames.keys())       
        for sourceCol in colmapNames.keys():
            pipeline.addColumn(colmapNames[sourceCol],newcols[sourceCol])

        # restore original display settings
        pipeline.colourFilter.setColour(dispColor)

    def OnSelectImgAndProcess(self, event):
        from PYME.DSView import dsviewer
        selection = selectWithDialog(dsviewer.openViewers.keys())
        if selection is not None:
            img = dsviewer.openViewers[selection].image
            self.OnSetDriftPars(None,img=img)
            self.OnGetIDsfromImage(None,img=img)
            self.OnMeasureTau(None)
            # areas last so that qpMeasures get the area info 
            self.OnAreaFromLabels(None,img=img)

    def OnSelectImgAndProcessMulticol(self, event):
        from PYME.DSView import dsviewer
        selection = selectWithDialog(dsviewer.openViewers.keys())
        if selection is not None:
            img = dsviewer.openViewers[selection].image
            self.OnSetDriftPars(None,img=img)
            self.OnGetIDsfromImage(None,img=img)
            self.OnTimedSpeciesFromImage(None,img=img)
            for chan in self.pipeline.colourFilter.getColourChans():
                self.OnChannelMeasureTau(None,chan=chan)
            # areas last so that qpMeasures get the area info 
            self.OnAreaFromLabels(None,img=img)

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


    def OnGeneralHist(self, event):
        from PYMEcs.recipes import localisations
        from PYME.recipes.base import ModuleCollection

        rec = ModuleCollection()

        HistByID = localisations.HistByID(rec)
        rec.namespace = {HistByID.inputName: self.pipeline}
        if HistByID.configure_traits(kind='modal'):
            # we call this with the pipeline to allow filtering etc
            HistByID.execute(rec.namespace)


    def OnQindexCalibrate(self, event):
        from PYMEcs.recipes import localisations
        from PYME.recipes.base import ModuleCollection

        rec = ModuleCollection()
        QIScaler = localisations.QindexScale(rec)
        
        # we call this with the pipeline to allow filtering etc
        rec.namespace = {QIScaler.inputName: self.pipeline}
        if QIScaler.configure_traits(kind='modal'):
            QIScaler.execute(rec.namespace)
            self.pipeline.addColumn(QIScaler.newKey, rec.namespace[QIScaler.outputName][QIScaler.newKey])

    def OnQindexRatio(self, event):
        from PYMEcs.recipes import localisations
        from PYME.recipes.base import ModuleCollection

        rec = ModuleCollection()
        QIRatio = localisations.QindexRatio(rec)
        
        # we call this with the pipeline to allow filtering etc
        rec.namespace = {QIRatio.inputName: self.pipeline}
        if QIRatio.configure_traits(kind='modal'):
            QIRatio.execute(rec.namespace)
            self.pipeline.addColumn(QIRatio.qIndexRatio, rec.namespace[QIRatio.outputName][QIRatio.qIndexRatio])

    def OnScatterByID(self, event):
        from PYMEcs.recipes import localisations
        from PYME.recipes.base import ModuleCollection

        rec = ModuleCollection()
        ScatterbyID = localisations.ScatterbyID(rec)
        rec.namespace = {ScatterbyID.inputName: self.pipeline}
        if ScatterbyID.configure_traits(kind='modal'):
            # we call this with the pipeline to allow filtering etc
            ScatterbyID.execute(rec.namespace)

    def OnSaveMeasurements(self,event):
        from PYME.recipes import runRecipe
        import os.path
        
        if len(self.qpMeasurements.keys()) > 0:
            filename = wx.FileSelector('Save all measurements (select basename)...',
                                       wildcard="CSV files (*.csv)|*.csv|Excell files (*.xlsx)|*.xlsx|HDF5 files (*.hdf)|*.hdf", 
                                       flags = wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT)
                                   
            if not filename == '':
                base, ext = os.path.splitext(filename)
                for chan in self.qpMeasurements.keys():
                    runRecipe.saveOutput(self.qpMeasurements[chan], base + '_' + chan + ext)

    def OnDarkT(self,event):
        import StringIO
        from PYMEcs.Analysis import fitDarkTimes
        
        visFr = self.visFr
        pipeline = visFr.pipeline
        mdh = pipeline.mdh

        NTMIN = 5
        maxPts = 1e4
        t = pipeline['t']
        if len(t) > maxPts:
            Warn(None,'aborting darktime analysis: too many events, current max is %d' % maxPts)
            return
        x = pipeline['x']
        y = pipeline['y']

        # determine darktime from gaps and reject zeros (no real gaps) 
        dts = t[1:]-t[0:-1]-1
        dtg = dts[dts>0]
        nts = dtg.shape[0]

        if nts > NTMIN:
            # now make a cumulative histogram from these
            cumux,cumuy = fitDarkTimes.cumuhist(dtg)
            binctrs,hc,binctrsg,hcg = fitDarkTimes.cumuhistBinned(dtg)
            
            bbx = (x.min(),x.max())
            bby = (y.min(),y.max())
            voxx = 1e3*mdh['voxelsize.x']
            voxy = 1e3*mdh['voxelsize.y']
            bbszx = bbx[1]-bbx[0]
            bbszy = bby[1]-bby[0]

        # fit theoretical distributions
        popth,pcovh,popt,pcov = (None,None,None,None)
        if nts > NTMIN:
            from scipy.optimize import curve_fit

            idx = (np.abs(hcg - 0.63)).argmin()
            tauesth = binctrsg[idx]
            popth,pcovh,infodicth,errmsgh,ierrh = curve_fit(fitDarkTimes.cumuexpfit,binctrsg,hcg, p0=(tauesth),full_output=True)
            chisqredh = ((hcg - infodicth['fvec'])**2).sum()/(hcg.shape[0]-1)
            idx = (np.abs(cumuy - 0.63)).argmin()
            tauest = cumux[idx]
            popt,pcov,infodict,errmsg,ierr = curve_fit(fitDarkTimes.cumuexpfit,cumux,cumuy, p0=(tauest),full_output=True)
            chisqred = ((cumuy - infodict['fvec'])**2).sum()/(nts-1)

            import matplotlib.pyplot as plt
            # plot data and fitted curves
            plt.figure()
            plt.subplot(211)
            plt.plot(cumux,cumuy,'o')
            plt.plot(cumux,fitDarkTimes.cumuexpfit(cumux,popt[0]))
            plt.plot(binctrs,hc/float(nts),'o')
            plt.plot(binctrs,fitDarkTimes.cumuexpfit(binctrs,popth[0]))
            plt.ylim(-0.2,1.2)
            plt.subplot(212)
            plt.semilogx(cumux,cumuy,'o')
            plt.semilogx(cumux,fitDarkTimes.cumuexpfit(cumux,popt[0]))
            plt.semilogx(binctrs,hc/float(nts),'o')
            plt.semilogx(binctrs,fitDarkTimes.cumuexpfit(binctrs,popth[0]))
            plt.ylim(-0.2,1.2)
            plt.show()
            
            outstr = StringIO.StringIO()

            analysis = {
                'Nevents' : t.shape[0],
                'Ndarktimes' : nts,
                'filterKeys' : pipeline.filterKeys.copy(),
                'darktimes' : (popt[0],popth[0]),
                'darktimeErrors' : (np.sqrt(pcov[0][0]),np.sqrt(pcovh[0][0]))
            }

            #if not hasattr(self.visFr,'analysisrecord'):
            #    self.visFr.analysisrecord = []
            #    self.visFr.analysisrecord.append(analysis)

            print >>outstr, "events: %d, dark times: %d" % (t.shape[0],nts)
            print >>outstr, "region: %d x %d nm (%d x %d pixel)" % (bbszx,bbszy,bbszx/voxx,bbszy/voxy)
            print >>outstr, "centered at %d,%d (%d,%d pixels)" % (x.mean(),y.mean(),x.mean()/voxx,y.mean()/voxy)
            print >>outstr, "darktime: %.1f+-%d (%.1f+-%d) frames - chisqr %.2f (%.2f)" % (popt[0],np.sqrt(pcov[0][0]),
                                                                                           popth[0],np.sqrt(pcovh[0][0]),
                                                                                           chisqred,chisqredh)
            print >>outstr, "darktime: starting estimates: %.1f (%.1f)" % (tauest,tauesth)
            print >>outstr, "qunits: %.2f (%.2f), eunits: %.2f" % (100.0/popt[0], 100.0/popth[0],t.shape[0]/500.0)

            labelstr = outstr.getvalue()
            plt.annotate(labelstr, xy=(0.5, 0.1), xycoords='axes fraction',
                         fontsize=10)
        else:
            Warn(None, 'not enough data points (<%d)' % NTMIN, 'Error')


def Plug(visFr):
    """Plugs this module into the gui"""
    visFr.qPAINTCalc = QPCalc(visFr)

