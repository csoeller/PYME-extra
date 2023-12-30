from __future__ import print_function

import wx
import numpy as np
import sys
from scipy import ndimage
from PYMEcs.misc.guiMsgBoxes import Warn
from PYME.recipes import tablefilters
from PYMEcs.Analysis import fitDarkTimes

import logging
logger = logging.getLogger(__file__)

from traits.api import HasTraits, Str, Int, CStr, List, Enum, Float, Bool
from traitsui.api import View, Item, Group
from traitsui.menu import OKButton, CancelButton, OKCancelButtons

class KeyChoice(HasTraits):
    clist = List([])
    Key = Enum(values='clist')

    traits_view = View(
        'Key',
        title = 'Select Measure',
        resizable = True,
        buttons = OKCancelButtons
    )

    def add_keys(self,chans):
        for chan in chans:
            if chan not in self.clist:
                self.clist.append(chan)


class myCChoice(HasTraits):
    clist = List([])
    RatioChannel1 = Enum(values='clist')
    RatioChannel2 = Enum(values='clist')
    Channel1Calibration = Float(1.0)
    Channel2Calibration = Float(1.0)

    traits_view = View(Group(Item(name = 'RatioChannel1'),
                             Item(name = 'Channel1Calibration'),
                             Item('_'),
                             Item(name = 'RatioChannel2'),
                             Item(name = 'Channel2Calibration'),
                             label = 'Select Channels and Calibration',
                             show_border = True),
                       buttons = OKCancelButtons)

    cal_view = View(Group(Item(name = 'Channel1Calibration'),
                          Item('_'),
                          Item(name = 'Channel2Calibration'),
                          label = 'Select Calibration',
                          show_border = True),
                    buttons = OKCancelButtons)

    one_view = View(Group(Item(name = 'RatioChannel1'),
                          Item(name = 'Channel1Calibration'),
                          label = 'Select Channel and Calibration',
                          show_border = True),
                    buttons = OKCancelButtons)

    def add_channels(self,chans):
        for chan in chans:
            if chan not in self.clist:
                self.clist.append(chan)

                
class TimedSpecies(HasTraits):
    Species1 = CStr()
    Species1FromTime = Float()
    Species1ToTime = Float()
    
    Species2 = CStr()
    Species2FromTime = Float()
    Species2ToTime = Float()

    Species3 = CStr()
    Species3FromTime = Float()
    Species3ToTime = Float()


    traits_view = View(Group(Item(name = 'Species1'),
                             Item(name = 'Species1FromTime'),
                             Item(name = 'Species1ToTime'),
                             Item('_'),
                             Item(name = 'Species2'),
                             Item(name = 'Species2FromTime'),
                             Item(name = 'Species2ToTime'),
                             Item('_'),
                             Item(name = 'Species3'),
                             Item(name = 'Species3FromTime'),
                             Item(name = 'Species3ToTime'),
                             label = 'Specify Timed Species',
                             show_border = True),
                       buttons = OKCancelButtons)

    def getSpeciesDescriptor(self):
        speclist = {}
        if self.Species1: # empty strings will be ignored
            speclist[self.Species1] = (self.Species1FromTime,
                                               self.Species1ToTime)
        if self.Species2: # empty strings will be ignored
            speclist[self.Species2] = (self.Species2FromTime,
                                               self.Species2ToTime)
        if self.Species3: # empty strings will be ignored
            speclist[self.Species3] = (self.Species3FromTime,
                                               self.Species3ToTime)

        logger.info('speclist is ' + repr(speclist))
        return speclist
                
def uniqueByID(ids,column):
    uids, idx = np.unique(ids, return_index=True)
    ucol = column[idx]
    valid = uids > 0
    return uids[valid], ucol[valid]

def selectWithDialog(choices, message='select image from list', caption='Selection'):
    dlg = wx.SingleChoiceDialog(None, message, caption, list(choices), wx.CHOICEDLG_STYLE)
    if dlg.ShowModal() == wx.ID_OK:
        item = dlg.GetStringSelection()
    else:
        item = None
    dlg.Destroy()
    return item

def unique_name(stem,names):
    if stem not in names:
        return stem
    for i in range(1,11):
        stem2 = "%s%d" % (stem,i)
        if stem2 not in names:
            return stem2

    return stem2 # here we just give up and accept a duplicate name
    

class QPCalc:
    """

    """
    def __init__(self, visFr):
        self.visFr = visFr
        self.pipeline = visFr.pipeline
        self.qpMeasurements = {}
        self.fitSettings = fitDarkTimes.FitSettings()
        self.useTau = 2 # we change to using the proper histogram for fitting
        if self.useTau == 1:
            self.tausrc = {
                'tau'    : 'tau1',
                'tauerr' : 'tau1err',
                'chisq'  : 'chisqr1'}
        elif self.useTau == 2:
            self.tausrc = {
                'tau'    : 'tau2',
                'tauerr' : 'tau2err',
                'chisq'  : 'chisqr2'}
        else:
            raise RuntimeError("Invalid useTau mode %d (must be 1 or 2)" % self.useTau)

        
        visFr.AddMenuItem('Experimental>qPAINT', "From Image - Set driftpars",self.OnSetDriftPars)
        visFr.AddMenuItem('Experimental>qPAINT', "From Image - Set ROI clipping",self.OnClipFromImage)
        visFr.AddMenuItem('Experimental>qPAINT', "From Image - Set objectIDs",self.OnGetIDsfromImage)
        visFr.AddMenuItem('Experimental>qPAINT', "From Image - Get Areas by ID",self.OnAreaFromLabels)
        visFr.AddMenuItem('Experimental>qPAINT', "Qindex - Measure object ID dark times",self.OnMeasureTau2)
        visFr.AddMenuItem('Experimental>qPAINT', "All in 1 go: select Image, set drift, IDs, measure qindex, areas",
                          self.OnSelectImgAndProcess)

        visFr.AddMenuItem('Experimental>qPAINT', itemType='separator') #--------------------------
        visFr.AddMenuItem('Experimental>qPAINT', "Multicolour - set timed species by Dialog",self.OnTimedSpecies)
        visFr.AddMenuItem('Experimental>qPAINT', "Multicolour - set timed species from image",self.OnTimedSpeciesFromImage)
        visFr.AddMenuItem('Experimental>qPAINT', "Multicolour - qIndex by channel",self.OnChannelMeasureTau)
        visFr.AddMenuItem('Experimental>qPAINT', "Multicolour - Merge Channel Measures", self.OnMergeChannelMeasures)
        visFr.AddMenuItem('Experimental>qPAINT', "Multicolour - Calculate Channel Ratios", self.OnChannelMeasurementRatios)
        visFr.AddMenuItem('Experimental>qPAINT', "Multicolour - calibrate channel qIndices (Pseudo test function)",self.OnChannelCalibrate)
        visFr.AddMenuItem('Experimental>qPAINT', "Multicolour - in 1 go: select Image, set drift, IDs, areas, species, qIndex & merge & ratio",
                          self.OnSelectImgAndProcessMulticol)
        
        visFr.AddMenuItem('Experimental>qPAINT', itemType='separator') #--------------------------
        visFr.AddMenuItem('Experimental>qPAINT', 'Points based: Set objectIDs by DBSCAN clumping', self.OnClumpObjects,
                          helpText='Calculate objectID using DBSCAN algorithm')
        visFr.AddMenuItem('Experimental>qPAINT', "Points based: Measure object ID volumes (area if 2D) by convex hull",self.OnMeasureVol)

        visFr.AddMenuItem('Experimental>qPAINT', itemType='separator') #--------------------------
        visFr.AddMenuItem('Experimental>qPAINT', "Darktime distribution of selected events\tCtrl+T",self.OnDarkT)
        visFr.AddMenuItem('Experimental>qPAINT', "Select Fit Settings",self.OnFitSettings)
        visFr.AddMenuItem('Experimental>qPAINT', "Calibrate a qIndex (or any column)",self.OnQindexCalibrate)
        visFr.AddMenuItem('Experimental>qPAINT', "Ratio two qIndices (or any two columns)",self.OnQindexRatio)

        visFr.AddMenuItem('Experimental>qPAINT', itemType='separator') #--------------------------
        visFr.AddMenuItem('Experimental>qPAINT', "Plot histogram of one data column",self.OnGeneralHist)
        visFr.AddMenuItem('Experimental>qPAINT', "Scatter plot by ID",self.OnScatterByID)
        visFr.AddMenuItem('Experimental>qPAINT', "Image of qPAINT measure from label image",self.OnLabelLookupByID)
        visFr.AddMenuItem('Experimental>qPAINT', "Load Measurements",self.OnLoadMeasurements)
        visFr.AddMenuItem('Experimental>qPAINT', "Save Measurements",self.OnSaveMeasurements)

        
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

        pipeline = self.pipeline

        if img is None:
            selection = selectWithDialog(dsviewer.openViewers.keys())
            if selection is not None:
                img = dsviewer.openViewers[selection].image

        if img is not None:
            with_ids = unique_name('with_ids',pipeline.dataSources.keys())
            valid_ids = unique_name('valid_ids',pipeline.dataSources.keys())
            recipe = pipeline.recipe
            mapp = tablefilters.Mapping(recipe,inputName=pipeline.selectedDataSourceKey,
                                        outputName=with_ids)
            recipe.add_module(mapp)
            recipe.execute()

            withIDs = recipe.namespace[with_ids]
           
            pixX = np.round((pipeline['x'] - img.imgBounds.x0 )/img.pixelSize).astype('i')
            pixY = np.round((pipeline['y'] - img.imgBounds.y0 )/img.pixelSize).astype('i')

            ind = (pixX < img.data.shape[0])*(pixY < img.data.shape[1])*(pixX >= 0)*(pixY >= 0)

            ids = np.zeros_like(pixX)
            #assume there is only one channel
            ids[ind] = img.data[:,:,:,0].squeeze()[pixX[ind], pixY[ind]].astype('i')

            numPerObject, b = np.histogram(ids, np.arange(ids.max() + 1.5) + .5)

            withIDs.addColumn('objectID', ids)
            withIDs.addColumn('NEvents', numPerObject[ids-1])

            recipe.add_module(tablefilters.FilterTable(recipe,inputName=with_ids, outputName=valid_ids,
                                                         filters={'objectID' : [.5, sys.maxsize]}))
            recipe.execute()

            pipeline.selectDataSource(valid_ids)
            pipeline.Rebuild()

            
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

    def OnTimedSpecies(self, event):
        # a reworked version that uses a simpler TraitsUI interface
        timedSpecies = TimedSpecies()
        if timedSpecies.configure_traits(kind='modal'):
            from PYME.LMVis import renderers
            speclist = timedSpecies.getSpeciesDescriptor()
            if len(speclist.keys())>0:
                # not sure if this should be Source.TimedSpecies or just TimedSpecies 
                #renderers.renderMetadataProviders.append(lambda mdh:
                #                                         mdh.setEntry('Source.TimedSpecies', speclist))
            
                pipeline = self.visFr.pipeline
                pipeline.mdh.setEntry('TimedSpecies', speclist)
                if pipeline.selectedDataSource is not None:
                    pipeline.selectedDataSource.setMapping('ColourNorm', '1.0 + 0*t')
                    for species in speclist.keys():
                        pipeline.selectedDataSource.setMapping('p_%s' % species,
                                                               '(t>= %d)*(t<%d)' % speclist[species])
                    pipeline.Rebuild()
                    # self.visFr.RegenFilter()
                    self.visFr.CreateFoldPanel()


    def OnClipFromImage(self, event, img=None):
        from PYME.DSView import dsviewer
        
        if img is None:
            selection = selectWithDialog(list(dsviewer.openViewers.keys()))
            if selection is not None:
                img = dsviewer.openViewers[selection].image

        if img is not None:
            if img.mdh.getOrDefault('Filter.Keys',None) is None:
                logger.debug('no Filter.Keys in image metadata')
                return
            try:
                self.pipeline.filterKeys['x'] = img.mdh['Filter.Keys']['x']
            except:
                logger.debug('cannot read or set filterKeys x-component')
                return
            try:
                self.pipeline.filterKeys['y'] = img.mdh['Filter.Keys']['y']
            except:
                logger.debug('cannot read or set filterKeys y-component')
                return
                
            self.pipeline.Rebuild()
            #self.visFr.RegenFilter()
            self.visFr.CreateFoldPanel()


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
                # there are a couple different formats for timedSpecies, so check which one it is
                try:
                    specs = timedSpecies.keys()
                    isdict = True
                except AttributeError:
                    isdict = False
                if isdict:
                    for species in timedSpecies.keys():
                        pipeline.selectedDataSource.setMapping('p_%s' % species,
                                                               '(t>= %d)*(t<%d)' % timedSpecies[species])
                else:
                    for entry in timedSpecies:
                        pipeline.selectedDataSource.setMapping('p_%s' % entry['name'],
                                                               '(t>= %d)*(t<%d)' % (entry['t_start'],entry['t_end']))
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

    def OnMeasureTau2(self, event):

        pipeline = self.pipeline
        #idCol = 'objectID'
        idCol = self.fitSettings.IDcolumn
        measure = fitDarkTimes.measureObjectsByID2(pipeline, idname=idCol, settings=self.fitSettings)
        self.qpMeasurements[pipeline.selectedDataSourceKey] = {'Everything' : measure}
        
        newcols =  fitDarkTimes.retrieveMeasuresForIDs2(measure,pipeline[idCol])       

        qp_measured = unique_name('qp_measured',pipeline.dataSources.keys())
        recipe = pipeline.recipe
        mapp = tablefilters.Mapping(recipe,inputName=pipeline.selectedDataSourceKey,
                                    outputName=qp_measured)
        recipe.add_module(mapp)
        recipe.execute()

        qpm = recipe.namespace[qp_measured]

        for sourceCol in newcols.keys():
            qpm.addColumn("qp_%s" % sourceCol,newcols[sourceCol])

    def OnMeasureTau(self, event):

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

    def OnChannelMeasureTau(self, event, chan=None, mapMeasuresToEvents=True):

        if chan is None:
            chan = selectWithDialog(self.pipeline.colourFilter.getColourChans(), message='select channel')
        if chan is None:
            return

        pipeline = self.pipeline
        dispColor = self.pipeline.colourFilter.currentColour

        # make sure we get our IDs when everything is present -
        # otherwise colourfiltering can create an issue with some IDs not present
        # in some channels
        pipeline.colourFilter.setColour('Everything')
        ids = np.unique(pipeline['objectID'].astype('int'))
        idsnz = ids[ids > 0]

        pipeline.colourFilter.setColour(chan)

        tausrc = self.tausrc
        self.qpMeasurements[chan] = fitDarkTimes.measureObjectsByID(pipeline, set(idsnz), sigDefocused=165.0)
        qmc = self.qpMeasurements[chan]
        fitDarkTimes.makeRatio(qmc, 'qIndex', 100.0, qmc[tausrc['tau']])
        fitDarkTimes.makeRatio(qmc, 'NTauDiff', np.abs(qmc['tau1']-qmc['tau2']), qmc['tau2'])
        
        if mapMeasuresToEvents:
            # switch back to all channels
            pipeline.colourFilter.setColour('Everything')
            colmapNames = {
                tausrc['tau']     : 'taudark_%s' % chan,
                tausrc['tauerr']  : 'taudark_error_%s' % chan,
                tausrc['chisq']   : 'taudark_chisq2_%s' % chan,
                'NDarktimes'      : 'NDarktimes_%s' % chan,
                'qIndex'          : 'qIndex_%s' % chan,
                'NEvents'         : 'NEvents_%s' % chan,
                'NDefocusedFrac'  : 'NDefocusedFrac_%s' % chan,
                'NTauDiff'        : 'NTauDiff_%s' % chan
            }

            newcols =  fitDarkTimes.retrieveMeasuresForIDs(self.qpMeasurements[chan],pipeline['objectID'],
                                                           columns=colmapNames.keys())       
            for sourceCol in colmapNames.keys():
                pipeline.addColumn(colmapNames[sourceCol],newcols[sourceCol])
            
        # restore original display settings
        pipeline.colourFilter.setColour(dispColor)


    def OnMergeChannelMeasures(self, event):
        if len(self.pipeline.colourFilter.getColourChans()) > 0:
            channels = self.pipeline.colourFilter.getColourChans()
            if np.all([chan in self.qpMeasurements.keys() for chan in channels]):
               mergedMeas = fitDarkTimes.mergeChannelMeasurements(channels,[self.qpMeasurements[chan] for chan in channels])
               self.qpMeasurements['Everything'] = mergedMeas


    def OnChannelMeasurementRatios(self, event=None, channels=None, cals=None):
        from PYME.recipes.traits import HasTraits, Enum, Float

        if not len(self.pipeline.colourFilter.getColourChans()) > 0:
            return

        if channels is None:
            chans = self.pipeline.colourFilter.getColourChans()
            cChoice = myCChoice()
            cChoice.add_channels(chans)
            if cChoice.configure_traits(kind='modal'):
                channels = (cChoice.Channel1, cChoice.Channel2)
                cals = (cChoice.Channel1Calibration,cChoice.Channel2Calibration)
                
        # here if cancel from configure_traits
        if channels is None:
            return
        
        if ('Everything' in self.qpMeasurements):
            fitDarkTimes.mergedMeasurementsRatios(self.qpMeasurements['Everything'],
                                                  channels[0], channels[1], cals[0], cals[1])            

        def chanName(key,chan):
            return '%s_%s' % (key,chan)
        
        pipeline = self.pipeline
        dispColor = pipeline.colourFilter.currentColour
        pipeline.colourFilter.setColour('Everything')
        
        tausrc = self.tausrc
        colmapNames = {
            tausrc['tau']:       'taudark',
            'NDarktimes': 'NDarktimes',
            'qIndex':     'qIndex',
            'qIndexC':    'qIndexC',
            tausrc['tauerr']:    'taudark',
            tausrc['chisq']:    'taudark_chisq2',
            'qDensity':   'qDensity',
            'qDensityC':   'qDensityC'}

        for chan in channels:
            cmapChan = { chanName(key,chan) : chanName(value,chan) for key, value in colmapNames.items()}
            newcols =  fitDarkTimes.retrieveMeasuresForIDs(self.qpMeasurements['Everything'],pipeline['objectID'],
                                                               columns=cmapChan.keys())       
            for sourceCol in cmapChan.keys():
                if cmapChan[sourceCol] not in pipeline.keys():
                    pipeline.addColumn(cmapChan[sourceCol],newcols[sourceCol])

        for ratioName in ['qRatio','qRatioC']:
            ratiokey = '%s_%svs%s' % (ratioName,channels[0],channels[1])
            newcol = fitDarkTimes.retrieveMeasuresForIDs(self.qpMeasurements['Everything'],pipeline['objectID'],
                                                         columns=[ratiokey])
            pipeline.addColumn(ratiokey, newcol[ratiokey])

        meas = self.qpMeasurements['Everything']
        chisq_chan1 = meas[chanName(self.tausrc['chisq'],channels[0])]
        chisq_chan2 = meas[chanName(self.tausrc['chisq'],channels[0])]
        ratioChisq = 'qRatio_chisq2'
        fitDarkTimes.makeSum(meas,ratioChisq,chisq_chan1,chisq_chan2)
        newcol = fitDarkTimes.retrieveMeasuresForIDs(self.qpMeasurements['Everything'],pipeline['objectID'],
                                                     columns=[ratioChisq])
        pipeline.addColumn(ratioChisq, newcol[ratioChisq])

        # restore original display settings
        pipeline.colourFilter.setColour(dispColor)

    def OnChannelCalibrate(self, event, channels=None):
        # currently just an interface test function
        
        chans = ['ryrtest','jphtest']
        if channels is None:
            cChoice = myCChoice()
            cChoice.add_channels(chans)
            if cChoice.configure_traits(view='one_view', kind='modal'):
                channels = (cChoice.RatioChannel1, cChoice.RatioChannel2)
                cals = (cChoice.Channel1Calibration,cChoice.Channel2Calibration)
 

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


    def OnLabelLookupByID(self, event):
        from PYME.DSView import dsviewer, ViewIm3D
        import PYME.IO.image as im
        
        if ('Everything' in self.qpMeasurements):
            meas = self.qpMeasurements['Everything']
            selection = selectWithDialog(dsviewer.openViewers.keys())
            if selection is None:
                return
            keyChoice = KeyChoice()
            keyChoice.add_keys(sorted(meas.keys()))
            if not keyChoice.configure_traits(kind='modal'):
                return

            labelimg = dsviewer.openViewers[selection].image
            labels = labelimg.data[:,:,:].squeeze()
            measures = meas.lookupByID(labels,keyChoice.Key)
            newimg = im.ImageStack(measures, titleStub = 'Measure %s' % (keyChoice.Key))
            newimg.mdh.copyEntriesFrom(labelimg.mdh)
            newimg.mdh['Parent'] = labelimg.filename
            newimg.mdh['Processing.qpMeasure'] = keyChoice.Key

            ViewIm3D(newimg, mode='visGUI', title='Measure %s' % (keyChoice.Key),
                     glCanvas=self.visFr.glCanvas, parent=self.visFr)

            
    def OnSelectImgAndProcessMulticol(self, event):
        from PYME.DSView import dsviewer

        selection = selectWithDialog(dsviewer.openViewers.keys())
        if selection is None:
            return

        img = dsviewer.openViewers[selection].image
        self.OnTimedSpeciesFromImage(None,img=img)

        # get channel ratio info
        if len(self.pipeline.colourFilter.getColourChans()) > 0:
            chans = self.pipeline.colourFilter.getColourChans()
            cChoice = myCChoice()
            cChoice.add_channels(chans)
            if cChoice.configure_traits(kind='modal'):
                channels = (cChoice.RatioChannel1, cChoice.RatioChannel2)
                cals = (cChoice.Channel1Calibration,cChoice.Channel2Calibration)
            else:
                channels = None
                cals = None

        prog = wx.ProgressDialog("Process Multicolour Data", "Setting Drift...", 100,
                                 style=wx.PD_ELAPSED_TIME | wx.PD_AUTO_HIDE)
        prog.Update(15)
        self.OnSetDriftPars(None,img=img)

        prog.Update(40,"Setting IDs...")        
        self.OnGetIDsfromImage(None,img=img)

        prog.Update(50,"Calculating qIndices...")        
        for chan in self.pipeline.colourFilter.getColourChans():
            self.OnChannelMeasureTau(None,chan=chan)

        prog.Update(85,"Calculating Areas and Ratios...")        
        # areas last so that qpMeasures get the area info 
        self.OnAreaFromLabels(None,img=img)
        self.OnMergeChannelMeasures(None)

        prog.Update(95,"Nearly Done")
        if channels is not None:
            self.OnChannelMeasurementRatios(None,channels=channels,cals=cals)
        prog.Update(100)
        prog.Destroy()

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
        from PYMEcs.IO.darkTimes import saveDarkTimesCSV
        import os.path
        import yaml
        
        if len(self.qpMeasurements.keys()) > 0:
            filename = wx.FileSelector('Save all measurements (select basename)...',
                                       wildcard="CSV files (*.csv)|*.csv|Excell files (*.xlsx)|*.xlsx|HDF5 files (*.hdf)|*.hdf", 
                                       flags = wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT)
                                   
            if not filename == '':
                base, ext = os.path.splitext(filename)
                for ds in self.qpMeasurements.keys():
                    for chan in self.qpMeasurements[ds]:
                        qpm = self.qpMeasurements[ds][chan]
                        runRecipe.saveOutput(qpm['measures'], base + '_' + ds + '_' + chan + ext)
                        saveDarkTimesCSV(qpm['measures']['objectID'],qpm['darkTimes'],
                                         base + '_' + ds + '_' + chan + '_' + 'dts.csv')
                        saveDarkTimesCSV(qpm['measures']['objectID'],qpm['times'],
                                         base + '_' + ds + '_' + chan + '_' + 'times.csv')
                        with open(base + '_' + ds + '_' + chan + '_' + 'state.yaml', 'w') as stream:
                            yaml.dump(qpm['state'], stream)

    def OnLoadMeasurements(self,event):
        import PYMEcs.IO.tabular as tb
        import os.path
        
        filename = wx.FileSelector('Load measurements (select basename)...',
                                   wildcard="CSV files (*.csv)|*.csv", 
                                   flags = wx.FD_OPEN | wx.FD_FILE_MUST_EXIST)
                                   
        if not filename == '':
            self.qpMeasurements['Everything'] = tb.tabularFromCSV(filename)


    def OnFitSettings(self,event):
        self.fitSettings.configure_traits(kind='modal')

    def OnDarkT(self,event):
        try:
            from StringIO import StringIO ## for Python 2
        except ImportError:
            from io import StringIO ## for Python 3

        
        visFr = self.visFr
        pipeline = visFr.pipeline
        mdh = pipeline.mdh

        NTMIN = 5
        maxPts = 1e5

        ## STEP 1: Extract times of relevant events, preprocess as needed and get darktimes from that (dtg)

        tc, dtg, nevents, nevents_corrected, usingClumpIndex, usingTminTmax = \
                fitDarkTimes.extractEventTimes(pipeline,
                                               useTminTmax = (self.fitSettings.coalescedProcessing == 'useTminTmaxIfAvailable'),
                                               return_modes=True)

        ## STEP 2: from the extracted times get the actual darktimes
        
        # determine darktime from gaps and reject zeros (no real gaps) 
        #dts = tc[1:]-tc[0:-1]-1
        #dtg = dts[dts>0]
        nts = dtg.size

        ## STEP 3: fit the dark time distribution
        
        if nts > NTMIN:
            # now make a cumulative histogram from these
            cumux,cumuy = fitDarkTimes.cumuhist(dtg)
            binctrs,hc,binctrsg,hcg = fitDarkTimes.cumuhistBinnedLog(dtg,dlog=0.05,return_good=True)
            binctrsCoarse,hcCoarse,binctrsgCoarse,hcgCoarse, hCoarse = fitDarkTimes.cumuhistBinnedLog(dtg,dlog=0.2,return_hist=True,return_good=True)
            
            # fit theoretical distributions
            from scipy.optimize import curve_fit

            try:
                tau1, tau1err, chisqr, tau1est = fitDarkTimes.fitTaus(binctrs,hc,return_tau1est=True)
                logger.info("fitTaus: tau1 %f, tau1err %f, chisq %.1g, tau1est %f" % (tau1, tau1err, chisqr, tau1est))
            except RuntimeError:
                logger.warn("fitting tau value from histogram failed")
                tau1ok = False
            else:
                tau1ok = True
                
            try:
                tau1c, tau1errc, chisqrc, tau1estc = fitDarkTimes.fitTaus(cumux,cumuy,return_tau1est=True)
                logger.info("fitTaus: tau1 %f, tau1err %f, chisq %.1g, tau1est %f" % (tau1c, tau1errc, chisqrc, tau1estc))
            except RuntimeError:
                logger.warn("fitting tau value from cumulative dist failed")
                tau1cok = False
            else:
                tau1cok = True
                
            try:
                tau1m, tau2m, atau1m, tau1errm, tau2errm, atau1errm, chisqrm = fitDarkTimes.fitTaus(binctrs,hc,
                                                                                                    fitTau2=not self.fitSettings.Tau2Constant,
                                                                                                    tau2const=self.fitSettings.Tau2FixedValue)
                logger.info("fitTaus: tau1 %f, tau2 %f, atau1 %f, tau1err %f, tau2err %f, atau1err %f, chisq %.1g" %
                            (tau1m, tau2m, atau1m, tau1errm, tau2errm, atau1errm, chisqrm))
            except RuntimeError:
                logger.warn("fitting two tau values failed")
                tau1mok = False
            else:
                tau1mok = True
            
            ## STEP 4: plotting and reporting

            import matplotlib.pyplot as plt

            def safepltshow():
                try:
                    plt.show()
                except RuntimeError:
                    pass
                
            # make a string with info
            outstr = StringIO()
            
            if usingClumpIndex:
                if usingTminTmax:
                    print("events: %d (raw clumpSize: %d), dark times: %d (using clumpIndices + Tmin & Tmax)"
                          % (nevents_corrected,nevents,nts),file=outstr)
                else:
                    print("events: %d (from clumpSize), dark times: %d (using clumpIndices)"
                          % (nevents,nts),file=outstr)
            else:    
                print("events: %d, dark times: %d" % (nevents,nts),file=outstr)

            print("darktime: %.1f+-%d (%.1f+-%d) frames - chisqr %.2g (%.2g)"
                  % (tau1c,tau1errc,
                     tau1,tau1err,
                     chisqrc,chisqr),file=outstr)
            if tau1mok:
                print("darktime: tau1 %.1f+-%d, tau2 %.1f+-%d"
                      % (tau1m, tau1errm, tau2m, tau2errm),file=outstr)
            if tau1mok:
                print("amplitudes: (a %.2f+-%.2f, b %.2f) - chisqr %.2g"
                      % (atau1m,atau1errm, 1-atau1m,chisqrm),file=outstr) 
            print("darktime: starting estimates: %.1f (%.1f)" % (tau1estc,tau1est),file=outstr)
            print("qunits: %.2f (%.2f)" % (100.0/tau1c, 100.0/tau1),file=outstr)

            labelstr = outstr.getvalue()

            # plot data and fitted curves
            plt.figure()
            plt.subplot(211)
            plt.plot(cumux,cumuy,'o')
            if tau1cok: plt.plot(cumux,fitDarkTimes.cumuexpfit(cumux,tau1c))
            
            plt.plot(binctrs,hc,'o')
            if tau1ok: plt.plot(binctrs,fitDarkTimes.cumuexpfit(binctrs,tau1))

            if tau1mok: plt.plot(binctrs,fitDarkTimes.cumumultiexpfit(binctrs,tau1m,tau2m,atau1m),'--')

            plt.ylim(-0.2,1.2)
            plt.subplot(212)
            plt.semilogx(cumux,cumuy,'o')
            if tau1cok: plt.semilogx(cumux,fitDarkTimes.cumuexpfit(cumux,tau1c))
            plt.semilogx(binctrs,hc,'o')
            if tau1ok: plt.semilogx(binctrs,fitDarkTimes.cumuexpfit(binctrs,tau1))

            if tau1mok: plt.semilogx(binctrs,fitDarkTimes.cumumultiexpfit(binctrs,tau1m,tau2m,atau1m),'--')
            
            plt.ylim(-0.2,1.2)
            safepltshow()
            
            plt.annotate(labelstr, xy=(0.5, 0.1), xycoords='axes fraction',
                         fontsize=10)

            plt.figure()
            plt.subplot(211)
            if True:
                plt.semilogx(cumux,cumuy,'o',color='darkblue')
                if tau1ok: plt.semilogx(binctrs,fitDarkTimes.cumuexpfit(binctrs,tau1))
                if tau1mok: plt.semilogx(binctrs,fitDarkTimes.cumumultiexpfit(binctrs,tau1m,tau2m,atau1m),'--')
                plt.ylim(-0.2,1.2)
                plt.annotate(labelstr, xy=(0.5, 0.1), xycoords='axes fraction',
                             fontsize=10)
            else:
                plt.semilogx(binctrsCoarse, hCoarse, 'o')

            plt.subplot(212)
            plt.semilogx(binctrs,hc,'go')
            plt.semilogx(binctrsg,hcg,'bo')
            if tau1ok: plt.semilogx(binctrs,fitDarkTimes.cumuexpfit(binctrs,tau1))
            if tau1mok: plt.semilogx(binctrs,fitDarkTimes.cumumultiexpfit(binctrs,tau1m,tau2m,atau1m),'--',color='orange')
            
            plt.ylim(-0.2,1.2)
            safepltshow()

        else:
            Warn(None, 'not enough data points, only found %d dark times (need at least  %d)' % (nts,NTMIN), 'Error')


def Plug(visFr):
    """Plugs this module into the gui"""
    visFr.qPAINTCalc = QPCalc(visFr)

