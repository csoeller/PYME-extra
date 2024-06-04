from PYME.warnings import warn
from PYMEcs.Analysis.NPC import estimate_nlabeled, npclabel_fit, plotcdf_npc3d
from PYME.recipes import tablefilters
import wx
from traits.api import HasTraits, Str, Int, CStr, List, Enum, Float, Bool
import numpy as np
import matplotlib.pyplot as plt
from PYMEcs.misc.utils import unique_name
import logging

logger = logging.getLogger(__name__)

def selectWithDialog(choices, message='select image from list', caption='Selection'):
    dlg = wx.SingleChoiceDialog(None, message, caption, list(choices), wx.CHOICEDLG_STYLE)
    if dlg.ShowModal() == wx.ID_OK:
        item = dlg.GetStringSelection()
    else:
        item = None
    dlg.Destroy()
    return item

def setNPCsfromImg(pipeline,img):
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

    ind = (pixX < img.data_xyztc.shape[0])*(pixY < img.data_xyztc.shape[1])*(pixX >= 0)*(pixY >= 0)
        
    ids = np.zeros_like(pixX)
    #assume there is only one channel
    ids[ind] = img.data_xyztc[:,:,0,0,0].squeeze()[pixX[ind], pixY[ind]].astype('i')
        
    numPerObject, b = np.histogram(ids, np.arange(ids.max() + 1.5) + .5)
        
    withIDs.addColumn('objectID', ids)
    withIDs.addColumn('NEvents', numPerObject[ids-1])
        
    recipe.add_module(tablefilters.FilterTable(recipe,inputName=with_ids, outputName=valid_ids,
                                                   filters={'objectID' : [.5, 1e5]}))
    recipe.execute()
        
    pipeline.selectDataSource(valid_ids)


class NPCsettings(HasTraits):
    SegmentThreshold_2D = Int(10)
    SegmentThreshold_3D = Int(1)
    SecondPass = Bool(False)
    FitMode = Enum(['abs','square'])

class NPCcalc():
    def __init__(self, visFr):
        self.visFr = visFr

        visFr.AddMenuItem('Experimental>NPCs', "Select NPCs by Mask", self.OnSelectNPCsByMask)
        visFr.AddMenuItem('Experimental>NPCs', "Analyse single 2D NPC\tCtrl+N", self.OnAnalyseSingleNPC)
        visFr.AddMenuItem('Experimental>NPCs', "Analyse 2D NPCs by ID", self.OnAnalyseNPCsByID)
        visFr.AddMenuItem('Experimental>NPCs', "Show 2D NPC labeling Statistics", self.OnNPCstats)
        visFr.AddMenuItem('Experimental>NPCs', "Select by mask, analyse and show stats (2D)", self.OnNPCcombinedFuncs)
        visFr.AddMenuItem('Experimental>NPCs', "Analyse 3D NPCs by ID", self.OnAnalyse3DNPCsByID)
        visFr.AddMenuItem('Experimental>NPCs', "Save 3D NPC Measurements",self.OnNPC3DSaveMeasurements)
        visFr.AddMenuItem('Experimental>NPCs', "Load and display saved 3D NPC Measurements",self.OnNPC3DLoadMeasurements)
        visFr.AddMenuItem('Experimental>NPCs', 'NPC Analysis settings', self.OnNPCsettings)

        self.NPCsettings = NPCsettings()


    def OnNPCsettings(self, event=None):
        if self.NPCsettings.configure_traits(kind='modal'):
            pass

    def OnAnalyseSingleNPC(self, event):
        pipeline = self.visFr.pipeline
        Nevents = pipeline['x'].size
        if Nevents > 500:
            warn('More than 500 events in pipeline - check if really a single NPC')
            return

        xExtent = pipeline['x'].max() - pipeline['x'].min()
        yExtent = pipeline['y'].max() - pipeline['y'].min()

        maxextent_nm = 300
        if (xExtent > maxextent_nm) or (yExtent > maxextent_nm):
            warn('x or y bounding box > %d nm (%d,%d) - check if single NPC' % (maxextent_nm,xExtent,yExtent))
            return

        estimate_nlabeled(pipeline['x'],pipeline['y'],nthresh=self.NPCsettings.SegmentThreshold_2D,
                          do_plot=True,secondpass=self.NPCsettings.SecondPass,fitmode=self.NPCsettings.FitMode)
        

    def OnAnalyseNPCsByID(self, event=None):
        from PYMEcs.recipes.localisations import NPCAnalysisByID
        pipeline = self.visFr.pipeline
        recipe = pipeline.recipe
        with_npclinfo = unique_name('with_npcinfo',pipeline.dataSources.keys())
        # for the NPCAnalysisByID module use the current NPCsettings
        npcanalysis = NPCAnalysisByID(inputName=pipeline.selectedDataSourceKey,outputName=with_npclinfo,
                                      SegmentThreshold=self.NPCsettings.SegmentThreshold_2D,
                                      SecondPass=self.NPCsettings.SecondPass,FitMode=self.NPCsettings.FitMode)
        if not npcanalysis.configure_traits(kind='modal'):
            return

        recipe.add_module(npcanalysis)
        recipe.execute()
        
        pipeline.selectDataSource(with_npclinfo)

    def OnAnalyse3DNPCsByID(self, event=None):
        from PYMEcs.Analysis.NPC import NPC3DSet
        pipeline = self.visFr.pipeline

        if 'npcs' in dir(pipeline):
            npcs = pipeline.npcs
            do_plot = False
        else:
            npcs = NPC3DSet()
            do_plot = True
            for oid in np.unique(pipeline['objectID']):
                npcs.addNPCfromPipeline(pipeline,oid)

        # for example use of ProgressDialog see also
        # https://github.com/Metallicow/wxPython-Sample-Apps-and-Demos/blob/master/101_Common_Dialogs/ProgressDialog/ProgressDialog_extended.py
        progress = wx.ProgressDialog("NPC analysis in progress", "please wait", maximum=len(npcs.npcs),
                                     parent=self.visFr,
                                     style=wx.PD_SMOOTH
                                     | wx.PD_AUTO_HIDE
                                     | wx.PD_CAN_ABORT
                                     | wx.PD_ESTIMATED_TIME
                                     | wx.PD_REMAINING_TIME)
        if do_plot:
            fig, axes=plt.subplots(2,3)
        cancelled = False
        npcs.measurements = []
        for i,npc in enumerate(npcs.npcs):
            if not npc.fitted:
                npc.fitbymll(npcs.llm,plot=True,printpars=False,axes=axes)
            nt,nb = npc.nlabeled(nthresh=self.NPCsettings.SegmentThreshold_3D,dr=20.0)
            npcs.measurements.append([nt,nb])
            (keepGoing, skip) = progress.Update(i+1)
            if not keepGoing:
                logger.info('OnAnalyse3DNPCsByID: progress cancelled, aborting NPC analysis')
                cancelled = True
                progress.Destroy()
                wx.Yield()
                # Cancelled by user.
                break
            wx.Yield()

        if not cancelled:
            pipeline.npcs = npcs # overwriting with the same object should be fine if pipeline.npcs already existed
            npcs.plot_labeleff()

    def OnNPC3DSaveMeasurements(self, event=None):
        pipeline = self.visFr.pipeline
        if 'npcs' not in dir(pipeline) or 'measurements' not in dir(pipeline.npcs):
            warn('no valid NPC measurements found, therefore cannot save...')
            return
        
        fdialog = wx.FileDialog(self.visFr, 'Save NPC measurements as ...',
                                wildcard='CSV (*.csv)|*.csv',
                                style=wx.FD_SAVE)
        if fdialog.ShowModal() != wx.ID_OK:
            return

        fpath = fdialog.GetPath()
        meas = np.array(pipeline.npcs.measurements, dtype='i')
        import pandas as pd
        df = pd.DataFrame({'Ntop_NPC3D': meas[:, 0], 'Nbot_NPC3D': meas[:, 1]})

        from pathlib import Path
        with open(fpath, 'w') as f:
            f.write('# threshold %d, source data file %s\n' %
                    (self.NPCsettings.SegmentThreshold_3D,Path(pipeline.filename).name))

        df.to_csv(fpath,index=False, mode='a')


    def OnNPC3DLoadMeasurements(self, event=None):
        fdialog = wx.FileDialog(self.visFr, 'Load NPC measurements from ...',
                                wildcard='CSV (*.csv)|*.csv',
                                style=wx.FD_OPEN)
        if fdialog.ShowModal() != wx.ID_OK:
            return
        import pandas as pd
        meas = pd.read_csv(fdialog.GetPath(),comment='#')
        # here we need some plotting code
        nlab = meas['Ntop_NPC3D'] + meas['Nbot_NPC3D']
        plt.figure()
        plotcdf_npc3d(nlab)
        
    def OnSelectNPCsByMask(self,event=None):
        from PYME.DSView import dsviewer

        pipeline = self.visFr.pipeline
            
        selection = selectWithDialog(dsviewer.openViewers.keys())
        if selection is not None:
            img = dsviewer.openViewers[selection].image

        if img is None:
            return

        setNPCsfromImg(pipeline,img)

    def OnNPCstats(self,event=None):
        pipeline = self.visFr.pipeline
        for prop in ['NPCnlabeled','objectID']:
            if prop not in pipeline.keys():
                warn("property '%s' not found in datasource" % prop)
                return
  
        uids,indices,idcounts = np.unique(pipeline['objectID'],return_index=True,return_counts=True)
        uidnz = uids > 0

        nlabel_uq = pipeline['NPCnlabeled'][indices[uidnz]]

        ds_mdh = pipeline.selectedDataSource.mdh
        plt.figure()
        nlab, bions, patches = plt.hist(nlabel_uq,bins=-0.5+np.arange(10,dtype='int'))
        plabel, nlabels, perr = npclabel_fit(nlab)
        plt.plot(np.arange(9),nlabels,'-o')
        plt.xlabel('Number of segments labeled')
        plt.ylabel('Frequency')
        title = 'NPC labeling statistics, pl_fit = %.2f+-%.3f' % (plabel,perr)
        if 'NPCAnalysis.EventThreshold' in ds_mdh:
            title += "\nEvent threshold = %d, mode = %s" % (int(ds_mdh['NPCAnalysis.EventThreshold']),
                                                              ds_mdh['NPCAnalysis.RotationAlgorithm'])
        plt.title(title)
        plt.tight_layout()

    def OnNPCcombinedFuncs(self,event):
        self.OnSelectNPCsByMask()
        self.OnAnalyseNPCsByID()
        self.OnNPCstats()
        
def Plug(visFr):
    """Plugs this module into the gui"""
    NPCcalc(visFr)
