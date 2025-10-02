import logging
import os

import matplotlib.pyplot as plt
import numpy as np
import wx
from PYME.recipes import tablefilters
from traits.api import Bool, Enum, Float, HasTraits, Int

from PYMEcs.Analysis.NPC import estimate_nlabeled, npclabel_fit, plotcdf_npc3d
from PYMEcs.IO.NPC import findNPCset
from PYMEcs.misc.utils import unique_name
from PYMEcs.pyme_warnings import warn

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
    SegmentThreshold_2D = Int(10,label='Threshold localisation count (2D)',
                              desc="a minimum number of localisations in a segment in 2D fitting to be counted as labeled; "+
                              "NOTE: definition changed from original code where 11 localisations where required with a threshold of 10 etc!")
    SegmentThreshold_3D = Int(1,label='Threshold localisation count (3D)',
                              desc="a minimum number of localisations in a segment in 3D fitting to be counted as labeled")
    SecondPass_2D = Bool(False,label='Second pass for NPC fitting (2D)',
                         desc="a second pass for 2D fitting should be run; we have experimented with a second pass rotation "+
                         "estimate and fitting hoping to improve on the first estimate; still experimental")
    StartHeight_3D = Float(50.0,label='Starting ring spacing for 3D fitting',
                           desc="starting ring spacing value for the 3D fit; note only considered when doing the initial full fit; "+
                           "not considered when re-evaluating existing fit")
    StartDiam_3D = Float(107.0,label='Starting ring diameter for 3D fitting',
                           desc="starting ring diameter value for the 3D fit; note only considered when doing the initial full fit; "+
                           "not considered when re-evaluating existing fit")
    TemplateMode_3D = Enum(['standard','detailed','twostage'],desc="standard or detailed NPC template")
    TemplateSigma_3D = Float(7.0,label='Sigma value for 3D template smoothing',desc="smoothing SD (sigma) to smooth the 3D fitting template, default is 7 nm")
    FitMode = Enum(['abs','square'],label='Fit mode for NPC rotation',
                   desc="fit mode for NPC rotation; in 2D and 3D estimation of the NPC lateral rotation a simple algorithm is used to find the start of the 'pizza pieces'; "+
                   "this mode refers to use of absolute or squared differences in the calculation; default should be ok")
    RotationLocked_3D = Bool(True,label='NPC rotation estimate locked (3D)',
                             desc="when estimating the NPC rotation (pizza slice boundaries), the top and bottom rings in 3D should be locked, "+
                             "i.e. have the same rotation from the underlying structure")
    RadiusUncertainty_3D = Float(20.0,label="Radius scatter in nm",
                                 desc="a radius scatter that determines how much localisations can deviate from the mean ring radius "+
                                 "and still be accepted as part of the NPC; allows for distortions of NPCs and localisation errors")
    # the two next things seem to set the same thing, so unify
    Zclip_3D = Float(75.0,label='Z-clip value from center of NPC',
                     desc='the used zrange from the (estimated) center of the NPC, from (-zclip..+zclip) in 3D fitting')
    OffsetMode_3D = Enum(['median','mean'],label='Method to estimate NPC center',
                         desc="Method to estimate the likely 3D center of an NPC; median seems more robust against outliners (say in z)")
    SkipEmptyTopOrBottom_3D = Bool(False,desc="if to skip NPCs with empty top or bottom ring")
    KnownNumber_3D = Int(-1, desc="if a known number of NPCs are present but some have 0 events use this to set; only considered if value >0")
    NPCRotationAngle = Enum(['positive','negative','zero'],desc="way to treat rotation for NPC gallery")
    NPCGalleryArrangement = Enum(['SingleAverageSBS','TopOverBottom','TopBesideBottom','SingleAverage'],desc="how to arrange 3D NPC parts in NPC gallery; SBS = SideBySide top and bottom")
    NoRotationForSingleNPCFit = Bool(False,desc="if rotation is disabled for single NPC fit")
    IncludeSegmentLinesWithGallery = Bool(True,desc="if 3D segement lines are generated when NPC 3D gallery is created")


# this should be a backwards compatible way to access the main filename associated with the pipeline/datasource
def pipeline_filename(pipeline):
    try:
        filename = pipeline.filename
    except (KeyError,AttributeError):
        # latest versions with session saving associated filenames with the actual datasource read in
        filename = pipeline.dataSources['FitResults'].filename
    return filename
        
class NPCcalc():
    def __init__(self, visFr):
        self.visFr = visFr

        visFr.AddMenuItem('Experimental>NPC2D', "Select NPCs by Mask", self.OnSelectNPCsByMask)
        visFr.AddMenuItem('Experimental>NPC3D', "Select NPCs by Mask", self.OnSelectNPCsByMask)
        visFr.AddMenuItem('Experimental>NPC2D', "Analyse single 2D NPC\tCtrl+N", self.OnAnalyseSingleNPC)
        visFr.AddMenuItem('Experimental>NPC2D', "Analyse 2D NPCs by ID", self.OnAnalyseNPCsByID)
        visFr.AddMenuItem('Experimental>NPC2D', "Show 2D NPC labeling Statistics", self.OnNPCstats)
        visFr.AddMenuItem('Experimental>NPC2D', "Select by mask, analyse and show stats", self.OnNPCcombinedFuncs)
        visFr.AddMenuItem('Experimental>NPC3D', "Analyse 3D NPCs by ID", self.OnAnalyse3DNPCsByID)
        visFr.AddMenuItem('Experimental>NPC3D', "Add 3D NPC templates", self.On3DNPCaddTemplates)
        visFr.AddMenuItem('Experimental>NPC3D', "Save NPC Set with full fit analysis", self.OnNPC3DSaveNPCSet)
        visFr.AddMenuItem('Experimental>NPC3D', "Load stored NPC Set with full fit analysis", self.OnNPC3DLoadNPCSet)
        visFr.AddMenuItem('Experimental>NPC3D', "Save Measurements Only (csv, no fit info saved)",self.OnNPC3DSaveMeasurements)
        visFr.AddMenuItem('Experimental>NPC3D', "Load and display saved Measurements (from csv)",self.OnNPC3DLoadMeasurements)
        visFr.AddMenuItem('Experimental>NPC3D', "Show NPC geometry statistics",self.OnNPC3DGeometryStats)
        visFr.AddMenuItem('Experimental>NPC3D', "Save NPC geometry statistics as CSV",self.OnNPC3DSaveGeometryStats)
        visFr.AddMenuItem('Experimental>NPC3D', "Show NPC template fit statistics",self.OnNPC3DTemplateFitStats)
        visFr.AddMenuItem('Experimental>NPC2D', 'NPC Analysis settings', self.OnNPCsettings)
        visFr.AddMenuItem('Experimental>NPC3D', 'NPC Analysis settings', self.OnNPCsettings)
        visFr.AddMenuItem('Experimental>NPC3D', 'Add NPC Gallery', self.On3DNPCaddGallery)
        visFr.AddMenuItem('Experimental>NPC3D', 'Plot NPC by-segment data', self.OnNPC3DPlotBySegments)
        visFr.AddMenuItem('Experimental>NPC3D', 'Save NPC by-segment data', self.OnNPC3DSaveBySegments)
        visFr.AddMenuItem('Experimental>NPC3D', 'Analysis 3D NPCs by ID and save all outputs', self.OnNPC3DRunAllActions) # Add combined MenuItem for all requested actions


        self._npcsettings = None
        self.gallery_layer = None
        self.segment_layer = None

        # --- Alex B addition for running several NPC calc action at once ---
    def OnNPC3DRunAllActions(self, event=None):
        """Performs all key NPC 3D analysis actions in sequence, prompting user for output folder."""
        # Prompt user for save directory
        with wx.DirDialog(self.visFr, "Select folder to save all NPC outputs", style=wx.DD_DEFAULT_STYLE | wx.DD_NEW_DIR_BUTTON) as dirdialog:
            if dirdialog.ShowModal() != wx.ID_OK:
                return  # User cancelled
            else:
                save_dir = dirdialog.GetPath()

        # Pass save_dir to all auto-save functions
        self.OnAnalyse3DNPCsByID_auto_save(save_dir=save_dir)
        self.OnNPC3DSaveNPCSet_auto_save(save_dir=save_dir)
        self.OnNPC3DSaveMeasurements_auto_save(save_dir=save_dir)
        self.OnNPC3DGeometryStats_auto_save(save_dir=save_dir)
        self.OnNPC3DTemplateFitStats_auto_save(save_dir=save_dir)
        self.OnNPC3DPlotBySegments_auto_save(save_dir=save_dir)
        self.OnNPC3DSaveBySegments_auto_save(save_dir=save_dir)
        self.OnNPC3DSaveGeometryStats_auto_save(save_dir=save_dir)
        # --- End of Alex B addition test for running several action at once ---

    @property
    def NPCsettings(self):
        if self._npcsettings is None:
            foreshortening=self.visFr.pipeline.mdh.get('MINFLUX.Foreshortening',1.0)
            if foreshortening < 1.0:
                self._npcsettings = NPCsettings(StartHeight_3D=70.0*foreshortening,Zclip_3D=75.0*foreshortening)
            else:
                self._npcsettings = NPCsettings(StartHeight_3D=50.0,Zclip_3D=55.0)
        return self._npcsettings
    
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

        if self.NPCsettings.NoRotationForSingleNPCFit:
            rot = 0
        else:
            rot = None
        estimate_nlabeled(pipeline['x'],pipeline['y'],nthresh=self.NPCsettings.SegmentThreshold_2D,rotation=rot,
                          do_plot=True,secondpass=self.NPCsettings.SecondPass_2D,fitmode=self.NPCsettings.FitMode)
        

    def OnAnalyseNPCsByID(self, event=None):
        from PYMEcs.recipes.localisations import NPCAnalysisByID
        pipeline = self.visFr.pipeline
        recipe = pipeline.recipe
        with_npclinfo = unique_name('with_npcinfo',pipeline.dataSources.keys())
        # for the NPCAnalysisByID module use the current NPCsettings
        npcanalysis = NPCAnalysisByID(inputName=pipeline.selectedDataSourceKey,outputName=with_npclinfo,
                                      SegmentThreshold=self.NPCsettings.SegmentThreshold_2D,
                                      SecondPass=self.NPCsettings.SecondPass_2D,FitMode=self.NPCsettings.FitMode)
        if not npcanalysis.configure_traits(kind='modal'):
            return

        recipe.add_module(npcanalysis)
        recipe.execute()
        
        pipeline.selectDataSource(with_npclinfo)

    def OnAnalyse3DNPCsByID(self, event=None):
        from PYMEcs.Analysis.NPC import NPC3DSet
        pipeline = self.visFr.pipeline

        if findNPCset(pipeline,warnings=False) is not None:
            npcs = findNPCset(pipeline)
            do_plot = False
        else:
            npcs = NPC3DSet(filename=pipeline_filename(pipeline),
                            zclip=self.NPCsettings.Zclip_3D,
                            offset_mode=self.NPCsettings.OffsetMode_3D,
                            NPCdiam=self.NPCsettings.StartDiam_3D,
                            NPCheight=self.NPCsettings.StartHeight_3D,
                            foreshortening=pipeline.mdh.get('MINFLUX.Foreshortening',1.0),
                            known_number=self.NPCsettings.KnownNumber_3D,
                            templatemode=self.NPCsettings.TemplateMode_3D,
                            sigma=self.NPCsettings.TemplateSigma_3D)
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
            if 'templatemode' in dir(npcs) and npcs.templatemode == 'twostage':
                figpre, axespre=plt.subplots(2,3,label='pre-llm')
        cancelled = False
        npcs.measurements = []
        if 'templatemode' in dir(npcs) and npcs.templatemode == 'detailed':
            rotation = 22.5 # this value may need adjustment
        else:
            rotation = None

        # keep track if any fits were performed
        anyfits = False
        for i,npc in enumerate(npcs.npcs):
            if not npc.fitted:
                if 'templatemode' in dir(npcs) and npcs.templatemode == 'twostage':
                    npc.fitbymll(npcs.llm,plot=True,printpars=False,axes=axes,preminimizer=npcs.llmpre,axespre=axespre)
                else:
                    npc.fitbymll(npcs.llm,plot=True,printpars=False,axes=axes)
                anyfits = True
            nt,nb = npc.nlabeled(nthresh=self.NPCsettings.SegmentThreshold_3D,
                                 dr=self.NPCsettings.RadiusUncertainty_3D,
                                 rotlocked=self.NPCsettings.RotationLocked_3D,
                                 zrange=self.NPCsettings.Zclip_3D,
                                 rotation=rotation)
            if self.NPCsettings.SkipEmptyTopOrBottom_3D and (nt == 0 or nb == 0):
                pass # we skip NPCs with empty rings in this case
            else:
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
        else:
            if anyfits:
                pipeline.npcs = npcs # we update the pipeline npcs attribute only if the for loop completed normally and we fitted

        if cancelled:
            return
        
        npcs.plot_labeleff(thresh=self.NPCsettings.SegmentThreshold_3D)

# --- Alex B addition --- 
# AIM: perform all actions 3D NPC actions and save output automatically

    def OnAnalyse3DNPCsByID_auto_save(self, event=None, save_dir=None):
        from PYMEcs.Analysis.NPC import NPC3DSet
        pipeline = self.visFr.pipeline

        # --- Alex B addition ---
        # We define a few variables used for automatic saving later
        if save_dir is None:
            save_dir = os.getcwd() # Default to current working directory
        MINFLUXts = pipeline.mdh.get('MINFLUX.TimeStamp') # Get the timestamp and use it for naming the file to save
        if MINFLUXts is not None:
            fitplot_filename = f"{MINFLUXts}-NPC_fit_plot.png"
            leplot_filename = f"{MINFLUXts}-LE_plot.png"
        else: # If no timestamp is found, use default filenames
            fitplot_filename = "NPC_fit_plot.png"
            leplot_filename = "LE_plot.png"
        fitplot_save_path = os.path.join(save_dir, fitplot_filename) # Save path for fit plot
        leplot_save_path = os.path.join(save_dir, leplot_filename) # Save path for LE plot
        # --- End Alex B addition ---

        if findNPCset(pipeline,warnings=False) is not None:
            npcs = findNPCset(pipeline)
            do_plot = False # If NPC set exists, do not re-analyze / re-plot
        else: # If NPC set does not exists, do analyze plot
            npcs = NPC3DSet(filename=pipeline_filename(pipeline),
                            zclip=self.NPCsettings.Zclip_3D,
                            offset_mode=self.NPCsettings.OffsetMode_3D,
                            NPCdiam=self.NPCsettings.StartDiam_3D,
                            NPCheight=self.NPCsettings.StartHeight_3D,
                            foreshortening=pipeline.mdh.get('MINFLUX.Foreshortening',1.0),
                            known_number=self.NPCsettings.KnownNumber_3D,
                            templatemode=self.NPCsettings.TemplateMode_3D,
                            sigma=self.NPCsettings.TemplateSigma_3D)
            do_plot = True # If NPC set does not exists, do analyze plot (lines 326)
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
        if do_plot: # IInitialize plots (only creates canvas and axes)
            fig, axes=plt.subplots(2,3)
            if 'templatemode' in dir(npcs) and npcs.templatemode == 'twostage':
                figpre, axespre=plt.subplots(2,3,label='pre-llm')
        cancelled = False
        npcs.measurements = []
        if 'templatemode' in dir(npcs) and npcs.templatemode == 'detailed':
            rotation = 22.5 # this value may need adjustment
        else:
            rotation = None


        # keep track if any fits were performed
        anyfits = False
        for i,npc in enumerate(npcs.npcs):
            if not npc.fitted:
                if 'templatemode' in dir(npcs) and npcs.templatemode == 'twostage':
                    npc.fitbymll(npcs.llm,plot=True,printpars=False,axes=axes,preminimizer=npcs.llmpre,axespre=axespre) # Plot generated here
                else:
                    npc.fitbymll(npcs.llm,plot=True,printpars=False,axes=axes) # Plot generated here
                anyfits = True
            nt,nb = npc.nlabeled(nthresh=self.NPCsettings.SegmentThreshold_3D,
                                 dr=self.NPCsettings.RadiusUncertainty_3D,
                                 rotlocked=self.NPCsettings.RotationLocked_3D,
                                 zrange=self.NPCsettings.Zclip_3D,
                                 rotation=rotation)
            if self.NPCsettings.SkipEmptyTopOrBottom_3D and (nt == 0 or nb == 0):
                pass # we skip NPCs with empty rings in this case
            else:
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
        else:
            if anyfits:
                pipeline.npcs = npcs # we update the pipeline npcs attribute only if the for loop completed normally and we fitted

        if cancelled:
            return

        #--- Alex B modif addition---
        
        # Save the main NPC fit plot 
        # Uses the variables (filename+path defined earlier)
            # Note: Save only the last fit plot 
        if do_plot and not cancelled:
            fig.savefig(fitplot_save_path)
            print(f"NPC fit plot automatically saved as: {fitplot_save_path}")
        
        # Create the labeling efficiency plot
        fig = npcs.plot_labeleff_for_auto_save(thresh=self.NPCsettings.SegmentThreshold_3D) # We are calling the alternative function: plot_labeleff_for_auto_save
                                                                                            # which returns the fig object, we then can use fig for automatic saving
        fig.savefig(leplot_save_path)
        print(f"Labeling efficiency plot automatically saved as: {leplot_save_path}")
        # --- End of Alex B addition ---


    def OnNPC3DSaveBySegments(self, event=None):
        pipeline = self.visFr.pipeline
        if findNPCset(pipeline) is not None:
            nbs = findNPCset(pipeline).n_bysegments()
            if nbs is None:
                warn("could not find npcs with by-segment fitting info, have you carried out fitting with recent code?")
                return
            with wx.FileDialog(self.visFr, 'Save NPC by-segment data as ...',
                                wildcard='CSV (*.csv)|*.csv',
                                style=wx.FD_SAVE) as fdialog:
                if fdialog.ShowModal() != wx.ID_OK:
                    return
                else:
                    fpath = fdialog.GetPath()
                    
            import pandas as pd
            df = pd.DataFrame.from_dict(dict(top=nbs['top'].flatten(),bottom=nbs['bottom'].flatten()))
            df.to_csv(fpath,index=False)
        else:
            warn("could not find valid NPC set, have you carried out fitting?")
    
    # --- Alex B addition ---
    
    def OnNPC3DSaveBySegments_auto_save(self, event=None, save_dir=None):
        pipeline = self.visFr.pipeline
        
                # --- Alex B addition ---
        # We define a few variables used for automatic saving later
        if save_dir is None:
            save_dir = os.getcwd()
        MINFLUXts = pipeline.mdh.get('MINFLUX.TimeStamp') # Get the timestamp and use it for naming the file to save
        if MINFLUXts is not None:
            NPC_segments_stats = f"{MINFLUXts}-NPC_segments.csv"
        else:
            NPC_segments_stats = "NPC_segments.csv"
        NPC_segments_stats_save_path = os.path.join(save_dir, NPC_segments_stats) # Save path for csv file
        # --- End of Alex B addition ---
        
        if findNPCset(pipeline) is not None:
            nbs = findNPCset(pipeline).n_bysegments()
            if nbs is None:
                warn("could not find npcs with by-segment fitting info, have you carried out fitting with recent code?")
                return
                    
            import pandas as pd
            df = pd.DataFrame.from_dict(dict(top=nbs['top'].flatten(),bottom=nbs['bottom'].flatten()))
            df.to_csv(NPC_segments_stats_save_path,index=False) # Alex B automatic saving
            print(f"NPC by-segment data automatically saved as: {NPC_segments_stats_save_path}")
            
        else:
            warn("could not find valid NPC set, have you carried out fitting?")
    
    # --- End of Alex B addition ---
        
    def OnNPC3DPlotBySegments(self, event=None):
        pipeline = self.visFr.pipeline
        if findNPCset(pipeline) is not None:
            nbs = findNPCset(pipeline).n_bysegments()
            if nbs is None:
                warn("could not find npcs with by-segment fitting info, have you carried out fitting with recent code?")
                return
            fig = plt.figure()
            plt.hist(nbs['bottom'].flatten(),bins=np.arange(nbs['bottom'].max()+2)-0.5,alpha=0.5,label='bottom',density=True,histtype='step')
            plt.plot([nbs['bottom'].mean(),nbs['bottom'].mean()],[0,0.2],'b--')
            plt.hist(nbs['top'].flatten(),bins=np.arange(nbs['top'].max()+2)-0.5,alpha=0.5,label='top',density=True,histtype='step')
            plt.plot([nbs['top'].mean(),nbs['top'].mean()],[0,0.2],'r--')
            plt.legend()
            nbsflat = nbs['bottom'].flatten()
            b_meanall = 0.5*nbsflat.mean() # 0.5 to get per site because we have two sites per segment
            b_meannz =  0.5*nbsflat[nbsflat>0].mean() # 0.5 to get per site because we have two sites per segment
            plt.text(0.65,0.5,'cytop per site: %.1f (%.1f per labeled)' % (b_meanall,b_meannz), horizontalalignment='center',
                     verticalalignment='center', color='b', transform=fig.transFigure)
            nbsflat = nbs['top'].flatten()
            t_meanall = 0.5*nbsflat.mean() # 0.5 to get per site because we have two sites per segment
            t_meannz =  0.5*nbsflat[nbsflat>0].mean() # 0.5 to get per site because we have two sites per segment
            plt.text(0.65,0.44,'nucleop per site: %.1f (%.1f per labeled)' % (t_meanall,t_meannz), horizontalalignment='center',
                     verticalalignment='center', color='r', transform=fig.transFigure)
        else:
            warn("could not find valid NPC set, have you carried out fitting?")
            
    # --- Alex B addition ---
    
    def OnNPC3DPlotBySegments_auto_save(self, event=None, save_dir=None):
        pipeline = self.visFr.pipeline
        
        # --- Alex B addition ---
        # We define a few variables used for automatic saving later
        
        if save_dir is None:
            save_dir = os.getcwd()
        MINFLUXts = pipeline.mdh.get('MINFLUX.TimeStamp') # Get the timestamp and use it for naming the file to save
        if MINFLUXts is not None:
            NPC_plot_segments = f"{MINFLUXts}-NPC_segments.png"
        else:
            NPC_plot_segments = "NPC_segments.png"
        NPC_plot_segments_save_path = os.path.join(save_dir, NPC_plot_segments) # Save path for csv file
        
        # --- End of Alex B addition ---
        
        if findNPCset(pipeline) is not None:
            nbs = findNPCset(pipeline).n_bysegments()
            if nbs is None:
                warn("could not find npcs with by-segment fitting info, have you carried out fitting with recent code?")
                return
            fig = plt.figure()
            plt.hist(nbs['bottom'].flatten(),bins=np.arange(nbs['bottom'].max()+2)-0.5,alpha=0.5,label='bottom',density=True,histtype='step')
            plt.plot([nbs['bottom'].mean(),nbs['bottom'].mean()],[0,0.2],'b--')
            plt.hist(nbs['top'].flatten(),bins=np.arange(nbs['top'].max()+2)-0.5,alpha=0.5,label='top',density=True,histtype='step')
            plt.plot([nbs['top'].mean(),nbs['top'].mean()],[0,0.2],'r--')
            plt.legend()
            nbsflat = nbs['bottom'].flatten()
            b_meanall = 0.5*nbsflat.mean() # 0.5 to get per site because we have two sites per segment
            b_meannz =  0.5*nbsflat[nbsflat>0].mean() # 0.5 to get per site because we have two sites per segment
            plt.text(0.65,0.5,'cytop per site: %.1f (%.1f per labeled)' % (b_meanall,b_meannz), horizontalalignment='center',
                        verticalalignment='center', color='b', transform=fig.transFigure)
            nbsflat = nbs['top'].flatten()
            t_meanall = 0.5*nbsflat.mean() # 0.5 to get per site because we have two sites per segment
            t_meannz =  0.5*nbsflat[nbsflat>0].mean() # 0.5 to get per site because we have two sites per segment
            plt.text(0.65,0.44,'nucleop per site: %.1f (%.1f per labeled)' % (t_meanall,t_meannz), horizontalalignment='center',
                        verticalalignment='center', color='r', transform=fig.transFigure)
            # --- Alex B addition ---
            # Save the NPC geometry stats plot
            fig.savefig(NPC_plot_segments_save_path)
            print(f"NPC geometry stats plot automatically saved as: {NPC_plot_segments_save_path}")
            # --- End of Alex B addition ---
        else:
            warn("could not find valid NPC set, have you carried out fitting?")
    
    # --- End of Alex B addition ---
        
    def On3DNPCaddGallery(self, event=None):
        pipeline = self.visFr.pipeline
        npcs = findNPCset(pipeline)
        if npcs is None or 'measurements' not in dir(npcs):
            warn('no valid NPC measurements found, therefore cannot add gallery...')
            return

        gallery_ds = None
        segment_ds = None
        npcmod = findNPCset(pipeline,return_mod=True)
        if npcmod is not None:
            gallery_ds = npcmod.outputGallery
            segment_ds = npcmod.outputSegments
        else:
            if 'npc_gallery' in pipeline.dataSources.keys():
                # this should not happen if we use a module, so only if somebody already manually made this ds
                warn("dataSource 'npc_gallery' already exists, currently we do not support recalculating a new gallery")
                return
            from PYMEcs.Analysis.NPC import mk_NPC_gallery
            gallery, segments = mk_NPC_gallery(npcs,self.NPCsettings.NPCGalleryArrangement,self.NPCsettings.Zclip_3D,self.NPCsettings.NPCRotationAngle)
            pipeline.addDataSource('npc_gallery',gallery,False) # should check if already exists!
            gallery_ds = 'npc_gallery'
            if self.NPCsettings.IncludeSegmentLinesWithGallery:
                pipeline.addDataSource('npc_segments',segments,False)
                segment_ds = 'npc_segments'
            pipeline.Rebuild() # check, is this the right way when add a new non-module based dataSource?

        if gallery_ds is not None and self.gallery_layer is None:
            from PYME.LMVis.layers.pointcloud import PointCloudRenderLayer
            self.gallery_layer = PointCloudRenderLayer(pipeline, dsname=gallery_ds, method='pointsprites', cmap='plasma', point_size=4.0, alpha=0.5, vertexColour='z')
            self.visFr.add_layer( self.gallery_layer)

        if segment_ds is not None and self.segment_layer is None:
            from PYME.LMVis.layers.tracks import TrackRenderLayer
            self.segment_layer = TrackRenderLayer(pipeline, dsname=segment_ds, method='tracks', clump_key='polyIndex', line_width=2.0,
                                     alpha=0.5,cmap='SolidWhite')
            self.visFr.add_layer(self.segment_layer)

    def On3DNPCaddTemplates(self, event=None):
        pipeline = self.visFr.pipeline

        npcs = findNPCset(pipeline)
        if npcs is None or 'measurements' not in dir(npcs):
            warn('no valid NPC measurements found, therefore cannot add templates...')
            return 

        ds_template_name = None
        npcmod = findNPCset(pipeline,return_mod=True)
        if npcmod is not None:
            ds_template_name = npcmod.outputTemplates
        else:
            warn("You are not using the NPCAnalysisInput module!\n\nWill attempt to create NPC templates, but session saving will not work")
            ds_template_name = 'NPCtemplates'
            if ds_template_name in pipeline.dataSources.keys():
                warn("dataSource '%s' already exists, we do not support recalculating a new template set, use the NPCAnalysisInput module instead" % ds_template_name)
                return
            from PYMEcs.Analysis.NPC import mk_npctemplates
            ds_template = mk_npctemplates(npcs)
            pipeline.addDataSource(ds_template_name,ds_template,False) # should check if already exists!
            pipeline.Rebuild() # check, is this the right way when add a new non-module based dataSource?

        # now we add a track layer to render our template polygons
        # TODO - we may need to check if this happened before or not!
        from PYME.LMVis.layers.tracks import (
            TrackRenderLayer,  # NOTE: we may rename the clumpIndex variable in this layer to polyIndex or similar
        )
        layer = TrackRenderLayer(pipeline, dsname=ds_template_name, method='tracks', clump_key='polyIndex', line_width=2.0, alpha=0.5)
        self.visFr.add_layer(layer)        

    def OnNPC3DSaveMeasurements(self, event=None):
        import pandas as pd
        pipeline = self.visFr.pipeline
        npcs = findNPCset(pipeline)
        if npcs is None or 'measurements' not in dir(npcs):
            warn('no valid NPC measurements found, therefore cannot save...')
            return 
        
        fdialog = wx.FileDialog(self.visFr, 'Save NPC measurements as ...',
                                wildcard='CSV (*.csv)|*.csv',
                                style=wx.FD_SAVE)
        if fdialog.ShowModal() != wx.ID_OK:
            return

        fpath = fdialog.GetPath()
        meas = np.array(npcs.measurements, dtype='i')

        df = pd.DataFrame({'Ntop_NPC3D': meas[:, 0], 'Nbot_NPC3D': meas[:, 1]})
        entries = len(np.unique(pipeline.objectID))
        MINFLUX_filename = [pipeline.mdh.getEntry('MINFLUX.Filename')]*entries
        NPC_threshold = [self.NPCsettings.SegmentThreshold_3D]*entries
        fshortening = [pipeline.mdh.getEntry('MINFLUX.Foreshortening')]*entries
        t_min = [np.round(np.min(pipeline.tim))]*entries
        t_max = [np.round(np.max(pipeline.tim))]*entries
        t_maxhr = [np.round(ti/3600) for ti in t_max]

        duration_hours = [np.round((np.max(pipeline.tim)-np.min(pipeline.tim))/3600,2)]*entries
        duration_hours_rounded = [np.round(duration) for duration in duration_hours]
        dim_z_nm = [np.round(np.max(pipeline.z)-np.min(pipeline.z),2)]*entries
        NPCLE = (meas[:,0]+meas[:,1])/16
        
        unique_ids, counts = np.unique(pipeline.objectID, return_counts=True)
        nEvents = counts.tolist()  
        
        df = pd.DataFrame({'Ntop_NPC3D': meas[:, 0], 'Nbot_NPC3D': meas[:, 1], "NPC_LE": NPCLE,
                            'objectID': np.unique(pipeline.objectID),
                            'diameters': np.round(pipeline.npcs.diam(), 4),
                            'heights': np.round(pipeline.npcs.height(), 4),
                            't_min': t_min,
                            't_max': t_max,
                            't_maxhr': t_maxhr,
                            'duration_hours': duration_hours,
                            'duration_hours_rounded': duration_hours_rounded,
                            'nEvents': nEvents,
                            'dim_z_nm': dim_z_nm,
                            'foreshortening': fshortening,
                            'NPC_threshold': NPC_threshold,                            
                            'filename': MINFLUX_filename})

        pipeline.selectDataSource('Localizations')
        dim_x = (np.max(pipeline.x)-np.min(pipeline.x))/1000
        dim_y = (np.max(pipeline.y)-np.min(pipeline.y))/1000
        area_xy = [np.round(dim_x*dim_y,2)]*entries

        df.insert(3, 'area_xy', area_xy)

        try:
            pipeline.selectDataSource('valid_npcs')
            print("Pipeline data source successfully changed")
        except:
            print("Unable to change the pipeline to 'valid_npcs'")

        df.to_csv(fpath,index=False, mode='a')

    # --- Alex B addition for auto-save of measurements ---

    def OnNPC3DSaveMeasurements_auto_save(self, event=None, save_dir=None):
        import pandas as pd
        pipeline = self.visFr.pipeline
        npcs = findNPCset(pipeline)
        if npcs is None or 'measurements' not in dir(npcs):
            warn('no valid NPC measurements found, therefore cannot save...')
            return 
        
        # --- Alex B addition ---
        # We define a few variables used for automatic saving later
        
        if save_dir is None:
            save_dir = os.getcwd()
        MINFLUXts = pipeline.mdh.get('MINFLUX.TimeStamp') # Get the timestamp and use it for naming the file to save
        if MINFLUXts is not None:
            csv_filename = f"{MINFLUXts}-LE_stats.csv"
        else:
            csv_filename = "LE_stats.csv"
        csv_save_path = os.path.join(save_dir, csv_filename) # Save path for csv file
        
        # --- End of Alex B addition ---


        meas = np.array(npcs.measurements, dtype='i')

        df = pd.DataFrame({'Ntop_NPC3D': meas[:, 0], 'Nbot_NPC3D': meas[:, 1]})
        entries = len(np.unique(pipeline.objectID))
        MINFLUX_filename = [pipeline.mdh.getEntry('MINFLUX.Filename')]*entries
        NPC_threshold = [self.NPCsettings.SegmentThreshold_3D]*entries
        fshortening = [pipeline.mdh.getEntry('MINFLUX.Foreshortening')]*entries
        t_min = [np.round(np.min(pipeline.tim))]*entries
        t_max = [np.round(np.max(pipeline.tim))]*entries
        t_maxhr = [np.round(ti/3600) for ti in t_max]

        duration_hours = [np.round((np.max(pipeline.tim)-np.min(pipeline.tim))/3600,2)]*entries
        duration_hours_rounded = [np.round(duration) for duration in duration_hours]
        dim_z_nm = [np.round(np.max(pipeline.z)-np.min(pipeline.z),2)]*entries
        NPCLE = (meas[:,0]+meas[:,1])/16
        
        unique_ids, counts = np.unique(pipeline.objectID, return_counts=True)
        nEvents = counts.tolist()  
        
        df = pd.DataFrame({'Ntop_NPC3D': meas[:, 0], 'Nbot_NPC3D': meas[:, 1], "NPC_LE": NPCLE,
                            'objectID': np.unique(pipeline.objectID),
                            'diameters': np.round(pipeline.npcs.diam(), 4),
                            'heights': np.round(pipeline.npcs.height(), 4),
                            't_min': t_min,
                            't_max': t_max,
                            't_maxhr': t_maxhr,
                            'duration_hours': duration_hours,
                            'duration_hours_rounded': duration_hours_rounded,
                            'nEvents': nEvents,
                            'dim_z_nm': dim_z_nm,
                            'foreshortening': fshortening,
                            'NPC_threshold': NPC_threshold,                            
                            'filename': MINFLUX_filename})

        pipeline.selectDataSource('Localizations')
        dim_x = (np.max(pipeline.x)-np.min(pipeline.x))/1000
        dim_y = (np.max(pipeline.y)-np.min(pipeline.y))/1000
        area_xy = [np.round(dim_x*dim_y,2)]*entries

        df.insert(3, 'area_xy', area_xy)

        try:
            pipeline.selectDataSource('valid_npcs')
            print("Pipeline data source successfully changed")
        except:
            print("Unable to change the pipeline to 'valid_npcs'")

        df.to_csv(csv_save_path,index=False, mode='a')

# --- End of Alex B addition for auto-save ---


    def OnNPC3DLoadMeasurements(self, event=None):
        from PYMEcs.misc.utils import get_timestamp_from_filename
        fdialog = wx.FileDialog(self.visFr, 'Load NPC measurements from ...',
                                wildcard='CSV (*.csv)|*.csv',
                                style=wx.FD_OPEN)
        if fdialog.ShowModal() != wx.ID_OK:
            return
        import pandas as pd
        fname = fdialog.GetPath()
        meas = pd.read_csv(fname,comment='#')
        # here we need some plotting code
        nlab = meas['Ntop_NPC3D'] + meas['Nbot_NPC3D']
        plt.figure()
        plotcdf_npc3d(nlab,timestamp=get_timestamp_from_filename(fname))

    def OnNPC3DLoadNPCSet(self, event=None):
        from PYMEcs.IO.NPC import load_NPC_set
        
        pipeline = self.visFr.pipeline
        fdialog = wx.FileDialog(self.visFr, 'Load NPC measurements from ...',
                                wildcard='Pickle (*.pickle)|*.pickle',
                                style=wx.FD_OPEN)
        if fdialog.ShowModal() != wx.ID_OK:
            return

        pipeline.npcs = load_NPC_set(fdialog.GetPath(),
                                     ts=pipeline.mdh.get('MINFLUX.TimeStamp'),
                                     foreshortening=pipeline.mdh.get('MINFLUX.Foreshortening',1.0))

    def OnNPC3DSaveNPCSet(self, event=None):
        from PYMEcs.IO.NPC import save_NPC_set
        pipeline = self.visFr.pipeline
        defaultFile = ''
        MINFLUXts = pipeline.mdh.get('MINFLUX.TimeStamp')
        if MINFLUXts is not None:
            defaultFile = "%s-NPCset.pickle" % MINFLUXts
        npcs = findNPCset(pipeline)
        if npcs is None:
            warn('no valid NPC Set found, therefore cannot save...')
            return
        fdialog = wx.FileDialog(self.visFr, 'Save NPC Set as ...',
                                wildcard='Pickle (*.pickle)|*.pickle',
                                defaultFile=defaultFile,
                                style=wx.FD_SAVE)
        if fdialog.ShowModal() != wx.ID_OK:
            return
        
        save_NPC_set(npcs,fdialog.GetPath())
        
# --- Alex B addition ---
# AIM: perform all actions 3D NPC actions and save output automatically

    def OnNPC3DSaveNPCSet_auto_save(self, event=None, save_dir=None): 
        """Automatically save the NPC set to a default file path without user dialog."""

        from PYMEcs.IO.NPC import save_NPC_set

        pipeline = self.visFr.pipeline

        # Get the MINFLUX timestamp from the pipeline metadata, if available
        MINFLUXts = pipeline.mdh.get('MINFLUX.TimeStamp')
        # Construct the default filename using the timestamp if present
        if MINFLUXts is not None:
            defaultFile = f"{MINFLUXts}-NPCset.pickle"
        else:
            defaultFile = "NPCset.pickle"
        # Save in the current directory (same as the session file
        if save_dir is None:
            save_dir = os.getcwd()
        save_path = os.path.join(save_dir, defaultFile)
        # Find the current NPC set in the pipeline
        npcs = findNPCset(pipeline)
        print(f"Attempting to automatically save NPC Set to: {save_path}.")
        # If no NPC set is found, warn and exit
        if npcs is None:
            warn('no valid NPC Set found, therefore cannot save...')
            return
        # Save the NPC set to the constructed path
        save_NPC_set(npcs, save_path)

# --- End of Alex B addition ---

    def OnNPC3DSaveGeometryStats(self,event=None):
        pipeline = self.visFr.pipeline
        npcs = findNPCset(pipeline)
        if npcs is None:
            warn('no valid NPC measurements found, thus no geometry info available...')
            return
        diams = np.asarray(npcs.diam())
        heights = np.asarray(npcs.height())
        import pandas as pd
        geo_df = pd.DataFrame.from_dict(dict(diameter=diams,height=heights))
        with wx.FileDialog(self.visFr, 'Save NPC measurements as ...',
                                wildcard='CSV (*.csv)|*.csv',
                                style=wx.FD_SAVE) as fdialog:
            if fdialog.ShowModal() != wx.ID_OK:
                return
            fpath = fdialog.GetPath()
 
        geo_df.to_csv(fpath,index=False)
        
# --- Alex B addition ---

    def OnNPC3DSaveGeometryStats_auto_save(self,event=None, save_dir=None):
        pipeline = self.visFr.pipeline
        npcs = findNPCset(pipeline)
        if npcs is None:
            warn('no valid NPC measurements found, thus no geometry info available...')
            return
        diams = np.asarray(npcs.diam())
        heights = np.asarray(npcs.height())
        import pandas as pd
        geo_df = pd.DataFrame.from_dict(dict(diameter=diams,height=heights))

        if save_dir is None:
            save_dir = os.getcwd()
        MINFLUXts = pipeline.mdh.get('MINFLUX.TimeStamp') # Get the timestamp and use it for naming the file to save
        if MINFLUXts is not None:
            geo_stats_csv = f"{MINFLUXts}-NPC_geometry_stats.csv"
        else:
            geo_stats_csv = "NPC_geometry_stats.csv"
        geo_df.to_csv(os.path.join(save_dir, geo_stats_csv), index=False)

# --- End of Alex B addition ---
        
     
    def OnNPC3DGeometryStats(self,event=None):
        pipeline = self.visFr.pipeline
        npcs = findNPCset(pipeline)
        if npcs is None:
            warn('no valid NPC measurements found, thus no geometry info available...')
            return
        import pandas as pd

        from PYMEcs.misc.matplotlib import boxswarmplot, figuredefaults
        diams = np.asarray(npcs.diam())
        heights = np.asarray(npcs.height())
        geo_df = pd.DataFrame.from_dict(dict(diameter=diams,height=heights))
        figuredefaults(fontsize=12)
        plt.figure()
        from scipy.stats import iqr
        iqrh = iqr(heights)
        sdh = np.std(heights)
        iqrd = iqr(diams)
        sdd = np.std(diams)
        bp = boxswarmplot(geo_df,width=0.35,annotate_medians=True,annotate_means=True,showmeans=True,swarmalpha=0.4,swarmsize=4)
        plt.text(0.0,50,"IQR %.1f\nSD %.1f" % (iqrd,sdd), horizontalalignment='center')
        plt.text(1.0,120,"IQR %.1f\nSD %.1f" % (iqrh,sdh), horizontalalignment='center')
        # res = plt.boxplot([diams,heights],showmeans=True,labels=['diameter','height'])
        plt.title("%d NPCs, mean diam %.0f nm, mean ring spacing %.0f nm" % (heights.size,diams.mean(),heights.mean()), fontsize=11)
        plt.ylim(0,150)

# --- Alex B addition ---

    def OnNPC3DGeometryStats_auto_save(self,event=None, save_dir=None):
        pipeline = self.visFr.pipeline
        npcs = findNPCset(pipeline)
        if npcs is None:
            warn('no valid NPC measurements found, thus no geometry info available...')
            return
        import pandas as pd

        # --- Alex B addition ---
        # We define a few variables used for automatic saving later
        
        if save_dir is None:
            save_dir = os.getcwd()
        MINFLUXts = pipeline.mdh.get('MINFLUX.TimeStamp') # Get the timestamp and use it for naming the file to save
        if MINFLUXts is not None:
            geom_stats_fig = f"{MINFLUXts}-Geom_stats.png"
        else:
            geom_stats_fig = "Geom_stats.png"
        geom_stats_fig_save_path = os.path.join(save_dir, geom_stats_fig) # Save path for csv file
        
        # --- End of Alex B addition ---

        from PYMEcs.misc.matplotlib import boxswarmplot, figuredefaults
        diams = np.asarray(npcs.diam())
        heights = np.asarray(npcs.height())
        geo_df = pd.DataFrame.from_dict(dict(diameter=diams,height=heights))
        figuredefaults(fontsize=12)
        plt.figure()
        from scipy.stats import iqr
        iqrh = iqr(heights)
        sdh = np.std(heights)
        iqrd = iqr(diams)
        sdd = np.std(diams)
        bp = boxswarmplot(geo_df,width=0.35,annotate_medians=True,annotate_means=True,showmeans=True,swarmalpha=0.4,swarmsize=4)
        plt.text(0.0,50,"IQR %.1f\nSD %.1f" % (iqrd,sdd), horizontalalignment='center')
        plt.text(1.0,120,"IQR %.1f\nSD %.1f" % (iqrh,sdh), horizontalalignment='center')
        # res = plt.boxplot([diams,heights],showmeans=True,labels=['diameter','height'])
        plt.title("%d NPCs, mean diam %.0f nm, mean ring spacing %.0f nm" % (heights.size,diams.mean(),heights.mean()), fontsize=11)
        plt.ylim(0,150)
        # --- Alex B addition ---
        # Save the geometry stats plot
        plt.savefig(geom_stats_fig_save_path)
        print(f"Geometry stats plot automatically saved as: {geom_stats_fig_save_path}")

# --- End of Alex B addition ---

    def OnNPC3DTemplateFitStats(self,event=None):
        pipeline = self.visFr.pipeline
        npcs = findNPCset(pipeline)
        if npcs is None:
            warn('no valid NPC measurements found, thus no geometry info available...')
            return
        import pandas as pd

        from PYMEcs.misc.matplotlib import boxswarmplot, figuredefaults
        id = [npc.objectID for npc in npcs.npcs] # not used right now
        llperloc = [npc.opt_result.fun/npc.npts.shape[0] for npc in npcs.npcs]
        ll_df = pd.DataFrame.from_dict(dict(llperloc=llperloc))
        figuredefaults(fontsize=12)
        plt.figure()
        bp = boxswarmplot(ll_df,width=0.35,annotate_medians=True,annotate_means=True,showmeans=True,swarmalpha=0.6,swarmsize=5)
        plt.title("NPC neg-log-likelihood per localization")

# --- Alex B addition ---
# Original function 'OnNPC3DTemplateFitStats', copied and modified for automatic saving.

    def OnNPC3DTemplateFitStats_auto_save(self,event=None, save_dir=None):
        pipeline = self.visFr.pipeline
        npcs = findNPCset(pipeline)
        if npcs is None:
            warn('no valid NPC measurements found, thus no geometry info available...')
            return
        import pandas as pd

        # --- Alex B addition ---
        # We define a few variables used for automatic saving later
        
        if save_dir is None:
            save_dir = os.getcwd()
        MINFLUXts = pipeline.mdh.get('MINFLUX.TimeStamp') # Get the timestamp and use it for naming the file to save
        if MINFLUXts is not None:
            template_fit_stats_fig = f"{MINFLUXts}-Template_fit_stats.png"
        else:
            template_fit_stats_fig = "Template_fit_stats.png"
        template_fit_stats_fig_save_path = os.path.join(save_dir, template_fit_stats_fig) # Save path for csv file

        # --- End of Alex B addition ---

        from PYMEcs.misc.matplotlib import boxswarmplot, figuredefaults
        id = [npc.objectID for npc in npcs.npcs] # not used right now
        llperloc = [npc.opt_result.fun/npc.npts.shape[0] for npc in npcs.npcs]
        ll_df = pd.DataFrame.from_dict(dict(llperloc=llperloc))
        figuredefaults(fontsize=12)
        plt.figure()
        bp = boxswarmplot(ll_df,width=0.35,annotate_medians=True,annotate_means=True,showmeans=True,swarmalpha=0.6,swarmsize=5)
        plt.title("NPC neg-log-likelihood per localization")

        # --- Alex B addition ---
        # Save the template fit stats plot
        plt.savefig(template_fit_stats_fig_save_path)
        print(f"Template fit stats plot automatically saved as: {template_fit_stats_fig_save_path}")

# --- End of Alex B addition ---


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
