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
    from PYMEcs.IO.MINFLUX import foreshortening
    SegmentThreshold_2D = Int(10,label='Threshold localisation count (2D)',
                              desc="a minimum number of localisations in a segment in 2D fitting to be counted as labeled; "+
                              "NOTE: definition changed from original code where 11 localisations where required with a threshold of 10 etc!")
    SegmentThreshold_3D = Int(1,label='Threshold localisation count (3D)',
                              desc="a minimum number of localisations in a segment in 3D fitting to be counted as labeled")
    SecondPass_2D = Bool(False,label='Second pass for NPC fitting (2D)',
                         desc="a second pass for 2D fitting should be run; we have experimented with a second pass rotation "+
                         "estimate and fitting hoping to improve on the first estimate; still experimental")
    StartHeight_3D = Float(70.0*foreshortening,label='Starting ring spacing for 3D fitting',
                           desc="starting ring spacing value for the 3D fit; note only considered when doing the initial full fit; "+
                           "not considered when re-evaluating existing fit")
    StartDiam_3D = Float(107.0,label='Starting ring diameter for 3D fitting',
                           desc="starting ring diameter value for the 3D fit; note only considered when doing the initial full fit; "+
                           "not considered when re-evaluating existing fit")
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
    Zclip_3D = Float(75.0*foreshortening,label='Z-clip value from center of NPC',
                     desc='the used zrange from the (estimated) center of the NPC, from (-zclip..+zclip) in 3D fitting')
    OffsetMode_3D = Enum(['median','mean'],label='Method to estimate NPC center',
                         desc="Method to estimate the likely 3D center of an NPC; median seems more robust against outliners (say in z)")


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
        visFr.AddMenuItem('Experimental>NPC2D', 'NPC Analysis settings', self.OnNPCsettings)
        visFr.AddMenuItem('Experimental>NPC3D', 'NPC Analysis settings', self.OnNPCsettings)

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

        if 'npcs' in dir(pipeline) and pipeline.npcs is not None:
            npcs = pipeline.npcs
            do_plot = False
        else:
            from PYMEcs.IO.MINFLUX import foreshortening
            npcs = NPC3DSet(filename=pipeline_filename(pipeline),
                            zclip=self.NPCsettings.Zclip_3D,
                            offset_mode=self.NPCsettings.OffsetMode_3D,
                            NPCdiam=self.NPCsettings.StartDiam_3D,
                            NPCheight=self.NPCsettings.StartHeight_3D,
                            foreshortening=foreshortening)
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
            nt,nb = npc.nlabeled(nthresh=self.NPCsettings.SegmentThreshold_3D,
                                 dr=self.NPCsettings.RadiusUncertainty_3D,
                                 rotlocked=self.NPCsettings.RotationLocked_3D,
                                 zrange=self.NPCsettings.Zclip_3D)
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

        if cancelled:
            return
        
        pipeline.npcs = npcs # overwriting with the same object should be fine if pipeline.npcs already existed
        npcs.plot_labeleff()


    def On3DNPCaddTemplates(self, event=None):
        pipeline = self.visFr.pipeline
        if 'npcs' not in dir(pipeline) or 'measurements' not in dir(pipeline.npcs):
            warn('no valid NPC measurements found, therefore cannot add templates...')
            return
        npcs = pipeline.npcs

        if 'NPCtemplates' in pipeline.dataSources.keys():
            # TODO: if practically this becomes an issue learn how to update an existing datasource
            warn("dataSource 'NPCtemplates' already exists, currently we do not support recalculating a new template set")
            return
        
        x = np.empty((0))
        y = np.empty((0))
        z = np.empty((0))
        polyIndex = np.empty((0),int)
        polySize = np.empty((0),int)
        objectID = np.empty((0),int)
        NtopLabelled = np.empty((0),int)
        NbotLabelled = np.empty((0),int)
        NLabelled = np.empty((0),int)
        diams = np.empty((0),float)
        heights = np.empty((0),float)
        ci = 1
        for npc in npcs.npcs:
            nt, nb = (npc.n_top,npc.n_bot)
            glyph = npc.get_glyph()
            pars = npc.opt_result.x
            diam = npc.get_glyph_diam() / (0.01*pars[5])
            height = npc.get_glyph_height() / (0.01*pars[6])
            for poly in ['circ_bot','circ_top','axis']:
                c3 = glyph[poly]
                xg = c3[:,0]
                yg = c3[:,1]
                zg = c3[:,2]
                x = np.append(x,xg)
                y = np.append(y,yg)
                z = np.append(z,zg)
                polyIndex = np.append(polyIndex,np.full_like(xg,ci,dtype=int))
                polySize = np.append(polySize,np.full_like(xg,xg.size,dtype=int))
                ci += 1
                objectID = np.append(objectID,np.full_like(xg,npc.objectID,dtype=int))
                NtopLabelled = np.append(NtopLabelled,np.full_like(xg,nt,dtype=int))
                NbotLabelled = np.append(NbotLabelled,np.full_like(xg,nb,dtype=int))
                NLabelled = np.append(NLabelled,np.full_like(xg,nt+nb,dtype=int))               
                diams = np.append(diams,np.full_like(xg,diam,dtype=float))
                heights = np.append(heights,np.full_like(xg,height,dtype=float))
        t = np.arange(x.size)
        A = np.full_like(x,10.0,dtype='f')
        error_x = np.full_like(x,1.0,dtype='f')
        error_y = np.full_like(x,1.0,dtype='f')
        error_z = np.full_like(x,1.0,dtype='f')
        
        dsdict = dict(x=x,y=y,z=z,polyIndex=polyIndex,polySize=polySize,
                      NtopLabelled=NtopLabelled,NbotLabelled=NbotLabelled,NLabelled=NLabelled,
                      objectID=objectID,t=t,A=A,
                      error_x=error_x,error_y=error_y,error_z=error_z,
                      npc_height=heights,npc_diam=diams)

        from PYME.IO.tabular import DictSource
        pipeline.addDataSource('NPCtemplates',DictSource(dsdict),False) # should check if already exists!
        pipeline.Rebuild() # check, is this the right way when add a new non-module based dataSource?

        # now we add a track layer to render our template polygons
        # TODO - we may need to check if this happened before or not!
        from PYME.LMVis.layers.tracks import TrackRenderLayer # NOTE: we may rename the clumpIndex variable in this layer to polyIndex or similar
        layer = TrackRenderLayer(pipeline, dsname='NPCtemplates', method='tracks', clump_key='polyIndex', line_width=2.0, alpha=0.5)
        self.visFr.add_layer(layer)        

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
                    (self.NPCsettings.SegmentThreshold_3D,Path(pipeline_filename(pipeline)).name))

        df.to_csv(fpath,index=False, mode='a')


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
        import pickle
        from pathlib import Path
        pipeline = self.visFr.pipeline
        fdialog = wx.FileDialog(self.visFr, 'Load NPC measurements from ...',
                                wildcard='Pickle (*.pickle)|*.pickle',
                                style=wx.FD_OPEN)
        if fdialog.ShowModal() != wx.ID_OK:
            return
        fname = Path(fdialog.GetPath())
        MINFLUXts = pipeline.mdh.get('MINFLUX.TimeStamp')
        if MINFLUXts is not None and MINFLUXts not in fname.stem:
            warn("MINFLUX time stamp not in NPC dataset filename; check correct filename chosen: %s vs %s" %
                 (MINFLUXts,fname.stem))
        with open(fname,'rb') as fi:
            pipeline.npcs=pickle.load(fi)

    def OnNPC3DSaveNPCSet(self, event=None):
        import pickle
        pipeline = self.visFr.pipeline
        defaultFile = None
        MINFLUXts = pipeline.mdh.get('MINFLUX.TimeStamp')
        if MINFLUXts is not None:
            defaultFile = "%s-NPCset.pickle" % MINFLUXts
        if 'npcs' not in dir(pipeline):
            warn('no valid NPC Set found, therefore cannot save...')
            return
        fdialog = wx.FileDialog(self.visFr, 'Save NPC Set as ...',
                                wildcard='Pickle (*.pickle)|*.pickle',
                                defaultFile=defaultFile,
                                style=wx.FD_SAVE)
        if fdialog.ShowModal() != wx.ID_OK:
            return

        fpath = fdialog.GetPath()
        with open(fpath, "wb") as file:
            pickle.dump(pipeline.npcs,file)

    def OnNPC3DGeometryStats(self,event=None):
        pipeline = self.visFr.pipeline
        if 'npcs' not in dir(pipeline) or pipeline.npcs is None:
            warn('no valid NPC measurements found, thus no geometry info available...')
            return
        diams = np.asarray(pipeline.npcs.diam())
        heights = np.asarray(pipeline.npcs.height())
        plt.figure()
        res = plt.boxplot([diams,heights],showmeans=True,labels=['diameter','height'])
        plt.title("NPC mean diam %.0f nm, mean ring spacing %.0f nm" % (diams.mean(),heights.mean()))
        plt.ylim(0,150)

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
