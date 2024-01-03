from PYME.warnings import warn
from PYMEcs.Analysis.NPC import estimate_nlabeled, npclabel_fit
from PYME.recipes import tablefilters
import wx
from traits.api import HasTraits, Str, Int, CStr, List, Enum, Float, Bool
import numpy as np
import matplotlib.pyplot as plt

def unique_name(stem,names):
    if stem not in names:
        return stem
    for i in range(1,11):
        stem2 = "%s%d" % (stem,i)
        if stem2 not in names:
            return stem2

    return stem2 # here we just give up and accept a duplicate name

def selectWithDialog(choices, message='select image from list', caption='Selection'):
    dlg = wx.SingleChoiceDialog(None, message, caption, list(choices), wx.CHOICEDLG_STYLE)
    if dlg.ShowModal() == wx.ID_OK:
        item = dlg.GetStringSelection()
    else:
        item = None
    dlg.Destroy()
    return item

class NPCsettings(HasTraits):
    SegmentThreshold = Int(10)
    SecondPass = Bool(False)
    FitMode = Enum(['abs','square'])

class NPCcalc():
    def __init__(self, visFr):
        self.visFr = visFr

        visFr.AddMenuItem('Experimental>NPCs', "Analyse single NPC", self.OnAnalyseSingleNPC)
        visFr.AddMenuItem('Experimental>NPCs', "Select NPCs by Mask", self.OnSelectNPCsByMask)
        visFr.AddMenuItem('Experimental>NPCs', "Analyse NPCs by ID", self.OnAnalyseNPCsByID)
        visFr.AddMenuItem('Experimental>NPCs', "Show NPC labeling Statistics", self.OnNPCstats)
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

        if (xExtent > 300) or (yExtent > 300):
            warn('x or y bounding box > 200 nm (%d,%d) - check if single NPC' % (xExtent,yExtent))
            return

        estimate_nlabeled(pipeline['x'],pipeline['y'],nthresh=self.NPCsettings.SegmentThreshold,
                          do_plot=True,secondpass=self.NPCsettings.SecondPass,fitmode=self.NPCsettings.FitMode)
        

    def OnAnalyseNPCsByID(self, event):
        from PYMEcs.recipes.localisations import NPCAnalysisByID
        pipeline = self.visFr.pipeline
        recipe = pipeline.recipe
        with_npclinfo = unique_name('with_npcinfo',pipeline.dataSources.keys())
        recipe.add_module(NPCAnalysisByID(inputName=pipeline.selectedDataSourceKey,outputName=with_npclinfo))
        recipe.execute()

        pipeline.selectDataSource(with_npclinfo)


    def OnSelectNPCsByMask(self,event):
        from PYME.DSView import dsviewer

        pipeline = self.visFr.pipeline
            
        selection = selectWithDialog(dsviewer.openViewers.keys())
        if selection is not None:
            img = dsviewer.openViewers[selection].image

        if img is None:
            return
        
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
                                                   filters={'objectID' : [.5, 1e5]}))
        recipe.execute()
        
        pipeline.selectDataSource(valid_ids)

    def OnNPCstats(self,event):
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

        
def Plug(visFr):
    """Plugs this module into the gui"""
    NPCcalc(visFr)
