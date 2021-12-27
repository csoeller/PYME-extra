from traitsui.api import View, Item, Group
from traitsui.menu import OKButton, CancelButton, OKCancelButtons
from traits.api import HasTraits, Str, Int, CStr, List, Enum, Float, Bool


#from PYME.DSView.dsviewer import View3D
#from PYME.IO.MetaDataHandler import NestedClassMDHandler

import numpy as np
import wx

def rendHistw(x,y, imageBounds, pixelSize, weights=None):
    X = np.arange(imageBounds.x0, imageBounds.x1 + 1.01*pixelSize, pixelSize)
    Y = np.arange(imageBounds.y0, imageBounds.y1 + 1.01*pixelSize, pixelSize)
    
    im, edx, edy = np.histogram2d(x,y, bins=(X,Y), weights=weights)

    return im


from PYME.LMVis.renderers import ColourRenderer
class BackgroundHistogramRenderer(ColourRenderer):
    """2D background histogram rendering"""

    name = 'Histogram Background'
    mode = 'histogram'

    def genIm(self, settings, imb, mdh):
        im_hist = rendHistw(self.colourFilter['x'],self.colourFilter['y'], imb, settings['pixelSize'])
        im_background = rendHistw(self.colourFilter['x'],self.colourFilter['y'], imb, settings['pixelSize'],
                                  weights=self.colourFilter['fitResults_background'])
        Inds = im_hist > 0
        im_bgmean = np.zeros_like(im_background)
        im_bgmean[Inds] = im_background[Inds]/im_hist[Inds]

        return im_bgmean

def voxelsize_mdh(pixelsize_nm):
    from PYME.IO.MetaDataHandler import NestedClassMDHandler
    
    mdh = NestedClassMDHandler()

    #voxelsize
    mdh.setEntry('voxelsize.x',1e-3*pixelsize_nm)
    mdh.setEntry('voxelsize.y',1e-3*pixelsize_nm)
    mdh.setEntry('voxelsize.z',0.2)
    mdh.setEntry('voxelsize.units','um')

    return mdh

class InterpolChoice(HasTraits):
    InterpolationMethod = Enum(['Splines','Radial Basis functions','Constraint Model'])

import wx
#import deClump
#from pylab import *
#import numpy as np

class simplerFilterDialog(wx.Dialog):
    def __init__(self, *args, **kwargs):
        wx.Dialog.__init__(self, *args, **kwargs)

        vsizer = wx.BoxSizer(wx.VERTICAL)

        hsizer = wx.BoxSizer(wx.HORIZONTAL)

        hsizer.Add(wx.StaticText(self, -1, 'Clump Radius: '), 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)

        self.tClumpRadMult = wx.TextCtrl(self, -1, '2.0', size=[30,-1])
        hsizer.Add(self.tClumpRadMult, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)

        hsizer.Add(wx.StaticText(self, -1, 'X'), 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)

        self.cClumpRadVar = wx.Choice(self, -1, choices=['1.0', 'error_x'])
        self.cClumpRadVar.SetSelection(1)
        hsizer.Add(self.cClumpRadVar,1, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
        
        vsizer.Add(hsizer, 0, wx.ALL, 5)
        
        hsizer = wx.BoxSizer(wx.HORIZONTAL)
        hsizer.Add(wx.StaticText(self, -1, 'Minimum Clump Size (post filter): '), 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)

        self.tClumpMin = wx.TextCtrl(self, -1, '3')
        hsizer.Add(self.tClumpMin,1, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)

        vsizer.Add(hsizer, 0, wx.ALL, 5)

        btSizer = wx.StdDialogButtonSizer()

        btn = wx.Button(self, wx.ID_OK)
        btn.SetDefault()

        btSizer.AddButton(btn)

        btn = wx.Button(self, wx.ID_CANCEL)

        btSizer.AddButton(btn)

        btSizer.Realize()

        vsizer.Add(btSizer, 0, wx.ALL, 5)

        self.SetSizerAndFit(vsizer)

    def GetClumpRadiusMultiplier(self):
        return float(self.tClumpRadMult.GetValue())

    def GetClumpRadiusVariable(self):
        return self.cClumpRadVar.GetStringSelection()

    def GetClumpMinSize(self):
        return int(self.tClumpMin.GetValue())


    
# things still to add
# - done -  "Filter events for SIMPLER": choose clump radius and min clumpsize interactively
# - new function to calculate N0 from data background
# - done - new function to insert SIMPLERzgenerator
# - add error estimates for SIMPLER (combination of formula and/or stddev) -> do in recipe module
class SIMPLER():
    def __init__(self, visFr):
        self.visFr = visFr
        self.intMeth = InterpolChoice()
        # pass visFr = None in constructor to avoid autogenerated menu entry
        self.bgrend = BackgroundHistogramRenderer(visFr=None,pipeline=visFr.pipeline,mainWind=visFr)
        # note we need to set the visFr posthoc so that rendering works as expected
        self.bgrend.visFr = visFr
        
        visFr.AddMenuItem('Experimental>SIMPLER', "Filter events for SIMPLER", self.OnFilterSIMPLER)
        visFr.AddMenuItem('Experimental>SIMPLER', "Generate z estimates with SIMPLER", self.OnGenerateZSIMPLER)
        
        visFr.AddMenuItem('Experimental>SIMPLER', "Cluster modes for SIMPLER N0 distribution", self.OnModesN0SIMPLER)
        visFr.AddMenuItem('Experimental>SIMPLER', "Interpolate N0 Map", self.OnInterpolateN0Field)
        visFr.AddMenuItem('Experimental>SIMPLER', "Generate N0 Image Map from Background", self.bgrend.GenerateGUI)
        visFr.AddMenuItem('Experimental>SIMPLER', "Set N0 from Image Map", self.OnN0FromImage)
        visFr.AddMenuItem('Experimental>SIMPLER', "Set N0 from Interpolation Map", self.OnN0FromIinterpolationField)
        
        visFr.AddMenuItem('Experimental>SIMPLER', "Show N0 Map from Interpolation", self.OnShowN0Map)
        
    def OnFilterSIMPLER(self, event):
        from PYME.recipes.tablefilters import FilterTable
        from PYME.recipes.tracking import FindClumps
        recipe = self.visFr.pipeline.recipe
        pipeline = self.visFr.pipeline

        dlg = simplerFilterDialog(None)
        ret = dlg.ShowModal()
        if ret == wx.ID_OK:
            minClumpSizePostFilt = dlg.GetClumpMinSize() + 2 # since 2 events are removed for edges
            recipe.add_modules_and_execute([FindClumps(recipe, inputName=pipeline.selectedDataSourceKey,
                                                       outputName='with_clumps',
                                                       timeWindow=0,
                                                       clumpRadiusVariable=dlg.GetClumpRadiusVariable(),
                                                       clumpRadiusScale=dlg.GetClumpRadiusMultiplier()),
                                            FilterTable(recipe, inputName='with_clumps',
                                                        outputName='simpler_filtered',
                                                        filters={'clumpEdge' : [-0.1,0.1],
                                                                 'clumpSize' : [minClumpSizePostFilt-0.5,1e5]})])
            self.visFr.pipeline.selectDataSource('simpler_filtered')

        dlg.Destroy()

    def OnGenerateZSIMPLER(self, event):
        from PYMEcs.recipes.simpler import SIMPLERzgenerator
        recipe = self.visFr.pipeline.recipe
        pipeline = self.visFr.pipeline
        
        szg = SIMPLERzgenerator(recipe,inputName='simpler_filtered',outputName='with_simpler_z')
        if szg.configure_traits(kind='modal'):
           recipe.add_modules_and_execute([szg])
           self.visFr.pipeline.selectDataSource(szg.outputName)

    def OnN0FromIinterpolationField(self, event):
        from PYMEcs.recipes.simpler import N0FromInterpolationMap
        recipe = self.visFr.pipeline.recipe
        pipeline = self.visFr.pipeline
        
        n0i = N0FromInterpolationMap(recipe,inputName=pipeline.selectedDataSourceKey,outputName='with_n0')
        if n0i.configure_traits(kind='modal'):
           recipe.add_modules_and_execute([n0i])
           self.visFr.pipeline.selectDataSource(n0i.outputName)     

    def OnModesN0SIMPLER(self, event):
        raise RuntimeError('not implemented yet!')

    def OnInterpolateN0Field(self, event):
        from PYME.Analysis.points import twoColour
        from PYME.DSView import View3D
        from scipy.interpolate import Rbf, SmoothBivariateSpline
        
        if not self.intMeth.configure_traits(kind='modal'):
            return # do nothing
        
        p = self.visFr.pipeline

        idu,idx = np.unique(p['dbscanClumpID'],return_index=True)
        ccx = p['clusterCentroid_x'][idx]
        ccy = p['clusterCentroid_y'][idx]
        cm = p['clusterMode'][idx]
        cme = p['clusterModeError'][idx]

        wonky = twoColour.findWonkyVectors(ccx, ccy, cm, cm, tol=2.0*cme.mean())
        good = wonky == 0
        cmn2 = 0.01*(cm-cm.mean())
        cme2 = 0.01*cme
        
        pixelsize = 100.0
        x0=ccx.min()
        x1=ccx.max()
        y0=ccy.min()
        y1=ccy.max()
        nx = int((x1-x0)/pixelsize)
        xn = x0+pixelsize*np.arange(nx+1)
        ny = int((y1-y0)/pixelsize)
        yn = y0+pixelsize*np.arange(ny+1)
        xv,yv = np.meshgrid(xn,yn,indexing='ij')

        if self.intMeth.InterpolationMethod == 'Splines':
            spm = SmoothBivariateSpline(ccx[good], ccy[good], cm[good], 1./cme[good],bbox=[x0,x1,y0,y1])
            #n0map = 100.0*spms(xv,yv,grid=False)+cm.mean()
            n0map = spm(xv,yv,grid=False)
        elif self.intMeth.InterpolationMethod == 'Radial Basis functions':
            spm = Rbf(ccx[good], ccy[good], cm[good], epsilon=1)
            n0map = spm(xv,yv)
        elif self.intMeth.InterpolationMethod == 'Constraint Model':
            spm = twoColour.lin3Model(ccx[good],ccy[good],cm[good],cme[good]**2)
            # n0map = 100.0*spm(xv,yv)+cm.mean()
            n0map = spm(xv,yv)

        from PYME.IO.MetaDataHandler import origin_nm
        origin = origin_nm(p.mdh)
        origin_img = [origin[0]+xn[0], origin[1]+yn[0], origin[2]]
 
        mdh=voxelsize_mdh(pixelsize)
        mdh['ImageBounds.x0'] = xn[0]
        mdh['ImageBounds.y0'] = yn[0]
        mdh['ImageBounds.x1'] = xn[-1] + pixelsize
        mdh['ImageBounds.y1'] = yn[-1] + pixelsize
        mdh['ImageBounds.z0'] = 0
        mdh['ImageBounds.z1'] = 0
        mdh['Origin.x'] = origin_img[0]
        mdh['Origin.y'] = origin_img[1]
        mdh['Origin.z'] = 0
            

        View3D(n0map, titleStub='Interpolated N0 - Method %s' % (self.intMeth.InterpolationMethod),
                 glCanvas=self.visFr.glCanvas, parent=self.visFr, mdh=mdh)
        
        from six.moves import cPickle
        bbox = [x0,x1,y0,y1]
        fdialog = wx.FileDialog(None, 'Save N0 map as ...',
                                wildcard='N0 map file (*.n0m)|*.n0m',
                                style=wx.FD_SAVE)
        succ = fdialog.ShowModal()
        if (succ == wx.ID_OK):
            fpath = fdialog.GetPath()
            #save as a pickle containing the data and voxelsize

            fid = open(fpath, 'wb')
            cPickle.dump((spm,bbox,origin_img), fid, 2)
            fid.close()

    def OnShowN0Map(self, event):
        from six.moves import cPickle
        fdialog = wx.FileDialog(None, 'Please select N0 map to use ...',
                    wildcard='N0 map file (*.n0m)|*.n0m', style=wx.FD_OPEN)
        succ = fdialog.ShowModal()
        if (succ == wx.ID_OK):
            nmFilename = fdialog.GetPath()
            with open(nmFilename, 'rb') as fid:
                n0m,bb,origin = cPickle.load(fid)
            pixelsize = 100.0
            x0,x1,y0,y1 = bb
            nx = int((x1-x0)/pixelsize)
            xn = x0+pixelsize*np.arange(nx+1)
            ny = int((y1-y0)/pixelsize)
            yn = y0+pixelsize*np.arange(ny+1)
            xv,yv = np.meshgrid(xn,yn,indexing='ij')

            mdh=voxelsize_mdh(pixelsize)
            mdh['ImageBounds.x0'] = xn[0]
            mdh['ImageBounds.y0'] = yn[0]
            mdh['ImageBounds.x1'] = xn[-1] + pixelsize
            mdh['ImageBounds.y1'] = yn[-1] + pixelsize
            mdh['ImageBounds.z0'] = 0
            mdh['ImageBounds.z1'] = 0
            mdh['Origin.x'] = origin[0]
            mdh['Origin.y'] = origin[1]
            mdh['Origin.z'] = 0
            
            from PYME.DSView import View3D
            try:
                n0map = n0m(xv,yv)
            except ValueError:
                n0map = n0m(xv,yv,grid=False)

            View3D(n0map, titleStub='N0 Map', mdh=mdh,
                   glCanvas=self.visFr.glCanvas, parent=self.visFr)

    def OnN0FromImage(self, event):
        """

        Function to propagate N0 values from an image (or stack of images) to localizations within the pipeline.
        Localizations in the same area (or volume) of image. ImageStack values will be given the N0 value from
        corresponding values. The ImageStack containing the values is selected through the GUI.

        Parameters
        ----------
        event: GUI event

        Returns
        -------
        Nothing, but adds N0 column to the pipeline
            N0: N0 values from image, mapped to each localization within that label
        """
        
        from PYME.IO import image
        from PYMEcs.recipes.simpler import N0FromImage
        #from PYME.Analysis.points import objectMeasure

        visFr = self.visFr
        pipeline = visFr.pipeline
        recipe = pipeline.recipe

        dlg = wx.SingleChoiceDialog(
                None, 'choose the image which contains N0 values', 'Set N0',
                list(image.openImages.keys()),
                wx.CHOICEDLG_STYLE
                )

        if dlg.ShowModal() == wx.ID_OK:
            img_name = dlg.GetStringSelection()
            img = image.openImages[img_name]
            
            recipe.namespace[img_name] = img
            
            output_name = recipe.new_output_name('with_N0')
            mod = N0FromImage(recipe, inputName=pipeline.selectedDataSourceKey,inputImage=img_name,
                                  outputName=output_name,label_key_name='N0')
            if mod.configure_traits(kind='modal'):
                recipe.add_modules_and_execute([mod,])
                pipeline.selectDataSource(mod.outputName)

        dlg.Destroy()

def Plug(visFr):
    return SIMPLER(visFr)
