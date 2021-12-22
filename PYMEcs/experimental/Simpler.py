from traitsui.api import View, Item, Group
from traitsui.menu import OKButton, CancelButton, OKCancelButtons
from traits.api import HasTraits, Str, Int, CStr, List, Enum, Float, Bool


#from PYME.DSView.dsviewer import View3D
#from PYME.IO.MetaDataHandler import NestedClassMDHandler

import numpy as np
import wx

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


# things still to add
# - "Filter events for SIMPLER": choose clump radius and min clumpsize interactively
# - new function to calculate N0 from data background
# - new function to insert SIMPLERzgenerator
# - add error estimates for SIMPLER (combination of formula and/or stddev) -> do in recipe module
class SIMPLER():
    def __init__(self, visFr):
        self.visFr = visFr
        self.intMeth = InterpolChoice()
        
        visFr.AddMenuItem('Experimental>SIMPLER', "Filter events for SIMPLER", self.OnFilterSIMPLER)
        visFr.AddMenuItem('Experimental>SIMPLER', "Cluster modes for SIMPLER N0 distribution", self.OnModesN0SIMPLER)
        visFr.AddMenuItem('Experimental>SIMPLER', "Interpolate N0 Map", self.OnCalculateN0Field)
        visFr.AddMenuItem('Experimental>SIMPLER', "Load N0 Map", self.OnLoadN0Map)
        
    def OnFilterSIMPLER(self, event):
        from PYME.recipes.tablefilters import FilterTable
        from PYME.recipes.tracking import FindClumps
        recipe = self.visFr.pipeline.recipe
        pipeline = self.visFr.pipeline
        recipe.add_modules_and_execute([FindClumps(recipe, inputName=pipeline.selectedDataSourceKey,
                                                   outputName='with_clumps',
                                                   timeWindow=0,
                                                   clumpRadiusVariable='error_x',
                                                   clumpRadiusScale=2.0),
                                        FilterTable(recipe, inputName='with_clumps',
                                                    outputName='simpler_filtered',
                                                    filters={'clumpEdge' : [-0.1,0.1]})])
        self.visFr.pipeline.selectDataSource('simpler_filtered')

    def OnModesN0SIMPLER(self, event):
        raise RuntimeError('not implemented yet!')

    def OnCalculateN0Field(self, event):
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

        View3D(n0map, titleStub='Interpolated N0 - Method %s' % (self.intMeth.InterpolationMethod),
                 glCanvas=self.visFr.glCanvas, parent=self.visFr, mdh=voxelsize_mdh(pixelsize))
        
        # def spmsrescale(x,y):
        #     return 100.0*spms(x,y,grid=False)+cm.mean()
    
        # def spmrescale(x,y):
        #     return 100.0*spm(x,y)+cm.mean()

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
            cPickle.dump((spm,bbox), fid, 2)
            fid.close()

    def OnLoadN0Map(self, event):
        from six.moves import cPickle
        fdialog = wx.FileDialog(None, 'Please select N0 map to use ...',
                    wildcard='N0 map file (*.n0m)|*.n0m', style=wx.FD_OPEN)
        succ = fdialog.ShowModal()
        if (succ == wx.ID_OK):
            nmFilename = fdialog.GetPath()
            with open(nmFilename, 'rb') as fid:
                n0m,bb = cPickle.load(fid)
            pixelsize = 100.0
            x0,x1,y0,y1 = bb
            nx = int((x1-x0)/pixelsize)
            xn = x0+pixelsize*np.arange(nx+1)
            ny = int((y1-y0)/pixelsize)
            yn = y0+pixelsize*np.arange(ny+1)
            xv,yv = np.meshgrid(xn,yn,indexing='ij')

            from PYME.DSView import View3D
            try:
                n0map = n0m(xv,yv)
            except ValueError:
                n0map = n0m(xv,yv,grid=False)

            View3D(n0map, titleStub='N0 Map', mdh=voxelsize_mdh(pixelsize),
                   glCanvas=self.visFr.glCanvas, parent=self.visFr)
       
def Plug(visFr):
    return SIMPLER(visFr)
