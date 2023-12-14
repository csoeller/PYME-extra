# -*- coding: utf-8 -*-
"""
Created on Tue Jun 30 10:24:45 2015

@author: Kenny
"""

# this is Kenny's code with minimal changes to run as standonlone plugin
# using the new PYME.config system

# mostly used as a comparison to the newer recipe based implementation
# will be dropped once the recipe based version is complete

import wx
import numpy as np
from scipy import ndimage
from collections import OrderedDict

import logging
logger = logging.getLogger(__name__)

def foffset(t,ft,navg=100):
    tu,idx = np.unique(t.astype('int'), return_index=True)
    fu = ft[idx]
    offs = fu[0:min(navg,fu.shape[0])].mean()
    return offs

def makeFilter(filtFunc):
    '''wrapper function for different filters'''
    def ffcn(t, data, scale):
        out = {}
        for k, v in data.items():
            r_v = v[~np.isnan(v)]
            r_t = t[~np.isnan(v)]
            out[k] = filtFunc(np.interp(t, r_t, r_v), scale)
        return out
    return ffcn
    
FILTER_FUNCS = {
    'Gaussian' : makeFilter(ndimage.gaussian_filter),
    'Uniform' : makeFilter(ndimage.uniform_filter),
    'Median' : makeFilter(ndimage.median_filter)
    } 

def _extractAverageTrajectory(pipeline, clumpRadiusVar = 'error_x', clumpRadiusMultiplier=5.0, 
                                  timeWindow=25, filter='Gaussian', filterScale=10.0):
                                      
    #import PYME.Analysis.trackUtils as trackUtils
    import PYME.Analysis.points.DeClump.deClump as deClump
    from scipy.optimize import fmin
    #track beads through frames
    if clumpRadiusVar == '1.0':
        delta_x = 0*pipeline['x'] + clumpRadiusMultiplier
    else:
        delta_x = clumpRadiusMultiplier*pipeline[clumpRadiusVar]
        
    t = pipeline['t'].astype('i')
    x = pipeline['x'].astype('f4')
    y = pipeline['y'].astype('f4')
    delta_x = delta_x.astype('f4')
    
    I = np.argsort(t)

    clumpIndex = np.zeros(len(x), dtype='i')
    clumpIndex[I] = deClump.findClumpsN(t[I], x[I], y[I], delta_x[I], timeWindow)
    #trackUtils.findTracks(pipeline, clumpRadiusVar,clumpRadiusMultiplier, timeWindow)
    
    #longTracks = pipeline['clumpSize'] > 50            
    
    #x = x[longTracks].copy()
    #y = pipeline['y_raw'][longTracks].copy()
    #t = pipeline['t'][longTracks].copy() #.astype('i')
    #clumpIndex = pipeline['clumpIndex'][longTracks].copy()
    
    tMax = t.max()
    
    clumpIndices = list(set(clumpIndex))

    x_f = []
    y_f = []
    clump_sizes = []
    
    t_f = np.arange(0, tMax + 1, dtype='i')
    
    #loop over all our clumps and extract trajectories
    for ci in clumpIndices:
        if ci > 0:
            clump_mask = (clumpIndex == ci)
            x_i = x[clump_mask]
            clump_size = len(x_i)
            
            if clump_size > 50:
                y_i = y[clump_mask]
                t_i = t[clump_mask].astype('i')
                
                x_i_f = np.NaN*np.ones_like(t_f)
                x_i_f[t_i]= x_i - x_i.mean()
                
                y_i_f = np.NaN*np.ones_like(t_f)
                y_i_f[t_i]= y_i - y_i.mean()
                
                #clumps.append((x_i_f, y_i_f))
                x_f.append(x_i_f)
                y_f.append(y_i_f)
                clump_sizes.append(len(x_i))
    
    #re-order to start with the largest clump
    clumpOrder = np.argsort(clump_sizes)[::-1]
    x_f = np.array(x_f)[clumpOrder,:]
    y_f = np.array(y_f)[clumpOrder,:]
    
    def _mf(p, meas):
        '''calculate the offset between trajectories'''
        m_adj = meas + np.hstack([[0], p])[:,None]
        
        return np.nansum(np.nanvar(m_adj, axis=0))
        
    #print x_f.shape, np.hstack([[0], np.random.randn(x_f.shape[0]-1)]).shape
        
    def _align(meas, tol=.1):
        n_iters = 0
        
        dm_old = 5e12
        dm = 4e12
        
        mm = np.nanmean(meas, 0)
        
        while ((dm_old - dm) > tol) and (n_iters < 50):  
            dm_old = dm
            mm = np.nanmean(meas, 0)        
            d = np.nanmean(meas - mm, 1)
            dm = sum(d**2)
            meas = meas - d[:,None]
            n_iters +=1
            print(n_iters, dm)
         
        mm = np.nanmean(meas, 0)
        print('Finished:', n_iters, dm)
        return mm
        
    x_corr = _align(x_f)
    y_corr = _align(y_f)
     
    filtered_corr = FILTER_FUNCS[filter](t_f, {'x' : x_corr, 'y':y_corr}, filterScale)
    
    return t_f, filtered_corr

class FiducialAnalyser:
    def __init__(self, visFr):
        self.visFr = visFr

        visFr.AddMenuItem('Experimental>Deprecated>FiducialsOld', "Estimate drift from Fiducials", self.FindBeadsAndTrack)
        visFr.AddMenuItem('Experimental>Deprecated>FiducialsOld', 'Apply fiducial correction', self.OnApplyFiducial,
                          helpText='Apply fiducial to x, y, z')
        visFr.AddMenuItem('Experimental>Deprecated>FiducialsOld', 'Revert fiducial correction', self.OnRevertFiducial,
                          helpText='Revert fiducial correction to x, y, z')

    def FindBeadsAndTrack(self, event):      
        dlg = ExtractTrajectoriesDialog(self.visFr)
        succ = dlg.ShowModal()
        if succ == wx.ID_OK:
            pipeline = self.visFr.pipeline
            
            beadTraj = _extractAverageTrajectory(pipeline, dlg.GetClumpRadiusVariable(),
                                                 dlg.GetClumpRadiusMultiplier(), dlg.GetClumpTimeWindow(),
                                                 filter=dlg.GetFilterMethod(),filterScale=dlg.GetFilterScale())

            self.ApplyCorrections([beadTraj,])
                
    def ApplyCorrections(self, interpInfos):
        """Averages drift from multiple files.
        Adds the drift as fiducial_? to pipeline.inputMapping.
        Overwrites the drift panel with ? + fiducial_?
        """
        pipeline = self.visFr.pipeline
        count = len(interpInfos)
        #master dim key list assuming first file is correct
        dims = interpInfos[0][1].keys()
        
        #loop over data sources and add drift info to each
        logger.debug('about to update datasources...')
        curds = pipeline.selectedDataSourceKey
        for dsname in pipeline.dataSources:
            logger.debug('updating datasource %s' % dsname)
            pipeline.selectDataSource(dsname)
            tCache = pipeline['t']
            fudical_multi = np.zeros((count, len(dims), len(tCache))) #num of files, num of dims, num of time points
            for i, j in enumerate(interpInfos):
                realTime, filtered = j
                for k, l in enumerate(dims):
                    print(l)
                    fudical_multi[i, k, :] = np.interp(tCache, realTime, filtered[l])
                
            fuducial_mean = fudical_multi.mean(0)
            for i, dim in enumerate(dims):
                fiducial = 'fiducial_%s' % dim
                pipeline.addColumn(fiducial, fuducial_mean[i,:]-foffset(tCache,fuducial_mean[i,:]))
                logger.debug('setting attribute %s' % fiducial)
            pipeline.Rebuild()
        pipeline.selectDataSource(curds)

    def OnApplyFiducial(self, event=None):
        pipeline = self.visFr.pipeline
        for dim in ('x','y','z'):
            fiducial = 'fiducial_%s' % dim
            if fiducial in pipeline.keys():
                pipeline.mapping.setMapping(dim,'%s - %s' % (dim, fiducial))
        pipeline.Rebuild()

    def OnRevertFiducial(self, event=None):
        pipeline = self.visFr.pipeline
        changed = False
        for dim in ('x','y','z'):
            if dim in pipeline.mapping.mappings:
                pipeline.mapping.mappings.pop(dim)
                changed = True
        if changed:
            pipeline.Rebuild()

    def Calculate(self, func, arg):
        """Returns tuple of 3 objects.
        realTime: arange from 0 to max time
        interp: interpolated x, y, z
        filtered: 'interp' filtered using the given filter function
        """
        pipeline = self.visFr.pipeline
        
        realTime = np.arange(0, pipeline['t'].max())
        
        interp = OrderedDict()
        filtered = OrderedDict()
        for dim in ['x', 'y', 'z']:
#            print(dim)
            if dim in pipeline.keys():
                interp[dim] = np.interp(realTime, pipeline['t'], pipeline[dim])
                filtered[dim] = func(interp[dim], arg)
        
        
        return (realTime, interp, filtered)
    
        
class ExtractTrajectoriesDialog(wx.Dialog):
    def __init__(self, *args, **kwargs):
        wx.Dialog.__init__(self, *args, **kwargs)

        vsizer = wx.BoxSizer(wx.VERTICAL)

        hsizer = wx.BoxSizer(wx.HORIZONTAL)

        hsizer.Add(wx.StaticText(self, -1, 'Clump Radius: '), 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)

        self.tClumpRadMult = wx.TextCtrl(self, -1, '5.0', size=[30,-1])
        hsizer.Add(self.tClumpRadMult, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)

        hsizer.Add(wx.StaticText(self, -1, 'X'), 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)

        self.cClumpRadVar = wx.Choice(self, -1, choices=['1.0', 'error_x'])
        self.cClumpRadVar.SetSelection(1)
        hsizer.Add(self.cClumpRadVar,1, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
        
        vsizer.Add(hsizer, 0, wx.ALL, 5)
        
        hsizer = wx.BoxSizer(wx.HORIZONTAL)
        hsizer.Add(wx.StaticText(self, -1, 'Time Window: '), 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)

        self.tClumpTime = wx.TextCtrl(self, -1, '25')
        hsizer.Add(self.tClumpTime,1, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)

        vsizer.Add(hsizer, 0, wx.ALL, 5)
        
        hsizer = wx.BoxSizer(wx.HORIZONTAL)
        hsizer.Add(wx.StaticText(self, wx.ID_ANY, 'filter:'), 1, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
        self.filter = wx.Choice(self, wx.ID_ANY, choices=['Gaussian', 'Uniform', 'Median'])
        self.filter.SetSelection(0)
        #self.filter.Bind(wx.EVT_CHOICE, self.OnFilterSelected)            
        hsizer.Add(self.filter, 1, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
        
        
        vsizer.Add(hsizer)
        
        hsizer = wx.BoxSizer(wx.HORIZONTAL)
        self.argText = wx.StaticText(self, wx.ID_ANY, 'Filter scale:')
        hsizer.Add(self.argText, 1, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
        self.filtScale = wx.TextCtrl(self, wx.ID_ANY, '10')
        hsizer.Add(self.filtScale, 1, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 5)
        vsizer.Add(hsizer)

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

    def GetClumpTimeWindow(self):
        return int(self.tClumpTime.GetValue())
    
    def GetFilterMethod(self):
        return self.filter.GetStringSelection()
    
    def GetFilterScale(self):
        return float(self.filtScale.GetValue())


def Plug(visFr):
    '''Plugs this module into the gui'''
    visFr.experimentalFiducialAnalyzer = FiducialAnalyser(visFr)
