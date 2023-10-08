from __future__ import print_function    # (at top of module)

import numpy as np
from scipy import ndimage
from collections import OrderedDict
import warnings

import logging
logger = logging.getLogger(__name__)

# this is just a slightly cleaned up version of the core code of
# Kenny's fiducial tracker from PYMEnf

# currently only does x and y dims

# in this version we enforce that the returned fiducial track
# starts out of zero

def foffset(t,ft,navg=50):
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
    # must be int(size), otherwise gets float type error!!!!!
    'Median' : makeFilter(lambda input, size: ndimage.median_filter(input,int(size)))
    } 

def extractTrajectoriesClump(ds, clumpRadiusVar = 'error_x', clumpRadiusMultiplier=5.0, 
                                  timeWindow=25, clumpMinSize=50, align=True):
                                      
    import PYME.Analysis.points.DeClump.deClump as deClump
    #track beads through frames
    if clumpRadiusVar == '1.0':
        delta_x = 0*ds['x'] + clumpRadiusMultiplier
    else:
        delta_x = clumpRadiusMultiplier*ds[clumpRadiusVar]
        
    t = ds['t'].astype('i')
    x = ds['x'].astype('f4')
    y = ds['y'].astype('f4')
    z = ds['z'].astype('f4')
    delta_x = delta_x.astype('f4')
    
    I = np.argsort(t)

    clumpIndex = np.zeros(len(x), dtype='i')
    isFiducial = np.zeros(len(x), dtype='i')
    clumpIndex[I] = deClump.findClumpsN(t[I], x[I], y[I], delta_x[I], timeWindow)
    
    tMax = t.max()
    
    clumpIndices = list(set(clumpIndex))

    x_f = []
    y_f = []
    z_f = []
    clump_sizes = []
    
    t_f = np.arange(0, tMax + 1, dtype='i')
    
    #loop over all our clumps and extract trajectories
    for ci in clumpIndices:
        if ci > 0:
            clump_mask = (clumpIndex == ci)
            x_i = x[clump_mask]
            clump_size = len(x_i)
            
            if clump_size > clumpMinSize:
                y_i = y[clump_mask]
                z_i = z[clump_mask]
                t_i = t[clump_mask].astype('i')
                isFiducial[clump_mask] = 1 # mark the event mask that this is a fiducial
                
                x_i_f = np.NaN*np.ones_like(t_f)
                if align:
                    x_i_f[t_i]= x_i - x_i.mean()
                else:
                    x_i_f[t_i]= x_i
                    
                y_i_f = np.NaN*np.ones_like(t_f)
                if align:
                    y_i_f[t_i]= y_i - y_i.mean()
                else:
                    y_i_f[t_i]= y_i
                    
                z_i_f = np.NaN*np.ones_like(t_f)
                if align:
                    z_i_f[t_i]= z_i - z_i.mean()
                else:
                    z_i_f[t_i]= z_i
                
                #clumps.append((x_i_f, y_i_f))
                x_f.append(x_i_f)
                y_f.append(y_i_f)
                z_f.append(z_i_f)
                clump_sizes.append(len(x_i))
    
    #re-order to start with the largest clump
    clumpOrder = np.argsort(clump_sizes)[::-1]
    x_f = np.array(x_f)[clumpOrder,:]
    y_f = np.array(y_f)[clumpOrder,:]
    z_f = np.array(z_f)[clumpOrder,:]

    return (t_f, x_f, y_f, z_f, isFiducial)

def AverageTrack(ds, tracks, filter='Gaussian', filterScale=10.0, align=True):

    t_f, x_f, y_f, z_f = tracks
    t = ds['t'].astype('i')
    
    # this function does not appear to be used anywhere
    # def _mf(p, meas):
    #     '''calculate the offset between trajectories'''
    #     m_adj = meas + np.hstack([[0], p])[:,None]
        
    #     return np.nansum(np.nanvar(m_adj, axis=0))
        
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

        # nanmean can generate warnings with all nan values
        # we catch them in this block only
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message='Mean of empty slice')
            mm = np.nanmean(meas, 0)
        
        print('Finished:', n_iters, dm)
        return mm

    def _simpleav(meas):
        # as above
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message='Mean of empty slice')
            mm = np.nanmean(meas, 0)
        return mm
        
    if align:
        x_corr = _align(x_f)
        y_corr = _align(y_f)
        z_corr = _align(z_f)
    else:
        x_corr = _simpleav(x_f)
        y_corr = _simpleav(y_f)
        z_corr = _simpleav(z_f)
        
    filtered_corr_woffs = FILTER_FUNCS[filter](t_f, {'x' : x_corr, 'y':y_corr, 'z':z_corr}, filterScale)
    
    dims = filtered_corr_woffs.keys()
    filtered_corr = {}
    for dim in dims:
        filtered_corr[dim] = filtered_corr_woffs[dim] - foffset(t_f,filtered_corr_woffs[dim])

    fcorr_columns = {}
    for dim in dims:
        fcorr_columns[dim] = np.interp(t, t_f, filtered_corr[dim]) 
    
    return fcorr_columns
