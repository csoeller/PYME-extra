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

def extractAverageTrajectory(pipeline, clumpRadiusVar = 'error_x', clumpRadiusMultiplier=5.0, 
                                  timeWindow=25, filter='Gaussian', filterScale=10.0):
                                      
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
            print n_iters, dm
         
        mm = np.nanmean(meas, 0)
        print 'Finished:', n_iters, dm
        return mm
        
    x_corr = _align(x_f)
    y_corr = _align(y_f)
     
    filtered_corr_woffs = FILTER_FUNCS[filter](t_f, {'x' : x_corr, 'y':y_corr}, filterScale)
    
    dims = filtered_corr_woffs.keys()
    filtered_corr = {}
    for dim in dims:
        filtered_corr[dim] = filtered_corr_woffs[dim] - foffset(t_f,filtered_corr_woffs[dim])

    fcorr_columns = {}
    for dim in dims:
        fcorr_columns[dim] = np.interp(t, t_f, filtered_corr[dim]) 
    
    return fcorr_columns
