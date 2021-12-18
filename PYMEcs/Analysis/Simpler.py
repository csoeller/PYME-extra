import numpy as np
from scipy.stats import iqr

import sys
if sys.version_info > (3,):
    xrange = range

# halfmode is apparently a robust measure of the mode:
# 1.Bickel, D. R. & Frühwirth, R. On a fast, robust estimator of the mode: Comparisons to other robust estimators with applications. Computational Statistics & Data Analysis 50, 3500–3530 (2006).

def halfmode(inputData, axis=None, dtype=None):
    """
    Robust estimator of the mode of a data set using the half-sample mode.
    
    .. versionadded: 1.0.3
    """
    
    if axis is not None:
        fnc = lambda x: mode(x, dtype=dtype)
        dataMode = np.apply_along_axis(fnc, axis, inputData)
    else:
        # Create the function that we can use for the half-sample mode
        def _hsm(data):
            if data.size == 1:
                return data[0]
            elif data.size == 2:
                return data.mean()
            elif data.size == 3:
                i1 = data[1] - data[0]
                i2 = data[2] - data[1]
                if i1 < i2:
                    return data[:2].mean()
                elif i2 > i1:
                    return data[1:].mean()
                else:
                    return data[1]
            else:
                wMin = data[-1] - data[0]
                N = int(data.size / 2) + data.size % 2 
                for i in xrange(0, int(N)):
                    w = data[i+N-1] - data[i] 
                    if w < wMin:
                        wMin = w
                        j = i
                return _hsm(data[j:j+N])
                
        data = inputData.ravel()
        if type(data).__name__ == "MaskedArray":
            data = data.compressed()
        if dtype is not None:
            data = data.astype(dtype)
            
        # The data need to be sorted for this to work
        data = np.sort(data)
        
        # Find the mode
        dataMode = _hsm(data)
        
    return dataMode


def clusterModes(x,y,eid,prop):
    clusterMode = np.zeros_like(x)
    clusterErr = np.zeros_like(x)
    clusterCentroid_x = np.zeros_like(x)
    clusterCentroid_y = np.zeros_like(x)
    
    uids = np.unique(eid)
    for j,i in enumerate(uids):
        if not i == 0: # cluster id 0 means not a cluster
            ind = eid == i
            xi = x[ind]
            yi = y[ind]
            pi = prop[ind]

            clusterMode[ind] = halfmode(pi)
            clusterErr[ind] = iqr(pi)
            clusterCentroid_x[ind] = np.sum(xi*pi)/pi.sum()
            clusterCentroid_y[ind] = np.sum(yi*pi)/pi.sum()

    return clusterMode, clusterErr, clusterCentroid_x, clusterCentroid_y
