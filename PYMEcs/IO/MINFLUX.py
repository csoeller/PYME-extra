# here we provide a few routines to translate MINFLUX provided data structures
# that are read from NPY files

# we translate the NPY based datastructure into a PYME compatible data structure that we hold
# in a pandas dataframe

# currently, for reading into PYME we provide the functionality to write out as a CSV
# from the pandas dataframe; PYME can parse the generated CSV pretty well upon reading

from scipy.stats import binned_statistic
import pandas as pd
import numpy as np

def get_stddev_property(ids, prop, statistic='std'):
    maxid = int(ids.max())
    edges = -0.5+np.arange(maxid+2)
    idrange = (0,maxid)
        
    propstd, bin_edge, binno = binned_statistic(ids, prop, statistic=statistic,
                                                bins=edges, range=idrange)
    propstd[np.isnan(propstd)] = 1000.0 # (mark as huge error)
    std_events = propstd[ids]
    return std_events


# here we check for size either 5 (2D) or 10 (3D); any other size raises an error
def minflux_npy_detect_3D(data):
    if data['itr'].shape[1] == 10:
        return True # 3D
    elif data['itr'].shape[1] == 5:
        return False # 2D
    else:
        raise RuntimeError('unknown size of itr array, neither 5 (2D) nor 10 (3D), is actually: %d' %
                            (data['itr'].shape[1]))

# this one should be able to deal both with 2d and 3D
def minflux_npy2pyme(fname,return_original_array=False):
    data = np.load(fname)
    
    if minflux_npy_detect_3D(data):
        is_3D = True
        iterno_loc = 9 # we pick up the most precise localisation from this iteration, also fbg
        iterno_other = 6 # we pick up cfr, efo from this iteration
    else:
        is_3D = False
        iterno_loc = 4
        iterno_other = 4

    posnm = 1e9*data['itr']['loc'][:,iterno_loc] # we keep all distances in units of nm
    if 'lnc' in data['itr'].dtype.fields:
        posnm_raw = 1e9*data['itr']['lnc'][:,iterno_loc]
        beamline_monitoring = True
    else:
        beamline_monitoring = False
    rawids = data['tid']
    # we replace the non-sequential trace ids from MINFLUX data with a set of sequential ids
    # this works better for clumpIndex assumptions in the end
    uids,revids = np.unique(rawids,return_inverse=True)
    newids = np.arange(1,uids.size+1,dtype='int32')[revids]
    stdx = get_stddev_property(newids,posnm[:,0])
    counts = get_stddev_property(newids,posnm[:,0],statistic='count')
    stdy = get_stddev_property(newids,posnm[:,1])
    pymedct =  {'x' : posnm[:,0],
                'y': posnm[:,1],
                # for t we use time to ms precision (without rounding); this is a reasonably close
                # correspondence to frame numbers as time coordinates in SMLM data
                't': (1e3*data['tim']).astype('i'),
                'clumpIndex': newids,
                'cfr':data['itr']['cfr'][:,iterno_other],
                'efo':data['itr']['efo'][:,iterno_other],
                'error_x' : stdx,
                'error_y' : stdy,
                'clumpSize' : counts,
                'fbg': data['itr']['fbg'][:,iterno_loc],
                # we assume for now the offset counts can be used to sum up
                # and get the total photons harvested
                # check with abberior
                'nPhotons' : data['itr']['eco'].sum(axis=1),
               }
    if is_3D:
        stdz = get_stddev_property(newids,posnm[:,2])
        pymedct.update({'z':posnm[:,2], 'error_z' : stdz})
    if beamline_monitoring:
        pymedct.update({'x_raw' : posnm_raw[:,0],
                        'y_raw' : posnm_raw[:,1]})
        if is_3D:
            pymedct.update({'z_raw' : posnm_raw[:,2]})
    pymepd = pd.DataFrame.from_dict(pymedct)
    if return_original_array:
        return (pymepd,data)
    else:
        return pymepd
                          
                           
def save_minflux_as_csv(pd_data, fname):
    pd_data.to_csv(fname,index=False)
