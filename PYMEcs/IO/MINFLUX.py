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
def minflux_npy2pyme(fname,return_original_array=False,make_clump_index=True,with_cfr_std=False):
    data = np.load(fname)
    
    if minflux_npy_detect_3D(data):
        is_3D = True
        iterno_loc = 9 # we pick up the most precise localisation from this iteration, also fbg
        iterno_other = 9 # we pick up dcr, efo from this iteration
        iterno_cfr = 6
    else:
        is_3D = False
        iterno_loc = 4
        iterno_other = 4
        iterno_cfr = 4

    posnm = 1e9*data['itr']['loc'][:,iterno_loc] # we keep all distances in units of nm
    if 'lnc' in data['itr'].dtype.fields:
        posnm_nc = 1e9*data['itr']['lnc'][:,iterno_loc]
        has_lnc = True
    else:
        has_lnc = False

    pymedct = {}
        
    # this way we ensure that the valid vs invalid portions of the same trace get separate ids
    #  it becomes important for calculation std_devs for traces
    rawids = 2*data['tid'] + data['vld']

    if make_clump_index:
        # we replace the non-sequential trace ids from MINFLUX data with a set of sequential ids
        # this works better for clumpIndex assumptions in the end
        uids,revids = np.unique(rawids,return_inverse=True)
        ids = np.arange(1,uids.size+1,dtype='int32')[revids]
        counts = get_stddev_property(ids,posnm[:,0],statistic='count')
        pymedct.update({'clumpIndex': ids,
                        'clumpSize' : counts,
                        })
    else:
        ids = rawids

    stdx = get_stddev_property(ids,posnm[:,0])
    # we expect this to only happen when clumpSize == 1, because then std dev comes back as 0
    stdx[stdx < 1e-3] = 100.0 # if error estimate is too small, replace with 100 as "large" flag
    stdy = get_stddev_property(ids,posnm[:,1])
    stdy[stdy < 1e-3] = 100.0
    if is_3D:
        stdz = get_stddev_property(ids,posnm[:,2])
        stdz[stdz < 1e-3] = 100.0
        pymedct.update({'z':posnm[:,2], 'error_z' : stdz})

    if with_cfr_std: # we also compute on request a cfr std dev across a trace ID (=clump in PYME)
        pymedct.update({'cfr_std':get_stddev_property(ids,data['itr']['cfr'][:,iterno_cfr])})
        
    pymedct.update({'x' : posnm[:,0],
                    'y': posnm[:,1],
                    # for t we use time to ms precision (without rounding); this is a reasonably close
                    # correspondence to frame numbers as time coordinates in SMLM data
                    't': (1e3*data['tim']).astype('i'),
                    'cfr':data['itr']['cfr'][:,iterno_cfr],
                    'efo':data['itr']['efo'][:,iterno_other],
                    'dcr':data['itr']['dcr'][:,iterno_other],
                    'error_x' : stdx,
                    'error_y' : stdy,
                    'fbg': data['itr']['fbg'][:,iterno_other],
                    # we assume for now the offset counts can be used to sum up
                    # and get the total photons harvested
                    # check with abberior
                    'nPhotons' : data['itr']['eco'].sum(axis=1),
                    'tim': data['tim'], # we also keep the original float time index, units are [s]                  
                    })

    if has_lnc:
        pymedct.update({'x_nc' : posnm_nc[:,0],
                        'y_nc' : posnm_nc[:,1]})
        if is_3D:
            pymedct.update({'z_nc' : posnm_nc[:,2]})

    # copy a few entries verbatim
    for key in ['tid','act','vld']:
        pymedct[key] = data[key].astype('i') # these are either integer types or should be converted to integer

    pymepd = pd.DataFrame.from_dict(pymedct)
    if return_original_array:
        return (pymepd,data)
    else:
        return pymepd


# convenience function                           
def save_minflux_as_csv(pd_data, fname):
    pd_data.to_csv(fname,index=False)

# we are monkeypatching pipeline and VisGUIFrame methods to sneak MINFLUX npy IO in;
# in future we will ask for a way to get this considered by David B for a proper hook
# in the IO code
def monkeypatch_npy_io(visFr):
    import types
    import logging
    import os
    from PYME.IO import MetaDataHandler

    logger = logging.getLogger(__name__)
    logger.info("MINFLUX monkeypatching IO")
    def _populate_open_args_npy(self, filename):
        # this is currently just the minmal functionality for .npy,
        # we should really check a few things before going any further
        # .mat and CSV files give examples...
        if os.path.splitext(filename)[1] == '.npy':
            return {}
        else:
            return self._populate_open_args_original(filename)

    visFr._populate_open_args_original = visFr._populate_open_args
    visFr._populate_open_args = types.MethodType(_populate_open_args_npy,visFr)

    def _ds_from_file_npy(self, filename, **kwargs):
        if os.path.splitext(filename)[1] == '.npy': # MINFLUX NPY file
            logger.info('.npy file, trying to load as MINFLUX npy ...')
            from PYMEcs.IO.tabular import MinfluxNpySource
            ds = MinfluxNpySource(filename)
            ds.mdh = MetaDataHandler.NestedClassMDHandler()
            return ds
        else:
            return self._ds_from_file_original(filename, **kwargs)

    visFr.pipeline._ds_from_file_original = visFr.pipeline._ds_from_file
    visFr.pipeline._ds_from_file = types.MethodType(_ds_from_file_npy,visFr.pipeline)
    logger.info("MINFLUX monkeypatching IO completed")

# below we make a class Pipeline that inherits from PYME.LMVis.pipeline.Pipeline
# and changes the relevant method in the subclass
#
# in your own code (e.g. Python notebook) use as
#
#    from PYMEcs.IO.MINFLUX import Pipeline # use this instead of PYME.LMVis.pipeline
#    data = Pipeline('my_minflux_file.npy')
#
from PYME.LMVis import pipeline
from PYME.IO import MetaDataHandler
import os
import logging
class Pipeline(pipeline.Pipeline):
    
    def _ds_from_file(self, filename, **kwargs):
        if os.path.splitext(filename)[1] == '.npy': # MINFLUX NPY file
            logging.getLogger(__name__).info('.npy file, trying to load as MINFLUX npy ...')
            from PYMEcs.IO.tabular import MinfluxNpySource
            ds = MinfluxNpySource(filename)
            ds.mdh = MetaDataHandler.NestedClassMDHandler()
            return ds
        else:
            return super()._ds_from_file(filename, **kwargs)
