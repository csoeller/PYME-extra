# here we provide a few routines to translate MINFLUX provided data structures
# that are read from NPY files

# we translate the NPY based datastructure into a PYME compatible data structure that we hold
# in a pandas dataframe

# currently, for reading into PYME we provide the functionality to write out as a CSV
# from the pandas dataframe; PYME can parse the generated CSV pretty well upon reading

from scipy.stats import binned_statistic
import pandas as pd
import numpy as np
import os

import PYME.config
# foreshortening factor estimate, see also
# Gwosch, K. C. et al. MINFLUX nanoscopy delivers 3D multicolor nanometer
# resolution in cells. Nature Methods 17, 217–224 (2020), who use 0.7.
foreshortening = PYME.config.get('MINFLUX_foreshortening',0.75)

warning_msg = ""

def get_stddev_property(ids, prop, statistic='std'):
    maxid = int(ids.max())
    edges = -0.5+np.arange(maxid+2)
    idrange = (0,maxid)
        
    propstd, bin_edge, binno = binned_statistic(ids, prop, statistic=statistic,
                                                bins=edges, range=idrange)
    propstd[np.isnan(propstd)] = 1000.0 # (mark as huge error)
    std_events = propstd[ids]
    return std_events

from PYME.warnings import warn
def npy_is_minflux_data(filename, warning=False, return_msg=False):
    data = np.load(filename)
    valid = True
    msg = None
    if data.dtype.fields is None:
        valid = False
        msg = 'no fields in NPY data, likely not a MINFLUX data set'
    else:
        for field in ['itr','tim','tid','vld']:
            if not field in data.dtype.fields:
                valid = False
                msg = 'no "%s" field in NPY data, likely not a MINFLUX data set' % field
                break

    if not valid and warning:
        if not msg is None:
                warn(msg)

    if return_msg:
        return (valid,msg)
    else:
        return valid

# here we check for size either 5 (2D) or 10 (3D); any other size raises an error
def minflux_npy_detect_3D(data):
    if data['itr'].shape[1] == 10 or data['itr'].shape[1] == 11:
        return True # 3D
    elif data['itr'].shape[1] == 5 or data['itr'].shape[1] == 6:
        return False # 2D
    else:
        raise RuntimeError('unknown size of itr array, neither 5 (2D) nor 10 (3D), is actually: %d' %
                            (data['itr'].shape[1]))


def minflux_npy_has_extra_iter(data):
    if data['itr'].shape[1] == 6 or data['itr'].shape[1] == 11:
        return True # has a spare empty starting position
    else:
        return False

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

    # NOTE CS 3/2024: latest data with MBM active seems to generate an "empty" iteration (at position 0)
    # that has NaNs or zeros in the relevant properties
    # we seem to be able to deal with this by just moving our pointers into the iteration just one position up
    # this is subject to confirmation
    if minflux_npy_has_extra_iter(data):
        has_extra_iter = True
        iterno_loc += 1
        iterno_other += 1
        iterno_cfr += 1
    else:
        has_extra_iter = False


    posnm = 1e9*data['itr']['loc'][:,iterno_loc] # we keep all distances in units of nm
    posnm[:,2] *= foreshortening
    if 'lnc' in data['itr'].dtype.fields:
        posnm_nc = 1e9*data['itr']['lnc'][:,iterno_loc]
        posnm_nc[:,2] *= foreshortening
        has_lnc = True
    else:
        has_lnc = False

    pymedct = {}
        
    # this way we ensure that the valid vs invalid portions of the same trace get separate ids
    #      it becomes important for calculating std_devs for traces which are otherwise contamined by NaNs
    #      from the invalid part of a trace
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
                    # NOTE CS 3/2024: there seems to be an extra iteration in the newer files with MBM
                    #  in some properties these are NAN, for eco this seems 0, so ok to still use sum along whole axis
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
        if key in data.dtype.fields:
            pymedct[key] = data[key].astype('i') # these are either integer types or should be converted to integer

    # TODO: think this through - we don't really need a dataframe here,
    # could return a record array, or at least make that optional
    pymepd = pd.DataFrame.from_dict(pymedct)
    if return_original_array:
        return (pymepd,data)
    else:
        return pymepd


# convenience function                           
def save_minflux_as_csv(pd_data, fname):
    pd_data.to_csv(fname,index=False)

# we are monkeypatching pipeline and VisGUIFrame methods to sneak MINFLUX npy IO in;
# this gets called from the MINFLUX plugin in the Plug routine;
# this way it can patch the relevant VisGUIFrame and Pipeline methods in the instances
# of these classes in the visGUI app
#
# in future we will ask for a way to get this considered by David B for a proper hook
# in the file loading code and possibly allow registering file load hooks for new formats
def monkeypatch_npy_io(visFr):
    import types
    import logging
    import os
    import wx
    from PYME.IO import MetaDataHandler
    from PYME.IO.FileUtils import nameUtils

    logger = logging.getLogger(__name__)
    logger.info("MINFLUX monkeypatching IO")
    def _populate_open_args_npy(self, filename):
        # this is currently just the minmal functionality for .npy,
        # we should really check a few things before going any further
        # .mat and CSV files give examples...
        if os.path.splitext(filename)[1] == '.npy':
            valid, warnmsg = npy_is_minflux_data(filename,warning=False,return_msg=True)
            if not valid:
                warn('file "%s" does not look like a valid MINFLUX NPY file:\n"%s"\n\nOPENING ABORTED'
                     % (os.path.basename(filename),warnmsg))
                return # this is not MINFLUX NPY data - we give up
            return {} # all good, just return empty args
        else:
            return self._populate_open_args_original(filename)

    visFr._populate_open_args_original = visFr._populate_open_args
    visFr._populate_open_args = types.MethodType(_populate_open_args_npy,visFr)

    def _load_ds_npy(filename):
        from PYMEcs.IO.tabular import MinfluxNpySource
        ds = MinfluxNpySource(filename)
        ds.filename = filename
            
        ds.mdh = MetaDataHandler.NestedClassMDHandler()
        data = np.load(filename)

        from pathlib import Path
        ds.mdh['MINFLUX.Is3D'] = minflux_npy_detect_3D(data)
        ds.mdh['MINFLUX.ExtraIteration'] = minflux_npy_has_extra_iter(data)
        ds.mdh['MINFLUX.Filename'] = Path(filename).name # the MINFLUX filename holds some metadata
        ds.mdh['MINFLUX.Foreshortening'] = foreshortening
        from PYMEcs.misc.utils import get_timestamp_from_filename, parse_timestamp_from_filename
        ts = get_timestamp_from_filename(filename)
        if ts is not None:
            ds.mdh['MINFLUX.TimeStamp'] = ts
            # we add the zero to defeat the regexp that checks for names ending with 'time$'
            # this falls foul of the comparison with an int (epoch time) in the metadata repr function
            # because our time stamp is a pandas time stamp and comparison with int fails
            ds.mdh['MINFLUX.StartTime0'] = parse_timestamp_from_filename(filename)
        return ds
    
    def _ds_from_file_npy(self, filename, **kwargs):
        if os.path.splitext(filename)[1] == '.npy': # MINFLUX NPY file
            logger.info('.npy file, trying to load as MINFLUX npy ...')
            return _load_ds_npy(filename)
        else:
            return self._ds_from_file_original(filename, **kwargs)

    visFr.pipeline._ds_from_file_original = visFr.pipeline._ds_from_file
    visFr.pipeline._ds_from_file = types.MethodType(_ds_from_file_npy,visFr.pipeline)


    ### we now also need to monkey_patch the _load_input method of the pipeline recipe
    ### this should allow session loading to succeed
    def _load_input_npy(self, filename, key='input', metadata_defaults={}, cache={}, default_to_image=True, args={}):
        """
        Load input data from a file and inject into namespace
        """
        from PYME.IO import unifiedIO
        import os

        if '?' in filename:
            self._load_input_original(filename,key=key,metadata_defaults=metadata_defaults,
                                      cache=cache,default_to_image=default_to_image,args=args)
        if os.path.splitext(filename)[1] == '.npy': # MINFLUX NPY file
            logger.info('.npy file, trying to load as MINFLUX npy ...')
            self.namespace[key] = _load_ds_npy(filename)
        else:
            self._load_input_original(filename,key=key,metadata_defaults=metadata_defaults,
                                      cache=cache,default_to_image=default_to_image,args=args)

    if '_load_input' in dir(visFr.pipeline.recipe):
        visFr.pipeline.recipe._load_input_original = visFr.pipeline.recipe._load_input
        visFr.pipeline.recipe._load_input = types.MethodType(_load_input_npy,visFr.pipeline.recipe)
 
    # we install this as new Menu item as File>Open is already assigned
    # however the new File>Open MINFLUX NPY entry can also open all other allowed file types
    def OnOpenFileNPY(self, event):
        filename = wx.FileSelector("Choose a file to open", 
                                   nameUtils.genResultDirectoryPath(), 
                                   wildcard='|'.join(['All supported formats|*.h5r;*.txt;*.mat;*.csv;*.hdf;*.3d;*.3dlp;*.npy',
                                                      'PYME Results Files (*.h5r)|*.h5r',
                                                      'Tab Formatted Text (*.txt)|*.txt',
                                                      'Matlab data (*.mat)|*.mat',
                                                      'Comma separated values (*.csv)|*.csv',
                                                      'HDF Tabular (*.hdf)|*.hdf',
                                                      'MINFLUX NPY (*.npy)|*.npy']))

        if not filename == '':
            self.OpenFile(filename)

    
    visFr.OnOpenFileNPY = types.MethodType(OnOpenFileNPY,visFr)
    visFr.AddMenuItem('File', "Open MINFLUX NPY", visFr.OnOpenFileNPY)
    
    logger.info("MINFLUX monkeypatching IO completed")

    # set option to make choosing filetype options available in FileDialogs on macOS
    # seems to be ok to be set on non-macOS systems, too
    wx.SystemOptions.SetOption(u"osx.openfiledialog.always-show-types", 1)

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
            if not npy_is_minflux_data(filename,warning=True):
                raise RuntimeError("can't read pipeline data from NPY file - not a MINFLUX data set")
            from PYMEcs.IO.tabular import MinfluxNpySource
            ds = MinfluxNpySource(filename)
            ds.mdh = MetaDataHandler.NestedClassMDHandler()
            data = np.load(filename)
            ds.mdh['MINFLUX.Is3D'] = minflux_npy_detect_3D(data)
            ds.mdh['MINFLUX.ExtraIteration'] = minflux_npy_has_extra_iter(data)
            return ds
        else:
            return super()._ds_from_file(filename, **kwargs)
