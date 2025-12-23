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
foreshortening = PYME.config.get('MINFLUX-foreshortening',0.72)

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

from PYMEcs.pyme_warnings import warn
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

def zip_is_minflux_zarr_data(filename, warning=False, return_msg=False): # currently just placeholder
    valid = True
    msg = None

    if not valid and warning:
        if not msg is None:
                warn(msg)

    if return_msg:
        return (valid,msg)
    else:
        return valid

def minflux_npy_new_format(data):
    return 'fnl' in data.dtype.fields

# wrapper around legacy vs new format IO
def minflux_npy2pyme(fname,return_original_array=False,make_clump_index=True,with_cfr_std=False):
    data = np.load(fname)
    
    if minflux_npy_new_format(data):
        pymedf = minflux_npy2pyme_new(data,
                                    make_clump_index=make_clump_index,with_cfr_std=with_cfr_std)
    else:
        pymedf = minflux_npy2pyme_legacy(data,
                                         make_clump_index=make_clump_index,with_cfr_std=with_cfr_std)

    pyme_recArray = pymedf.to_records(index=False) # convert into NUMPY recarray
    if return_original_array:
        return (pyme_recArray,data)
    else:
        return pyme_recArray

def minflux_zarr2pyme(archz,return_original_array=False,make_clump_index=True,with_cfr_std=False):
    # make data array
    mfx = archz['mfx']
    mfxv = mfx[:][mfx['vld'] == 1]
    seqidsm, incseqs = mk_seqids_maxpos(mfxv)
    data = mfxv[np.logical_not(np.isin(seqidsm,incseqs))] # remove any incomplete sequences
    pymedf = minflux_npy2pyme_new(data,
                                  make_clump_index=make_clump_index,with_cfr_std=with_cfr_std)

    pyme_recArray = pymedf.to_records(index=False) # convert into NUMPY recarray
    if return_original_array:
        return (pyme_recArray,data)
    else:
        return pyme_recArray


###############################
### MINFLUX property checks ###
###############################

# here we check for size either 5 (2D) or 10 (3D); any other size raises an error
def minflux_npy_detect_3D_legacy(data):
    if data['itr'].shape[1] == 10 or data['itr'].shape[1] == 11:
        return True # 3D
    elif data['itr'].shape[1] == 5 or data['itr'].shape[1] == 6:
        return False # 2D
    else:
        raise RuntimeError('unknown size of itr array, neither 5 (2D) nor 10 (3D), is actually: %d' %
                            (data['itr'].shape[1]))

def minflux_check_poperties(data): # this is aiming at becoming a single stop to check MINFLUX file/dataset properties
    props = {}
    props['Is3D'] = minflux_npy_detect_3D_new(data)
    props['Tracking'] = minflux_npy_detect_2Dtracking_new(data)
    if minflux_npy_new_format(data):
        props['Format'] = 'RevAutumn2024'
    else:
        props['Format'] = 'Legacy'
    return props
    
def minflux_npy_detect_3D_new(data):
    dfin = data[data['fnl'] == True]
    if dfin['itr'][0] == 9:
        if not np.all(dfin['itr'] == 9):
            raise RuntimeError('3D detected but some "last iterations" have an index different from 9, giving up')
        return True # 3D
    elif dfin['itr'][0] == 4:
        if not np.all(dfin['itr'] == 4):
            raise RuntimeError('2D detected but some "last iterations" have an index different from 4, giving up')
        return False # 2D
    elif dfin['itr'][0] == 3: # 2D tracking
        if not np.all(dfin['itr'] == 3):
            raise RuntimeError('2D tracking detected but some "last iterations" have an index different from 3, giving up')
        return False # 2D
    else:
        raise RuntimeError('unknown number of final iteration, neither 3, (2D tracking), 4 (2D) nor 9 (3D), is actually: %d' %
                            (dfin['itr'][0]))

def minflux_npy_detect_2Dtracking_new(data):
    dfin = data[data['fnl'] == True]
    if np.all(dfin['itr'] == 3):
        return True
    else:
        return False

def minflux_npy_has_extra_iter_legacy(data):
    if data['itr'].shape[1] == 6 or data['itr'].shape[1] == 11:
        return True # has a spare empty starting position
    else:
        return False


##################
### LEGACY IO ####
##################

    
# this one should be able to deal both with 2d and 3D
def minflux_npy2pyme_legacy(data,make_clump_index=True,with_cfr_std=False):

    if minflux_npy_detect_3D_legacy(data):
        is_3D = True
        iterno_loc = 9 # we pick up the most precise localisation from this iteration, also fbg
        iterno_other = 9 # we pick up dcr, efo from this iteration
        iterno_cfr = 6
    else:
        is_3D = False
        iterno_loc = 4
        iterno_other = 4
        iterno_cfr = 3

    # NOTE CS 3/2024: latest data with MBM active seems to generate an "empty" iteration (at position 0)
    # that has NaNs or zeros in the relevant properties
    # we seem to be able to deal with this by just moving our pointers into the iteration just one position up
    # this is subject to confirmation
    if minflux_npy_has_extra_iter_legacy(data):
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
        posinid = mk_posinid(ids)
        pymedct.update({'clumpIndex': ids,
                        'clumpSize' : counts,
                        'posInClump': posinid,
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
    return pymepd

#########################
### RevAutumn2024 IO ####
#########################

# below is code to generate sequence IDs for all sequences present in the mfx data
# goal is to use only fast vectorized expressions

# this one uses "final iteration" as end of sequence marker
# we noticed later that this can lead to issues with "incomplete sequences",
#  i.e. sequences that are not terminated by a valid final localisation
def mk_seqids(data):
    indexlast = data['fnl'] == True
    seq_uid = np.arange(1,indexlast.sum()+1,dtype='i')
    seqidwnans = np.full(data.shape[0],np.nan)
    seqidwnans[indexlast] = seq_uid
    dfidnan = pd.DataFrame({'seqid':seqidwnans})
    seqid = dfidnan.bfill().to_numpy(dtype='i').squeeze() # we use pandas fast backfill to mark the other events that are part of this sequence
    return seqid

# this one uses "change to a lower (or equal) sequence number" as end of sequence marker
# this seems safer than looking for an iteration with the 'fnl' marker as there can be incomplete sequences, see below
# note we now also look for "<=" in the idxmax computation which should only happen if a valid itr 0 is followed directly by another valid itr 0
# we also return a list (actually numpy array) of seqids of incomplete sequences
def mk_seqids_maxpos(data):
    idxmax = np.nonzero((data['itr'][1:]-data['itr'][0:-1]) <= 0)[0]
    seq_uid = np.arange(1,idxmax.size+1,dtype='i')
    seqidwnans = np.full(data.shape[0],np.nan)
    seqidwnans[idxmax] = seq_uid
    if np.isnan(seqidwnans[-1]):
        seqidwnans[-1] = seq_uid.max()+1 # we may need to marke the last event with a unique id
    dfidnan = pd.DataFrame({'seqid':seqidwnans})
    seqid = dfidnan.bfill().to_numpy(dtype='i').squeeze() # we use pandas fast backfill to mark the other events that are part of this sequence
    # also mark incomplete sequences for weeding out
    idxincp = idxmax[data['fnl'][idxmax] != 1] # incomplete sequences end with an event that is not marked as 'fnl'
    incomplete_seqs = seqid[idxincp]
    if data['fnl'][-1] != 1:
        incomplete_seqs = np.append(incomplete_seqs,seqid[-1])
    return seqid, incomplete_seqs

# number the position within clumps from 0 to clumpSize-1
# here we assume that the data is already strictly orderd by time of occurence
# this should generally be the case!
# the implementation is currently not as fast as would ideally be the case (we iterate over all ids)
# ideally a full vector expression would be used - but need to figure out how
# however, not yet timed if this computation is rate-limiting for the import, it may not be
#  in which case no further optimization would be currently needed
def mk_posinid(ids):
    posinid = np.zeros_like(ids)
    for curid in np.unique(ids):
        isid = ids == curid
        posinid[isid] = np.arange(int(np.sum(isid)))
    return posinid

# this one should be able to deal both with 2d and 3D
def minflux_npy2pyme_new(data,make_clump_index=True,with_cfr_std=False):
    lastits = data['fnl'] == True
    wherelast = np.nonzero(lastits)[0]
    dfin = data[lastits]

    props = minflux_check_poperties(data)

    if props['Is3D']:
        wherecfr = wherelast - 3
        if not np.all(data[wherecfr]['itr'] == 6):
            raise RuntimeError('CFR check_3D: 3D detected but some "cfr iterations" have an index different from 6, giving up')
    else:
        if props['Tracking']:
            wherecfr = wherelast # this is bogus for now; we really need to get CFR from previous itr==2 that belongs to the same trace
        else:
            wherecfr = wherelast - 1 # in 2D we do use the last but one iteration (iteration 3)
            if not np.all(data[wherecfr]['itr'] == 3):
                raise RuntimeError('CFR check_2D: 2D detected but some "cfr iterations" have an index different from 3, giving up')

    posnm = 1e9*dfin['loc'] # we keep all distances in units of nm
    posnm[:,2] *= foreshortening
    if 'lnc' in data.dtype.fields:
        posnm_nc = 1e9*dfin['lnc']
        posnm_nc[:,2] *= foreshortening
        has_lnc = True
    else:
        has_lnc = False

    pymedct = {}
        
    # this way we ensure that the valid vs invalid portions of the same trace get separate ids
    #      it becomes important for calculating std_devs for traces which are otherwise contamined by NaNs
    #      from the invalid part of a trace
    rawids = 2*dfin['tid'] + dfin['vld']

    if make_clump_index:
        # we replace the non-sequential trace ids from MINFLUX data with a set of sequential ids
        # this works better for clumpIndex assumptions in the end
        uids,revids = np.unique(rawids,return_inverse=True)
        ids = np.arange(1,uids.size+1,dtype='int32')[revids]
        posinid = mk_posinid(ids)
        counts = get_stddev_property(ids,posnm[:,0],statistic='count')
        pymedct.update({'clumpIndex': ids,
                        'clumpSize' : counts,
                        'posInClump': posinid,
                        })
    else:
        ids = rawids

    # we are currently not using the info on incomplete sequences
    seqid,incomplete_seqid  = mk_seqids_maxpos(data) # we give every sequence a unique id to allow summing up the photons
    # we assume for now the counts at offset can be used to sum up
    # and get the total photons harvested in a sequence
    nphotons_all = get_stddev_property(seqid,data['eco'],statistic='sum')
    niterations_all = get_stddev_property(seqid,data['eco'],statistic='count') # we also count how many iterations were done, to see complete vs partial sequences
    
    if with_cfr_std: # we also compute on request a cfr std dev across a trace ID (=clump in PYME)
        pymedct.update({'cfr_std':get_stddev_property(ids,data[wherecfr]['cfr'])})

    pymedct.update({'x' : posnm[:,0],
                    'y': posnm[:,1],
                    # for t we use time to ms precision (without rounding); this is a reasonably close
                    # correspondence to frame numbers as time coordinates in SMLM data
                    't': (1e3*dfin['tim']).astype('i'),
                    'cfr':data[wherecfr]['cfr'],
                    'efo':dfin['efo'],
                    'fbg': dfin['fbg'],
                    # check with abberior
                    # NOTE CS 3/2024: there seems to be an extra iteration in the newer files with MBM
                    #  in some properties these are NAN, for eco this seems 0, so ok to still use sum along whole axis
                    'tim': dfin['tim'], # we also keep the original float time index, units are [s]
                    'nPhotons': nphotons_all[wherelast],
                    'nIters': niterations_all[wherelast],
                    'itr': dfin['itr']
                    })
    # copy a few entries verbatim
    for key in ['tid','act','vld','sta','sqi','thi','gri']:
        if key in data.dtype.fields:
            pymedct[key] = dfin[key].astype('i') # these are either integer types or should be converted to integer

    # spectral colour info
    pymedct.update({'dcr':dfin['dcr'][:,0]})
    if dfin['dcr'].shape[1] > 1: # first element is ch1/(ch1 + ch2), second is ch2/(ch1 + ch2) if present
        pymedct.update({'dcr2':dfin['dcr'][:,1]})

    stdx = get_stddev_property(ids,posnm[:,0])
    # we expect this to only happen when clumpSize == 1, because then std dev comes back as 0
    stdx[stdx < 1e-3] = 100.0 # if error estimate is too small, replace with 100 as "large" flag
    stdy = get_stddev_property(ids,posnm[:,1])
    stdy[stdy < 1e-3] = 100.0
    if props['Is3D']:
        stdz = get_stddev_property(ids,posnm[:,2])
        stdz[stdz < 1e-3] = 100.0

    if props['Tracking']: # NOTE: for now 2D only, must fix in future for 3D!

        # estimating the experimental localization precision σ for each track by calculating the
        # standard deviation (SD) of coordinate difference between consecutive localizations
        # from supplement in Deguchi, T. et al. Direct observation of motor protein stepping in
        #                             living cells using MINFLUX. Science 379, 1010–1015 (2023).
        def diffstd(data):
            # take differential and then look at std_dev of that
            # 1/sqrt(2) to account for variance increase on differences
            return np.diff(data).std()/1.41
        
        track_stdx = stdx
        track_stdy = stdy
        #LOCERR_MAX = 15.0
        #stdx = np.clip(stdx,None,LOCERR_MAX) # current workaround, need better loc err estimation
        #stdy = np.clip(stdy,None,LOCERR_MAX) # current workaround, need better loc err estimation
        stdx = get_stddev_property(ids,posnm[:,0],statistic=diffstd)
        stdy = get_stddev_property(ids,posnm[:,1],statistic=diffstd)
        track_tmin = get_stddev_property(ids,dfin['tim'],'min')
        track_tms = 1e3*(dfin['tim']-track_tmin)
        track_lims = np.zeros_like(ids)
        track_lims[np.diff(ids,prepend=0) > 0] = 1 # mark beginning of tracks with 1
        track_lims[np.diff(ids,append=ids.max()+1) > 0] = 2 # mark end of tracks with 2
        pymedct.update({'track_stdx':track_stdx, 'track_stdy':track_stdy, 'track_tms':track_tms,
                        # we return track_err[xy] in addition to error_x, error_y since it avoids
                        # special treatment on coalescing and therefore allows comparison between
                        # track_stdx and track_errx etc on a per track basis
                        'track_errx':stdx.copy(), 'track_erry':stdy.copy(),
                        'track_lims':track_lims,
                        })

    pymedct.update({'error_x' : stdx,'error_y' : stdy})
    if props['Is3D']:
        pymedct.update({'z':posnm[:,2], 'error_z' : stdz})

    if has_lnc:
        pymedct.update({'x_nc' : posnm_nc[:,0],
                        'y_nc' : posnm_nc[:,1]})
        if props['Is3D']:
            pymedct.update({'z_nc' : posnm_nc[:,2]})

    # TODO: think this through - we don't really need a dataframe here,
    # could return a record array, or at least make that optional
    pymepd = pd.DataFrame.from_dict(pymedct)
    return pymepd

#################################
### MBM utility functionality ###
#################################

# we try to find an MBM collection attached
# to an MBMcorrection module generated data source
# returns None if unsuccesful
def findmbm(pipeline,warnings=True,return_mod=False):
    from PYMEcs.recipes.localisations import MBMcorrection
    dsname = None
    # search/check for instance
    for mod in pipeline.recipe.modules:
        if isinstance(mod,MBMcorrection):
            dsname = mod.output
            break
    if dsname is None:
        if warnings:
            warn("we rely on MBM info present in a datasource generated by the MBMcorrection module.\n\n" +
                 "Can't find such a datasource, please add MBMcorrection module to your recipe.\n\nAborting...")
        return None
    mbm = pipeline.dataSources[dsname].mdh.get('Processing.MBMcorrection.mbm')
    if mbm is None:
        if warnings:
            warn(("found no MBM collection in metadata of datasource '%s' generated by MBMcorrection module.\n\n" % dsname )+
                 "Have you loaded valid MBM data into module yet?\n\nAborting..." )
        return None
    if return_mod:
        return mod
    else:
        return mbm


#####################################
### metadata utility functions ######
#####################################

def _get_basic_MINFLUX_metadata(filename):
    from pathlib import Path
    mdh = MetaDataHandler.NestedClassMDHandler()

    mdh['MINFLUX.Filename'] = Path(filename).name # the MINFLUX filename often holds some metadata
    mdh['MINFLUX.Foreshortening'] = foreshortening
    from PYMEcs.misc.utils import get_timestamp_from_filename, parse_timestamp_from_filename
    ts = get_timestamp_from_filename(filename)
    if ts is not None:
        mdh['MINFLUX.TimeStamp'] = ts
        # we add the zero to defeat the regexp that checks for names ending with 'time$'
        # this falls foul of the comparison with an int (epoch time) in the metadata repr function
        # because our time stamp is a pandas time stamp and comparison with int fails
        mdh['MINFLUX.StartTime0'] = parse_timestamp_from_filename(filename).strftime("%Y-%m-%d %H:%M:%S")

    return mdh
    
def _get_mdh(data,filename):
    mdh = _get_basic_MINFLUX_metadata(filename)
    if minflux_npy_new_format(data):
        props = minflux_check_poperties(data)
        mdh['MINFLUX.Format'] = props['Format']
        mdh['MINFLUX.Is3D'] = props['Is3D']
        mdh['MINFLUX.Tracking'] = props['Tracking']
    else:
        mdh['MINFLUX.Format'] = 'Legacy'
        mdh['MINFLUX.Is3D'] = minflux_npy_detect_3D_legacy(data)
        mdh['MINFLUX.ExtraIteration'] = minflux_npy_has_extra_iter_legacy(data)
        mdh['MINFLUX.Tracking'] = False # for now we do not support tracking with legacy data

    return mdh

def _get_mdh_zarr(filename,arch):
    mdh = _get_basic_MINFLUX_metadata(filename)
    mfx_attrs = arch['mfx'].attrs.asdict()
    if not '_legacy' in mfx_attrs:
        mdh['MINFLUX.Format'] = 'RevAutumn2024'
        mdh['MINFLUX.AcquisitionDate'] = mfx_attrs['acquisition_date']
        mdh['MINFLUX.DataID'] = mfx_attrs['did']
        mdh['MINFLUX.Is3D'] = mfx_attrs['measurement']['dimensionality'] > 2
        # now do some checks of acquisitiondate vs any filename derived info
        from PYMEcs.misc.utils import get_timestamp_from_mdh_acqdate, compare_timestamps_s
        ts = get_timestamp_from_mdh_acqdate(mdh)
        if ts is not None:
            mts = mdh.get('MINFLUX.TimeStamp')
            if mts is not None:
                if mts != ts:
                    delta_s = compare_timestamps_s(mts,ts)
                    if delta_s > 5: # there can be rounding errors from the different TS sources, we tolerate up to 5s difference
                        warn("acq time stamp (%s) not equal to filename time stamp (%s), delta in s is %d" % (ts,mts,delta_s))
            else:
                mdh['MINFLUX.TimeStamp'] = ts

        md_by_itrs,mfx_global_par = get_metadata_from_mfx_attrs(mfx_attrs)
        for par in mfx_global_par:
            mdh['MINFLUX.Globals.%s' % par] = mfx_global_par[par]
        for pars in md_by_itrs:
            # make sure we convert to list; otherwise we cannot easily convert to JSON as JSON does not like ndarray
            mdh['MINFLUX.ByItrs.%s' % pars] = md_by_itrs[pars].to_numpy().tolist()
        import re
        mdh['MINFLUX.Tracking'] = re.search('tracking', mfx_global_par['ID'], re.IGNORECASE) is not None
    else:
        mdh['MINFLUX.Format'] = 'LegacyZarrConversion'
        mdh['MINFLUX.Is3D'] = mfx_attrs['_legacy']['_seqs'][0]['Itr'][0]['Mode']['dim'] > 2
            
    return mdh

def get_metadata_from_mfx_attrs(mfx_attrs):
    mfx_itrs = mfx_attrs['measurement']['threads'][0]['sequences'][0]['Itr']
    mfx_globals = mfx_attrs['measurement']['threads'][0]['sequences'][0]
    
    md_by_itrs = pd.DataFrame(columns=['IterationNumber','PinholeAU','ActivationLaser', 'ExcitationLaserAbbrev',
                                       'ExcitationWavelength_nm', 'ExcitationPower_percent', 'ExcitationDAC',
                                       'DetectionChannel01','DetectionChannel02','BackgroundThreshold',
                                       'PhotonLimit', 'CCRLimit', 'DwellTime_ms',
                                       'PatternGeoFactor','PatternRepeat', 'PatternGeometryAbbrev','Strategy'],
                              index=range(len(mfx_itrs)))
    for i, itr in enumerate(mfx_itrs):
        md_by_itrs.loc[i].IterationNumber = i
        md_by_itrs.loc[i].PinholeAU = itr['Mode']['phDiaAU']
        md_by_itrs.loc[i].ActivationLaser = itr['_activation']['laser'] if itr['_activation']['laser'] != '' else 'NA'
        md_by_itrs.loc[i].ExcitationLaserAbbrev = itr['_excitation']['laser'].replace('MINFLUX','M')
        md_by_itrs.loc[i].ExcitationWavelength_nm = np.rint(1e9*itr['_excitation']['wavelength'])
        md_by_itrs.loc[i].ExcitationPower_percent = itr['_excitation']['power']
        md_by_itrs.loc[i].ExcitationDAC = itr['_excitation']['dac']
        md_by_itrs.loc[i].DetectionChannel01 = itr['_detection']['channels'][0]
        md_by_itrs.loc[i].DetectionChannel02 = itr['_detection']['channels'][1] if len(itr['_detection']['channels']) >1 else 'NA'
        md_by_itrs.loc[i].BackgroundThreshold = itr['bgcThreshold']
        md_by_itrs.loc[i].PhotonLimit = itr['phtLimit']
        md_by_itrs.loc[i].CCRLimit = itr['ccrLimit']
        md_by_itrs.loc[i].DwellTime_ms = 1e3*itr['patDwellTime']
        md_by_itrs.loc[i].PatternGeoFactor = itr['patGeoFactor']
        md_by_itrs.loc[i].PatternRepeat = itr['patRepeat']
        md_by_itrs.loc[i].PatternGeometryAbbrev = itr['Mode']['pattern'].replace('hexagon','hex').replace('zline','zl').replace('square','sq')
        md_by_itrs.loc[i].Strategy = itr['Mode']['strategy']

    mfx_global_pars = {}
    
    mfx_global_pars['BgcSense'] = mfx_globals['bgcSense']
    mfx_global_pars['CtrDwellFactor'] = mfx_globals['ctrDwellFactor']
    mfx_global_pars['Damping'] = mfx_globals['damping']
    mfx_global_pars['Headstart'] = mfx_globals['headstart']
    mfx_global_pars['ID'] = mfx_globals['id']
    mfx_global_pars['Liveview'] = mfx_globals['liveview']['show']
    mfx_global_pars['LocLimit'] = mfx_globals['locLimit']
    mfx_global_pars['Stickiness'] = mfx_globals['stickiness']
    mfx_global_pars['FieldAlgorithm'] = mfx_globals['field']['algo']
    mfx_global_pars['FieldGeoFactor'] = mfx_globals['field']['fldGeoFactor']
    mfx_global_pars['FieldStride'] = mfx_globals['field']['stride']

    return (md_by_itrs,mfx_global_pars)


##############################
### tabular classes  #########
##############################

from PYME.IO.tabular import TabularBase

# closely modeled on RecArraySource
class MinfluxNpySource(TabularBase):
    _name = "MINFLUX NPY File Source"
    def __init__(self, filename):
        """ Input filter for use with NPY data exported from MINFLUX data (typically residing in MSR files)."""

        self.res = minflux_npy2pyme(filename)

        # check for invalid localisations:
        # possible TODO - is this needed/helpful, or should we propagate missing values further?
        # FIXED - minflux_npy2pyme should now also work properly when invalid data is present
        #         so returning just the valid events to PYME should be ok
        if np.any(self.res['vld'] < 1):
            self.res = self.res[self.res['vld'] >= 1]

        self._keys = list(self.res.dtype.names)
       

    def keys(self):
        return self._keys

    def __getitem__(self, keys):
        key, sl = self._getKeySlice(keys)
        
        if not key in self._keys:
            raise KeyError('Key (%s) not found' % key)

       
        return self.res[key][sl]

    
    def getInfo(self):
        return 'MINFLUX NPY Data Source\n\n %d points' % len(self.res['x'])

class MinfluxZarrSource(MinfluxNpySource):
    _name = "MINFLUX zarr File Source"
    def __init__(self, filename):
        """ Input filter for use with ZARR data exported from MINFLUX data (originally residing in MSR files)."""
        import zarr
        archz = zarr.open(filename)
        self.zarr = archz
        self._own_file = True # is this necessary? Normally only used by HDF to close HFD on destroy, zarr does not need "closing"

        # NOTE: no further 'locations valid' check should be necessary - we filter already in the conversion function
        self.res = minflux_zarr2pyme(archz)
        
        self._keys = list(self.res.dtype.names)

        # note: aparently, closing an open zarr archive is not required; accordingly no delete and close methods necessary
        self._paraflux_analysis = None
    
##############################
### Register IO with PYME ####
##############################

# we are monkeypatching pipeline and VisGUIFrame methods to sneak MINFLUX npy IO in;
# this gets called from the MINFLUX plugin in the Plug routine;
# this way it can patch the relevant VisGUIFrame and Pipeline methods in the instances
# of these classes in the visGUI app
#
# in future we will ask for a way to get this considered by David B for a proper hook
# in the file loading code and possibly allow registering file load hooks for new formats
def monkeypatch_npyorzarr_io(visFr):
    import types
    import logging
    import os
    import wx
    from PYME.IO import MetaDataHandler
    from PYME.IO.FileUtils import nameUtils

    logger = logging.getLogger(__name__)
    logger.info("MINFLUX monkeypatching IO")
    def _populate_open_args_npyorzarr(self, filename):
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
        elif  os.path.splitext(filename)[1] == '.zip':
            valid, warnmsg = zip_is_minflux_zarr_data(filename,warning=False,return_msg=True)
            if not valid:
                warn('file "%s" does not look like a valid MINFLUX zarr file:\n"%s"\n\nOPENING ABORTED'
                     % (os.path.basename(filename),warnmsg))
                return # this is not MINFLUX zarr data - we give up
            return {} # all good, just return empty args
        else:
            return self._populate_open_args_original(filename)

    visFr._populate_open_args_original = visFr._populate_open_args
    visFr._populate_open_args = types.MethodType(_populate_open_args_npyorzarr,visFr)

    def _load_ds_npy(filename):
        ds = MinfluxNpySource(filename)
        ds.filename = filename
        
        data = np.load(filename)
        ds.mdh = _get_mdh(data,filename)

        return ds

    def _load_ds_zarrzip(filename):
        ds = MinfluxZarrSource(filename)
        ds.filename = filename
        
        ds.mdh = _get_mdh_zarr(filename,ds.zarr)

        return ds
    
    def _ds_from_file_npyorzarr(self, filename, **kwargs):
        if os.path.splitext(filename)[1] == '.npy': # MINFLUX NPY file
            logger.info('.npy file, trying to load as MINFLUX npy ...')
            return _load_ds_npy(filename)
        elif os.path.splitext(filename)[1] == '.zip': # MINFLUX ZARR file in zip format
            logger.info('.zip file, trying to load as MINFLUX zarr ...')
            return _load_ds_zarrzip(filename)
        else:
            return self._ds_from_file_original(filename, **kwargs)

    visFr.pipeline._ds_from_file_original = visFr.pipeline._ds_from_file
    visFr.pipeline._ds_from_file = types.MethodType(_ds_from_file_npyorzarr,visFr.pipeline)

    from PYMEcs.IO.NPC import findNPCset
    visFr.pipeline.get_npcs = types.MethodType(findNPCset,visFr.pipeline) # we make this a method for pipeline to make access easier

    ### we now also need to monkey_patch the _load_input method of the pipeline recipe
    ### this should allow session loading to succeed
    def _load_input_npyorzarr(self, filename, key='input', metadata_defaults={}, cache={}, default_to_image=True, args={}):
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
        elif os.path.splitext(filename)[1] == '.zip': # MINFLUX NPY file
            logger.info('.npy file, trying to load as MINFLUX zarr ...')
            self.namespace[key] = _load_ds_zarrzip(filename)
        else:
            self._load_input_original(filename,key=key,metadata_defaults=metadata_defaults,
                                      cache=cache,default_to_image=default_to_image,args=args)

    if '_load_input' in dir(visFr.pipeline.recipe):
        visFr.pipeline.recipe._load_input_original = visFr.pipeline.recipe._load_input
        visFr.pipeline.recipe._load_input = types.MethodType(_load_input_npyorzarr,visFr.pipeline.recipe)
 
    # we install this as new Menu item as File>Open is already assigned
    # however the new File>Open MINFLUX NPY entry can also open all other allowed file types
    def OnOpenFileNPYorZARR(self, event):
        filename = wx.FileSelector("Choose a file to open", 
                                   nameUtils.genResultDirectoryPath(), 
                                   wildcard='|'.join(['All supported formats|*.h5r;*.txt;*.mat;*.csv;*.hdf;*.3d;*.3dlp;*.npy;*.zip;*.pvs',
                                                      'PYME Results Files (*.h5r)|*.h5r',
                                                      'Tab Formatted Text (*.txt)|*.txt',
                                                      'Matlab data (*.mat)|*.mat',
                                                      'Comma separated values (*.csv)|*.csv',
                                                      'HDF Tabular (*.hdf)|*.hdf',
                                                      'MINFLUX NPY (*.npy)|*.npy',
                                                      'MINFLUX ZARR (*.zip)|*.zip',
                                                      'Session files (*.pvs)|*.pvs',]))

        if not filename == '':
            self.OpenFile(filename)

    
    visFr.OnOpenFileNPYorZARR = types.MethodType(OnOpenFileNPYorZARR,visFr)
    visFr.AddMenuItem('File', "Open MINFLUX NPY, zarr or session", visFr.OnOpenFileNPYorZARR)
    
    logger.info("MINFLUX monkeypatching IO completed")

    # set option to make choosing filetype options available in FileDialogs on macOS
    # seems to be ok to be set on non-macOS systems, too
    wx.SystemOptions.SetOption(u"osx.openfiledialog.always-show-types", 1)

    def _get_session_datasources_whook(self): # with hook for saving lowess cache
        # try to save an mbm lowess cache if present
        mod = findmbm(visFr.pipeline,warnings=False,return_mod=True)
        mbm = findmbm(visFr.pipeline,warnings=False,return_mod=False)
        if mod is not None and mbm is not None:
            if mod.MBM_lowess_fraction > 1e-5:
                if not mod.lowess_cachefilepath().exists():
                    mod.lowess_cachesave()

        return self._get_session_datasources_original()

    visFr.pipeline._get_session_datasources_original = visFr.pipeline._get_session_datasources
    visFr.pipeline._get_session_datasources = types.MethodType(_get_session_datasources_whook,visFr.pipeline)

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
                raise RuntimeError("can't read pipeline data from NPY file - not a MINFLUX data set?")
            ds = MinfluxNpySource(filename)
            data = np.load(filename)
            ds.mdh =  _get_mdh(data,filename)
            return ds
        elif os.path.splitext(filename)[1] == '.zip': # MINFLUX zarr file
            logging.getLogger(__name__).info('.zip file, trying to load as MINFLUX zarr ...')
            if not zip_is_minflux_zarr_data(filename,warning=True):
                raise RuntimeError("can't read pipeline data from MINFLUX zarr file - not a MINFLUX data set?")
            ds = MinfluxZarrSource(filename)
            ds.mdh = _get_mdh_zarr(filename,ds.zarr)
            return ds
        else:
            return super()._ds_from_file(filename, **kwargs)
