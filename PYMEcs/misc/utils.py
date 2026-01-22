import logging
logger = logging.getLogger(__file__)

# filtering for a few sources of messages that we can generally blend out in notenooks

def pyme_logging_filter(loglevel=logging.WARN,
                        filterTextFileSource=True,
                        filterDictMDHandlerWarning=True,
                        filterTrackUtilsWarning=True):
    if filterTextFileSource:
        from PYME.IO.tabular import logger as tabular_logger
        # this will filter our logging message as long as it uses logging
        def textFileSource_filter(record):
            if 'TextFileSource-use_pandas' in record.msg:
                return False
            return True
        tabular_logger.addFilter(textFileSource_filter)
        
    if filterDictMDHandlerWarning:
        import warnings
        # supress warnings from the DictMDHandler about inability to handle localisations
        warnings.filterwarnings("ignore",message=r'DictMDHandler')
    if filterTrackUtilsWarning:
        import warnings
        # supress warnings from trackutils about lacking mpld3 (which we do not really need)
        warnings.filterwarnings("ignore",message=r'Could not import mpld3')
        
    logging.basicConfig()
    logging.getLogger().setLevel(loglevel)

# get unique name for recipe output
def unique_name(stem,names):
    if stem not in names:
        return stem
    for i in range(1,11):
        stem2 = "%s_%d" % (stem,i)
        if stem2 not in names:
            return stem2

    return stem2 # here we just give up and accept a duplicate name

import pandas as pd


def setTempDataFolder(modename,configvar,parent):
    import wx
    import PYME.config as config
    curfolder = config.get(configvar)
    if curfolder is None:
        warn("currently the %s temperature file folder is not set. Set the path in the following file dialog" % modename)
    else:
        warn("%s temperature file folder currently set to '%s'.\n\nAlter in the following dialog if needed or cancel that dialog if want to leave as is" %
             (modename,curfolder))
    with wx.DirDialog(parent, "Choose folder containing temperature CSV files") as dialog:
        if dialog.ShowModal() == wx.ID_CANCEL:
            return
        folder = dialog.GetPath()
    if config.get(configvar) == folder:
        warn("config option '%s' already set to %s, leaving as is" % (configvar,folder))
        return # already set to this value, return

    config.update_config({configvar: folder},
                         config='user', create_backup=True)


def read_room_temperature_csv():
    import PYME.config as config
    from glob import glob
    import os
    import pandas as pd

    def read_temp_csv(file):
        trec = pd.read_csv(file,sep=';',index_col='Timestamp',parse_dates=['Timestamp'],dayfirst=True)
        return trec

    configvar = 'MINFLUX-room-temperature_folder'
    folder = config.get(configvar)
    if folder is None:
        warn("Need to set room temperature file location first by setting config variable %s" % configvar)
        return
    elif not os.path.isdir(folder):
        warn(("Config variable %s is not set to a folder;\n" % (configvar)) +
             ("needs to be a **folder** location, currently set to %s" % (folder)))
        return
    
    tempdata = []
    reps = []
    for year in [2025,2026,2027,2028,2029,2030]:
        if glob(os.path.join(folder, 'Report_*_%d.csv' % year)):
            reps += sorted(glob(os.path.join(folder, 'Report_*_%d.csv' % year)))
    if len(reps) > 0:
        for rep in reps:
            tempdata.append(read_temp_csv(rep))
        df = pd.concat(tempdata)
        df.index.names = ['datetime']
        return df.reset_index() # move datetime into a value and replace with numeric index
    else:
        return None

# makes the reading a little more flexible
# contributed by Alex B
def read_temperature_csv(filename, timeformat=['%d.%m.%Y %H:%M:%S', # Newest format
                                        '%d/%m/%Y %H:%M:%S' # Original format
                                        ]):
    import re
    def remap_names(name):
        if re.search(r'\bRack\b', name, re.IGNORECASE):
            return 'Rack'
        elif re.search(r'\bBox\b', name, re.IGNORECASE):
            return 'Box'
        elif re.search(r'\bStativ\b', name, re.IGNORECASE):
            return 'Stand'
        elif re.search(r'\bTime\b', name, re.IGNORECASE):
            return 'Time'
        else:
            return name
        
    trec = pd.read_csv(filename, encoding="ISO-8859-1")
    trec.columns = [remap_names(col) for col in trec.columns]
    
    # Ensure timeformat is a list (even if only one format is provided)
    if isinstance(timeformat, str):
        timeformat = [timeformat]
        
    # Try all provided time formats
    for fmt in timeformat:
        try:
            trec['datetime'] = pd.to_datetime(trec['Time'], format=fmt)
            break
        except ValueError:
            continue
    else: # we get here if the for loop terminates without breaking implying no format matched
        raise ValueError("None of the provided time formats matched the 'Time' column.")

    return trec

def set_diff(trec,t0):
    trec['tdiff'] = trec['datetime'] - t0
    trec['tdiff_s'] = trec['tdiff'].dt.total_seconds().astype('f')

from PYMEcs.pyme_warnings import warn
def get_timestamp_from_filename(fname):
    from pathlib import Path
    import re
    
    basename = Path(fname).name
    match = re.search(r'[2-3][0-9]\d{4}-\d{6}',basename) # we match with data from 2020 to 2039
    if match:
        timestamp = match.group()
        return timestamp
    else:
        warn("no timestamp match found in %s" % basename)
        return None

def get_timestamp_from_mdh_acqdate(mdh):
    from datetime import datetime
    acqdate = mdh.get('MINFLUX.AcquisitionDate')
    if acqdate is not None:
        ti = datetime.strptime(acqdate,'%Y-%m-%dT%H:%M:%S%z')
        return ti.strftime('%y%m%d-%H%M%S')
    else:
        return None

def compare_timestamps_s(ts1,ts2):
    t1 = timestamp_to_datetime(ts1)
    t2 = timestamp_to_datetime(ts2)
    if t1 > t2:
        delta = t1-t2
    else:
        delta = t2-t1
    return delta.seconds

def timestamp_to_datetime(ts):
    t0 = pd.to_datetime(ts,format="%y%m%d-%H%M%S")
    return t0

def parse_timestamp_from_filename(fname):   
    timestamp = get_timestamp_from_filename(fname)
    if timestamp is None:
        return None
    t0 = timestamp_to_datetime(timestamp)
    return t0

def recipe_from_mdh(mdh):
    import re
    separator = '|'.join([ # we "or"-combine the following regexs
        '(?<=:) (?= )', # a space preceeded by a colon AND also followed by another space ; this is therefore not a "key: value" type YAML line
        '(?<![ :-]) '   # a space NOT preceded by a colon, dash or another space
    ])
    recipe = mdh.get('Pipeline.Recipe')
    if recipe is not None:
        return('\n'.join(re.split(separator,recipe)))
    else:
        warn("could not retrieve Pipeline.Recipe")
        return None

def load_sessionfile(filename,substitute=True):
    import yaml
    from PYME.LMVis.sessionpaths import substitute_sessiondir
        
    with open(filename, 'r') as f:
        session_txt = f.read()

    if substitute:
        # replace any possibly present SESSIONDIR_TOKEN
        session = yaml.safe_load(substitute_sessiondir(session_txt, filename))
    else:
        session = yaml.safe_load(session_txt)
        
    return session

from pathlib import Path
def zarrtozipstore(zarr_root,archive_name,verbose=False):
    zarr_root = Path(zarr_root)
    archive_name = Path(archive_name)

    from shutil import get_archive_formats
    if 'zip' not in dict(get_archive_formats()):
        raise RuntimeError('shutil.make_archive does not support zip format, aborting')
        
    if not (zarr_root.exists() and zarr_root.is_dir()):
        raise RuntimeError('path "%s" does not exist or is not a directory' % (zarr_root))
    if not (zarr_root / '.zgroup').exists():
        warn("did not find .zgroup file in directory, this may not be a zarr directory")

    if verbose:
        warn("zarr file archive at\n'%s'\n, zipping to dir\n'%s'\n with name '%s'" %
             (zarr_root,archive_name.parent,archive_name.name))
    
    from shutil import make_archive
    created = make_archive(archive_name,
                           'zip',
                           root_dir=zarr_root)
    return created

def get_ds_path(pipeline,ds='FitResults'):
    if 'filename' in dir(pipeline):
        return pipeline.filename()
    try:
        fpath = pipeline.get_session()['datasources'][ds]
    except AttributeError:
        fpath = None

    return fpath

def fname_from_timestamp(datapath,mdh,stemsuffix,ext='.csv'):
    from pathlib import Path
    storagepath = Path(datapath).parent
    tstamp = mdh.get('MINFLUX.TimeStamp','timestamp_unknown')

    fname = (storagepath / (tstamp + stemsuffix)).with_suffix(ext)
    return fname
  
def autosave_csv(df,datapath,mdh,suffix):
    from pathlib import Path
    fname = fname_from_timestamp(datapath,mdh,suffix,ext='.csv')
    logger.debug(f"autosaving file {fname}...")
    df.to_csv(fname, index=False, header=True)

def autosave_check():
    import PYME.config
    return PYME.config.get('MINFLUX-autosave',False)

