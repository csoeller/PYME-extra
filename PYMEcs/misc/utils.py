import logging
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

def read_temp_csv(filename,timeformat='%d/%m/%Y %H:%M'):
    def remap_names(name): # for slightly more robust comlumn renaming
        if 'Rack' in name:
            return 'Rack'
        elif 'Box' in name:
            return 'Box'
        elif 'Stativ' in name:
            return 'Stand'
        else:
            return name
    trec = pd.read_csv(filename,encoding = "ISO-8859-1")
    trec['datetime'] = pd.to_datetime(trec['Time'],format=timeformat)
    return trec.rename(columns=remap_names)

def set_diff(trec,t0):
    trec['tdiff'] = trec['datetime'] - t0
    trec['tdiff_s'] = trec['tdiff'].dt.total_seconds().astype('f')

from PYME.warnings import warn
def get_timestamp_from_filename(fname):
    from pathlib import Path
    import re
    
    basename = Path(fname).name
    match = re.search(r'2[3-5]\d{4}-\d{6}',basename)
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
