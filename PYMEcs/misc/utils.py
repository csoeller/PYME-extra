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
