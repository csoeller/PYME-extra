from PYME.localization.remFitBuf import CameraInfoManager
import PYME.Analysis.gen_sCMOS_maps as gmaps
from PYME.IO.MetaDataHandler import NestedClassMDHandler
from PYME.IO.image import ImageStack
from PYME.DSView import ViewIm3D
from PYME.IO.FileUtils import nameUtils
import numpy as np

import logging
logger = logging.getLogger(__name__)

defaultCalibrationDir = nameUtils.getCalibrationDir('',create=False)

def defaultMapName(source, createPath=False, calibrationDir=defaultCalibrationDir):
    resname = source.mdh.getOrDefault('Analysis.resultname',None)
    if resname is None:
        return None

    if resname == 'mean':
        if source.mdh.getOrDefault('Analysis.FlatField',False):
            maptype = 'flatfield'
        else:
            maptype = 'dark'
    else:
        maptype = 'variance'

    mapname = gmaps.mkDefaultPath(maptype, source.mdh, create=createPath, calibrationDir=calibrationDir)
    return mapname

# return a list of existing camera directories that have tif files inside (assuming these are camera maps)
def installedCams(calibrationDir=defaultCalibrationDir):
    from glob import glob
    import os

    camdirs = []
    for (dirpath, dirnames, filenames) in os.walk(calibrationDir):
        camdirs.extend(dirnames)
        break
    fulldirs = [os.path.join(calibrationDir,cdir) for cdir in camdirs]
    validdirs = [cdir for cdir in fulldirs if glob(os.path.join(cdir, '*.tif'))]

    return validdirs


# source passed as PYME ImageStack
def install_map(source, calibrationDir=defaultCalibrationDir):
    """Installs a map file to a calibration directory. By default uses the system claibration directory."""

    import os
    if source.mdh.getOrDefault('Analysis.name', '') != 'mean-variance':
        msg = 'map %s, Analysis.name is not equal to "mean-variance" - probably not a map' % source.filename
        return msg

    validROIHeight = source.mdh.getOrDefault('Analysis.valid.ROIHeight',
                                             source.mdh['Camera.ROIHeight'])
    validROIWidth = source.mdh.getOrDefault('Analysis.valid.ROIWidth',
                                             source.mdh['Camera.ROIWidth'])
    if not (validROIHeight == source.mdh['Camera.ROIHeight']
            and validROIWidth == source.mdh['Camera.ROIWidth']):
        msg = 'Partial (ROI based) maps cannot be installed to a calibration directory'
        return msg

    if source.mdh.getOrDefault('Analysis.isuniform', False):
        msg = 'Uniform maps cannot be installed to ba calibration directory'
        return msg

    if source.mdh['Analysis.resultname'] == 'mean':
        if source.mdh.getOrDefault('Analysis.FlatField',False):
            maptype = 'flatfield'
        else:
            maptype = 'dark'
    else:
        maptype = 'variance'

    mapname = gmaps.mkDefaultPath(maptype, source.mdh, create=True, calibrationDir=calibrationDir)
    if os.path.isfile(mapname):
        msg = 'map %s exists, not overwriting' % mapname
        return msg
    
    source.Save(filename=mapname)
    return None
    
# attempt to install a map in a calibration dir but check a few things
# 1) does it have the .tif extension?
# 2) can it be opened as a tiffstack
# 3) if all ok so far pass on to install_map which does additional checks; note that
#    this function does not overwrite existing maps at the destination
def checkAndInstallMap(mapf, calibrationDir=defaultCalibrationDir):
    import os
    inst = 0
    
    ext = os.path.splitext(mapf)[-1].lower()
    if ext != ".tif":
        msg = 'asked to install %s, not a tif extension' % mapf
        return (inst, msg)
    try:
        source = ImageStack(filename=mapf)
    except:
        msg = 'asked to install %s, could not open as PYME ImageStack, not a map?' % mapf
        return (inst, msg)

    msg = install_map(source, calibrationDir=calibrationDir)
    if msg is None:
        msg = 'installed map %s in default location' % mapf
        msg += "\n\t-> %s" % defaultMapName(source, calibrationDir=calibrationDir)
        inst = 1

    return (inst, msg)


# install maps, potentially several
# if fromfile is a directory attempt to install all maps below that directory
# if fromfile is a file, attempt to install that single file
def installMapsFrom(fromfile, calibrationDir=defaultCalibrationDir):
    from glob import glob
    import os
    
    msgs = []
    ntotal = 0
    if os.path.isdir(fromfile):
        # this is a directory, walk it
        msgs.append("Installing maps from directory %s -> %s Folder\n" % (fromfile,calibrationDir))
        mapfiles = [y for x in os.walk(fromfile) for y in glob(os.path.join(x[0], '*.tif'))]
        for mapf in mapfiles:
            ninst, msg = checkAndInstallMap(mapf, calibrationDir=calibrationDir)
            ntotal += ninst
            msgs.append(msg)
    else:
        ninst, msg = checkAndInstallMap(fromfile, calibrationDir=calibrationDir)
        ntotal = 1
        msgs.append(msg)
    msgs.append("\ninstalled %d maps" % ntotal)
    return ntotal, "\n".join(msgs)

    
def getInstalledMapList():
    from glob import glob
    import os
    rootdir = defaultCalibrationDir
    result = [y for x in os.walk(rootdir) for y in glob(os.path.join(x[0], '*.tif'))]

    return result


def check_mapexists(mdh, type = 'dark'):
    import os
    if type == 'dark':
        id = 'Camera.DarkMapID'
    elif type == 'variance':
        id = 'Camera.VarianceMapID'
    elif type == 'flatfield':
        id = 'Camera.FlatfieldMapID'
    else:
        raise RuntimeError('unknown map type %s' % type)
        

    mPath = mdh.getOrDefault(id,'')
    defPath = gmaps.mkDefaultPath(type,mdh,create=False)
    if os.path.exists(mPath):
        return mPath
    elif os.path.exists(defPath):
        mdh[id] = defPath
        return defPath
    else:
        mdh[id] = ''
        return None

def mk_compositeMap(sourcemap):
    # make new composite map using sourcemap to populate
    # check source is a map
    # composite map has first channel with map data, second channel is mask where map is valid
    # returns composite map (identified by its metadata)
    pass

def addMap2composite(map,compMap):
    # insert map according to valid ROI into composite map
    # update the valid channel accordingly
    # possibly perform a number of checks
    pass

def export_mapFromComposite(compMap):
    # export a normal map from the composite map
    # in the exported map we set the valid ROI to the whole chip area
    # return as image?
    pass

import os

def _getDefaultMap(ci,mdh,maptype = 'dark', return_loadedmdh = False):
    # logic: try map in default calibration dir first
    #        if not present just return default value
    mapname = gmaps.mkDefaultPath(maptype, mdh, create=False)
    mdh2 = NestedClassMDHandler(mdh)
    logger.debug("checking map file %s" % mapname)
    
    if maptype == 'dark':
        id = 'Camera.DarkMapID'
    elif maptype == 'variance':
        id = 'Camera.VarianceMapID'
    elif maptype == 'flatfield':
        id = 'Camera.FlatfieldMapID'
    else:
        raise RuntimeError('unknown map type %s' % type)

    if os.path.isfile(mapname):
        mdh2[id] = mapname
    else:
        mdh2[id] = ''
        logger.debug("map file %s does not exist" % mapname)

    if maptype == 'dark':
        theMap = ci.getDarkMap(mdh2)
    elif maptype == 'variance':
        theMap = ci.getVarianceMap(mdh2)
    elif maptype == 'flatfield':
        theMap = ci.getFlatfieldMap(mdh2)

    if return_loadedmdh:
        return (theMap,mdh2)
    else:
        return theMap

def get_dark_default(ci,mdh):
    return _getDefaultMap(ci,mdh,maptype = 'dark')

def get_variance_default(ci,mdh):
    return _getDefaultMap(ci,mdh,maptype = 'variance')
    
def get_flatfield_default(ci,mdh):
    return _getDefaultMap(ci,mdh,maptype = 'flatfield')
