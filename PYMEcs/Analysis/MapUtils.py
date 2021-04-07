import numpy as np
from PYME.IO.MetaDataHandler import NestedClassMDHandler
from PYME.IO.image import ImageStack


def defaultmapname(mdh):
    if mdh['CameraMap.Type'] == 'mean':
        prefix = 'dark'
    else:
        prefix = 'variance'

    itime = int(1000*mdh['Camera.IntegrationTime'])

    return "%s_%dms.%s" % (prefix, itime, 'tif')

def mkdestarr(img):
    sensorSize = [2048,2048]
    sensorSize[0] = int(img.mdh['Camera.SensorWidth'])
    sensorSize[1] = int(img.mdh['Camera.SensorHeight'])

    destmap = np.zeros(sensorSize,dtype='float64')


    if img.mdh['CameraMap.Type'] == 'mean':
        maptype = 'dark'
        destmap.fill(img.mdh['Camera.ADOffset'])
    else:
        maptype = 'variance'
        destmap.fill(img.mdh['Camera.ReadNoise']**2)

    return destmap

def insertvmap(sourceim, destarr,validMap):
    smdh = sourceim.mdh
    px = int(smdh['CameraMap.ValidROI.ROIOriginX']) # zero based index
    py = int(smdh['CameraMap.ValidROI.ROIOriginY']) # zero based index
    wx = int(smdh['CameraMap.ValidROI.ROIWidth'])
    wy = int(smdh['CameraMap.ValidROI.ROIHeight'])

    destarr[px:px+wx,py:py+wy] = sourceim.data[px:px+wx,py:py+wy,:].squeeze()
    validMap[px:px+wx,py:py+wy] = 1

# not implemented yet
def checkMapCompat(img,mdh):
    pass

# we need to add some sanity checking so that the right maps are combined
# and not incorrect ones are combined by accident
def combine_maps(maps, return_validMap=False):
    destarr = None
    mapimgs = []
    for map in maps:
        mapimg = ImageStack(filename=map)
        mapimgs.append(mapimg)
        if destarr is None:
            mdh = NestedClassMDHandler(mapimg.mdh)
            destarr = mkdestarr(mapimg)
            validMap = np.zeros_like(destarr,dtype='int')
        else:
            checkMapCompat(mapimg,mdh)
        insertvmap(mapimg, destarr, validMap)

    mdh.setEntry('CameraMap.combinedFromMaps', maps)
    mdh.setEntry('CameraMap.ValidROI.ROIHeight', mapimgs[0].mdh['Camera.SensorHeight'])
    mdh.setEntry('CameraMap.ValidROI.ROIWidth', mapimgs[0].mdh['Camera.SensorWidth'])
    mdh.setEntry('CameraMap.ValidROI.ROIOriginX', 0)
    mdh.setEntry('CameraMap.ValidROI.ROIOriginY', 0)

    combinedMap = ImageStack(destarr, mdh=mdh)
    if return_validMap:
        vmdh = NestedClassMDHandler(mdh)
        vmdh.setEntry('CameraMap.ValidMask', True)
        return (combinedMap,ImageStack(validMap, mdh=vmdh))
    else:
        return combinedMap


