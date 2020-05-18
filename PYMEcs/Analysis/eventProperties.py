import numpy as np

# in the below p is supposed to be a pipeline object
def getarea(p):
    if 'x' in p.filterKeys:
        width = p.filterKeys['x'][1]-p.filterKeys['x'][0]
    else:
        width = p.mdh.Camera.ROIWidth*p.mdh.voxelsize_nm.x
    if 'y' in p.filterKeys:
        height = p.filterKeys['y'][1]-p.filterKeys['y'][0]
    else:
        height = p.mdh.Camera.ROIHeight*p.mdh.voxelsize_nm.y
    area = 1e-6*width*height
    return area # area in um^2

def evtDensity(p):
    area = getarea(p)
    nEvents = p['x'].size

    if area > 1e-6:
        return nEvents / area
    else:
        return None
