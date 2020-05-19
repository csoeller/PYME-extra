import numpy as np


# in the below p is supposed to be a pipeline object
# some heuristics to get the area covered by the pipeline data
def getarea(p):
    if 'x' in p.filterKeys:
        width = p.filterKeys['x'][1]-p.filterKeys['x'][0]
    else:
        width = p.mdh.Camera.ROIWidth*p.mdh.voxelsize_nm.x
    if 'y' in p.filterKeys:
        height = p.filterKeys['y'][1]-p.filterKeys['y'][0]
    else:
        height = p.mdh.Camera.ROIHeight*p.mdh.voxelsize_nm.y
    area1 = 1e-6*width*height

    nEvents = p['x'].size
    if nEvents > 1:
        xrange = [p['x'].min(),p['x'].max()]
        yrange = [p['y'].min(),p['y'].max()]

        width = xrange[1]-xrange[0]
        height = yrange[1]-yrange[0]

        area2 = 1e-6*width*height

        if abs(area1-area2)/area1 > 0.2: # we have a > 20 % difference
            area = area2
        else:
            area = area1
    else:
        area = area1
        
    return area # area in um^2

def evtDensity(p):
    area = getarea(p)
    trange = p['t'].max()-p['t'].min()+1
    nEvents = p['x'].size

    if area > 1e-6:
        dens = nEvents / area
        return (dens, dens / trange * 5e3, trange)
    else:
        return None