# calculate object volumes for given set of points
# objects are characterised by ID
# ID vector must match first dim of points
import numpy as np
from scipy.spatial import ConvexHull

def objectVolumes(points, ids):
    idi = ids.astype('int')
    volumes = np.zeros_like(idi,dtype='float')
    idu = np.unique(idi)
    for i in range(idu.shape[0]):
        objectid = idu[i]
        thisID = (idi == objectid)
        if thisID.sum() > 2:
            hull = ConvexHull(points[thisID,:])
            volumes[thisID] = hull.volume
        else:
            volumes[thisID] = 0

    return volumes
