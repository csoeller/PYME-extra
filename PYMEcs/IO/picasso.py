import numpy as np
import pandas as pd
import yaml
import math
import os

def picasso2pyme(locs,mdh):
    try:
        nmperpix = mdh['Camera.Pixelsize']
    except KeyError:
        nmperpix = mdh['Pixelsize']
        
    pyme = {}
    pyme['x'] = locs['x']*nmperpix
    pyme['y'] = locs['y']*nmperpix
    pyme['sig'] = locs['sx']*nmperpix
    pyme['nPhotons'] = locs['photons']
    pyme['error_x'] = locs['lpx']*nmperpix
    pyme['error_y'] = locs['lpy']*nmperpix
    pyme['t'] = locs['frame']
    pyme['A'] = locs['photons']/(2*math.pi)/(locs['sx']**2) # this is from the PYME kinmodels formula - check
    if 'z' in locs.keys():
        pyme['z'] = locs['z'] # from astigmatoc data z should already be scaled to nm?
    return pyme

def struc_from_mdh(mdh):
    if mdh is None or 'Structure.HandleStruct' not in mdh.keys():
        return None
    struc = {}
    nmperpix = float(mdh['Camera.Pixelsize'])
    npixels = float(mdh['Camera.Image Size'])
    for key in ('HandleStruct','HandleEx'):
        struc[key] = np.array([float(i) for i in mdh["Structure.%s" % key].split(',')],dtype='i')
    offshalfpix = 0.5 # half a pixel offset seems to be required (origin in middle of pixel?)
    for key in ['HandleX']:
        struc[key] = np.array([(float(i)-offshalfpix)*nmperpix for i in mdh["Structure.%s" % key].split(',')])
    # for a reason not quite clear to me yet the Y coordinates need mirroring
    for key in ['HandleY']:
        struc[key] = np.array([(npixels-float(i)-offshalfpix)*nmperpix for i in mdh["Structure.%s" % key].split(',')])
    return struc

def parse_dicts(dicts):
    import re
    mdh = None
    strucmdh = None
    for d in dicts:
        if any (re.match('pixelsize',key,re.I) for key in d.keys()):
            mdh = d
    for d in dicts:
        if any (re.match('Structure',key,re.I) for key in d.keys()):
            strucmdh = d    
    return (mdh,strucmdh)

def read_picasso_hdf(name):
    locs = pd.read_hdf(name,'locs')
    basename, file_extension = os.path.splitext(name)
    with open(basename+'.yaml') as file:
        finf = list(yaml.load_all(file, Loader=yaml.FullLoader))
    mdh, strucmdh = parse_dicts(finf)
    pymedf = pd.DataFrame(picasso2pyme(locs,mdh))
    struc = struc_from_mdh(strucmdh)
    return pymedf, mdh, struc

def pymedf2csv(pymedf,filename):
    cols = ['A','x','y','t','sig','error_x','error_y','nPhotons']
    pymedf.to_csv(filename,columns=cols,index=False)

def pymedf2ds(pymedf):
    from PYME.IO.tabular import DictSource
    pymedict={}
    for key in pymedf.keys():
        pymedict[key] = pymedf[key].values
    ds = DictSource(pymedict)
    return ds

import skimage.morphology as morph
from PYME.IO.image import ImageStack
#from PYME.IO.MetaDataHandler import NestedClassMDHandler

def picasso_structure_mask(inputim,struc,dilations=2,
                           dilationselem=np.ones((5,5),dtype='float32')):
    ox = inputim.mdh['Origin.x']
    oy = inputim.mdh['Origin.y']
    vxsz = 1e3*inputim.mdh['voxelsize.x']

    mask0 = np.zeros(inputim.data.shape[0:2],dtype='float32')
    labels = []
    nsites = []
    cxa = []
    cya = []
    
    for label in range(0,int(struc['HandleStruct'].max()+1)):
        strucn = struc['HandleStruct'] == label
        newim = np.zeros(inputim.data.shape[0:2],dtype='float32')
        strucx = ((struc['HandleX'][strucn]) - ox) / vxsz
        strucy = ((struc['HandleY'][strucn]) - oy) / vxsz
        cx = struc['HandleX'][strucn].mean()
        cy = struc['HandleY'][strucn].mean()
        cxi = np.rint(strucx.mean()).astype('int') # integer centroid x
        cyi = np.rint(strucy.mean()).astype('int') # integer centroid y
        labels.append(label+1)
        nsites.append(strucx.size)
        cxa.append(cx)
        cya.append(cy)
        ind =  (strucx < newim.shape[0])*(strucy < newim.shape[1])*(strucx >= 0)*(strucy >= 0)
        if np.any(ind):
            newim[strucx[ind].astype('i'),strucy[ind].astype('i')] = 1.0
            newim2 = morph.convex_hull_image(newim)
            for i in range(dilations):
                newim2 = morph.binary_dilation(newim2,selem=dilationselem)
            mask0[newim2 > 0.5] = label+1

    sitesdf = pd.DataFrame(list(zip(labels, nsites, cxa, cya)), 
                           columns =['Label', 'NSites', 'CentroidX', 'CentroidY'])
    return ImageStack(mask0, mdh=inputim.mdh), sitesdf

# same as the routine above but everything is done in a small ROI that is then inserted into the
# full size image; this speeds everything up considerably, mostly as hoped
def picasso_structure_mask_roi(inputim,struc,dilations=2,
                               dilationselem=np.ones((5,5),dtype='float32'),roisize=20):
    ox = inputim.mdh['Origin.x']
    oy = inputim.mdh['Origin.y']
    vxsz = 1e3*inputim.mdh['voxelsize.x']

    mask0 = np.zeros(inputim.data.shape[0:2],dtype='float32')
    labels = []
    nsites = []
    cxa = []
    cya = []
    
    rszh = roisize
    xroi = np.arange(0,2*rszh+1,dtype='int')
    yroi = np.arange(0,2*rszh+1,dtype='int')
    xroi2 = np.outer(xroi,np.ones(2*rszh+1,dtype='i'))
    yroi2 = np.outer(np.ones(2*rszh+1,dtype='i'),yroi)

    for label in range(0,int(struc['HandleStruct'].max()+1)):
        strucn = struc['HandleStruct'] == label
        strucx = ((struc['HandleX'][strucn]) - ox) / vxsz
        strucy = ((struc['HandleY'][strucn]) - oy) / vxsz
        cx = struc['HandleX'][strucn].mean()
        cy = struc['HandleY'][strucn].mean()
        cxi = np.rint(strucx.mean()).astype('int') # integer centroid x
        cyi = np.rint(strucy.mean()).astype('int') # integer centroid y
        labels.append(label+1)
        nsites.append(strucx.size)
        cxa.append(cx)
        cya.append(cy)
        # make empty small ROI
        roinew = np.zeros((2*rszh+1,2*rszh+1),dtype='f')
        # transform coordinates to the ROI coordinate system
        # the ROI is centered on the structure centroid (cxi,cyi)
        roisx = strucx - (cxi-rszh)
        roisy = strucy - (cyi-rszh)
        # set pixels corresponding to structure to 1
        roinew[roisx.astype('i'),roisy.astype('i')] = 1.0
        # convex hull of structure + dilations
        roi2 = morph.convex_hull_image(roinew)
        for i in range(dilations):
            roi2 = morph.binary_dilation(roi2,selem=dilationselem)
        # transform ROI coordinates to full image coordiante system 
        x2 = xroi2 + (cxi-rszh)
        y2 = yroi2 + (cyi-rszh)
        # restrict to non-zero pixels
        mask = roi2 > 0.5
        x2m = x2[mask]
        y2m = y2[mask]
        # check which of these non-zero pixels are within the full image bounds
        valid = (x2m < mask0.shape[0])*(y2m < mask0.shape[1])*(x2m >= 0)*(y2m >= 0)
        if np.any(valid):
            mask0[x2[mask][valid],y2[mask][valid]] = label+1 # set to label of that structure

    sitesdf = pd.DataFrame(list(zip(labels, nsites, cxa, cya)),
                           columns =['Label', 'NSites', 'CentroidX', 'CentroidY'])
    return ImageStack(mask0, mdh=inputim.mdh), sitesdf

# import skimage.morphology as morph

# dilations=2
# dilationselem=np.ones((5,5),dtype='i')

# strucx = np.array([200,210,220,200,210])
# strucy = np.array([100,100,100,110,110])

# cx = np.round(strucx.mean()).astype('i')
# cy = np.round(strucy.mean()).astype('i')

# mask0 = np.zeros((256,256),dtype='f')

# rszh = 20

# # initialisation
# xroi = np.arange(0,2*rszh+1,dtype='int')
# yroi = np.arange(0,2*rszh+1,dtype='int')
# xroi2 = np.outer(xroi,np.ones(2*rszh+1,dtype='i'))
# yroi2 = np.outer(np.ones(2*rszh+1,dtype='i'),yroi)

# # this is done in each loop
# roinew = np.zeros((2*rszh+1,2*rszh+1),dtype='f')

# roisx = strucx - (cx-rszh)
# roisy = strucy - (cy-rszh)

# roinew[roisx.astype('i'),roisy.astype('i')] = 1.0
# roi2 = morph.convex_hull_image(roinew)
# for i in range(dilations):
#     roi2 = morph.binary_dilation(roi2,selem=dilationselem)

# x2 = xroi2 + (cx-rszh)
# y2 = yroi2 + (cy-rszh)

# mask = roi2 > 0.5

# x2m = x2[mask]
# y2m = y2[mask]

# valid = (x2m < mask0.shape[0])*(y2m < mask0.shape[1])*(x2m >= 0)*(y2m >= 0)

# mask0[x2[mask][valid],y2[mask][valid]] = roi2[xroi2[mask][valid],yroi2[mask][valid]]
