# code for coordinate algebra, specifically rotations given by the `rotation` transform

import numpy as np
import matplotlib.pyplot as plt

def vlen(v):
    return np.sqrt(v[0]**2+v[1]**2)

def rmat_ang(angle):
    theta = np.radians(-angle)
    # note: tests suggest that for the svg angle we use the *negative* angle for this matrix
    c, s = np.cos(theta), np.sin(theta)
    R = np.array(((c, -s), (s, c)))
    return R

def vec(x,y):
    return np.array((x,y))

def transform_rot(x,y,angle):
    return vec(x,y).dot(rmat_ang(angle))

# the code below parses out the rectangles from an SVG
# all other elements are ignored
# for the identified rectangles the code obtains the revelant parameters by retrieveal
# and string conversion to float

import re

# NOTE: the argument 'svg' is assumed to be a parsed svg file
# obtained with something like
#
# import xml.etree.ElementTree as ET
# svg = ET.parse('RyR-pattern-simulation-Ashagri14-Fig1.svg')
def find_svgrects(svg):
    root = svg.getroot()
    group = root.find('{http://www.w3.org/2000/svg}g')
    rects = group.findall('{http://www.w3.org/2000/svg}rect')
    return rects

def parse_rect(rect):
    # ['style', 'id', 'width', 'height', 'x', 'y', 'transform']
    prect = {}
    for key in ['width', 'height', 'x', 'y']:
        prect[key] = float(rect.get(key))
    for key in ['style', 'id']:
        prect[key] = rect.get(key)
    if rect.get('transform') is not None:
        m = re.search('^rotate\(([-+]*[\d.]+)\)',rect.get('transform'))
        if m:
            prect['rot_angle'] = float(m.group(1))
    return prect
    
def parse_rects(rects):
    return [parse_rect(rect) for rect in rects]

# code to display parsed rectangles via matplotlib
# we need this to check that things work as anticipated

from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle

def transformed_rect(prect):
    angle = prect.get('rot_angle',0)
    r_origin = list(transform_rot(prect['x'],prect['y'],angle))
    return Rectangle(r_origin, prect['width'], prect['height'], angle=angle)   

def make_matplotrects(prects):
    mprects = [transformed_rect(prect) for prect in prects]
    return mprects

def plot_rects(ax,matplotrects,alpha=1.0):
    # Create patch collection with specified colour/alpha
    pc = PatchCollection(matplotrects, facecolor='r', alpha=alpha,
                         edgecolor='black',linewidth=1.0)

    # Add collection to axes
    ax.add_collection(pc,autolim=True)
    
    # try to fix axes limits
    # this does not seem to work as unticipated
    # fine to ignore for now (because it does nothing)...
    ax._unstale_viewLim()
    datalim = pc.get_datalim(ax.transData)
    print(datalim.get_points())
    points = datalim.get_points()
    if not np.isinf(datalim.minpos).all():
        # By definition, if minpos (minimum positive value) is set
        # (i.e., non-inf), then min(points) <= minpos <= max(points),
        # and minpos would be superfluous. However, we add minpos to
        # the call so that self.dataLim will update its own minpos.
        # This ensures that log scales see the correct minimum.
        points = np.concatenate([points, [datalim.minpos]])
    ax.update_datalim(points)
    return points

# function to scale everything to be accurate in nm
# below we assume that rectangles should be 27 nm x 27 nm
def scale(prect,factor):
    srect = prect.copy()
    for key in ['width', 'height', 'x', 'y']:
        srect[key] = prect[key]*factor
    return srect

# code to get label points from scaled prect
def points_from_sprect(sprect,mode='default'):
    angle = sprect.get('rot_angle',0)
    origin = transform_rot(sprect['x'],sprect['y'],angle)
    if mode == 'default':
        corner_points = np.array([[3.5,3.5], [23.5,3.5], [23.5, 23.5], [3.5, 23.5]]) # anti-clock wise
    elif mode == 'RyR-T1366':
        corner_points = np.array([[6,5], [21.6,6.6], [20.2,22.5], [4.8,20.5]])
    elif mode == 'RyR-T2023':
        corner_points = np.array([[4.6,12.0], [15.5,4.8], [22.5,15.2], [11.9,22.5]])
    elif mode == 'RyR-D4365':
        corner_points = np.array([[2,13.2],[13.9,2.2],[25.3,14.0],[13.3,25.2]])

    else:
        raise RuntimeError('unknow RyR corner mode %s requested' % (mode))
    # first rotate corner points
    R = rmat_ang(sprect.get('rot_angle',0))
    cpsr = corner_points[:,None,:].dot(R).squeeze()
    # now shift corner points to the origin of this rect
    scpsr = cpsr + origin[None,:]
    return scpsr

# functions to pick points and rotate point sets for generating clusters etc

from numpy.random import default_rng
rng = default_rng()

def random_pick(size,p):
    return rng.random(size) < p

def rotpts(pts,angle):
    rotpts = []
    for ptnum in range(pts.shape[0]):
        rx,ry = transform_rot(pts[ptnum,0],pts[ptnum,1],angle)
        rotpts.append([rx,ry])
    return np.array(rotpts)

def ctrorigin(pts):
    mxy = pts.mean(axis=0)
    return pts - mxy[None,:]

def shftpts(pts,x,y):
    return pts + np.array([x,y])[None,:]

def pickpts(pts,prob):
    picked = random_pick(pts.shape[0],prob)

    return pts[picked,:]

def mergeClusterlist(clist):
    clusters = clist[0].copy()
    for c in clist[1:]:
        clusters = np.append(clusters,c,axis = 0)
    return clusters

import pandas as pd
def pts2df(pts):
    pymedct =  {'x' : pts[:,0], 'y' : pts[:,1], 'z' : np.zeros((pts.shape[0])), 't' : np.ones((pts.shape[0]))}
    pymedf = pd.DataFrame.from_dict(pymedct)
    return pymedf

def pymedf_add_err(pymedf,locerr=2.0): # error in nm
    errdf = pymedf.copy()
    nrows = errdf.shape[0]
    errdf['x'] += locerr*np.random.normal(size=nrows)
    errdf['y'] += locerr*np.random.normal(size=nrows)
    errdf['z'] += locerr*np.random.normal(size=nrows)
    errdf['error_x'] = locerr*np.ones((nrows))
    errdf['error_y'] = locerr*np.ones((nrows))
    errdf['error_z'] = locerr*np.ones((nrows))

    return errdf

def datasourcePick2df(ds,prob):
    picked = random_pick(ds['x'].size,prob)
    pymedct = {}
    for key in ['x','y','z','error_x','error_y','error_z','t']:
        pymedct[key] = ds[key][picked]
    return pd.DataFrame.from_dict(pymedct)
