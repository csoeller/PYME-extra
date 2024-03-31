###############################################
### A class definition for basic MBM processing
### and analysis
###############################################
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from warnings import warn

def interp_bead(tnew, bead, customdict=None, extrapisnan=False):
    ibead = {}
    if customdict is None:
        for i,axis in enumerate(['x','y','z']):
            if extrapisnan:
                ibead[axis] = np.interp(tnew, bead['tim'], 1e9*bead['pos'][:,i], right=np.nan) # everything in nm
            else:
                ibead[axis] = np.interp(tnew, bead['tim'], 1e9*bead['pos'][:,i]) # everything in nm
            #ibead[axis][tnew>bead['tim'].max()] = 0
    else:
        for key,value in customdict.items():
            if extrapisnan:
                ibead[key] = np.interp(tnew, bead['tim'], bead[value], right=np.nan)
            else:
                ibead[key] = np.interp(tnew, bead['tim'], bead[value])
            #ibead[key][tnew>bead['tim'].max()] = 0

    ibead['t'] = tnew
    return ibead

def stdev_bead(bead,samplewindow=9):
    sbead = {}
    for i,axis in enumerate(['x','y','z']):
        sbead["std_%s" % axis] = pd.Series(1e9*bead['pos'][:,i]).rolling(window=samplewindow).std() # everything in nm
    sbead['std'] = np.sqrt(sbead['std_x']**2 + sbead['std_y']**2 + sbead['std_z']**2)
    sbead['tim'] = bead['tim']
    return sbead

def stdev_beads(beads,samplewindow=9):
    sbeads = {}
    for bead in beads:
        sbeads[bead] = stdev_bead(beads[bead],samplewindow=samplewindow)
    return sbeads

def interp_sbeads(sbeads,extrapisnan=False):
    return interp_beads(sbeads,customdict=dict(std='std',std_x='std_x',
                                               std_y='std_y',std_z='std_z'),
                        extrapisnan=extrapisnan)

def interp_beads(beads,customdict=None,extrapisnan=False):
    mint = 1e6
    for bead in beads:
        mincur = beads[bead]['tim'].min()
        if mincur < mint:
            mint = mincur

    maxt = 0
    for bead in beads:
        maxcur = beads[bead]['tim'].max()
        if maxcur > maxt:
            maxt = maxcur

    # here we may need some checks if some bead tracks are a lot shorter than others (does this occur)?
    # this could lead to issues with interpolation unless these go to zero
    # so watch out for cases like that and consider code tweaks if needed
    tnew = np.arange(np.round(mint),np.round(maxt)+1)
    ibeads = {}

    for bead in beads:
        ibeads[bead] = interp_bead(tnew,beads[bead],customdict=customdict,extrapisnan=extrapisnan)

    return ibeads

import pandas as pd
def df_from_interp_beads(beads,customdict=None):
    ibeads = interp_beads(beads,customdict=customdict,extrapisnan=True)
    dictbeads = {}
    for axis in ['x','y','z','std_x','std_y','std_z','std']:
        dictbeads[axis] = {}
        
    for bead in ibeads:
        for axis in ['x','y','z']:
            dictbeads[axis][bead] = ibeads[bead][axis]
        t = ibeads[bead]['t'] # this is actually always the same t

    dfbeads = {}
    for axis in ['x','y','z']:
        dfbeads[axis] = pd.DataFrame(dictbeads[axis],index=t)

    sbeads = stdev_beads(beads)
    sibeads = interp_sbeads(sbeads,extrapisnan=True)
    for bead in sibeads:
        for axis in ['std_x','std_y','std_z','std']:
            dictbeads[axis][bead] = sibeads[bead][axis]
    for axis in ['std_x','std_y','std_z','std']:
        dfbeads[axis] = pd.DataFrame(dictbeads[axis],index=t)
    
    return dfbeads

def get_mbm(ds):
    mbm = {}
    mbm['t'] = 1e-3*ds['t']
    mbm['x'] = ds['x']-ds['x_nc']
    mbm['y'] = ds['y']-ds['y_nc']
    if 'z_nc' in ds.keys():
        mbm['z'] = ds['z']-ds['z_nc']
    return mbm

# minimal recipe to coalesce events
COALESCE_RECIPE = """
- localisations.AddPipelineDerivedVars:
    inputEvents: ''
    inputFitResults: FitResults
    outputLocalizations: Localizations
- localisations.MergeClumps:
    discardTrivial: true
    inputName: Localizations
    outputName: coalesced_nz
"""

import hashlib
import json
# we use this function to generate a unique hash from a dictionary
# see also https://stackoverflow.com/questions/16092594/how-to-create-a-unique-key-for-a-dictionary-in-python
# this will be used further below to check if our cached value of the mean is still usable
def hashdict(dict):
    hashkey = hashlib.sha1(json.dumps(dict, sort_keys=True).encode()).hexdigest()
    return hashkey
    
class MBMCollection(object):
    def __init__(self,name=None,filename=None,variance_window = 9):
        self.mbms = {}
        self.beadisgood = {}
        self.offsets = {}
        self.is3D = False
        self._mean = None
        self._hashkey = ''
        self._offsets_valid = False
        self.t = None
        self.tperiod = None
        self._trange= (None,None)
        self.variance_window = variance_window # by default use last 9 localisations for variance/std calculation
        
        if filename is not None:
            # this is a MBM bead file with raw bead tracks
            self.name=filename
            self._raw_beads = np.load(filename)
            ibeads = interp_beads(self._raw_beads)
            self.add_beads(ibeads)
            # now add info on std deviation
            sbeads = stdev_beads(self._raw_beads)
            sibeads = interp_sbeads(sbeads)
            for bead in sibeads:
                if not np.allclose(self.t,sibeads[bead]['t'],1e-3):
                    raise RuntimeError("time vector for new bead variance differs by more than 0.1 %")
                for property in sibeads[bead].keys():
                    if property != 't':
                        self.mbms[bead][property] = sibeads[bead][property]
        else:
            self.name = name

    @property    
    def beads(self):
        return self.mbms.keys()

    def _validaxis(self,axis):
        if self.is3D:
            axes = ['x','y','z','std_x','std_y','std_z','std']
        else:
            axes = ['x','y','std_x','std_y','std']
        return axis in axes
    
    def beadtrack(self,bead,axis,unaligned=False):
        if not bead in self.beads:
            raise RuntimeError("asking for non existing bead track for bead %s" % bead)
        if not self._validaxis(axis):
            raise RuntimeError("asking for invalid axis %s" % axis)
        if self._offsets_valid and not unaligned:
            return self.mbms[bead][axis] - self.offsets[bead][axis]
        else:
            return self.mbms[bead][axis]

    def mean(self):
        if self._mean is not None and self._hashkey == hashdict(self.beadisgood):
            # print("cache hit!")
            return self._mean
        else:
            self._mean = {}
            for axis in ['x','y']:
                self._mean[axis] = np.mean([self.beadtrack(bead,axis) for bead in self.beads if self.beadisgood[bead]], axis=0)
            if self.is3D:
                self._mean['z'] = np.mean([self.beadtrack(bead,'z') for bead in self.beads if self.beadisgood[bead]], axis=0)
            self._hashkey = hashdict(self.beadisgood)
            return self._mean

    def add_bead(self,bead,mbm):
        if self.t is None:
            self.t = mbm['t']
            self.is3D = 'z' in mbm.keys()
        else:
            if self.t.size != mbm['t'].size:
                raise RuntimeError("register bead: size of time vectors do not match\nold size %s, new size %s, bead %s"
                                   % (self.t.size,mbm['t'].size,bead) )
            if not np.allclose(self.t,mbm['t'],1e-3):
                raise RuntimeError("time vector for new bead differs by more than 0.1 %")
            if not 'z' in mbm.keys() and self.is3D:
                raise RuntimeError('adding bead lacking z info to existing 3D MBM collection')
        self.mbms[bead] = mbm # note we may need copies of vectors, possibly at leas
        self.markasgood(bead)
        self._mean = None # invalidate cache
        self._offsets_valid = False

    def add_beads(self,beads):
        for bead in beads:
            self.add_bead(bead,beads[bead])

    def align_beads(self,tmin=None,tmax=None):
        if tmin is None:
            if self._trange[0] is None:
                tmin = self.t.min()
            else:
                tmin=self._trange[0]
        if tmax is None:
            if self._trange[1] is None:
                tmax = self.t.max()
            else:
                tmax=self._trange[0]
        self.tperiod = (self.t > tmin)*(self.t < tmax)
        if np.all(self.tperiod == False):
            raise RuntimeError("empty range, tmin: %d, tmax: %d" % (tmin,tmax))
        self._trange = (tmin,tmax)
        for bead in self.beads:
            self.offsets[bead] = {}            
            self.offsets[bead]['x'] = np.mean(self.mbms[bead]['x'][self.tperiod])
            self.offsets[bead]['y'] = np.mean(self.mbms[bead]['y'][self.tperiod])
            if self.is3D:
                self.offsets[bead]['z'] = np.mean(self.mbms[bead]['z'][self.tperiod])
        self._offsets_valid = True
        self._mean = None # invalidate cache
    
    def markasbad(self,*beads): # mark a bead as bad
        for bead in beads:
            if bead in self.beads:
                self.beadisgood[bead] = False

    def markasgood(self,*beads): # if currently bad, mark as good
        for bead in beads:
            if bead in self.beads:
                self.beadisgood[bead] = True

    def plot_tracks(self,axis,unaligned=False,use_tperiod=False,legend=True,tmin=None,tmax=None,plot_mean=True):
        if tmin is None:
            tmin=self._trange[0]
        if tmax is None:
            tmax=self._trange[1]
        if tmin != self._trange[0] or tmax != self._trange[1]:
            self._offsets_valid = False

        if axis.startswith('std'):
            unaligned = True # not sensible to align the std devs
            plot_mean = False # the mean also does not make much sense
            
        # we may need to check the alignment logic below
        if not unaligned:
            if not self._offsets_valid:
                self.align_beads(tmin=tmin,tmax=tmax)

        for bead in self.beads:
            if self.beadisgood[bead]:
                plt.plot(self.t,self.beadtrack(bead,axis,unaligned=unaligned),label=bead)
        if plot_mean:
            plt.plot(self.t,self.mean()[axis],'--',label='mean')
        if legend:
            plt.legend()
        if use_tperiod:
            plt.xlim(self._trange[0],self._trange[1])
        
    def plot_deviation_from_mean(self,axis,align=True,legend=True):
        if axis.startswith('std'):
            raise RuntimeError("this method is not suitable for standard deviation trajectories")
        if align and not self._offsets_valid:
            self.align_beads()
        for bead in self.beads:
            if self.beadisgood[bead]:
                plt.plot(self.t,self.beadtrack(bead,axis)-self.mean()[axis],label=bead)        
        if legend:
            plt.legend()

try:
    import plotly.express as px
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
except ImportError:
    warn("can't import plotly modules, new style bead plotting using MBMCollectionDF will not work")

# we use this function to generate a unique hash from a dataframe
# need to check if this is necessary or if it is ok to make the the filter settimngs into a unique hash for caching
def hashdf(df):
    return hashlib.sha1(pd.util.hash_pandas_object(df).values).hexdigest()

class MBMCollectionDF(object): # collection based on dataframe objects
    def __init__(self,name=None,filename=None,variance_window = 9):
        self.mbms = {}
        self.beadisgood = {}
        self.t = None
        self.tperiod = None
        self._trange= (None,None)
        self.variance_window = variance_window # by default use last 9 localisations for variance/std calculation
        self.median_window = 0 # 0 means not active
        
        if filename is not None:
            self.populate_df_from_npz(filename)
        else:
            self.name = name

    def populate_df_from_npz(self,filename):
        # this is a MBM bead file with raw bead tracks
        self.name=filename
        self._raw_beads = np.load(filename)
        self.beads = df_from_interp_beads(self._raw_beads)
        self.t = self.beads['x'].index
        
        for bead in self.beads['x']:
            self.beadisgood[bead] = True

    def markasbad(self,*beads): # mark a bead as bad
        for bead in beads:
            if bead in self.beads['x']:
                self.beadisgood[bead] = False

    def markasgood(self,*beads): # if currently bad, mark as good
        for bead in beads:
            if bead in self.beads['x']:
                self.beadisgood[bead] = True

    def plot_tracks(self,axis,unaligned=False,tmin=None,tmax=None):
        if tmin is None:
            tmin=self._trange[0]
        if tmax is None:
            tmax=self._trange[1]

        if tmin is None:
            tmin = self.t.min()
        if tmax is None:
            tmax = self.t.max()

        if axis.startswith('std'):
            unaligned = True # not sensible to align the std devs

        if self.median_window > 0:
            startdf = self.beads[axis].rolling(self.median_window).median()
        else:
            startdf = self.beads[axis]
        if not unaligned:
            dfplot = startdf-startdf.loc[tmin:tmax].mean(axis=0)
            dfplotg = dfplot[[bead for bead in self.beadisgood if self.beadisgood[bead]]]
            emptybeads = dfplotg.columns[dfplotg.isnull().all(axis=0)]
            if len(emptybeads)>0:
                warn('removing beads with no valid info after alignment %s...' % emptybeads)
                dfplotg = dfplotg[dfplotg.columns[~dfplotg.isnull().all(axis=0)]]
                
            fig1 = px.line(dfplotg)
            fig1.add_trace(go.Scatter(x=self.t, y=dfplotg.mean(axis=1), name='Mean',
                                     line=dict(color='firebrick', dash='dash')))
            fig2 = px.line(dfplotg.sub(dfplotg.mean(axis=1),axis=0))

            fig = make_subplots(rows=2, cols=1)

            for d in fig1.data:
                fig.add_trace((go.Scatter(x=d['x'], y=d['y'], name = d['name'])), row=1, col=1)
            for d in fig2.data:
                fig.add_trace((go.Scatter(x=d['x'], y=d['y'],  name = d['name'])), row=2, col=1)

            fig.update_layout(autosize=False, width=1000, height=700,title_text="aligned MBM tracks along %s" % axis)
            # Update axes properties
            fig.update_xaxes(title_text="time (s)", row=1, col=1)
            fig.update_xaxes(title_text="time (s)", row=2, col=1)
            fig.update_yaxes(title_text="drift (nm)", range=[np.min([-15.0,dfplotg.min().min()]),np.max([15.0,dfplotg.max().max()])], row=1, col=1)
            fig.update_yaxes(title_text="deviation (nm)", range=[-10,10], row=2, col=1)
            
            fig.show()

        else:
            if axis.startswith('std'):
                yaxis_title = "std dev (nm)"
                title = 'MBM localisation precisions (%s)' % axis
            else:
                title = 'tracks along %s, not aligned' % axis
                yaxis_title = "distance (nm)"
            dfplot = startdf
            dfplotg = dfplot[[bead for bead in self.beadisgood if self.beadisgood[bead]]]
            fig = px.line(dfplotg)
            fig.update_layout(xaxis_title="time (s)", yaxis_title=yaxis_title, title_text=title)
            if axis.startswith('std'):
                fig.update_yaxes(range = (0,np.max([10.0,dfplotg.max().max()])))
            fig.show()
        
