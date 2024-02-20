###############################################
### A class definition for basic MBM processing
### and analysis
###############################################
import numpy as np
import matplotlib.pyplot as plt

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
    
class mbmcollection(object):
    def __init__(self,name):
        self.name = name
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

    @property    
    def beads(self):
        return self.mbms.keys()

    def _validaxis(self,axis):
        if self.is3D:
            axes = ['x','y','z']
        else:
            axes = ['x','y']
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
        if align and not self._offsets_valid:
            self.align_beads()
        for bead in self.beads:
            if self.beadisgood[bead]:
                plt.plot(self.t,self.beadtrack(bead,axis)-self.mean()[axis],label=bead)        
        if legend:
            plt.legend()
