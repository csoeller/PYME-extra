import numpy as np

# this is code to obtain dark times from tabular event columns
#
# it works but likely needs a code overhaul
# this includes both the fitting but also
# how the input pipeline is handled - conceivably it
# should just received 1D vectors and be completely
# ignorant of any pipeline/datasource origin

import PYME.IO.tabular as tabular
# quick tabular class that wraps a recarray
# and allows adding new columns
# and inherits tabular I/O
class TabularRecArrayWrap(tabular.TabularBase):
    def __init__(self, recarray):
        self._recarray = recarray
        self.new_columns = {}
        
    def keys(self):
        return list(set(list(self._recarray.dtype.fields.keys() + self.new_columns.keys())))

    def __getitem__(self, keys):
        key, sl = self._getKeySlice(keys)
        if key in self._recarray.dtype.fields.keys():
            return self._recarray[key][sl]
        else:
            return self.new_columns[key][sl]

    def addColumn(self, name, values):
        """
        Adds a column of values to the tabular measure.

        Parameters
        ----------
        name : str
            The new column name
        values : array-like
            The values. This should be the same length as the existing columns.

        """

        #force to be an array
        values = np.array(values)

        if not len(values) == self._recarray.shape[0]:
            raise RuntimeError('New column does not match the length of existing columns')

        #insert into our __dict__ object (for backwards compatibility - TODO change me to something less hacky)
        #setattr(self, name, values)

        self.new_columns[name] = values

    def getZeroColumn(self, dtype='float64'):
        return np.zeros(self._recarray.shape, dtype=dtype)

    def addNewColumnByID(self, fromids, colname, valsByID):
        if not np.all(np.in1d(fromids,self['objectID'])):
            raise RuntimeError('some ids not present in measurements')
        # this expression finds the lookup table to locate
        # ids in fromids in self['objectID']
        # i.e. we should have self['objectID'][idx] == fromids
        idx = np.nonzero((fromids[None,:] == self['objectID'][:,None]).T)[1]
        if not np.all(self['objectID'][idx] == fromids):
            raise RuntimeError('Lookup error - this should not happen')
    
        newcol = self.getZeroColumn(dtype='float64')
        newcol[idx] = valsByID
        self.addColumn(colname,newcol)


def mergeChannelMeasurements(channels,measures):
    master = measures[0]['objectID']
    for chan,meas in zip(channels,measures):
        if meas['objectID'].size != master.size:
            raise RuntimeError('channel %s does not have same size as channel %s' % (chan, channels[0]))
        if not np.all(meas['objectID'] == master):
            raise RuntimeError('channel %s object IDs do not match channel %s object IDs' % (chan, channels[0]))
    mergedIDs = np.zeros(master.size, dtype=[('objectID','i4')])
    mergedIDs[:] = master
    mergedMeas = TabularRecArrayWrap(mergedIDs)

    for chan,meas in zip(channels,measures):
        for key in meas.keys():
            if key != 'objectID':
                mergedMeas.addColumn('%s_%s' % (key,chan), meas[key])

    return mergedMeas

from scipy.optimize import curve_fit
def cumuexpfit(t,tau):
    return 1-np.exp(-t/tau)

def notimes(ndarktimes):
    analysis = {
        'NDarktimes' : ndarktimes,
        'tau1' : [None,None,None,None],
        'tau2' : [None,None,None,None]
    }
    return analysis

def cumuhist(timeintervals):
    ti = timeintervals
    nIntervals = ti.shape[0]
    cumux = np.sort(ti+0.01*np.random.random(nIntervals)) # hack: adding random noise helps us ensure uniqueness of x values
    cumuy = (1.0+np.arange(nIntervals))/np.float(nIntervals)
    return (cumux,cumuy)
    
def cumuhistBinned(timeintervals):
    binedges = 0.5+np.arange(0,timeintervals.max())
    binctrs = 0.5*(binedges[0:-1]+binedges[1:])
    h,be2 = np.histogram(timeintervals,bins=binedges)
    hc = np.cumsum(h)
    hcg = hc[h>0]/float(timeintervals.shape[0]) # only nonzero bins and normalise
    binctrsg = binctrs[h>0]

    return (binctrs, hc, binctrsg, hcg)

def fitDarktimes(t):
    # determine darktime from gaps and reject zeros (no real gaps) 
    nts = 0 # initialise to safe default
    NTMIN = 5
    
    if t.size > NTMIN:
        dts = t[1:]-t[0:-1]-1
        dtg = dts[dts>0]
        nts = dtg.shape[0]

    if nts > NTMIN:
        # now make a cumulative histogram from these
        cumux, cumuy = cumuhist(dtg)
        tauEst = cumux[(np.abs(cumuy - 0.63)).argmin()]
        # generate alternative histogram with binning
        binctrs, hc, binctrsg, hcg = cumuhistBinned(dtg)
        tauEstH = binctrsg[(np.abs(hcg - 0.63)).argmin()]
        
        success = True
        # fit theoretical distributions
        try:
            popth,pcovh,infodicth,errmsgh,ierrh  = curve_fit(cumuexpfit,binctrsg,hcg, p0=(tauEstH),full_output=True)
        except:
            success = False
        else:
            chisqredh = ((hcg - infodicth['fvec'])**2).sum()/(hcg.shape[0]-1)
        try:
            popt,pcov,infodict,errmsg,ierr = curve_fit(cumuexpfit,cumux,cumuy, p0=(tauEst),full_output=True)
        except:
            success = False
        else:
            chisqred = ((cumuy - infodict['fvec'])**2).sum()/(nts-1)

        if success:
            analysis = {
                'NDarktimes' : nts,
                'tau1' : [popt[0],np.sqrt(pcov[0][0]),chisqred,tauEst], # cumuhist based
                'tau2' : [popth[0],np.sqrt(pcovh[0][0]),chisqredh,tauEstH] # cumuhistBinned based
            }
        else:
            analysis = notimes(nts)
    else:
        analysis = notimes(nts)

    return analysis

measureDType = [('objectID', 'i4'), ('t', 'i4'), ('x', 'f4'), ('y', 'f4'),
                ('NEvents', 'i4'), ('NDarktimes', 'i4'), ('tau1', 'f4'),
                ('tau2', 'f4'), ('tau1err', 'f4'), ('tau2err', 'f4'),
                ('chisqr1', 'f4'), ('chisqr2', 'f4'), ('tau1est', 'f4'), ('tau2est', 'f4')]


def measure(object, measurements = np.zeros(1, dtype=measureDType)):
    #measurements = {}

    measurements['NEvents'] = object['t'].shape[0]
    measurements['t'] = np.median(object['t'])
    measurements['x'] = object['x'].mean()
    measurements['y'] = object['y'].mean()

    t = object['t']

    darkanalysis = fitDarktimes(t)
    measurements['tau1'] = darkanalysis['tau1'][0]
    measurements['tau2'] = darkanalysis['tau2'][0]
    measurements['tau1err'] = darkanalysis['tau1'][1]
    measurements['tau2err'] = darkanalysis['tau2'][1]
    measurements['chisqr1'] = darkanalysis['tau1'][2]
    measurements['chisqr2'] = darkanalysis['tau2'][2]
    measurements['tau1est'] = darkanalysis['tau1'][3]
    measurements['tau2est'] = darkanalysis['tau2'][3]
    
    measurements['NDarktimes'] = darkanalysis['NDarktimes']
    
    return measurements

def measureObjectsByID(filter, ids):
    x = filter['x'] #+ 0.1*random.randn(filter['x'].size)
    y = filter['y'] #+ 0.1*random.randn(x.size)
    id = filter['objectID'].astype('i')
    t = filter['t']

    measurements = np.zeros(len(ids), dtype=measureDType)

    for j,i in enumerate(ids):
        if not i == 0:
            ind = id == i
            obj = {'x': x[ind], 'y': y[ind], 't': t[ind]}
            #print obj.shape
            measure(obj, measurements[j])
            measurements[j]['objectID'] = i

    # wrap recarray in tabular that allows us to
    # easily add columns and save using tabular methods
    return TabularRecArrayWrap(measurements)

def retrieveMeasuresForIDs(measurements,idcolumn,columns=['tau1','NDarktimes','qindex1']):
    newcols = {key: np.zeros_like(idcolumn, dtype = 'float64') for key in columns}

    for j,i in enumerate(measurements['objectID']):
        if not i == 0:
            ind = idcolumn == i
            for col in newcols.keys():
                if col.startswith('tau'):
                    if not np.isnan(measurements[col][j]) and measurements[col][j] is not None:
                        newcols[col][ind] = measurements[col][j]
                if col.startswith('qindex'):
                    source = col.replace('qindex','tau')
                    if measurements[source][j] > 0:
                        newcols[col][ind] = 100.0/measurements[source][j]
                else:
                    newcols[col][ind] = measurements[col][j]

    return newcols
