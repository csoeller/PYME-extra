import numpy as np

# this is code to obtain dark times from tabular event columns
#
# it works but likely needs a code overhaul
# this includes both the fitting but also
# how the input pipeline is handled - conceivably it
# should just received 1D vectors and be completely
# ignorant of any pipeline/datasource origin

import logging
logger = logging.getLogger(__name__)

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
            logger.warn('some ids not present in measurements')
        # limit everything below to those IDs present in the events
        fromids1 = fromids[np.in1d(fromids,self['objectID'])]
        # this expression finds the lookup table to locate
        # ids in fromids in self['objectID']
        # i.e. we should have self['objectID'][idx] == fromids
        idx = np.nonzero((fromids1[None,:] == self['objectID'][:,None]).T)[1]
        if not np.all(self['objectID'][idx] == fromids1):
            raise RuntimeError('Lookup error - this should not happen')
    
        newcol = self.getZeroColumn(dtype='float64')
        # make sure we match fromids1 shape in assignment
        newcol[idx] = valsByID[np.in1d(fromids,self['objectID'])]
        self.addColumn(colname,newcol)

    def lookupByID(self, ids, key):
        idi = ids.astype('int')
        uids = np.unique(idi[idi.nonzero()])
        uids_avail = uids[np.in1d(uids,self['objectID'])]
        idx = np.nonzero((uids_avail[None,:] == self['objectID'][:,None]).T)[1]
        valsByID = self[key][idx] # these are the values from column 'key' matching uids_avail
        
        lut = np.zeros(uids_avail.max()+1,dtype='float64')
        lut[uids_avail] = valsByID
        
        lu = np.zeros_like(idi,dtype='float64')
        idiflat = idi.ravel()
        luflat = lu.ravel()
        luflat[np.in1d(idiflat,uids_avail)] = lut[idiflat[np.in1d(idiflat,uids_avail)]]
        return lu


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


def safeRatio(mmeas, div11, div22):
    mzeros = mmeas.getZeroColumn(dtype='float')
    div1 = mzeros+div11 # this converts scalars if needed
    div2 = mzeros+div22
    ratio = np.zeros_like(div1)
    d1good = (np.logical_not(np.isnan(div1)))
    d2good = div2 > 0
    allgood = d1good * d2good
    ratio[allgood] = div1[allgood] / div2[allgood]
    return ratio


def makeRatio(meas, key, div1, div2):
    meas.addColumn(key,safeRatio(meas, div1, div2))

    
def makeSum(meas, key, add11, add22):
    mzeros = meas.getZeroColumn(dtype='float')
    add1 = mzeros+add11
    add2 = mzeros+add22
    msum = np.zeros_like(add1)
    a1good = (np.logical_not(np.isnan(add1)))
    a2good = (np.logical_not(np.isnan(add2)))
    allgood = a1good*a2good
    msum[allgood] = add1[allgood] + add2[allgood]
    meas.addColumn(key,msum)

    
def channelName(key, chan):
    return '%s_%s' % (key,chan)


def channelColumn(meas,key,chan):
    fullkey = channelName(key,chan)
    return meas[fullkey]


def mergedMeasurementsRatios(mmeas, chan1, chan2, cal1, cal2):
    for chan, cal in zip([chan1,chan2],[cal1,cal2]):
#        if  channelName('qIndex',chan) not in mmeas.keys():
#            makeRatio(mmeas, channelName('qIndex',chan), 100.0, channelColumn(mmeas,'tau1',chan))
        if  channelName('qIndexC',chan) not in mmeas.keys():
            makeRatio(mmeas, channelName('qIndexC',chan), channelColumn(mmeas,'qIndex',chan), cal)
        if  (channelName('qDensity',chan) not in mmeas.keys()) and (channelName('area',chan) in mmeas.keys()):
            makeRatio(mmeas, channelName('qDensity',chan), channelColumn(mmeas,'qIndex',chan),
                      channelColumn(mmeas,'area',chan)) 
        if  (channelName('qDensityC',chan) not in mmeas.keys()) and (channelName('area',chan) in mmeas.keys()):
            makeRatio(mmeas, channelName('qDensityC',chan), channelColumn(mmeas,'qIndexC',chan),
                      channelColumn(mmeas,'area',chan)) 
    makeRatio(mmeas, channelName('qRatio','%svs%s' % (chan1,chan2)),
              channelColumn(mmeas,'qIndex',chan1),
              channelColumn(mmeas,'qIndex',chan2))
    makeRatio(mmeas, channelName('qRatioC','%svs%s' % (chan1,chan2)),
              channelColumn(mmeas,'qIndexC',chan1),
              channelColumn(mmeas,'qIndexC',chan2))



# darktime fitting section
from scipy.optimize import curve_fit
def cumuexpfit(t,tau):
    return 1-np.exp(-t/tau)

def cumumultiexpfit(t,tau1,tau2,a):
    return a*(1-np.exp(-t/tau1))+(1-a)*(1-np.exp(-t/tau2))

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
    hc = np.cumsum(h)/float(timeintervals.shape[0]) # normalise
    hcg = hc[h>0] # only nonzero bins
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
        try:
            tauEst = cumux[(np.abs(cumuy - 0.63)).argmin()]
        except ValueError:
            tauEst = 100.0
        # generate alternative histogram with binning
        binctrs, hc, binctrsg, hcg = cumuhistBinned(dtg)
        try:
            tauEstH = binctrsg[(np.abs(hcg - 0.63)).argmin()]
        except ValueError:
            tauEstH = 100.0

        success = True
        # fit theoretical distributions
        try:
            popth,pcovh,infodicth,errmsgh,ierrh  = curve_fit(cumuexpfit,binctrs,hc, p0=(tauEstH),full_output=True)
        except:
            success = False
        else:
            if hc.shape[0] > 1:
                chisqredh = ((hc - infodicth['fvec'])**2).sum()/(hc.shape[0]-1)
            else:
                chisqredh = 0
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
                ('chisqr1', 'f4'), ('chisqr2', 'f4'), ('tau1est', 'f4'), ('tau2est', 'f4'),
                ('NDefocused', 'i4'), ('NDefocusedFrac', 'f4')]


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


def measureObjectsByID(filter, ids, sigDefocused = None):
    # IMPORTANT: repeated filter access is extremely costly!
    # need to cache any filter access here first
    x = filter['x'] #+ 0.1*random.randn(filter['x'].size)
    y = filter['y'] #+ 0.1*random.randn(x.size)
    id = filter['objectID'].astype('i')
    t = filter['t']
    sig = filter['sig'] # we must do our own caching!

    measurements = np.zeros(len(ids), dtype=measureDType)

    for j,i in enumerate(ids):
        if not i == 0:
            if np.all(np.in1d(i, id)): # check if this ID is present in data
                ind = id == i
                obj = {'x': x[ind], 'y': y[ind], 't': t[ind]}
                #print obj.shape
                measure(obj, measurements[j])
                # here we measure the fraction of defocused localisations to give us an idea how 3D something is
                if sigDefocused is not None:
                    measurements[j]['NDefocused'] = np.sum(sig[ind] > sigDefocused)
                    measurements[j]['NDefocusedFrac'] = float(measurements[j]['NDefocused'])/measurements[j]['NEvents'] 
            else:
                for key in  measurements[j].dtype.fields.keys():
                    measurements[j][key]=0
            measurements[j]['objectID'] = i

    # wrap recarray in tabular that allows us to
    # easily add columns and save using tabular methods
    return TabularRecArrayWrap(measurements)


def retrieveMeasuresForIDs(measurements,idcolumn,columns=['tau1','NDarktimes','qIndex']):
    newcols = {key: np.zeros_like(idcolumn, dtype = 'float64') for key in columns}

    for j,i in enumerate(measurements['objectID']):
        if not i == 0:
            ind = idcolumn == i
            for col in newcols.keys():
                if not np.isnan(measurements[col][j]):
                    newcols[col][ind] = measurements[col][j]

    return newcols
