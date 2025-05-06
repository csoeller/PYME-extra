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
    def __init__(self, recarray, validCols = None):
        self._recarray = recarray
        self.new_columns = {}
        if validCols is not None:
            self._recarrayKeys = validCols
        else:
            self._recarrayKeys = self._recarray.dtype.fields.keys()

    def keys(self):
        return list(set(list(self._recarrayKeys + self.new_columns.keys())))

    def __getitem__(self, keys):
        key, sl = self._getKeySlice(keys)
        if key in self._recarrayKeys:
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

def mkcmexpfit(tau2):
    def fitfunc(t,tau1,a):
        return a*(1-np.exp(-t/tau1))+(1-a)*(1-np.exp(-t/tau2))
    return fitfunc

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
    cumuy = (1.0+np.arange(nIntervals))/float(nIntervals)
    return (cumux,cumuy)


def cumuhistBinned(timeintervals):
    binedges = 0.5+np.arange(0,timeintervals.max())
    binctrs = 0.5*(binedges[0:-1]+binedges[1:])
    h,be2 = np.histogram(timeintervals,bins=binedges)
    hc = np.cumsum(h)/float(timeintervals.shape[0]) # normalise
    hcg = hc[h>0] # only nonzero bins
    binctrsg = binctrs[h>0]

    return (binctrs, hc, binctrsg, hcg)

import math
def cumuhistBinnedLog(timeintervals,dlog=0.1,return_hist=False, return_good=False):
    binmax = int((math.log10(timeintervals.max())-1.0-dlog)/dlog+2.0)
    binedges = np.append(0.5+np.arange(10), 10.0**(1.0+dlog*(np.arange(binmax)+1.0)))
    binctrs = 0.5*(binedges[0:-1]+binedges[1:])
    h,be2 = np.histogram(timeintervals,bins=binedges)
    hc = np.cumsum(h)/float(timeintervals.shape[0]) # normalise
    hcg = hc[h>0] # only nonzero bins
    binctrsg = binctrs[h>0]

    retvals = [binctrs, hc]
    if return_good:
       retvals = retvals + [binctrsg, hcg]
    if return_hist:
        retvals = retvals + [h/float(timeintervals.shape[0])]

    return retvals

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


from traits.api import HasTraits, Str, Int, CStr, List, Enum, Float, Bool
class FitSettings(HasTraits):
    coalescedProcessing = Enum(['useClumpIndexOnly','useTminTmaxIfAvailable'])
    cumulativeDistribution = Enum(['binned','empirical'])
    fitMode = Enum(['SingleMode','TwoModes'])
    Tau2Constant = Bool(False)
    Tau2FixedValue = Float(2.0)
    IDcolumn = CStr('objectID')

    
measureDType2 = [('objectID', 'i4'),
                 ('NEvents', 'i4'), ('NEventsCorr', 'i4'), ('NDarktimes', 'i4'),
                 ('tau1', 'f4'), ('tau2', 'f4'), ('tau1err', 'f4'), ('tau2err', 'f4'),
                 ('amp1','f4'), ('amp1err','f4'),
                 ('chisqr', 'f4'), ('tau1est', 'f4'), ('qindex', 'f4')]

def measureObjectsByID2(datasource, idname='objectID', settings=FitSettings()):
    # Probably an old statement - check: IMPORTANT: repeated filter access is extremely costly!
    ds = datasource
    idDs = ds[idname].astype('i')

    idall = np.unique(ds[idname].astype('int'))
    ids = idall[idall > 0] # only accept nonzero IDs
    
    meas = np.zeros(ids.size, dtype=measureDType2)
    darkTimes = []
    times = []

    # check which type of time processing we are going to use
    if ('clumpIndex' in ds.keys()) and not ('fitError_x0' in ds.keys()): # heuristic to only do on coalesced data
        usingClumpIndex = True
        if (settings.coalescedProcessing == 'useTminTmaxIfAvailable') and  ('tmin' in ds.keys()) and ('tmax' in ds.keys()):
            usingTminTmax = True
        else:
            usingTminTmax = False
    else:
        usingClumpIndex = False
        usingTminTmax = False

    fields = ['NEvents','NDarktimes', 'qindex']
    if usingTminTmax:
        fields.append('NEventsCorr')

    if settings.Tau2FixedValue:
        tau2const = settings.Tau2Constant
    else:
        tau2const = 0.0
        
    if settings.fitMode == 'SingleMode':
        #  retv = [tau1, tau1err, chisqr, tau1est]
        mfields = ['tau1','tau1err','chisqr','tau1est']
    else:
        # retv = [tau1, tau2, atau1, tau1err, tau2err, atau1err, chisqr, tauest]
        mfields = ['tau1','tau2','amp1','tau1err','tau2err','amp1err','chisqr','tau1est']
        
    validCols = fields + mfields

    # loop over object IDs
    ndtmin = 5
    for j,this_id in enumerate(ids):
        if not this_id == 0:
            if np.all(np.in1d(this_id, idDs)): # check if this ID is present in data
                idx = idDs == this_id
                # stuff to be done in the innermost loop
                te, dte, nevents, nevtscorr = extractEventTimes(ds,idx,useTminTmax=usingTminTmax)
                meas[j]['NEvents'] = nevents
                meas[j]['NDarktimes'] = dte.size
                darkTimes.append(dte)
                times.append(te)
                if usingTminTmax:
                    meas[j]['NEventsCorr'] = nevtscorr
                if dte.size >= ndtmin:
                    if settings.cumulativeDistribution == 'binned':
                        xc, yc = cumuhistBinnedLog(dte,dlog=0.05)
                    else:
                        xc, yc = cumuhist(dte)

                    try:
                        retv = fitTaus(xc,yc,fitTau2 = (settings.fitMode == 'TwoModes'), tau2const=tau2const, return_tau1est=True)
                    except RuntimeError:
                        pass # we got a convergence error
                    else:
                        for i, field in enumerate(mfields):
                            meas[j][field] = retv[i]
                        meas[j]['qindex'] = 100.0/meas[j]['tau1']
            else:
                for key in  meas[j].dtype.fields.keys():
                    meas[j][key]=0
            meas[j]['objectID'] = this_id

    # wrap recarray in tabular that allows us to
    # easily add columns and save using tabular methods
    rmeas = TabularRecArrayWrap(meas, validCols=validCols+['objectID'])
    
    return {'measures': rmeas,
            'darkTimes' : darkTimes,
            'times' : times,
            'validColumns': validCols,
            'state' : {
                'usingClumpIndex': usingClumpIndex,
                'usingTminTmax': usingTminTmax,
                'IDcolumn': settings.IDcolumn,
                'coalescedProcessing': settings.coalescedProcessing,
                'Tau2Constant': settings.Tau2Constant,
                'Tau2FixedValue': settings.Tau2FixedValue,
                'FitMode' : settings.fitMode
            }}

def retrieveMeasuresForIDs2(measurements,idcol):
    validCols = measurements['validColumns']
    measures = measurements['measures']
    
    newcols = {key: np.zeros_like(idcol, dtype = measures._recarray.dtype.fields[key][0]) for key in validCols}

    for j,id in enumerate(measures['objectID']):
        if not id == 0:
            ind = idcol == id
            for col in newcols.keys():
                if not np.isnan(measures[col][j]):
                    newcols[col][ind] = measures[col][j]

    return newcols

def extractEventTimes(datasource, idx = None, useTminTmax = True, return_modes = False):
    d = datasource
    t = d['t']
    if idx is None:
        idx = np.ones_like(t,dtype='bool')
    ti = t[idx]
    
    nevents_corrected = None
    # if we have coalesced events use this info
    if ('clumpIndex' in d.keys()) and not ('fitError_x0' in d.keys()): # heuristic to only do on coalesced data
        usingClumpIndex = True
        csz = d['clumpSize'][idx]
        nevents = csz.sum()
        if useTminTmax and  ('tmin' in d.keys()) and ('tmax' in d.keys()):
            tmin = d['tmin'][idx]
            tmax = d['tmax'][idx]
            tc = np.arange(tmin[0],tmax[0]+1)
            for i in range(1,tmin.size):
                tc = np.append(tc,np.arange(tmin[i],tmax[i]+1))
            tc.sort()
            usingTminTmax = True
            nevents_corrected = tc.size
        else:
            tc = np.arange(int(ti[0]-csz[0]/2),int(ti[0]+csz[0]/2))
            for i in range(1,ti.size):
                tc = np.append(tc,np.arange(int(t[i]-csz[i]/2),int(t[i]+csz[i]/2)))
            tc.sort()
            usingTminTmax = False
    else:
        tc = ti
        usingTminTmax = False
        usingClumpIndex = False
        nevents = tc.size

    dts = tc[1:]-tc[0:-1]-1
    dtg = dts[dts>0]
    
    if return_modes:
        return (tc, dtg, nevents, nevents_corrected, usingClumpIndex, usingTminTmax)
    else:
        return (tc, dtg, nevents, nevents_corrected)

   
def fitTaus(x_t, y_h, fitTau2 = False, tau2const = 0.0, return_tau1est = False, tau2max=8.0):

    # could be refined by subtracting off the histogram for values around 9 frames or so
    # and then ask for reaching 63% off the remaining difference to 1
    idx = (np.abs(y_h - 0.63)).argmin()
    tau1est = x_t[idx]
    
    # further possibilities:
    # use tau2 but keep it fixed
    # add bounds on atau1 (between 0..1) and tau2 (0..8)

    if fitTau2:
        popt,pcov = curve_fit(cumumultiexpfit,x_t,y_h, p0=(tau1est,2.0,0.8),bounds=(0, (np.inf, tau2max, 1.0)))
        (tau1, tau2, atau1) = popt
        (tau1err, tau2err, atau1err) = np.sqrt(np.diag(pcov))
        chisqr = ((y_h - cumumultiexpfit(x_t,*popt))**2).sum()/(x_t.size-1)
        results = [tau1, tau2, atau1, tau1err, tau2err, atau1err, chisqr]
    else:
        if tau2const < 1e-4:
            popt,pcov = curve_fit(cumuexpfit,x_t,y_h, p0=(tau1est))
            (tau1,tau1err) = (popt[0],np.sqrt(pcov[0][0]))
            chisqr = ((y_h - cumuexpfit(x_t,*popt))**2).sum()/(x_t.size-1)
            results = [tau1, tau1err, chisqr]
        else:
            popt,pcov = curve_fit(mkcmexpfit(tau2const),x_t,y_h, p0=(tau1est,0.8),bounds=(0, (np.inf, 1.0)))
            (tau1, atau1) = popt
            (tau1err, atau1err) = np.sqrt(np.diag(pcov))
            (tau2,tau2err) = (tau2const,0)
            chisqr = ((y_h - cumumultiexpfit(x_t, tau1, tau2, atau1))**2).sum()/(x_t.size-1)
            results = [tau1, tau2, atau1, tau1err, tau2err, atau1err, chisqr]
            
    if return_tau1est:
        results.append(tau1est)

    return results

