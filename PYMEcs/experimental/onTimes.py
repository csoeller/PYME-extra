import wx
import numpy as np
import sys
from scipy import ndimage
from PYMEcs.misc.guiMsgBoxes import Warn, Info

from traits.api import HasTraits, Str, Int, CStr, List, Enum, Float
from traitsui.api import View, Item, Group
from traitsui.menu import OKButton, CancelButton, OKCancelButtons


import logging
logger = logging.getLogger(__file__)

class propertyChoice(HasTraits):
    clist = List([])
    EventProperty = Enum(values='clist')

    traits_view = View(Group(Item(name = 'EventProperty'),
                             show_border = True),
                       buttons = OKCancelButtons)

    def add_channels(self,chans):
        for chan in chans:
            if chan not in self.clist:
                self.clist.append(chan)


class onTimer:
    """

    """
    def __init__(self, visFr):
        self.visFr = visFr
        self.pipeline = visFr.pipeline

        visFr.AddMenuItem('Experimental>Event Processing',
                          "onTimes from selected coalesced events",
                          self.OnTimes,
                          helpText='analyse the on time distribution of events in a region - needs the coalesced data source with nFrames event property')
        visFr.AddMenuItem('Experimental>Event Processing',
                          "plot time series of clumps",
                          self.OnPlotClumps,
                          helpText='plots the time course of the selected events in a "single channel trace" - needs the clumpIndex event property')
        visFr.AddMenuItem('Experimental>Event Processing',
                          "plot time series of event property",
                          self.OnPlotProperty,
                          helpText='plots a time series of the selected events in a "single channel trace" - using a user chosen event property')
        visFr.AddMenuItem('Experimental>Event Processing',
                          "plot time gating from selected events\tCtrl+G",
                          self.OnPlotTser,
                          helpText='plots the time series of gating of detected molecules from the selected group of events')
        visFr.AddMenuItem('Experimental>Event Processing',
                          "Event density in ROI\tCtrl+D",
                          self.OnEvtDensity,
                          helpText='calculate the event density (events/unit area) in the current ROI/image area')


    def OnEvtDensity(self,event):
        from PYMEcs.Analysis.eventProperties import evtDensity, getarea
        visFr = self.visFr
        pipeline = visFr.pipeline

        dens, intens1, intens2, trange = evtDensity(pipeline)
        area = getarea(pipeline)

        if dens is None:
            Warn(visFr,'area too small (%.2f um^2)' % (area))
        else:
            infostr = \
"""
Event density:        %.1f events/um^2
Norm. intens.: %.1f events/um^2/5K-fr
Norm. intens.: %.1f events/(20um)^2/fr
Number of Events:     %d events
Area probed:          %.2f um^2
Time range:           %.2fK frames""" % (dens, intens1, intens2, dens*area, area, trange/1e3)
            Info(visFr,infostr)


    def OnPlotTser(self,event):
        from PYMEcs.misc.shellutils import plotserpipeline

        visFr = self.visFr
        pipeline = visFr.pipeline
        mdh = pipeline.mdh

        t = pipeline['t']
        maxPts = 1e4
        if len(t) > maxPts:
            Warn(None,'aborting darktime analysis: too many events, current max is %d' % maxPts)
            return

        p = pipeline

        # if we have coalesced events use this info
        if ('tmin' in pipeline.keys()) and ('tmax' in pipeline.keys()):
            tmin = pipeline['tmin']
            tmax = pipeline['tmax']
            tc = np.arange(tmin[0],tmax[0]+1)
            for i in range(1,tmin.shape[0]):
                tc = np.append(tc,np.arange(tmin[i],tmax[i]+1))
            tc.sort()
        else:
            tc = t
        
        startn = pipeline.selectedDataSource['t'].min()
        endn = pipeline.selectedDataSource['t'].max()

        tt, v = plotserpipeline(tc, np.ones_like(tc))

        # add points to start from beginning of series
        tt2 = np.append([startn,tt[0],tt[0]],tt)
        v2 = np.append([0,0,1],v)
        # add points to go to end of series
        tt2 = np.append(tt2,[tt[-1],tt[-1],endn])
        v2 = np.append(v2,[1,0,0])
        
        import matplotlib.pyplot as plt
        # plot data and fitted curves
        plt.figure()
        plt.plot(tt2,v2)
        plt.show()
        plt.xlim(startn,endn)
        
    def OnTimes(self,event):
        import StringIO
        from PYMEcs.Analysis import fitDarkTimes
        
        visFr = self.visFr
        pipeline = visFr.pipeline
        mdh = pipeline.mdh

        NTMIN = 5
        maxPts = 1e5
        t = pipeline['t']
        if len(t) > maxPts:
            Warn(None,'aborting analysis: too many events, current max is %d' % maxPts)
            return
        x = pipeline['x']
        y = pipeline['y']

        # determine darktime from gaps and reject zeros (no real gaps) 

        try:
            dtg = pipeline['nFrames']
        except KeyError:
            Warn(None,'aborting analysis: could not find "nFrames" property, needs "Coalesced" Data Source','Error')
            return

        nts = dtg.shape[0]

        if nts > NTMIN:
            # now make a cumulative histogram from these
            cumux,cumuy = fitDarkTimes.cumuhist(dtg)
            binctrs,hc,binctrsg,hcg = fitDarkTimes.cumuhistBinned(dtg)
            
            bbx = (x.min(),x.max())
            bby = (y.min(),y.max())
            voxx = 1e3*mdh['voxelsize.x']
            voxy = 1e3*mdh['voxelsize.y']
            bbszx = bbx[1]-bbx[0]
            bbszy = bby[1]-bby[0]

        # fit theoretical distributions
        popth,pcovh,popt,pcov = (None,None,None,None)
        if nts > NTMIN:
            from scipy.optimize import curve_fit

            idx = (np.abs(hcg - 0.63)).argmin()
            tauesth = binctrsg[idx]
            popth,pcovh,infodicth,errmsgh,ierrh = curve_fit(fitDarkTimes.cumuexpfit,binctrs,hc, p0=(tauesth),full_output=True)
            chisqredh = ((hc - infodicth['fvec'])**2).sum()/(hc.shape[0]-1)
            idx = (np.abs(cumuy - 0.63)).argmin()
            tauest = cumux[idx]
            popt,pcov,infodict,errmsg,ierr = curve_fit(fitDarkTimes.cumuexpfit,cumux,cumuy, p0=(tauest),full_output=True)
            chisqred = ((cumuy - infodict['fvec'])**2).sum()/(nts-1)

            import matplotlib.pyplot as plt
            # plot data and fitted curves
            plt.figure()
            plt.subplot(211)
            plt.plot(cumux,cumuy,'o')
            plt.plot(cumux,fitDarkTimes.cumuexpfit(cumux,popt[0]))
            plt.plot(binctrs,hc,'o')
            plt.plot(binctrs,fitDarkTimes.cumuexpfit(binctrs,popth[0]))
            plt.ylim(-0.2,1.2)
            plt.subplot(212)
            plt.semilogx(cumux,cumuy,'o')
            plt.semilogx(cumux,fitDarkTimes.cumuexpfit(cumux,popt[0]))
            plt.semilogx(binctrs,hc,'o')
            plt.semilogx(binctrs,fitDarkTimes.cumuexpfit(binctrs,popth[0]))
            plt.ylim(-0.2,1.2)
            plt.show()
            
            outstr = StringIO.StringIO()

            analysis = {
                'Nevents' : t.shape[0],
                'Ndarktimes' : nts,
                'filterKeys' : pipeline.filterKeys.copy(),
                'darktimes' : (popt[0],popth[0]),
                'darktimeErrors' : (np.sqrt(pcov[0][0]),np.sqrt(pcovh[0][0]))
            }

            #if not hasattr(self.visFr,'analysisrecord'):
            #    self.visFr.analysisrecord = []
            #    self.visFr.analysisrecord.append(analysis)

            print >>outstr, "events: %d, on times: %d" % (t.shape[0],nts)
            print >>outstr, "region: %d x %d nm (%d x %d pixel)" % (bbszx,bbszy,bbszx/voxx,bbszy/voxy)
            print >>outstr, "centered at %d,%d (%d,%d pixels)" % (x.mean(),y.mean(),x.mean()/voxx,y.mean()/voxy)
            print >>outstr, "ontime: %.1f+-%d (%.1f+-%d) frames - chisqr %.2f (%.2f)" % (popt[0],np.sqrt(pcov[0][0]),
                                                                                           popth[0],np.sqrt(pcovh[0][0]),
                                                                                           chisqred,chisqredh)
            print >>outstr, "ontime: starting estimates: %.1f (%.1f)" % (tauest,tauesth)
            print >>outstr, "qunits: %.2f (%.2f), mean on time: %.2f" % (100.0/popt[0], 100.0/popth[0], dtg.mean())

            labelstr = outstr.getvalue()
            plt.annotate(labelstr, xy=(0.5, 0.1), xycoords='axes fraction',
                         fontsize=10)
        else:
            Warn(None, 'not enough data points, only found %d on times (need at least  %d)' % (nts,NTMIN), 'Error')

    def OnPlotClumps(self,event):
        import PYMEcs.Analysis.timeSeries as ts
        try:
            ci = self.pipeline['clumpIndex']
        except KeyError:
            Warn(None,'aborting analysis: could not find "clumpIndex" property','Error')
            return
        t = self.pipeline['t']

        maxPts = 1e4
        if len(t) > maxPts:
            Warn(None,'aborting analysis: too many events, current max is %d' % maxPts)
            return

        ts.plotClumpSeries(t, ci)

    def tserFromPipeline(self,key,base=0):
    # t is on integer times assumed (from pipeline)
    # val is the corresponding value
        t = self.pipeline['t']
        val = self.pipeline[key]
    
        tmin = t.min()
        tmax = t.max()
        tstart = tmin-1
        tend = tmax + 1
        tlen = tend-tstart+1

        to = t-tstart
        tidx = np.argsort(to)
        tos = to[tidx]
        vs = val[tidx]

        tvalid = np.zeros((tlen))
        tvalid[tos] = 1
    
        vals = np.zeros((tlen))
        vals[tos] = vs
    
        tup = tos[1:][(tos[1:]-tos[:-1] > 1)]
        tvalid[tup-1] = 1
        vals[tup-1] = base

        tdown = tos[:-1][(tos[1:]-tos[:-1] > 1)]
        tvalid[tdown+1] = 1
        vals[tdown+1] = base

        tplot = tvalid.nonzero()[0]
        vplot = vals[tplot]
        
        return (tplot+tstart,vplot)


    def OnPlotProperty(self,event):
        pChoice = propertyChoice()
        pChoice.add_channels(sorted(self.pipeline.keys()))
        if pChoice.configure_traits(kind='modal'):
            t,prop = self.tserFromPipeline(pChoice.EventProperty)
            if t.shape[0] > 2e4:
                 Warn(None,'aborting analysis: too many events, current max is %d' % 2e4)
                 return
            import matplotlib.pyplot as plt
            plt.figure()
            plt.plot(t,prop)
            plt.xlabel('Frame Number')
            plt.ylabel(pChoice.EventProperty)
            plt.show()
            
def Plug(visFr):
    """Plugs this module into the gui"""
    visFr.onTimer = onTimer(visFr)
