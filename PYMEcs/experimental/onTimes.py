import wx
import numpy as np
import sys
from scipy import ndimage

import logging
logger = logging.getLogger(__file__)

def Warn(parent, message, caption = 'Warning!'):
    dlg = wx.MessageDialog(parent, message, caption, wx.OK | wx.ICON_WARNING)
    dlg.ShowModal()
    dlg.Destroy()

class onTimer:
    """

    """
    def __init__(self, visFr):
        self.visFr = visFr
        self.pipeline = visFr.pipeline

        visFr.AddMenuItem('Experimental>Event Processing', "onTimes from selected coalesced events",
                          self.OnTimes)


        
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

def Plug(visFr):
    """Plugs this module into the gui"""
    visFr.onTimer = onTimer(visFr)
