import matplotlib.pyplot as plt
from PYMEcs.pyme_warnings import warn

class EventProcessing:
    """
    plugins to conduct some event processing from events in h5r data
    currently mostly using events from Piezobased tracking
    """
    def __init__(self, visFr):
        self.visFr = visFr
        self.pipeline = visFr.pipeline

        visFr.AddMenuItem('Experimental',
                          'Display Piezo events in h5r file',
                          self.OnDisplayEvents,
                          helpText='display recorded events (in the PYMEAcquire sense) from an h5r file')

    def OnDisplayEvents(self, event=None):
        from PYME.Analysis import piecewiseMapping
        p = self.pipeline

        if p.events is None:
            warn('No events in pipeline')
            return
        
        offupd = piecewiseMapping.GeneratePMFromEventList(p.events, p.mdh, p.mdh.getEntry('StartTime'), 0, b'PiezoOffsetUpdate',0)
        tminutes = offupd.xvals * p.mdh['Camera.CycleTime'] / 60.0

        offsets = piecewiseMapping.GeneratePMFromEventList(p.events, p.mdh, p.mdh.getEntry('StartTime'), 0, b'PiezoOffset',0)
        correlamp = piecewiseMapping.GeneratePMFromEventList(p.events, p.mdh, p.mdh.getEntry('StartTime'), 0, b'CorrelationAmplitude',0)
        
        has_offsets = offsets.yvals.size > 0
        has_offupd = offupd.yvals.size > 0 # not using this one at the moment
        has_correlamp = correlamp.yvals.size > 0
        has_drift = 'driftx' in p.keys()

        plot_rows = 0
        if has_offsets:
            plot_rows += 1
        if has_correlamp:
            plot_rows += 1
        if has_drift:
            plot_rows += 3

        row = 1
        fig, axs = plt.subplots(nrows=plot_rows,figsize=(6.4,6.4), num='OffsetPiezo Event Analysis')
        if has_drift:
            plt.subplot(plot_rows,1,row)
            plt.plot(p['t'],p['driftx'])
            plt.title('Drift in x (nm)')
            plt.ylabel('Drift')
            plt.subplot(plot_rows,1,row+1)
            plt.plot(p['t'],p['drifty'])
            plt.title('Drift in y (nm)')
            plt.ylabel('Drift')
            plt.subplot(plot_rows,1,row+2)
            plt.plot(p['t'],1e3*p['driftz']) # is driftz in um?
            plt.title('Drift in z (nm)')
            plt.ylabel('Drift')
            row += 3

        if has_offsets:
            offsVSt = offsets(p['t']-0.01)
            plt.subplot(plot_rows,1,row)
            plt.plot(p['t'],offsVSt)
            plt.xlabel('time (frame number)')
            plt.ylabel('offset (um)')
            plt.title('OffsetPiezo offsets from PiezoOffset events')
            row += 1

        if has_correlamp:
            campVSt = correlamp(p['t']-0.01)
            plt.subplot(plot_rows,1,row)
            plt.plot(p['t'],campVSt)
            plt.xlabel('time (frame number)')
            plt.ylabel('amp')
            plt.title('normalised correlation amplitude')
            plt.ylim(0.3,1.2)
            axs[-1].set_yticks([0.75],minor=True)
            plt.grid(which='both',axis='y')
            row += 1

        plt.tight_layout()

            #plt.step(tminutes,offupd.yvals,where='post')
            #plt.xlabel('time (minutes)')
            #plt.ylabel('OffsetPiezo offset (um)')
            #plt.title('OffsetPiezo offsets from PiezoOffsetUpdate events')


def Plug(visFr):
    """Plugs this module into the gui"""
    EventProcessing(visFr)
