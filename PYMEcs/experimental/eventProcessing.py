import matplotlib.pyplot as plt
from PYME.warnings import warn

class EventProcessing:
    """
    plugins to conduct some event processing in h5r data
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

        if offupd.yvals.size > 0:
            plt.figure()
            plt.step(tminutes,offupd.yvals,where='post')
            plt.xlabel('time (minutes)')
            plt.ylabel('OffsetPiezo offset (um)')
            plt.title('OffsetPiezo offsets from PiezoOffsetUpdate events')

        offsets = piecewiseMapping.GeneratePMFromEventList(p.events, p.mdh, p.mdh.getEntry('StartTime'), 0, b'PiezoOffset',0)
        if offsets.yvals.size > 0:
            offsVSt = offsets(p['t']-0.01)
            plt.figure()
            plt.plot(p['t'],offsVSt)
            plt.xlabel('time (frame number)')
            plt.ylabel('OffsetPiezo offset (um)')
            plt.title('OffsetPiezo offsets from PiezoOffset events')

        if 'driftx' in p.keys():
            plt.figure()
            plt.subplot(311)
            plt.plot(p['t'],p['driftx'])
            plt.title('Drift in x (nm)')
            plt.subplot(312)
            plt.plot(p['t'],p['drifty'])
            plt.title('Drift in y (nm)')
            plt.subplot(313)
            plt.plot(p['t'],1e3*p['driftz']) # is driftz in um?
            plt.title('Drift in z (nm)')
            plt.tight_layout()

def Plug(visFr):
    """Plugs this module into the gui"""
    EventProcessing(visFr)
